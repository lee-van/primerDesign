# designPrimers.py
# v1.2

import sys
from Bio import SeqIO
import re
import argparse
import subprocess
import os
import shutil
import datetime
import re
from distutils import spawn


PRIMER3=spawn.find_executable("/Data05/evlee/packages/primer3-2.3.7/src/primer3_core") #change once installed globally
if PRIMER3 is None:
   print("***ERROR: primer3_core is not found")
   sys.exit("Please install primer3_core or make sure it is in the PATH")

PYTHON2=spawn.find_executable("python2.7")
if PYTHON2 is None:
   print("***ERROR: python2.7 is not found")
   sys.exit("Please install / python2.7 or make sure it is in the PATH")

MFEprimer="/Data05/evlee/packages/MFEprimer-v2.0/MFEprimer.py" #change once installed globally

# command line arguments
parser = argparse.ArgumentParser(description="Takes a bam file that has been sorted with redundant reads removed and generates a HAMR predicted_mods.txt output")
parser.add_argument('inFasta',help='Per-transcript or per-gene fasta')
parser.add_argument('inMFEprimerDB',help='MFEprimer database prefix')
parser.add_argument('output',help='tabular output file of designed primers')
parser.add_argument('--inBED', '-b', action='store', dest='inBED', nargs='?', help='bed file of target transcripts or genes')
parser.add_argument('--outDir', '-o', action='store', dest='outDir', nargs='?', help='retain intermediatre files to an output folder')
parser.add_argument('--silent','-s',action='store_true',help='Suppress output messages')
parser.add_argument('--timeout','-t',action='store', dest='timeout', type=int, default=60, help='(Seconds) Skip primer3 or MFEprimer run if hanging longer than this time')

args=parser.parse_args()

# make tmp directory
if args.outDir:
	tmpDIR = args.outDir
else:
	now = datetime.datetime.now()
	datelist = [str(now.year),str(now.month),str(now.day),str(now.hour),str(now.minute),str(now.second),str(now.microsecond)]
	rightnow= "_".join(datelist)
	tmpDIR= "temp."+rightnow
subprocess.check_call(['mkdir', '-p', tmpDIR])


#read through gene list bed
geneList=[]
if args.inBED:
	openBED = open(args.inBED,'r')
	for i in openBED:
		j=i.rstrip().split("\t")
		geneList.append(j[3])
	openBED.close()

#write primer3 input files (boulder I/O) for each target
Primer3InputFiles = []
for record in SeqIO.parse(args.inFasta, "fasta"):
	if record.id in geneList or not args.inBED:
		gene = re.sub(r'\.\d+$', '', record.id)
		Primer3InputFile = tmpDIR + "/" + record.id + ".io"
		Primer3InputFiles.append(Primer3InputFile)
		Primer3InputFileFH = open(Primer3InputFile, "w")
		Primer3InputFileFH.write("SEQUENCE_ID="+str(record.id)+"\n")
		Primer3InputFileFH.write("SEQUENCE_TEMPLATE="+str(record.seq)+"\n")
		Primer3InputFileFH.write("PRIMER_TASK=pick_pcr_primers\n")
		Primer3InputFileFH.write("PRIMER_PICK_LEFT_PRIMER=1\n")
		Primer3InputFileFH.write("PRIMER_PICK_INTERNAL_OLIGO=0\n")
		Primer3InputFileFH.write("PRIMER_PICK_RIGHT_PRIMER=1\n")
		Primer3InputFileFH.write("PRIMER_OPT_SIZE=20\n")
		Primer3InputFileFH.write("PRIMER_MIN_SIZE=18\n")
		Primer3InputFileFH.write("PRIMER_MAX_SIZE=25\n")
		Primer3InputFileFH.write("PRIMER_MAX_NS_ACCEPTED=0\n")
		Primer3InputFileFH.write("PRIMER_PRODUCT_SIZE_RANGE=75-300\n")
		Primer3InputFileFH.write("P3_FILE_FLAG=0\n")
		Primer3InputFileFH.write("PRIMER_EXPLAIN_FLAG=1\n")
		Primer3InputFileFH.write("PRIMER_THERMODYNAMIC_PARAMETERS_PATH=/Data05/evlee/packages/primer3-2.3.7/src/primer3_config/\n")
		Primer3InputFileFH.write("=")
		Primer3InputFileFH.close()

#run primer3
Primer3ResultsFiles = []
Primer3TimedOut = []
for Primer3InputFile in Primer3InputFiles:
	if not args.silent:
		target = Primer3InputFile.replace(tmpDIR + "/", "").replace(".io", "")
		print ("Designing primers for target: "+target)
	Primer3ResultsFile = Primer3InputFile.replace(".io", ".p3")
	Primer3ResultsFileFH = open(Primer3ResultsFile, "w")
	try:
		subprocess.run([PRIMER3, Primer3InputFile], stdout=Primer3ResultsFileFH, timeout=args.timeout)
	except subprocess.TimeoutExpired:
		Primer3TimedOut.append(target)
		continue
	Primer3ResultsFiles.append(Primer3ResultsFile)
	Primer3ResultsFileFH.close()

#parse primer3 output. For each primer pair, run MFEprimer and parse output to tabular file
outputFH = open(args.output, 'w')
outputFH.write("target" + "\t" + "pair" + "\t" + "left" + "\t" + "right" + "\t" + "left.Tm" + "\t" + "right.Tm" + "\t" + "ampliconSize" + "\t" + "predictedTargets" + "\t" + "dimerTm" + "\t" + "left.selfTm" + "\t" + "right.selfTm" + "\t" + "left.hairpinTm" + "\t" + "right.hairpintTm" + "\n")

MFEprimerTimedOut = []
for Primer3ResultsFile in Primer3ResultsFiles:
	Primer3ResultsFileFH = open(Primer3ResultsFile,'r')
	rightPrimer = []
	leftPrimer = []
	penalty = []
	leftTm = []
	rightTm = []
	leftSelf = []
	rightSelf = []
	leftHairpin = []
	rightHairin = []
	complementarity = []
	size = []
	target = ""
	current_primer_pair = ""
	for line in Primer3ResultsFileFH:
		m = re.search(r'SEQUENCE_ID=(.+?)$', line)
		if m:
			target = m.group(1)
		m = re.search(r'PRIMER_PAIR_(\d+)_PENALTY=(.+?)$', line)
		if m:
			current_primer_pair = m.group(1)
			penalty.append(m.group(2))
		m = re.search(r'PRIMER_LEFT_(\d+)_SEQUENCE=(.+?)$', line)
		if m:
			primer_pair = m.group(1)
			if primer_pair != current_primer_pair:
				sys.exit("Primer 3 output not properly sorted")
			leftPrimer.append(m.group(2))
		m = re.search(r'PRIMER_RIGHT_(\d+)_SEQUENCE=(.+?)$', line)
		if m:
			primer_pair = m.group(1)
			if primer_pair != current_primer_pair:
				sys.exit("Primer 3 output not properly sorted")		
			rightPrimer.append(m.group(2))
		# Tm
		m = re.search(r'PRIMER_LEFT_(\d+)_TM=(.+?)$', line)
		if m:
			primer_pair = m.group(1)
			if primer_pair != current_primer_pair:
				sys.exit("Primer 3 output not properly sorted")		
			leftTm.append(m.group(2))
		m = re.search(r'PRIMER_RIGHT_(\d+)_TM=(.+?)$', line)
		if m:
			primer_pair = m.group(1)
			if primer_pair != current_primer_pair:
				sys.exit("Primer 3 output not properly sorted")		
			rightTm.append(m.group(2))
		# QC thermodynamics
		m = re.search(r'PRIMER_LEFT_(\d+)_SELF_ANY_TH=(.+?)$', line)
		if m:
			primer_pair = m.group(1)
			if primer_pair != current_primer_pair:
				sys.exit("Primer 3 output not properly sorted")		
			leftSelf.append(m.group(2))
		m = re.search(r'PRIMER_RIGHT_(\d+)_SELF_ANY_TH=(.+?)$', line)
		if m:
			primer_pair = m.group(1)
			if primer_pair != current_primer_pair:
				sys.exit("Primer 3 output not properly sorted")		
			rightSelf.append(m.group(2))
		m = re.search(r'PRIMER_LEFT_(\d+)_HAIRPIN_TH=(.+?)$', line)
		if m:
			primer_pair = m.group(1)
			if primer_pair != current_primer_pair:
				sys.exit("Primer 3 output not properly sorted")		
			leftHairpin.append(m.group(2))
		m = re.search(r'PRIMER_RIGHT_(\d+)_HAIRPIN_TH=(.+?)$', line)
		if m:
			primer_pair = m.group(1)
			if primer_pair != current_primer_pair:
				sys.exit("Primer 3 output not properly sorted")		
			rightHairin.append(m.group(2))
		m = re.search(r'PRIMER_PAIR_(\d+)_COMPL_ANY_TH=(.+?)$', line)
		if m:
			primer_pair = m.group(1)
			if primer_pair != current_primer_pair:
				sys.exit("Primer 3 output not properly sorted")		
			complementarity.append(m.group(2))
		#size
		m = re.search(r'PRIMER_PAIR_(\d+)_PRODUCT_SIZE=(.+?)$', line)
		if m:
			primer_pair = m.group(1)
			if primer_pair != current_primer_pair:
				sys.exit("Primer 3 output not properly sorted")		
			size.append(m.group(2))

	Primer3ResultsFileFH.close()

	#check that all attributes captured
	if len(set([len(rightPrimer), len(leftPrimer), len(penalty), len(leftTm), len(rightTm), len(size)])) != 1:
		sys.exit("Failed to capture all attributes for each primer pair")

	#For each primer pair, run MFEprimer. Parse output to tabular file
	if not args.silent:
		print("analyzing potential off-targets for: "+target)
	PrimerPairsFile = Primer3ResultsFile.replace(".p3", ".primerPairs.txt")
	PrimerPairsFileFH = open(PrimerPairsFile, "w")
	PrimerPairsFileFH.write("target" + "\t" + "pair" + "\t" + "left" + "\t" + "right" + "\t" + "left.Tm" + "\t" + "right.Tm" + "\t" + "ampliconSize" + "\t" + "predictedTargets" + "\t" + "dimerPenalty" + "\t" + "left.selfTm" + "\t" + "right.selfTm" + "\t" + "left.hairpinTm" + "\t" + "right.hairpintTm" + "\n")
	for index in range(0, len(rightPrimer)):
		pairFasta = tmpDIR + "/" + target + ".pair" + str(index) + ".fa"
		pairFastaFH = open(pairFasta, "w")
		pairFastaFH.write(">"+target+" Penalty="+penalty[index]+" Size="+size[index]+" Tm="+rightTm[index]+"\n")
		pairFastaFH.write(rightPrimer[index]+"\n")
		pairFastaFH.write(">"+target+" Penalty="+penalty[index]+" Size="+size[index]+" Tm="+leftTm[index]+"\n")
		pairFastaFH.write(leftPrimer[index]+"\n")
		pairFastaFH.close()
		MFEprimerResultsFile = tmpDIR + "/" + target + ".pair" + str(index) + ".MFEprimer.txt"
		MFEprimerResultsFileFH = open(MFEprimerResultsFile, 'w')
		try:
			subprocess.run([PYTHON2, MFEprimer, "-i", pairFasta, "-d", args.inMFEprimerDB, "--tab"], stdout=MFEprimerResultsFileFH, timeout=args.timeout)
		except subprocess.TimeoutExpired:
			MFEprimerTimedOut.append(target + ".pair" + str(index))
			continue
		MFEprimerResultsFileFH.close()
		MFEprimerResultsFileFH = open(MFEprimerResultsFile, 'r')
		targets = []
		header = MFEprimerResultsFileFH.readline()
		for line in MFEprimerResultsFileFH:
			(AmpID, FpID, RpID, HitID, PPC, Size, AmpGC, FpTm, RpTm, FpDg, RpDg, BindingStart, BindingStop, AmpSeq) = line.rstrip().split("\t")
			targets.append(HitID + " (" + PPC + ")")
		targetsOutput = ";".join(targets)
		PrimerPairsFileFH.write(target + "\t" + str(index) + "\t" + rightPrimer[index] + "\t" + leftPrimer
			[index] + "\t" + rightTm[index] + "\t" + leftTm[index] + "\t" + size[index] + "\t" + targetsOutput + "\n")
		outputFH.write(target + "\t" + str(index) + "\t" + leftPrimer[index] + "\t" + rightPrimer[index] + "\t" + leftTm[index] + "\t" + rightTm[index] + "\t" + size[index] + "\t" + targetsOutput + "\t" + complementarity[index] + "\t" + leftSelf[index] + "\t" + rightSelf[index] + "\t" + leftHairpin[index] + "\t" + rightHairin[index]+ "\n")

if not args.outDir:
	shutil.rmtree(tmpDIR)

## write TimedOut processes to log file
print (Primer3TimedOut)
logFH = open(args.output+".log", 'w')
logFH.write("###unable to run Primer3 in "+str(args.timeout)+" seconds for the the following targets:\n")
for target in Primer3TimedOut:
	logFH.write(target+"\n")
logFH.write("###unable to run MFEprimer in "+str(args.timeout)+" seconds for the the following primer pairs:\n")
for target in MFEprimerTimedOut:
	logFH.write(target+"\n")
logFH.close()








