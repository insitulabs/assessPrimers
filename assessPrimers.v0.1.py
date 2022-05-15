#!/usr/bin/env python
import os
import sys
import argparse
import math
import subprocess

global MAKEBLASTDB
global TBLASTN
global CDHIT
global MAFFT

ALPHABET = {'A','C','G','T'}
DEG_CODE = {
	'A': 'A',
	'C': 'C',
	'G': 'G',
	'T': 'T',
	'AG': 'R',
	'CT': 'Y',
	'CG': 'S',
	'AT': 'W',
	'GT': 'K',
	'AC': 'M',
	'CGT': 'B',
	'AGT': 'D',
	'ACT': 'H',
	'ACG': 'V',
	'ACGT': 'N',
	'R': 'AG',
	'Y': 'CT',
	'S': 'CG',
	'W': 'AT',
	'K': 'GT',
	'M': 'AC',
	'B': 'CGT',
	'D': 'AGT',
	'H': 'ACT',
	'V': 'ACG',
	'N': 'ACGT',
	'I': 'ACGT'
}

RC_CODE = {
	'A': 'T',
	'C': 'G',
	'T': 'A',
	'G': 'C',
	'R': 'Y',
	'Y': 'R',
	'S': 'S',
	'W': 'W',
	'K': 'M',
	'M': 'K',
	'B': 'V',
	'D': 'H',
	'H': 'D',
	'V': 'B',
	'N': 'N'
}

##TODO: better handling of Ns?	
def shannonEntropy(seq):
	seq_nogap = seq.replace('-','')
	
	H = 0
	for letter in ALPHABET:
		p = seq_nogap.count(letter)	* 1.0 / len(seq_nogap)
		#print p
		if p > 0:
			H = H + (p * math.log(p, 4))
	if H != 0: H = H * -1

	return(H)


def scanEntropy(faFile, outPrefix, output, PRIMER_NAMES, kmer):
	seqNamesAll = []
	seqNamesNoPrimers = []
	seqs = []
	cols = ""

	curSeq = ""
	
	fa = open(faFile, 'r')
	##read and organize data
	while True:

		line = fa.readline().strip()

		if line == "": 
			if curSeq != "":
				curSeq = curSeq.upper() 
				##add previous sequence
				seqs.append(curSeq)

				if seqNamesAll[-1] not in PRIMER_NAMES:
					if len(cols) == 0: cols = list(curSeq)
					else:
						for i in range(len(curSeq)):
							cols[i] = cols[i] + curSeq[i]
			break

		if line.startswith(">"):
			if curSeq != "":
				curSeq = curSeq.upper() 
				##add previous sequence
				seqs.append(curSeq)

				if seqNamesAll[-1] not in PRIMER_NAMES:
					if len(cols) == 0: cols = list(curSeq)
					else:
						for i in range(len(curSeq)):
							cols[i] = cols[i] + curSeq[i]

			curSeq = ""

			
			name = line[1:]
			
			seqNamesAll.append(name)
			if name not in PRIMER_NAMES: seqNamesNoPrimers.append(name)

		else:
			curSeq = curSeq + line
	
	seqStats = open(outPrefix+".05.sequence_statistics.txt", 'w')
	seqStats.write("#sequence\tprimer\tnumMismatches\n")	
	##get primer statistics
	coordToPrimer = {}
	primerToCoord = {}
	for i in range(len(seqNamesAll)):
		if seqNamesAll[i] in PRIMER_NAMES:
			primerName = seqNamesAll[i]
			primerSeq = seqs[i]
			for j in range(len(primerSeq)):
				if primerSeq[j] != "-":
					primerToCoord[primerName] = j

					if j in coordToPrimer:
						coordToPrimer[j] = coordToPrimer[j]+", "+primerName
					else:
						coordToPrimer[j] = primerName
					break
			
			seqToMismatch = {}
			
			entropy = 0.0
			count = 1
			conSeq = ""
			primerLength = 0
			for j in range(len(primerSeq)):
				if primerSeq[j] != "-":
					primerLength += 1
					##some primer statistics##
					entropy = entropy + shannonEntropy(cols[j])
					count = count * len(DEG_CODE[primerSeq[j]])
					conSeq = conSeq + primerSeq[j]
						
					for x in range(len(cols[j])):
						okayLetters  = DEG_CODE[primerSeq[j]]
						if cols[j][x] != okayLetters and cols[j][x] not in okayLetters:
							if seqNamesNoPrimers[x] in seqToMismatch:
								seqToMismatch[seqNamesNoPrimers[x]] = seqToMismatch[seqNamesNoPrimers[x]] + 1
							else:
								seqToMismatch[seqNamesNoPrimers[x]] = 1

			print "%s\tcoordinate: %d\tentropy: %0.2f\tnumPrimers:%d\t%s" % (primerName, (primerToCoord[primerName]+1), entropy/primerLength, count, conSeq)
			
			print "%s captures:" % primerName

			histMismatches = {}
			maxMismatch = 1
			for seq, mismatch in seqToMismatch.iteritems():
				seqStats.write("%s\t%s\t%d\n" % (seq, seqNamesAll[i], mismatch))
				#print seq, mismatch
				if mismatch > maxMismatch: maxMismatch = mismatch

				if mismatch in histMismatches:
					histMismatches[mismatch] = histMismatches[mismatch] + 1
				else:
					histMismatches[mismatch] = 1
			
			for seq in seqNamesNoPrimers:
				if seq not in seqToMismatch:
					seqStats.write("%s\t%s\t0\n" % (seq, seqNamesAll[i]))

			
			#print histMismatches
					
			print "\t%d sequences with 0 mismatches" % (len(seqNamesNoPrimers) - len(seqToMismatch))

			for n in range(1,maxMismatch+1):
				if n in histMismatches:
					print "\t%d sequences with %d mismatch(es)" % (histMismatches[n], n)
				else:
					print "\t0 sequences with %d mismatch(es)" %n
			
			print ""				
	

	seqStats.close()
		
	##get column statistics
	entropy = []
	gap = []
	consensus = []
	count = []

	for j in range(len(cols)):
		entropy.append(shannonEntropy(cols[j]))
		gap.append(cols[j].count("-"))
		
		noGaps = cols[j].replace("-", "")
		count.append(len(set(noGaps)))
		
			
		uniq = "".join(sorted(set(noGaps)))

		##get rid of any degenerate bases
		uniqNodeg = ""
		for base in uniq:
			if base not in {"A", "C", "T", "G"}:
				uniqNodeg += DEG_CODE[base]
			else:
				uniqNodeg += base
		
		uniqFinal = "".join(sorted(set(uniqNodeg)))
		 
		consensus.append(DEG_CODE[uniqFinal])


	outfile = open(output, 'w')
	outfile.write("#coordinate\tentropy\tnumGaps\tnumPrimers\tconsensusSequence\tprimerName\n")
	for j in range(len(cols)-kmer):
		cumEntropy = 0.0
		cumGaps = 0
		cumCount = 1
		conSeq = ""

		for k in range(kmer):
			cumEntropy += entropy[j+k]
			cumGaps += gap[j+k]
			cumCount = cumCount * count[j+k]
			conSeq += consensus[j+k]	
		if j in coordToPrimer:
			outfile.write("%d\t%.2f\t%d\t%d\t%s\t%s\n" %(j+1, cumEntropy/kmer, cumGaps, cumCount, conSeq, coordToPrimer[j]))
		else:	
			outfile.write("%d\t%.2f\t%d\t%s\t%s\n" % (j+1, cumEntropy/kmer, cumGaps, cumCount, conSeq))

def filterVirusDb(virusDb, refProtein, outPrefix, refProtOnly):
	alignStats = {}
		
	#find length of refProtein 
	refProt = open(refProtein, 'r')
	numEntries = 0
	refLength = 0
	while True:
		line = refProt.readline()
		if line == "": break

		if line.startswith(">"):
			numEntries += 1
		else:
			refLength += len(line.strip())
		if numEntries > 1:
			exit("Error: refProtein contains more than one reference protein sequence")
	
	CMD = "%s -dbtype nucl -in %s" % (MAKEBLASTDB, virusDb)
	print CMD
	os.system(CMD)

	CMD = "\n%s -query %s -db %s -outfmt 6 -out %s.01.tblastn.txt -max_target_seqs 9999999" % (TBLASTN, refProtein, virusDb, outPrefix)
	print CMD
	os.system(CMD)
	
	results = open("%s.01.tblastn.txt" % outPrefix)
	keep = set()
	while True:
		line = results.readline()
		if line == "": break
		line = line.split()
		refName = line[0].strip()
		seqName = line[1].strip()
		length = int(line[3].strip())
		if float(line[10].strip()) > 0.05: continue

		if length > (0.8 * refLength):
			keep.add(seqName)
		
		
		##if refProtOnly flag is on, keep track of alignment statistics
		start = int(line[8].strip())
		end = int(line[9].strip())
		
		if seqName in alignStats:
			curStats = alignStats[seqName]
			if start < curStats[0]: curStats[0] = start
			if end > curStats[1]: curStats[1] = end
		else:
			curStats = [start,end]
		alignStats[seqName] = curStats
	#print alignStats
	##Make filtered fasta

	virusFile = open(virusDb, 'r')
	
	virusDbCount = 0
	filtVirusDb = {}
	
	seq = ""
	for line in virusFile.readlines():
		if line.startswith(">"):
			virusDbCount += 1
			##add previous sequence to virusDb
			if seq != "" and name in keep: filtVirusDb[name] = seq

			name = line[1:].split()[0].strip()
			seq = ""
		else:
			seq += line.strip()

	##add last entry
	if name in keep: filtVirusDb[name] = seq

	
	filtVirusFile = open(outPrefix+".01.homologous.fa", 'w')
	for name,seq in filtVirusDb.iteritems():
		filtVirusFile.write(">"+name+"\n")
		if refProtOnly:
			curStats = alignStats[name]
			start = curStats[0]-201
			if start < 0: start = 0
			end = curStats[1]+200
			if end > len(seq): end = len(seq)
			
			if end < start:
				tmp = start
				start = end
				end = tmp

			filtVirusFile.write(seq[start:end]+"\n")
		else:
			filtVirusFile.write(seq+"\n")
	filtVirusFile.close()


	#print filtVirusDb
	print("\nKeeping %d sequences with homology to %s out of total of %d sequences. Results written to %s.01.homologous.fa\n" % (len(keep), refName, virusDbCount, outPrefix))
	return(outPrefix+".01.homologous.fa")

def filterN(fasta, outPrefix, filtN):
	fastaFile = open(fasta, 'r')
	statsFile = open(outPrefix+".02.filtN.stats.txt", 'w')
	filtFile = open(outPrefix+".02.filtN.fa", 'w')
	stats = {}		
	length = 0
	numN = 0
	percN = 0
	name = ""
	seq = ""

	numTotal = 0
	numKeep = 0

	for line in fastaFile.readlines():
		if line.startswith(">"):
			numTotal += 1

			##keep track of entry before
			if name != "": 
				percN = (numN * 1.0) / length
				stats[name] = [length, numN, percN, seq]
				statsFile.write("%s\t%d\t%d\t%0.2f\n" % (name.split()[0].strip(), length, numN, percN))
				
				if percN <= filtN:
					filtFile.write(">%s\n%s\n" % (name, seq))
					numKeep += 1
			##reset variables
			name = line[1:].strip()
			length = 0
			numN = 0
			seq = ""
		else:
			length += len(line.strip())
			numN += line.upper().count('N')
			seq += line.strip()

	#add last entry
	percN = (numN * 1.0) / length
	stats[name] = [length, numN, percN, seq]
	statsFile.write("%s\t%d\t%d\t%0.2f\n" % (name, length, numN, percN))

	if percN <= filtN:
		filtFile.write(">%s\n%s\n" % (name, seq))
		numKeep += 1

	filtFile.close()
	statsFile.close()

	print "\tKeeping %d out of %d records with N content <= %0.2f" % (numKeep, numTotal, filtN)

def rc(seq):

	newSeq = ""
	for i in range(len(seq),0,-1):
		newSeq += RC_CODE[seq[(i-1)]]
	#print newSeq

	return newSeq

def main():
	parser = argparse.ArgumentParser(description='assessPrimers')
	parser.add_argument('virusDb', type=str, help='fasta file of all viral sequences (NUCLEOTIDE)')
	parser.add_argument('primers', type=str, help='fasta file of primer sequences')
	#parser.add_argument('r_primers', type=str, help='fasta file of reverse primer sequences (5\'->3\' NUCLEOTIDE)')
	parser.add_argument('outPrefix', type=str, help='prefix for output files')
	parser.add_argument('--refProtein', type=str, help='fasta file of reference protein sequence (AMINO ACID)')
	parser.add_argument('--overwrite', action='store_true', help='by default, this software will use intermediate files if they exist rather than re-run steps. set this flag to overwrite existing intermediate files.')
	parser.add_argument('--blastDir', type=str, help='path to folder with BLAST tools. Specifically need makeblastdb and tblastn.')
	parser.add_argument('--cdhit', type=str, help='path to cd-hit-est')
	parser.add_argument('--mafft', type=str, help='path to mafft')
	parser.add_argument('--filtN', type=float, default=0.05, help='threshold for filtering out sequences with N content higher than filtN. Default is 0.05')
	parser.add_argument('--refProtOnly', action='store_true', help='If a refProtein is provided, setting this flag will extract only the portion of sequences that align to the refProtein.')
	parser.add_argument('--idCutoff', type=float, default=0.75, help='sequences with >idCutoff similarity will be collapsed into one representative sequence. default = 0.9')
	parser.add_argument('--singleMSA', action='store_true', help='flag to align all primers together in a single multiple sequence alignment. default is to align each primer separately.')
	parser.add_argument('--kmer', type=int, default = 20, help="length of kmer to calculate entropy for")
	args = parser.parse_args()
	
	global MAKEBLASTDB
	global TBLASTN
	global CDHIT
	global MAFFT
	#TODO: check that all software is installed		
	MAKEBLASTDB = 'makeblastdb'
	TBLASTN = 'tblastn'
	CDHIT = 'cd-hit'
	MAFFT = 'mafft'
	if args.blastDir is not None: 
		MAKEBLASTDB = args.blastDir+"/makeblastdb"
		TBLASTN = args.blastDir+"/tblastn"
	if args.cdhit is not None: CDHIT = args.cdhit
	if args.mafft is not None: MAFFT = args.mafft
	

	F_PRIMERS = []
	#R_PRIMERS = []
	PRIMER_NAMES = []

	#TODO: Make sure primers have no illegal basepairs
	#STORE PRIMER INFORMATION
	f_primers = open(args.primers, 'r')
	name = ""
	seq = ""
	for line in f_primers.readlines():
		if line[0] == ">": 
			if name != "": 
				F_PRIMERS.append([name,seq]) 
				name = ""
				seq = ""
			name = line[1:].strip()
			PRIMER_NAMES.append(name)
			
		else:
			seq += line.strip()
		#add last entry
	F_PRIMERS.append([name, seq])
				
	f_primers.close()

	#print F_PRIMERS	
	
	'''
	#STORE PRIMER INFORMATION
	r_primers = open(args.r_primers, 'r')
	name = ""
	seq = ""
	for line in r_primers.readlines():
		#print line
		if line[0] == ">":
			if name != "":  
				R_PRIMERS.append([name, rc(seq)])
				name = ""
				seq = ""

			name = line[1:].strip()
			PRIMER_NAMES.append(name)
		else:
			seq += line.strip()
			#print "seq: "+seq

		#add last entry
	R_PRIMERS.append([name,rc(seq)])
	r_primers.close()
	'''
	
	print "\n***************************************************"
	print "STEP 1: FILTER FOR SEQUENCES WITH HOMOLOGY TO REFERENCE PROTEIN\n"

	if args.refProtein is not None:
		if args.overwrite or not os.path.exists(args.outPrefix+".01.homologous.fa"):
			virusDbFilt = filterVirusDb(args.virusDb, args.refProtein, args.outPrefix, args.refProtOnly)
		else: 
			print "%s.01.homologous.fa found! Skipping step..." % args.outPrefix
			virusDbFilt = args.outPrefix+".01.homologous.fa"
	else:
		print "--refProtein not used. Skipping step..."
		virusDbFilt = args.virusDb
	
	print "\n***************************************************"
	print "STEP 2: REMOVE SEQUENCES WITH HIGH N CONTENT\n"
	if args.overwrite or not os.path.exists(args.outPrefix+".02.filtN.fa"):
		filterN(virusDbFilt, args.outPrefix, args.filtN)

	else:
		print "%s.02.filtN.fa found! Skipping step..." % args.outPrefix

	#check filtN file is not empty
	filtN = open(args.outPrefix+".02.filtN.fa", 'r').readlines()
	if len(filtN) == 0:
		print "\nNo sequences left after removing sequences with high N content! Considering increasing --filtN argument. Exiting..."
		sys.exit()

	print "\n***************************************************"
	print "STEP 3: CLUSTER SEQUENCES WITH >%0.2f SIMILARITY. This step may take a while.\n" % args.idCutoff
	
	if args.overwrite or not os.path.exists(args.outPrefix+".03.cd-hit.fa"):
		CMD = "%s -i %s -o %s -c %0.2f" % (CDHIT, args.outPrefix+".02.filtN.fa", args.outPrefix+".03.cd-hit.fa", args.idCutoff)
		print CMD
		os.system(CMD)

	else:
		print "%s.03.cd-hit.fa found! Skipping step..." % args.outPrefix
	

	print "\n****************************************************"
	print "STEP 4: MAKE MULTIPLE SEQUENCE ALIGNMENT."
	
	if args.singleMSA:
		print "Aligning all primers in a single multiple sequence alignment.\n"
		if args.overwrite or not os.path.exists(args.outPrefix+".04.mafft.all_primers.fa"):
			os.system("cp "+args.outPrefix+".03.cd-hit.fa "+args.outPrefix+".03.cd-hit.all_primers.fa")
				
			cdhit_primers = open(args.outPrefix+".03.cd-hit.all_primers.fa", 'a')

			for i in range(len(F_PRIMERS)):
				cdhit_primers.write(">%s\n%s\n" % (F_PRIMERS[i][0], F_PRIMERS[i][1]))
				#cdhit_primers.write(">%s\n%s\n" % (R_PRIMERS[i][0], R_PRIMERS[i][1]))

			cdhit_primers.close()
				

			CMD = "%s --maxiterate 1000 --localpair %s > %s" % (MAFFT, args.outPrefix+".03.cd-hit.all_primers.fa", args.outPrefix+".04.mafft.all_primers.fa")
			print CMD
			
			os.system(CMD)
		else:
			print "%s.04.mafft.all_primers.fa found! Skipping step..." % args.outPrefix

	else:
		print "This step may take a while. If you are testing multiple primers, you can speed up this step by aligning all primers together with the --singleMSA option.\n"
		for i in range(len(F_PRIMERS)):
			print "Current primer: "+F_PRIMERS[i][0]
			curNo = str((i+1))
			if len(curNo) == 1: curNo = "0"+curNo
			curPrimer = "primer"+curNo	
			
			if args.overwrite or not os.path.exists(args.outPrefix+".04.mafft."+curPrimer+".fa"):	
					
				##add primer sequences to cd-hit fasta
				#print "cp "+args.outPrefix+".02.cd-hit.fa "+args.outPrefix+".02.cd-hit."+curPrimer+".fa"
				os.system("cp "+args.outPrefix+".03.cd-hit.fa "+args.outPrefix+".03.cd-hit."+curPrimer+".fa")
				
				cdhit_primers = open(args.outPrefix+".03.cd-hit."+curPrimer+".fa", 'a')
				cdhit_primers.write(">%s\n%s\n" % (F_PRIMERS[i][0], F_PRIMERS[i][1]))
				#cdhit_primers.write(">%s\n%s\n" % (R_PRIMERS[i][0], R_PRIMERS[i][1]))

				cdhit_primers.close()
				

				CMD = "%s --maxiterate 1000 --localpair %s > %s" % (MAFFT, args.outPrefix+".03.cd-hit."+curPrimer+".fa", args.outPrefix+".04.mafft."+curPrimer+".fa")
				print CMD
				os.system(CMD)
				
				##cleanup
				os.system("rm "+args.outPrefix+".03.cd-hit."+curPrimer+".fa")
			else:
				print "%s.04.mafft.%s.fa found! Skipping %s..." % (args.outPrefix, curPrimer, curPrimer)
	
	print "\n****************************************************"
	print "STEP 5: CALCULATE PRIMER STATISTICS\n"
	if args.singleMSA:
			scanEntropy(args.outPrefix+".04.mafft.all_primers.fa", args.outPrefix, args.outPrefix+".05.all_kmer_statistics.all_primers.txt", PRIMER_NAMES, args.kmer)
		#do something
	else:	
		for i in range(len(F_PRIMERS)):	
			curNo = str((i+1))
			if len(curNo) == 1: curNo = "0"+curNo
			curPrimer = "primer"+curNo	
	
			scanEntropy(args.outPrefix+".04.mafft."+curPrimer+".fa", args.outPrefix, args.outPrefix+".05.all_kmer_statistics."+curPrimer+".txt", PRIMER_NAMES, args.kmer)


if __name__ == "__main__":
	main()
