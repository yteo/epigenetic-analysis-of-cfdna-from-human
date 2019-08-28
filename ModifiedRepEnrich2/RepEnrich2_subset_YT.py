#!/usr/bin/env python
import argparse
import csv
import numpy
import os
import shlex
import shutil
import subprocess
import sys

parser = argparse.ArgumentParser(description='Prepartion for downstream running of RepEnrich2 by subsetting reads to uniquely mapped and multi-mapped. For this script to run properly samtools and bedtools must be loaded.  You will need 1) An alignment file. 2) A defined MAPQ threshold (30 by default)   command-line usage EXAMPLE: python RepEnrich_subset.py /example/path/data/alignment.sam 30', prog='RepEnrich2_subset.py')
parser.add_argument('--version', action='version', version='%(prog)s 0.1')
parser.add_argument('alignment_file', action= 'store', metavar='alignment_file', help='bam file of aligned reads. Example /example/path/alignment.bam')
parser.add_argument('MAPQ', action= 'store', metavar='MAPQ', default='30', help='MAPQ threshold for subsetting uniquely mapping and multi-mapping reads.  This will be aligner-dependent as many aligners treat the MAPQ field differently (Default = 30) ')
parser.add_argument('sample_name', action='store', metavar='sample_name', help='The name of your sample (used to generate the name of the output files)')
#parser.add_argument('repeattype', action='store', metavar='repeattype', help='Either ALUY or L1HS')
parser.add_argument('--pairedend', action= 'store', dest='pairedend', default= 'FALSE', help='Designate this option for paired-end data. Default FALSE change to TRUE')


#parser.add_argument('--repeatbed', action= 'store', dest='repeatbed', default= '/users/yteo/data/yteo/Repeats/L1HS.bed', help='Path to repeats bed file in the genome')
#parser.add_argument('--conversionbed', action= 'store', dest='conversionbed', default= '/users/yteo/data/yteo/Repeats/hs37d5L1HS_consensus_cigar1bp_sorted.bed', help='Path to genome to consensus conversion bed file')
# conversion bed AluY is at /users/yteo/data/yteo/Repeats/Alu/hs37d5AluY_consensus_cigar1bp_sorted.bed
# repeat bed AluY is at /users/yteo/data/yteo/Repeats/Alu/AluY_dimeric_noAluYlinearge_chronly.bed
# bed files chromosome names are without the suffix "chr"
parser.add_argument('--counting', action= 'store', dest='counting', default= 'Count.sh', help='Path to Count.sh file')
# counting folder

# extend reads to cover the whole fragment and count the coverage (Extend = true), or false (only count coverage at mapped reads)
parser.add_argument('--extend', action= 'store', dest='extend', default= 'FALSE', help='Extend reads to cover the whole fragment and count the coverage. Can only be used with paired end reads. Default: FALSE')

#######################YT ADDED THE REPEATBED ARGUMENT
args = parser.parse_args()


# parameters specified in args_parse
alignment_file = args.alignment_file
MAPQ = args.MAPQ
sample_name = args.sample_name
pairedend = args.pairedend
extend = args.extend
#repeatbed = args.repeatbed
#conversionbed = args.conversionbed
counting = args.counting
#repeattype = args.repeattype


################################################################################
# check for loaded dependencies
try:
	subprocess.call(shlex.split("samtools --version"), stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
except OSError:
	print "Error: Samtools not loaded"
	raise
try:
	subprocess.call(shlex.split("bedtools --version"), stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
except OSError:
	print "Error: Bedtools not loaded"
	raise
################################################################################


# #output .bam (unique) and .fastq/s (multimapping) alignments

if pairedend == "FALSE":
	print "Subsetting uniquely mapping reads above MAPQ threshold..."
	fileout= sample_name + '_unique.bam'
	with open(fileout, 'w') as stdout:
		#filter for mapped reads, with MAPQ above the specified value
		command = shlex.split("samtools view -F 4 -bq " + MAPQ + " " + alignment_file)
		p = subprocess.Popen(command, stdout=stdout)
	p.communicate()
	stdout.close()
	print "Done."


if pairedend == "TRUE"	:
	print "Subsetting uniquely mapping paired-end reads in proper pair above MAPQ threshold..."
	fileout= sample_name + '_unique.bam'
	with open(fileout, 'w') as stdout:
		#filter for paired reads with read mapped in proper pair, with MAPQ above the specified value
		command = shlex.split("samtools view -f 3 -bq " + MAPQ + " " + alignment_file)
		p = subprocess.Popen(command, stdout=stdout)
		p.communicate()
		stdout.close()
	print "Done."







##############################################################
if pairedend == "FALSE":
	print "Subsetting multi-mapping reads..."
	fileout2= sample_name + '_map.bam'
	fileout3= sample_name + '_multimap_filtered.bam'
	fileout4= sample_name + '_multimap.fastq'
	with open(fileout2, 'w') as stdout:
		#filter for mapped reads only
		command = shlex.split("samtools view -h -F 4 -b " + alignment_file)
		p = subprocess.Popen(command, stdout=stdout)

		p.communicate()
		stdout.close()

	#filter from the mapped reads for everything that does not have a MAPQ meeting the threshold; everything should already map from the last filter
	#-U option requires samtools 1.3.0+
	command = shlex.split("samtools view -U " + fileout3 + " -bq " + MAPQ + " " + fileout2)
	#output to fileout3 and silence stdout
	p = subprocess.Popen(command, stdout=open(os.devnull, 'wb'))
	p.communicate()

	#convert bamfile to fastq for downstream use
	command = shlex.split("bedtools bamtofastq -i " + fileout3 + " -fq " + fileout4)
	p=subprocess.Popen(command)
	p.communicate()
	print "Done."

	print "Performing cleanup..."
	#cleanup intermediate files
	command = shlex.split("rm " + fileout2)
	p=subprocess.Popen(command)
	p.communicate()
	stdout.close()
	command = shlex.split("rm " + fileout3)
	p.communicate()
	stdout.close()
	print "Done."

##############################
if pairedend == "TRUE":
	print "Subsetting paired-end multi-mapping reads..."
	fileout2= sample_name + '_map.bam'
	fileout3= sample_name + '_multimap_filtered.bam'
	fileout4= sample_name + '_multimap_sorted.bam'
	fileout5= sample_name + '_multimap_R1.fastq'
	fileout6= sample_name + '_multimap_R2.fastq'

	with open(fileout2, 'w') as stdout:
    #filter for paired / mapped in proper pair only
		command = shlex.split("samtools view -f 3 -b " + alignment_file)
		p = subprocess.Popen(command, stdout=stdout)
		p.communicate()
		stdout.close()

	#filter from the mapped reads for everything that does not have a MAPQ meeting the threshold for being uniquely mapped
	command = shlex.split("samtools view -U " + fileout3 + " -bq " + MAPQ + " " + fileout2)
	#output to fileout3 and silence stdout
	p = subprocess.Popen(command, stdout=open(os.devnull, 'wb'))
	p.communicate()
	stdout.close()

	print "Performing sort on bamfile..."
	#sort bamfile by name for use with bamtofastq (sorting only needed for paired end reads)
	command = shlex.split("samtools sort -n " + fileout3 + " -o " + fileout4)
	p = subprocess.Popen(command)
	p.communicate()
	

	print "Converting to fastq..."
	#convert bamfile to fastq for downstream use
	command = shlex.split("bedtools bamtofastq -i " + fileout4 + " -fq " + fileout5 + " -fq2 " + fileout6)
	p=subprocess.Popen(command)
	p.communicate()
	

	print "Performing cleanup..."
	#cleanup intermediate files
	command = shlex.split("rm " + fileout2)
	p=subprocess.Popen(command)
	p.communicate()
	command = shlex.split("rm " + fileout3)
	p=subprocess.Popen(command)
	p.communicate()
	command = shlex.split("rm " + fileout4)
	p=subprocess.Popen(command)
	p.communicate()

	print "Done."



######################################################################
#######################YT ADDED THIS CHUNK
# extracting unique reads mapping to L1HS of Aluy
Repeattypes = ["L1HS","ALUY"]
for repeattype in Repeattypes:
	if repeattype == "L1HS": 
		repeatbed="./TE_Data/L1HS.bed"
		conversionbed="./TE_Data/hs37d5L1HS_consensus_cigar1bp_sorted.bed"
	elif repeattype =="ALUY":
		repeatbed="./TE_Data/Alu/AluY_dimeric_noAluYlinearge_chronly.bed"
		conversionbed="./TE_Data/hs37d5AluY_consensus_cigar1bp_sorted.bed"

	fileout= sample_name + '_unique.bam'
	fileout_rep= sample_name + '_unique_'+repeattype+'.bam'
	print "Extracting uniquely mapped reads mapping to L1HS or AluY..."
	with open(fileout_rep, 'w') as stdout:
		command = shlex.split("samtools view -b -L " + repeatbed + " " +fileout)
		p = subprocess.Popen(command, stdout=stdout)
		p.communicate()
		stdout.close()
	print "Done."



	# convert bam to bed
	fileout_rep_bed= sample_name + '_unique_'+repeattype+'.bed'
	print "Converting repeat mapping bam file to bed format..."
	with open(fileout_rep_bed, 'w') as stdout:
		command = shlex.split("bedtools bamtobed -i " + fileout_rep)
		p = subprocess.Popen(command, stdout=stdout)
		p.communicate()
		stdout.close()
	print "Done."


	fileout_rep_bed_sorted= sample_name + '_unique_'+repeattype+'_sorted.bed'
	with open(fileout_rep_bed_sorted, 'w') as stdout:
		command = shlex.split("sort -k1,1 -k2,2n " + fileout_rep_bed)
		p = subprocess.Popen(command, stdout=stdout)
		p.communicate()
		stdout.close()


	# counting uniquely mapped reads mapping to each L1HS or AluY consensus position
	print "Counting uniquely mapped reads mapping to each L1HS or AluY consensus position..."
	if extend == "TRUE" and pairedend == "TRUE" :
		##########
		# extend read pairs to find the coverage of the fragment
		# store the mapping coordinates of a read into read {} dictionary
		read={}
		chromosome={}
		fin=open(fileout_rep_bed_sorted,"r")
		for line in fin:
			# list to store the coordinates of mapping reads
			#readcoor=[]
			line=line.strip("\r").strip("\n").split("\t")
			readid=line[3].split("/")[0]
			#readcoor.extend((line[1],line[2]))
			coor1=line[1]
			coor2=line[2]
			# storing the chrososome that the read is mapped to in a dictionary, make sure that the paired reads are on the same chromosome
			if readid in chromosome:
				if chromosome[readid] != line[0]:
					sys.exit("Error with the chromosome mapping")
			else:
				chromosome[readid]=line[0]
			if readid in read:
				read[readid].append(coor1)
				read[readid].append(coor2)		
				#read[readid].append(readcoor)
			else:
				read[readid]=[coor1]
				read[readid].append(coor2)

		fin.close()
		extended=sample_name + '_unique_'+repeattype+'_extend.bed'
		fout=open(extended,"w")
		for readid, coord in read.items():
			fout.write(str(chromosome[readid])+"\t"+ str(min(coord))+"\t"+str(max(coord))+"\t"+str(readid)+"\n")
		fout.close()




		extended_sorted= sample_name + '_unique_'+repeattype+'_extend_sorted.bed'
		with open(extended_sorted, 'w') as stdout:
			command = shlex.split("sort -k1,1 -k2,2n " + extended)
			p = subprocess.Popen(command, stdout=stdout)
			p.communicate()
			stdout.close()

		##############
		fileout_con= sample_name + '_'+repeattype+'_consensuscoord.bed'
		print "Intersecting consensus position and mapping reads..."
		with open(fileout_con, 'w') as stdout:
			command = shlex.split("intersectBed -a "+ conversionbed + " -b " + extended_sorted + " -sorted -wa -wb ")
			p = subprocess.Popen(command, stdout=stdout)
			p.communicate()
			stdout.close()



		subprocess.call(shlex.split(counting+" "+fileout_con ))
	else:
		##############
		fileout_con= sample_name + '_'+repeattype+'_consensuscoord.bed'
		print "Intersecting consensus position and mapping reads..."
		with open(fileout_con, 'w') as stdout:
			command = shlex.split("intersectBed -a "+ conversionbed + " -b " + fileout_rep_bed_sorted + " -sorted -wa -wb ")
			p = subprocess.Popen(command, stdout=stdout)
			p.communicate()
			stdout.close()
		print "count the occurence of reads mapping to each consensus location. Collapse paired end reads"
		subprocess.call(shlex.split(counting+" "+fileout_con ))


	if os.path.exists(sample_name + '_'+repeattype+'_consensuscoord.bed'):
		os.remove(sample_name + '_'+repeattype+'_consensuscoord.bed')
	if os.path.exists(sample_name + '_unique_'+repeattype+'.bed'):
		os.remove(sample_name + '_unique_'+repeattype+'.bed')
	if os.path.exists(sample_name + '_unique_'+repeattype+'.bam'):
		os.remove(sample_name + '_unique_'+repeattype+'.bam')
	if os.path.exists(sample_name + '_unique_'+repeattype+'_sorted.bed'):
		os.remove(sample_name + '_unique_'+repeattype+'_sorted.bed')