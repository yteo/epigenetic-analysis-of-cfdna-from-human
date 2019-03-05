#!/usr/bin/env python
import argparse
import csv
import numpy
import os
import shlex
import shutil
import subprocess
import sys

parser = argparse.ArgumentParser(description='Part II: Conducting the alignments to the psuedogenomes.  Before doing this step you will require 1) a bamfile of the unique alignments with index 2) a fastq file of the reads mapping to more than one location.  These files can be obtained using the following bowtie options [EXAMPLE: bowtie -S -m 1 --max multimap.fastq mm9 mate1_reads.fastq]  Once you have the unique alignment bamfile and the reads mapping to more than one location in a fastq file you can run this step.  EXAMPLE: python master_output.py /users/nneretti/data/annotation/hg19/hg19_repeatmasker.txt /users/nneretti/datasets/repeatmapping/POL3/Pol3_human/HeLa_InputChIPseq_Rep1 HeLa_InputChIPseq_Rep1 /users/nneretti/data/annotation/hg19/setup_folder HeLa_InputChIPseq_Rep1_multimap.fastq HeLa_InputChIPseq_Rep1.bam')
parser.add_argument('--version', action='version', version='%(prog)s 0.1')
parser.add_argument('annotation_file', action= 'store', metavar='annotation_file', help='List RepeatMasker.org annotation file for your organism.  The file may be downloaded from the RepeatMasker.org website.  Example: /data/annotation/hg19/hg19_repeatmasker.txt')
parser.add_argument('outputfolder', action= 'store', metavar='outputfolder', help='List folder to contain results. Example: /outputfolder')
parser.add_argument('outputprefix', action= 'store', metavar='outputprefix', help='Enter prefix name for data. Example: HeLa_InputChIPseq_Rep1')
parser.add_argument('setup_folder', action= 'store', metavar='setup_folder', help='List folder that contains the repeat element psuedogenomes. Example /data/annotation/hg19/setup_folder')
parser.add_argument('fastqfile', action= 'store', metavar='fastqfile', help='Enter file for the fastq reads that map to multiple locations. Example /data/multimap.fastq')
parser.add_argument('alignment_bam', action= 'store', metavar='alignment_bam', help='Enter bamfile output for reads that map uniquely. Example /bamfiles/old.bam')
parser.add_argument('--pairedend', action= 'store', dest='pairedend', default= 'FALSE', help='Designate this option for paired-end sequencing. Default FALSE change to TRUE')
parser.add_argument('--collapserepeat', action= 'store', dest='collapserepeat', metavar='collapserepeat', default= 'Simple_repeat', help='Designate this option to generate a collapsed repeat type.  Uncollapsed output is generated in addition to collapsed repeat type.  Simple_repeat is default to simplify downstream analysis.  You can change the default to another repeat name to collapse a seperate specific repeat instead or if the name of Simple_repeat is different for your organism.  Default Simple_repeat')
parser.add_argument('--fastqfile2', action= 'store', dest='fastqfile2', metavar='fastqfile2', default= 'none', help='Enter fastqfile2 when using paired-end option.  Default none')
parser.add_argument('--cpus', action= 'store', dest='cpus', metavar='cpus', default= "1", type=int, help='Enter available cpus per node.  The more cpus the faster RepEnrich performs. RepEnrich is designed to only work on one node. Default: "1"')
parser.add_argument('--allcountmethod', action= 'store', dest='allcountmethod', metavar='allcountmethod', default= "FALSE", help='By default the pipeline only outputs the fraction count method.  Consdidered to be the best way to count multimapped reads.  Changing this option will include the unique count method, a conservative count, and the total count method, a liberal counting strategy. Our evaluation of simulated data indicated fraction counting is best. Default = FALSE, change to TRUE')
parser.add_argument('--is_bed', action= 'store', dest='is_bed', metavar='is_bed', default= 'FALSE', help='Is the annotation file a bed file.  This is also a compatible format.  The file needs to be a tab seperated bed with optional fields.  Ex. format chr\tstart\tend\tName_element\tclass\tfamily.  The class and family should identical to name_element if not applicable.  Default FALSE change to TRUE')
parser.add_argument('--multi', action= 'store', dest='multi', default= '/users/yteo/data/yteo/Software/RepEnrich2/YT/Multimapped.sh', help='Path to Multimapped.sh file')

args = parser.parse_args()
L1HS_consensus_bt="/gpfs/data/nneretti/yteo/Blood_Microbiome/Reads/Clean/annotation_retrofind/L1HS_Alu/L1HS"
ALUY_consensus_bt="/gpfs/data/nneretti/yteo/Blood_Microbiome/Reads/Clean/annotation_retrofind/L1HS_Alu/ALUY"

L1HS_1stepbed="/users/yteo/data/yteo/Software/RepEnrich2/YT/L1HS_1bpstep.bed"
ALUY_1stepbed="/users/yteo/data/yteo/Software/RepEnrich2/YT/AluY_1bpstep.bed"
# parameters
annotation_file = args.annotation_file
outputfolder = args.outputfolder
outputfile_prefix = args.outputprefix
setup_folder = args.setup_folder
repeat_bed = setup_folder + os.path.sep + 'repnames.bed' 
unique_mapper_bam = args.alignment_bam
fastqfile_1 = args.fastqfile
fastqfile_2 = args.fastqfile2
cpus = args.cpus
b_opt = "-k 1 -p " +str(1) +" --quiet"
simple_repeat = args.collapserepeat
paired_end = args.pairedend
allcountmethod = args.allcountmethod
is_bed = args.is_bed
multi = args.multi
################################################################################
# check that the programs we need are available
try:
	subprocess.call(shlex.split("coverageBed -h"), stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
	subprocess.call(shlex.split("bowtie2 --version"), stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
except OSError:
	print "Error: Bowtie2 or BEDTools not loaded"
	raise

# check version of bedtools and store numbers in variables for determining syntax of bedCoverage
btv1 = subprocess.Popen('bedtools --version | cut -d" " -f2 | cut -d"v" -f2 | cut -d"." -f1', shell=True, stdout=subprocess.PIPE)
v1 = btv1.stdout.read()
btv2 = subprocess.Popen('bedtools --version | cut -d" " -f2 | cut -d"v" -f2 | cut -d"." -f2', shell=True, stdout=subprocess.PIPE)
v2 = btv2.stdout.read()

################################################################################
# define a csv reader that reads space deliminated files
print 'Preparing for analysis using RepEnrich2...'
csv.field_size_limit(sys.maxsize)
def import_text(filename, separator):
	for line in csv.reader(open(filename), delimiter=separator, 
						   skipinitialspace=True):
		if line:
			yield line
			
################################################################################
# build dictionaries to convert repclass and rep families'
if is_bed == "FALSE":
	repeatclass = {}
	repeatfamily = {}
	fin = import_text(annotation_file, ' ')
	x = 0
	for line in fin:
		if x>2:
			classfamily =[]
			classfamily = line[10].split(os.path.sep)
			line9 = line[9].replace("(","_").replace(")","_").replace("/","_")
			repeatclass[line9] = classfamily[0]
			if len(classfamily) == 2:
				repeatfamily[line9] = classfamily[1]
			else:
				repeatfamily[line9] = classfamily[0]
		x +=1
if is_bed == "TRUE":
	repeatclass = {}
	repeatfamily = {}
	fin = open(annotation_file, 'r')
	for line in fin:
		line=line.strip('\n')
		line=line.split('\t')
		theclass =line[4]
		thefamily = line[5]
		line3 = line[3].replace("(","_").replace(")","_").replace("/","_")
		repeatclass[line3] = theclass 
		repeatfamily[line3] = thefamily
fin.close()

################################################################################
# build list of repeats initializing dictionaries for downstream analysis'
fin = import_text(setup_folder + os.path.sep + 'repgenomes_key.txt', '\t')
repeat_key ={}
rev_repeat_key ={}
repeat_list = []
reptotalcounts = {}
classfractionalcounts = {}
familyfractionalcounts = {}
classtotalcounts = {}
familytotalcounts = {}
reptotalcounts_simple = {}
fractionalcounts = {}
i = 0
for line in fin:
	reptotalcounts[line[0]] = 0
	fractionalcounts[line[0]] = 0
	if repeatclass.has_key(line[0]):
		classtotalcounts[repeatclass[line[0]]] = 0
		classfractionalcounts[repeatclass[line[0]]] = 0
	if repeatfamily.has_key(line[0]):
		familytotalcounts[repeatfamily[line[0]]] = 0
		familyfractionalcounts[repeatfamily[line[0]]] = 0
	if repeatfamily.has_key(line[0]):
		if repeatfamily[line[0]] == simple_repeat:
			reptotalcounts_simple[simple_repeat] = 0
	else:
		reptotalcounts_simple[line[0]] = 0
	repeat_list.append(line[0])
	repeat_key[line[0]] = int(line[1])
	rev_repeat_key[int(line[1])] = line[0]
fin.close()
################################################################################
# map the repeats to the psuedogenomes:
if not os.path.exists(outputfolder):
	os.mkdir(outputfolder)
################################################################################
#Conduct the regions sorting
print 'Conducting region sorting on unique mapping reads....'
fileout= outputfolder + os.path.sep + outputfile_prefix + '_regionsorter.txt'
with open(fileout, 'w') as stdout:
	#check version of bedtools to determine which syntax to use for coverageBed
	if (int(v1) == 2 and int(v2) < 24):
		command = shlex.split("coverageBed -abam " +unique_mapper_bam+" -b " +setup_folder + os.path.sep + 'repnames.bed')
	else:
		command = shlex.split("coverageBed -a " + setup_folder + os.path.sep + 'repnames.bed' + " -b " + unique_mapper_bam)
	p = subprocess.Popen(command, stdout=stdout)
p.communicate()
stdout.close()
filein = open(outputfolder + os.path.sep + outputfile_prefix + '_regionsorter.txt','r')
counts = {}
sumofrepeatreads=0
for line in filein:
	line= line.split('\t')
	if not counts.has_key(str(repeat_key[line[3]])):
		counts[str(repeat_key[line[3]])]=0
	counts[str(repeat_key[line[3]])]+=int(line[4])
	sumofrepeatreads+=int(line[4])
print 'Identified ' + str(sumofrepeatreads) + ' unique reads that mapped to repeats.'
################################################################################
if paired_end == 'TRUE':
	if not os.path.exists(outputfolder + os.path.sep + 'pair1_'):
		os.mkdir(outputfolder + os.path.sep + 'pair1_')
	if not os.path.exists(outputfolder + os.path.sep + 'pair2_'):
		os.mkdir(outputfolder + os.path.sep + 'pair2_')
	folder_pair1 = outputfolder + os.path.sep + 'pair1_'
	folder_pair2 = outputfolder + os.path.sep + 'pair2_'
################################################################################
	print "Processing repeat psuedogenomes..."
	ps = []
	psb= []
	ticker= 0
	for metagenome in repeat_list:
		metagenomepath = setup_folder + os.path.sep + metagenome
		file1=folder_pair1 + os.path.sep + metagenome + '.txt'
		file2 =folder_pair2 + os.path.sep + metagenome + '.txt'
		with open(file1, 'w') as stdout:
			p = subprocess.Popen('bowtie2 ' + b_opt + ' ' + '-x ' + metagenomepath + ' ' + fastqfile_1 + ' | grep \"repname\" -', shell=True, stdout=stdout)
			
		with open(file2, 'w') as stdout:
			pp = subprocess.Popen('bowtie2 ' + b_opt + ' ' + '-x ' + metagenomepath + ' ' + fastqfile_2 + ' | grep \"repname\" -', shell=True, stdout=stdout)
		


		ps.append(p)
		ticker +=1
		psb.append(pp)
		ticker +=1
		if ticker == cpus:
			for p in ps:
				p.communicate()
			for p in psb:
				p.communicate()
			ticker = 0
			psb =[]
			ps = []
	if len(ps) > 0:
		for p in ps:
			p.communicate()
	stdout.close()



	
################################################################################
# combine the output from both read pairs:
	print 'sorting and combining the output for both read pairs...'
	if not os.path.exists(outputfolder + os.path.sep + 'sorted_'):
		os.mkdir(outputfolder + os.path.sep + 'sorted_')
	sorted_ = outputfolder + os.path.sep + 'sorted_'
	for metagenome in repeat_list:
		file1 = folder_pair1 + os.path.sep + metagenome + '.txt'
		file2 = folder_pair2 + os.path.sep + metagenome + '.txt'
		fileout= sorted_ + os.path.sep + metagenome + '.txt'
		with open(fileout, 'w') as stdout:
			p1 = subprocess.Popen(['cat',file1,file2], stdout = subprocess.PIPE)
			p2 = subprocess.Popen(['cut', '-f1',"-d "], stdin = p1.stdout, stdout = subprocess.PIPE)
			p3 = subprocess.Popen(['cut', '-f1', "-d/"], stdin = p2.stdout, stdout = subprocess.PIPE)
			p4 = subprocess.Popen(['sort'], stdin=p3.stdout, stdout = subprocess.PIPE)
			p5 = subprocess.Popen(['uniq'], stdin=p4.stdout, stdout = stdout)
			p5.communicate()
		stdout.close()
	print 'completed ...'


################# yt added this chunk
# FOR MULTI MAPPED READS
### count how many repeats a read has mapped to

### count how many repeats a read has mapped to
	repeatdict1 = {}
	repeatdict = {}
	sorted_ = outputfolder + os.path.sep + 'sorted_'
	for rep in repeat_list:
		for data in import_text(sorted_ + os.path.sep + rep + '.txt', '\t'):
			repeatdict1[data[0]] = []

	for rep in repeat_list:
		for data in import_text(sorted_ + os.path.sep + rep + '.txt', '\t'):
			repeatdict1[data[0]].append(str(repeat_key[rep]))

	for readid,rep in repeatdict1.items():
		repeatdict[readid]=len(rep)

	file_l1hs1=outputfolder + os.path.sep + outputfile_prefix + '_L1HS1.bam'
	file_l1hs2=outputfolder + os.path.sep + outputfile_prefix + '_L1HS2.bam'
	file_aluy1=outputfolder + os.path.sep + outputfile_prefix + '_ALUY1.bam'
	file_aluy2=outputfolder + os.path.sep + outputfile_prefix + '_ALUY2.bam'
	for metagenome in repeat_list:
		metagenomepath = setup_folder + os.path.sep + metagenome
		if metagenome =="L1HS":
			with open(file_l1hs1, 'w') as stdout:
				command = shlex.split("bowtie2 " + b_opt + " " + "-x " + metagenomepath + " " + fastqfile_1)
				p = subprocess.Popen(command, stdout=stdout)
				p.communicate()
			stdout.close()
			with open(file_l1hs2, 'w') as stdout:
				command = shlex.split("bowtie2 " + b_opt + " " + "-x " + metagenomepath + " " + fastqfile_2)
				pp = subprocess.Popen(command, stdout=stdout)
				pp.communicate()
			stdout.close()
		if metagenome =="AluY":
			with open(file_aluy1, 'w') as stdout:
				command = shlex.split("bowtie2 " + b_opt + " " + "-x " + metagenomepath + " " + fastqfile_1)
				p = subprocess.Popen(command, stdout=stdout)
				p.communicate()
			stdout.close()
			with open(file_aluy2, 'w') as stdout:
				command = shlex.split("bowtie2 " + b_opt + " " + "-x " + metagenomepath + " " + fastqfile_2)
				pp = subprocess.Popen(command, stdout=stdout)
				pp.communicate()
			stdout.close()
		
	file_l1hs1_mapped=outputfolder + os.path.sep + outputfile_prefix + '_L1HS1_mapped.bam'
	with open(file_l1hs1_mapped, 'w') as stdout:
		#filter for mapped reads only
		command = shlex.split("samtools view -h -F 4 -b " + file_l1hs1)
		p = subprocess.Popen(command, stdout=stdout)
		p.communicate()
		stdout.close()
	file_l1hs2_mapped=outputfolder + os.path.sep + outputfile_prefix + '_L1HS2_mapped.bam'
	with open(file_l1hs2_mapped, 'w') as stdout:
		#filter for mapped reads only
		command = shlex.split("samtools view -h -F 4 -b " + file_l1hs2)
		p = subprocess.Popen(command, stdout=stdout)
		p.communicate()
		stdout.close()
	file_aluy1_mapped=outputfolder + os.path.sep + outputfile_prefix + '_ALUY1_mapped.bam'
	with open(file_aluy1_mapped, 'w') as stdout:
		#filter for mapped reads only
		command = shlex.split("samtools view -h -F 4 -b " + file_aluy1)
		p = subprocess.Popen(command, stdout=stdout)
		p.communicate()
		stdout.close()
	file_aluy2_mapped=outputfolder + os.path.sep + outputfile_prefix + '_ALUY2_mapped.bam'
	with open(file_aluy2_mapped, 'w') as stdout:
		#filter for mapped reads only
		command = shlex.split("samtools view -h -F 4 -b " + file_aluy2)
		p = subprocess.Popen(command, stdout=stdout)
		p.communicate()
		stdout.close()
		
	if os.path.exists(outputfolder + os.path.sep + outputfile_prefix + '_L1HS1.bam'):
		os.remove(outputfolder + os.path.sep + outputfile_prefix + '_L1HS1.bam')
	if os.path.exists(outputfolder + os.path.sep + outputfile_prefix + '_L1HS2.bam'):
		os.remove(outputfolder + os.path.sep + outputfile_prefix + '_L1HS2.bam')
	if os.path.exists(outputfolder + os.path.sep + outputfile_prefix + '_ALUY1.bam'):
		os.remove(outputfolder + os.path.sep + outputfile_prefix + '_ALUY1.bam')
	if os.path.exists(outputfolder + os.path.sep + outputfile_prefix + '_ALUY2.bam'):
		os.remove(outputfolder + os.path.sep + outputfile_prefix + '_ALUY2.bam')
# convert aluy or l1hs mapped reads to the human genome , to fastq (to be  mapped to the consensus)
	file_l1hsf1=outputfolder + os.path.sep + outputfile_prefix + '_L1HS1.fastq'
	file_l1hsf2=outputfolder + os.path.sep + outputfile_prefix + '_L1HS2.fastq'
	file_aluyf1=outputfolder + os.path.sep + outputfile_prefix + '_ALUY1.fastq'
	file_aluyf2=outputfolder + os.path.sep + outputfile_prefix + '_ALUY2.fastq'
	command = shlex.split("bedtools bamtofastq -i " +  file_l1hs1_mapped + " -fq " + file_l1hsf1)
	p=subprocess.Popen(command)
	p.communicate()
	
	command = shlex.split('bedtools bamtofastq -i '+ file_l1hs2_mapped + ' -fq ' + file_l1hsf2)
	p = subprocess.Popen(command)
	p.communicate()



	command = shlex.split('bedtools bamtofastq -i '+ file_aluy1_mapped + ' -fq ' + file_aluyf1)
	p = subprocess.Popen(command)
	p.communicate()



	command = shlex.split('bedtools bamtofastq -i '+ file_aluy2_mapped + ' -fq ' + file_aluyf2)
	p = subprocess.Popen(command)
	p.communicate()

# mapping reads that mapped to the pseudogenomes , to the consensus
# aligning to L1HS consensus, read 1 and read 2 separately
	file_l1hsCon1=outputfolder + os.path.sep + outputfile_prefix + '_L1HS_con1.sam'
	file_l1hsCon2=outputfolder + os.path.sep + outputfile_prefix + '_L1HS_con2.sam'
	file_aluyCon1=outputfolder + os.path.sep + outputfile_prefix + '_ALUY_con1.sam'
	file_aluyCon2=outputfolder + os.path.sep + outputfile_prefix + '_ALUY_con2.sam'	
	with open(file_l1hsCon1, 'w') as stdout:
		command = shlex.split("bowtie2 --no-unal -p 16 -N 1 --local -x " + L1HS_consensus_bt+ " -U " + file_l1hsf1)
		p = subprocess.Popen(command, stdout=stdout )
		p.communicate()
		stdout.close()
	with open(file_l1hsCon2, 'w') as stdout:
		command = shlex.split("bowtie2 --no-unal -p 16 -N 1 --local -x " + L1HS_consensus_bt+ " -U " + file_l1hsf2)
		p = subprocess.Popen(command, stdout=stdout )
		p.communicate()
		stdout.close()
# aligning to AluY consensus
	with open(file_aluyCon1, 'w') as stdout:
		command = shlex.split("bowtie2 --no-unal -p 16 -N 1 --local -x " + ALUY_consensus_bt+ " -U " + file_aluyf1)
		p = subprocess.Popen(command, stdout=stdout )
		p.communicate()
		stdout.close()
	with open(file_aluyCon2, 'w') as stdout:
		command = shlex.split("bowtie2 --no-unal -p 16 -N 1 --local -x " + ALUY_consensus_bt+ " -U " + file_aluyf2)
		p = subprocess.Popen(command, stdout=stdout )
		p.communicate()
		stdout.close()

# convert consensus mapped bam to bed
# read 1 L1HS
	file_l1hsCon1_mapped=outputfolder + os.path.sep + outputfile_prefix + '_L1HS_con1_mapped.bam'
	with open(file_l1hsCon1_mapped, 'w') as stdout:
		command=shlex.split("samtools view -F4 -b "+ file_l1hsCon1)
		p = subprocess.Popen(command, stdout=stdout )
		p.communicate()
		stdout.close()
	file_l1hsCon1_mapped2=outputfolder + os.path.sep + outputfile_prefix + '_L1HS_con1_mapped.bed'
	with open(file_l1hsCon1_mapped2, 'w') as stdout:
		command=shlex.split("bedtools bamtobed -i " +file_l1hsCon1_mapped)
		p = subprocess.Popen(command, stdout=stdout )
		p.communicate()
		stdout.close()
# read 2 L1HS
	file_l1hsCon2_mapped=outputfolder + os.path.sep + outputfile_prefix + '_L1HS_con2_mapped.bam'
	with open(file_l1hsCon2_mapped, 'w') as stdout:
		command=shlex.split("samtools view -F4 -b "+ file_l1hsCon2)
		p = subprocess.Popen(command, stdout=stdout )
		p.communicate()
		stdout.close()
	file_l1hsCon2_mapped2=outputfolder + os.path.sep + outputfile_prefix + '_L1HS_con2_mapped.bed'
	with open(file_l1hsCon2_mapped2, 'w') as stdout:
		command=shlex.split("bedtools bamtobed -i " +file_l1hsCon2_mapped)
		p = subprocess.Popen(command, stdout=stdout )
		p.communicate()
		stdout.close()
# alu read 1

	file_aluyCon1_mapped=outputfolder + os.path.sep + outputfile_prefix + '_ALUY_con1_mapped.bam'
	with open(file_aluyCon1_mapped, 'w') as stdout:
		command=shlex.split("samtools view -F4 -b "+ file_aluyCon1)
		p = subprocess.Popen(command, stdout=stdout )
		p.communicate()
		stdout.close()
	file_aluyCon1_mapped2=outputfolder + os.path.sep + outputfile_prefix + '_ALUY_con1_mapped.bed'
	with open(file_aluyCon1_mapped2, 'w') as stdout:
		command=shlex.split("bedtools bamtobed -i " +file_aluyCon1_mapped)
		p = subprocess.Popen(command, stdout=stdout )
		p.communicate()
		stdout.close()

# alu read 2

	file_aluyCon2_mapped=outputfolder + os.path.sep + outputfile_prefix + '_ALUY_con2_mapped.bam'
	with open(file_aluyCon2_mapped, 'w') as stdout:
		command=shlex.split("samtools view -F4 -b "+ file_aluyCon2)
		p = subprocess.Popen(command, stdout=stdout )
		p.communicate()
		stdout.close()
	file_aluyCon2_mapped2=outputfolder + os.path.sep + outputfile_prefix + '_ALUY_con2_mapped.bed'
	with open(file_aluyCon2_mapped2, 'w') as stdout:
		command=shlex.split("bedtools bamtobed -i " +file_aluyCon2_mapped)
		p = subprocess.Popen(command, stdout=stdout )
		p.communicate()
		stdout.close()
	L1HS_con_1bp=outputfolder + os.path.sep + outputfile_prefix + '_L1HS_consensus1bp.bed'
	with open(L1HS_con_1bp, 'w') as stdout:
		command=shlex.split("intersectBed -a " + file_l1hsCon1_mapped2 + " -b " + L1HS_1stepbed)
		p = subprocess.Popen(command, stdout=stdout )
		p.communicate()
		command=shlex.split("intersectBed -a " + file_l1hsCon2_mapped2 + " -b " + L1HS_1stepbed)
		p = subprocess.Popen(command, stdout=stdout )
		p.communicate()
	stdout.close()
	ALUY_con_1bp=outputfolder + os.path.sep + outputfile_prefix + '_ALUY_consensus1bp.bed'
	with open(ALUY_con_1bp, 'w') as stdout:
		command=shlex.split("intersectBed -a " + file_aluyCon1_mapped2 + " -b " + ALUY_1stepbed)
		p = subprocess.Popen(command, stdout=stdout )
		p.communicate()
		command=shlex.split("intersectBed -a " + file_aluyCon2_mapped2 + " -b " + ALUY_1stepbed)
		p = subprocess.Popen(command, stdout=stdout )
		p.communicate()
	stdout.close()
# combine read 1 and 2 mapping to the consensus
# L1HS
	fileout= outputfolder + os.path.sep + outputfile_prefix + '_L1HS_consensus1bp_collapsed.bed'
	with open(fileout, 'w') as stdout:
		p1 = subprocess.Popen(['cut', '-f1', "-d/",L1HS_con_1bp], stdout = subprocess.PIPE)
		p2 = subprocess.Popen(['sort', '-u'], stdin = p1.stdout, stdout = stdout)
		p2.communicate()
	stdout.close()


##########
# storing the fraction count into each key (position on the consensus) L1HS
	frac_consensus={}
	fin=open(fileout,"r")
	for line in fin:
		line=line.strip("\r").strip("\n").split("\t")
		# line[2] is the position on the consensus
		if line[2] not in frac_consensus:
			frac_consensus[line[2]]=float(1/float(repeatdict[line[3]]))
		else:
			frac_consensus[line[2]]+=(float(1/float(repeatdict[line[3]])))
	fin.close()

	fileout= outputfolder + os.path.sep + outputfile_prefix + '_L1HS_multi_frac.txt'
	fout=open(fileout,'w')
	for pos in frac_consensus:
		#frac_consensus[pos]=sum(frac_consensus[pos])
		fout.write("L1HS"+"\t"+str(pos)+"\t"+str(frac_consensus[pos])+"\n")
	fout.close()



##########
# storing the fraction count into each key (position on the consensus) ALUY
	fileout= outputfolder + os.path.sep + outputfile_prefix + '_ALUY_consensus1bp_collapsed.bed'
# fileout looks like this
#AluY	208	 209	 E00492:134:HFCJMALXX:6:1101:1742:36715
	with open(fileout, 'w') as stdout:
		p1 = subprocess.Popen(['cut', '-f1', "-d/",ALUY_con_1bp], stdout = subprocess.PIPE)
		p2 = subprocess.Popen(['sort', '-u'], stdin = p1.stdout, stdout = stdout)
		p2.communicate()
	stdout.close()	

	frac_consensus={}
	fin=open(fileout,"r")
	for line in fin:
		line=line.strip("\r").strip("\n").split("\t")
		# line[2] is the position on the consensus
		if line[2] not in frac_consensus:
			frac_consensus[line[2]]=float(1/float(repeatdict[line[3]]))
		else:
			frac_consensus[line[2]]+=(float(1/float(repeatdict[line[3]])))
	fin.close()
	
	fileout= outputfolder + os.path.sep + outputfile_prefix + '_ALUY_multi_frac.txt'
	fout=open(fileout,'w')
	for pos in frac_consensus:
		#frac_consensus[pos]=sum(frac_consensus[pos])
		fout.write("AluY"+"\t"+str(pos)+"\t"+str(frac_consensus[pos])+"\n")
	fout.close()


############## END OF YT CHUNK FOR MULTI MAPPED READS
############################################
################################################################################
if paired_end == 'FALSE':
	if not os.path.exists(outputfolder + os.path.sep + 'pair1_'):
		os.mkdir(outputfolder + os.path.sep + 'pair1_')
	folder_pair1 = outputfolder + os.path.sep + 'pair1_'
################################################################################
	ps = []
	ticker= 0
	print "Processing repeat psuedogenomes..."
	for metagenome in repeat_list:
		metagenomepath = setup_folder + os.path.sep + metagenome
		file1=folder_pair1 + os.path.sep + metagenome + '.txt'
		with open(file1, 'w') as stdout:
			p = subprocess.Popen('bowtie2 ' + b_opt + ' ' + '-x ' + metagenomepath + ' ' + fastqfile_1 + ' | grep \"repname\" -', shell=True, stdout=stdout)
		ps.append(p)
		ticker +=1
		if ticker == cpus:
			for p in ps:
				p.communicate()
			ticker = 0
			ps = []
	if len(ps) > 0:
		for p in ps:
			p.communicate()
	stdout.close()
	
################################################################################
# combine the output from both read pairs:
	print 'Sorting and combining the output for both read pairs....'
	if not os.path.exists(outputfolder + os.path.sep + 'sorted_'):
		os.mkdir(outputfolder + os.path.sep + 'sorted_')
	sorted_ = outputfolder + os.path.sep + 'sorted_'
	for metagenome in repeat_list:
		file1 = folder_pair1 + os.path.sep + metagenome + '.txt'
		fileout= sorted_ + os.path.sep + metagenome + '.txt'
		with open(fileout, 'w') as stdout:
			p1 = subprocess.Popen(['cat',file1], stdout = subprocess.PIPE)
			p2 = subprocess.Popen(['cut', '-f1'], stdin = p1.stdout, stdout = subprocess.PIPE)
			p3 = subprocess.Popen(['cut', '-f1', "-d/"], stdin = p2.stdout, stdout = subprocess.PIPE)
			p4 = subprocess.Popen(['sort'], stdin = p3.stdout,stdout = subprocess.PIPE)
			p5 = subprocess.Popen(['uniq'], stdin = p4.stdout,stdout = stdout)
			p5.communicate()
		stdout.close()
	print 'completed ...'
	

####### YT added this chunk for single end 
################# yt added this chunk
# FOR MULTI MAPPED READS
### count how many repeats a read has mapped to
## fileout2 contains all readid compiled from alignments to all repeats. have to collapse to count the number of repeats it maps to

	# fileout2= outputfolder+ os.path.sep + 'allreadsid_repeat.txt'
	# print 'Counting how many repeats a read has mapped to...'
	# with open(fileout2, 'w') as stdout:
	# 	for metagenome in repeat_list:
	# 		file1 = folder_pair1 + os.path.sep + metagenome + '.txt'
	# 		p1 = subprocess.Popen(['cat',file1], stdout = subprocess.PIPE)
	# 		p2 = subprocess.Popen(['cut', '-f1',"-d "], stdin = p1.stdout, stdout = subprocess.PIPE)
	# 		p3 = subprocess.Popen(['cut', '-f1', "-d/"], stdin = p2.stdout, stdout = subprocess.PIPE)
	# 		p4 = subprocess.Popen(['sort'], stdin=p3.stdout, stdout = subprocess.PIPE)
	# 		p5 = subprocess.Popen(['uniq'], stdin=p4.stdout, stdout = subprocess.PIPE)
	# 		p6 = subprocess.Popen(['grep', '-v', "^@ "], stdin=p5.stdout, stdout = stdout)
	# 		p6.communicate()
	# stdout.close()
	# print 'completed ...'
	# subprocess.call(shlex.split(multi+" "+fileout2+" "+ outputfolder))
	# # fileout3 contains the number of repeats each read has mapped to
	# fileout3=outputfolder+os.path.sep+ "allreads_repeat.txt"

	# # store the number of times a read mapped to a repeat into a dictionary repeatdict
	# repeatdict={}
	# fin=open(fileout3,'r')
	# for line in fin:
	# 	line=line.strip("\r").strip("\n").split("\t")
	# 	if not line[0].startswith('@'):
	# 		if line[0] not in repeatdict:
	# 			repeatdict[line[0]] = int(line[1])
	# fin.close()

	# if os.path.exists(outputfolder+ os.path.sep + 'allreadsid_repeat.txt'):
	# 	os.remove(outputfolder+ os.path.sep + 'allreadsid_repeat.txt')

	# if os.path.exists(outputfolder+os.path.sep+ "allreads_repeat.txt"):
	# 	os.remove(outputfolder+os.path.sep+ "allreads_repeat.txt")
### count how many repeats a read has mapped to
	repeatdict = {}
	repeatdict1 = {}
	sorted_ = outputfolder + os.path.sep + 'sorted_'
	for rep in repeat_list:
		for data in import_text(sorted_ + os.path.sep + rep + '.txt', '\t'):
			repeatdict1[data[0]] = []

	for rep in repeat_list:
		for data in import_text(sorted_ + os.path.sep + rep + '.txt', '\t'):
			repeatdict1[data[0]].append(str(repeat_key[rep]))

	for readid,rep in repeatdict1.items():
		repeatdict[readid]=len(rep)

	file_l1hs1=outputfolder + os.path.sep + outputfile_prefix + '_L1HS1.bam'

	file_aluy1=outputfolder + os.path.sep + outputfile_prefix + '_ALUY1.bam'

	for metagenome in repeat_list:
		metagenomepath = setup_folder + os.path.sep + metagenome
		if metagenome =="L1HS":
			with open(file_l1hs1, 'w') as stdout:
				command = shlex.split("bowtie2 " + b_opt + " " + "-x " + metagenomepath + " " + fastqfile_1)
				p = subprocess.Popen(command, stdout=stdout)
				p.communicate()
			stdout.close()

		if metagenome =="AluY":
			with open(file_aluy1, 'w') as stdout:
				command = shlex.split("bowtie2 " + b_opt + " " + "-x " + metagenomepath + " " + fastqfile_1)
				p = subprocess.Popen(command, stdout=stdout)
				p.communicate()
			stdout.close()

	file_l1hs1_mapped=outputfolder + os.path.sep + outputfile_prefix + '_L1HS1_mapped.bam'
	with open(file_l1hs1_mapped, 'w') as stdout:
		#filter for mapped reads only
		command = shlex.split("samtools view -h -F 4 -b " + file_l1hs1)
		p = subprocess.Popen(command, stdout=stdout)
		p.communicate()
		stdout.close()

	file_aluy1_mapped=outputfolder + os.path.sep + outputfile_prefix + '_ALUY1_mapped.bam'
	with open(file_aluy1_mapped, 'w') as stdout:
		#filter for mapped reads only
		command = shlex.split("samtools view -h -F 4 -b " + file_aluy1)
		p = subprocess.Popen(command, stdout=stdout)
		p.communicate()
		stdout.close()

		
	if os.path.exists(outputfolder + os.path.sep + outputfile_prefix + '_L1HS1.bam'):
		os.remove(outputfolder + os.path.sep + outputfile_prefix + '_L1HS1.bam')

	if os.path.exists(outputfolder + os.path.sep + outputfile_prefix + '_ALUY1.bam'):
		os.remove(outputfolder + os.path.sep + outputfile_prefix + '_ALUY1.bam')

# convert aluy or l1hs mapped reads to the human genome , to fastq (to be  mapped to the consensus)
	file_l1hsf1=outputfolder + os.path.sep + outputfile_prefix + '_L1HS1.fastq'

	file_aluyf1=outputfolder + os.path.sep + outputfile_prefix + '_ALUY1.fastq'

	command = shlex.split("bedtools bamtofastq -i " +  file_l1hs1_mapped + " -fq " + file_l1hsf1)
	p=subprocess.Popen(command)
	p.communicate()
	



	command = shlex.split('bedtools bamtofastq -i '+ file_aluy1_mapped + ' -fq ' + file_aluyf1)
	p = subprocess.Popen(command)
	p.communicate()





# mapping reads that mapped to the pseudogenomes , to the consensus
# aligning to L1HS consensus, read 1 and read 2 separately
	file_l1hsCon1=outputfolder + os.path.sep + outputfile_prefix + '_L1HS_con1.sam'

	file_aluyCon1=outputfolder + os.path.sep + outputfile_prefix + '_ALUY_con1.sam'

	with open(file_l1hsCon1, 'w') as stdout:
		command = shlex.split("bowtie2 --no-unal -p 16 -N 1 --local -x " + L1HS_consensus_bt+ " -U " + file_l1hsf1)
		p = subprocess.Popen(command, stdout=stdout )
		p.communicate()
		stdout.close()

# aligning to AluY consensus
	with open(file_aluyCon1, 'w') as stdout:
		command = shlex.split("bowtie2 --no-unal -p 16 -N 1 --local -x " + ALUY_consensus_bt+ " -U " + file_aluyf1)
		p = subprocess.Popen(command, stdout=stdout )
		p.communicate()
		stdout.close()

# convert consensus mapped bam to bed
# read 1 L1HS
	file_l1hsCon1_mapped=outputfolder + os.path.sep + outputfile_prefix + '_L1HS_con1_mapped.bam'
	with open(file_l1hsCon1_mapped, 'w') as stdout:
		command=shlex.split("samtools view -F4 -b "+ file_l1hsCon1)
		p = subprocess.Popen(command, stdout=stdout )
		p.communicate()
		stdout.close()
	file_l1hsCon1_mapped2=outputfolder + os.path.sep + outputfile_prefix + '_L1HS_con1_mapped.bed'
	with open(file_l1hsCon1_mapped2, 'w') as stdout:
		command=shlex.split("bedtools bamtobed -i " +file_l1hsCon1_mapped)
		p = subprocess.Popen(command, stdout=stdout )
		p.communicate()
		stdout.close()

# alu read 1

	file_aluyCon1_mapped=outputfolder + os.path.sep + outputfile_prefix + '_ALUY_con1_mapped.bam'
	with open(file_aluyCon1_mapped, 'w') as stdout:
		command=shlex.split("samtools view -F4 -b "+ file_aluyCon1)
		p = subprocess.Popen(command, stdout=stdout )
		p.communicate()
		stdout.close()
	file_aluyCon1_mapped2=outputfolder + os.path.sep + outputfile_prefix + '_ALUY_con1_mapped.bed'
	with open(file_aluyCon1_mapped2, 'w') as stdout:
		command=shlex.split("bedtools bamtobed -i " +file_aluyCon1_mapped)
		p = subprocess.Popen(command, stdout=stdout )
		p.communicate()
		stdout.close()


	L1HS_con_1bp=outputfolder + os.path.sep + outputfile_prefix + '_L1HS_consensus1bp.bed'
	with open(L1HS_con_1bp, 'w') as stdout:
		command=shlex.split("intersectBed -a " + file_l1hsCon1_mapped2 + " -b " + L1HS_1stepbed)
		p = subprocess.Popen(command, stdout=stdout )
		p.communicate()
	stdout.close()
	ALUY_con_1bp=outputfolder + os.path.sep + outputfile_prefix + '_ALUY_consensus1bp.bed'
	with open(ALUY_con_1bp, 'w') as stdout:
		command=shlex.split("intersectBed -a " + file_aluyCon1_mapped2 + " -b " + ALUY_1stepbed)
		p = subprocess.Popen(command, stdout=stdout )
		p.communicate()
	stdout.close()
# combine read 1 and 2 mapping to the consensus
# L1HS
	fileout= outputfolder + os.path.sep + outputfile_prefix + '_L1HS_consensus1bp_collapsed.bed'
	with open(fileout, 'w') as stdout:
		p1 = subprocess.Popen(['cut', '-f1', "-d/",L1HS_con_1bp], stdout = subprocess.PIPE)
		p2 = subprocess.Popen(['sort', '-u'], stdin = p1.stdout, stdout = stdout)
		p2.communicate()
	stdout.close()


##########
# storing the fraction count into each key (position on the consensus) L1HS
	frac_consensus={}
	fin=open(fileout,"r")
	for line in fin:
		line=line.strip("\r").strip("\n").split("\t")
		# line[2] is the position on the consensus
		if line[2] not in frac_consensus:
			frac_consensus[line[2]]=float(1/float(repeatdict[line[3]]))
		else:
			frac_consensus[line[2]]+=(float(1/float(repeatdict[line[3]])))
	fin.close()

	fileout= outputfolder + os.path.sep + outputfile_prefix + '_L1HS_multi_frac.txt'
	fout=open(fileout,'w')
	for pos in frac_consensus:
		#frac_consensus[pos]=sum(frac_consensus[pos])
		fout.write("L1HS"+"\t"+str(pos)+"\t"+str(frac_consensus[pos])+"\n")
	fout.close()



##########
# storing the fraction count into each key (position on the consensus) ALUY
	fileout= outputfolder + os.path.sep + outputfile_prefix + '_ALUY_consensus1bp_collapsed.bed'
# fileout looks like this
#AluY	208	 209	 E00492:134:HFCJMALXX:6:1101:1742:36715
	with open(fileout, 'w') as stdout:
		p1 = subprocess.Popen(['cut', '-f1', "-d/",ALUY_con_1bp], stdout = subprocess.PIPE)
		p2 = subprocess.Popen(['sort', '-u'], stdin = p1.stdout, stdout = stdout)
		p2.communicate()
	stdout.close()	

	frac_consensus={}
	fin=open(fileout,"r")
	for line in fin:
		line=line.strip("\r").strip("\n").split("\t")
		# line[2] is the position on the consensus
		if line[2] not in frac_consensus:
			frac_consensus[line[2]]=float(1/float(repeatdict[line[3]]))
		else:
			frac_consensus[line[2]]+=(float(1/float(repeatdict[line[3]])))
	fin.close()
	
	fileout= outputfolder + os.path.sep + outputfile_prefix + '_ALUY_multi_frac.txt'
	fout=open(fileout,'w')
	for pos in frac_consensus:
		#frac_consensus[pos]=sum(frac_consensus[pos])
		fout.write("AluY"+"\t"+str(pos)+"\t"+str(frac_consensus[pos])+"\n")
	fout.close()


############## END OF YT CHUNK FOR MULTI MAPPED READS SINGLE END


################################################################################
# build a file of repeat keys for all reads
print 'Writing and processing intermediate files...'
sorted_ = outputfolder + os.path.sep + 'sorted_'
readid = {}
sumofrepeatreads=0
for rep in repeat_list: 
	for data in import_text(sorted_ + os.path.sep + rep + '.txt', '\t'):
		readid[data[0]] = ''
for rep in repeat_list: 
	for data in import_text(sorted_ + os.path.sep + rep + '.txt', '\t'):	
		readid[data[0]]+=str(repeat_key[rep]) + str(',')
for subfamilies in readid.values():
	if not counts.has_key(subfamilies):
		counts[subfamilies]=0
	counts[subfamilies] +=1
	sumofrepeatreads+=1
del readid
print 'Identified ' + str(sumofrepeatreads) + ' reads that mapped to repeats for unique and multimappers.'

################################################################################
print "Conducting final calculations..."
# build a converter to numeric label for repeat and yield a combined list of repnames seperated by backslash
def convert(x):
	x = x.strip(',')
	x = x.split(',')
	global repname
	repname = ""
	for i in x:
		repname = repname + os.path.sep + rev_repeat_key[int(i)] 
# building the total counts for repeat element enrichment...
for x in counts.keys():
	count= counts[x]
	x = x.strip(',')
	x = x.split(',')
	for i in x:
		reptotalcounts[rev_repeat_key[int(i)]] += int(count)
# building the fractional counts for repeat element enrichment...

for x in counts.keys():
	count= counts[x]
	x = x.strip(',')
	x = x.split(',')
	splits = len(x)
	for i in x:
		fractionalcounts[rev_repeat_key[int(i)]] += float(numpy.divide(float(count),float(splits)))




# building categorized table of repeat element enrichment... 
repcounts = {}
repcounts['other'] = 0
for key in counts.keys():
	convert(key)
	repcounts[repname] = counts[key]
# building the total counts for class enrichment...
for key in reptotalcounts.keys():
	classtotalcounts[repeatclass[key]] += reptotalcounts[key]
# building total counts for family enrichment...
for key in reptotalcounts.keys():
	familytotalcounts[repeatfamily[key]] += reptotalcounts[key]
# building unique counts table'
repcounts2 = {}
for rep in repeat_list:
	if repcounts.has_key("/" +rep):
		repcounts2[rep] = repcounts["/" +rep]
	else:
		repcounts2[rep] = 0
# building the fractionalcounts counts for class enrichment...
for key in fractionalcounts.keys():
	classfractionalcounts[repeatclass[key]] += fractionalcounts[key]
# building fractional counts for family enrichment...
for key in fractionalcounts.keys():
	familyfractionalcounts[repeatfamily[key]] += fractionalcounts[key]   
	
################################################################################
print 'Writing final output and removing intermediate files...'
# print output to file of the categorized counts and total overlapping counts:
if allcountmethod  == "TRUE":
	fout1 = open(outputfolder + os.path.sep + outputfile_prefix + '_total_counts.txt' , 'w')
	for key in reptotalcounts.keys():
		print >> fout1, str(key) + '\t' + repeatclass[key] + '\t' + repeatfamily[key] + '\t' + str(reptotalcounts[key])
	fout2 = open(outputfolder + os.path.sep + outputfile_prefix + '_class_total_counts.txt' , 'w')
	for key in classtotalcounts.keys():
		print >> fout2, str(key) + '\t' + str(classtotalcounts[key])  
	fout3 = open(outputfolder + os.path.sep + outputfile_prefix + '_family_total_counts.txt' , 'w')
	for key in familytotalcounts.keys():
		print >> fout3, str(key) + '\t' + str(familytotalcounts[key])  
	fout4 = open(outputfolder + os.path.sep + outputfile_prefix + '_unique_counts.txt' , 'w')
	for key in repcounts2.keys():
		print >> fout4, str(key) + '\t' + repeatclass[key] + '\t' + repeatfamily[key] + '\t' + str(repcounts2[key])	 
		fout5 = open(outputfolder + os.path.sep + outputfile_prefix + '_class_fraction_counts.txt' , 'w')
	for key in classfractionalcounts.keys():
		print >> fout5, str(key) + '\t' + str(classfractionalcounts[key])  
	fout6 = open(outputfolder + os.path.sep + outputfile_prefix + '_family_fraction_counts.txt' , 'w')
	for key in familyfractionalcounts.keys():
		print >> fout6, str(key) + '\t' + str(familyfractionalcounts[key])
	fout7 = open(outputfolder + os.path.sep + outputfile_prefix + '_fraction_counts.txt' , 'w')
	for key in fractionalcounts.keys():
		print >> fout7, str(key) + '\t' + repeatclass[key] + '\t' + repeatfamily[key] + '\t' + str(int(fractionalcounts[key]))
		fout1.close()
	fout2.close()
	fout3.close()
	fout4.close()
	fout5.close()
	fout6.close()
	fout7.close()
else:
	fout1 = open(outputfolder + os.path.sep + outputfile_prefix + '_class_fraction_counts.txt' , 'w')
	for key in classfractionalcounts.keys():
		print >> fout1, str(key) + '\t' + str(classfractionalcounts[key])  
	fout2 = open(outputfolder + os.path.sep + outputfile_prefix + '_family_fraction_counts.txt' , 'w')
	for key in familyfractionalcounts.keys():
		print >> fout2, str(key) + '\t' + str(familyfractionalcounts[key])
	fout3 = open(outputfolder + os.path.sep + outputfile_prefix + '_fraction_counts.txt' , 'w')
	for key in fractionalcounts.keys():
		print >> fout3, str(key) + '\t' + repeatclass[key] + '\t' + repeatfamily[key] + '\t' + str(int(fractionalcounts[key]))
	fout1.close()
	fout2.close()
	fout3.close()
	
################################################################################
# Remove Large intermediate files
if os.path.exists(outputfolder + os.path.sep + outputfile_prefix + '_regionsorter.txt'):
	os.remove(outputfolder + os.path.sep + outputfile_prefix + '_regionsorter.txt')
if os.path.exists(outputfolder + os.path.sep + 'pair1_'):
	shutil.rmtree(outputfolder + os.path.sep + 'pair1_')
if os.path.exists(outputfolder + os.path.sep + 'pair2_'):
	shutil.rmtree(outputfolder + os.path.sep + 'pair2_')
if os.path.exists(outputfolder + os.path.sep + 'sorted_'):
	shutil.rmtree(outputfolder + os.path.sep + 'sorted_')

if os.path.exists(outputfolder + os.path.sep + outputfile_prefix + '_L1HS1_mapped.bam'):
	os.remove(outputfolder + os.path.sep + outputfile_prefix + '_L1HS1_mapped.bam')
if os.path.exists(outputfolder + os.path.sep + outputfile_prefix + '_L1HS2_mapped.bam'):
	os.remove(outputfolder + os.path.sep + outputfile_prefix + '_L1HS2_mapped.bam')
if os.path.exists(outputfolder + os.path.sep + outputfile_prefix + '_ALUY1_mapped.bam'):
	os.remove(outputfolder + os.path.sep + outputfile_prefix + '_ALUY1_mapped.bam')
if os.path.exists(outputfolder + os.path.sep + outputfile_prefix + '_ALUY2_mapped.bam'):
	os.remove(outputfolder + os.path.sep + outputfile_prefix + '_ALUY2_mapped.bam')


if os.path.exists(outputfolder + os.path.sep + outputfile_prefix + '_L1HS1.bam'):
	os.remove(outputfolder + os.path.sep + outputfile_prefix + '_L1HS1.bam')
if os.path.exists(outputfolder + os.path.sep + outputfile_prefix + '_L1HS2.bam'):
	os.remove(outputfolder + os.path.sep + outputfile_prefix + '_L1HS2.bam')
if os.path.exists(outputfolder + os.path.sep + outputfile_prefix + '_ALUY1.bam'):
	os.remove(outputfolder + os.path.sep + outputfile_prefix + '_ALUY1.bam')
if os.path.exists(outputfolder + os.path.sep + outputfile_prefix + '_ALUY2.bam'):
	os.remove(outputfolder + os.path.sep + outputfile_prefix + '_ALUY2.bam')

if os.path.exists(outputfolder + os.path.sep + outputfile_prefix + '_L1HS1.fastq'):
	os.remove(outputfolder + os.path.sep + outputfile_prefix + '_L1HS1.fastq')
if os.path.exists(outputfolder + os.path.sep + outputfile_prefix + '_L1HS2.fastq'):
	os.remove(outputfolder + os.path.sep + outputfile_prefix + '_L1HS2.fastq')
if os.path.exists(outputfolder + os.path.sep + outputfile_prefix + '_ALUY1.fastq'):
	os.remove(outputfolder + os.path.sep + outputfile_prefix + '_ALUY1.fastq')
if os.path.exists(outputfolder + os.path.sep + outputfile_prefix + '_ALUY2.fastq'):
	os.remove(outputfolder + os.path.sep + outputfile_prefix + '_ALUY2.fastq')


if os.path.exists(outputfolder + os.path.sep + outputfile_prefix + '_L1HS_con1.sam'):
	os.remove(outputfolder + os.path.sep + outputfile_prefix + '_L1HS_con1.sam')
if os.path.exists(outputfolder + os.path.sep + outputfile_prefix + '_L1HS_con2.sam'):
	os.remove(outputfolder + os.path.sep + outputfile_prefix + '_L1HS_con2.sam')
if os.path.exists(outputfolder + os.path.sep + outputfile_prefix + '_ALUY_con1.sam'):
	os.remove(outputfolder + os.path.sep + outputfile_prefix + '_ALUY_con1.sam')
if os.path.exists(outputfolder + os.path.sep + outputfile_prefix + '_ALUY_con2.sam'):
	os.remove(outputfolder + os.path.sep + outputfile_prefix + '_ALUY_con2.sam')


if os.path.exists(outputfolder + os.path.sep + outputfile_prefix + '_L1HS_con1_mapped.bam'):
	os.remove(outputfolder + os.path.sep + outputfile_prefix + '_L1HS_con1_mapped.bam')
if os.path.exists(outputfolder + os.path.sep + outputfile_prefix + '_L1HS_con1_mapped.bed'):
	os.remove(outputfolder + os.path.sep + outputfile_prefix + '_L1HS_con1_mapped.bed')
if os.path.exists(outputfolder + os.path.sep + outputfile_prefix + '_L1HS_con2_mapped.bam'):
	os.remove(outputfolder + os.path.sep + outputfile_prefix + '_L1HS_con2_mapped.bam')
if os.path.exists(outputfolder + os.path.sep + outputfile_prefix + '_L1HS_con2_mapped.bed'):
	os.remove(outputfolder + os.path.sep + outputfile_prefix + '_L1HS_con2_mapped.bed')
if os.path.exists(outputfolder + os.path.sep + outputfile_prefix + '_ALUY_con1_mapped.bam'):
	os.remove(outputfolder + os.path.sep + outputfile_prefix + '_ALUY_con1_mapped.bam')
if os.path.exists(outputfolder + os.path.sep + outputfile_prefix + '_ALUY_con1_mapped.bed'):
	os.remove(outputfolder + os.path.sep + outputfile_prefix + '_ALUY_con1_mapped.bed')
if os.path.exists(outputfolder + os.path.sep + outputfile_prefix + '_ALUY_con2_mapped.bam'):
	os.remove(outputfolder + os.path.sep + outputfile_prefix + '_ALUY_con2_mapped.bam')
if os.path.exists(outputfolder + os.path.sep + outputfile_prefix + '_ALUY_con2_mapped.bed'):
	os.remove(outputfolder + os.path.sep + outputfile_prefix + '_ALUY_con2_mapped.bed')

if os.path.exists(outputfolder + os.path.sep + outputfile_prefix + '_L1HS_consensus1bp.bed'):
	os.remove(outputfolder + os.path.sep + outputfile_prefix + '_L1HS_consensus1bp.bed')

if os.path.exists(outputfolder + os.path.sep + outputfile_prefix + '_ALUY_consensus1bp.bed'):
	os.remove(outputfolder + os.path.sep + outputfile_prefix + '_ALUY_consensus1bp.bed')

if os.path.exists(outputfolder + os.path.sep + outputfile_prefix + '_L1HS_consensus1bp_collapsed.bed'):
	os.remove(outputfolder + os.path.sep + outputfile_prefix + '_L1HS_consensus1bp_collapsed.bed')
if os.path.exists(outputfolder + os.path.sep + outputfile_prefix + '_ALUY_consensus1bp_collapsed.bed'):
	os.remove(outputfolder + os.path.sep + outputfile_prefix + '_ALUY_consensus1bp_collapsed.bed')

if os.path.exists(outputfolder+ os.path.sep + 'allreadsid_repeat.txt'):
	os.remove(outputfolder+ os.path.sep + 'allreadsid_repeat.txt')

if os.path.exists(outputfolder+os.path.sep+ "allreads_repeat.txt"):
	os.remove(outputfolder+os.path.sep+ "allreads_repeat.txt")
if os.path.exists(outputfolder + os.path.sep + outputfile_prefix + '_ytcount.txt'):
	os.remove(outputfolder + os.path.sep + outputfile_prefix + '_ytcount.txt')
print "... Done"