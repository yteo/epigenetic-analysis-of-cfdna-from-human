# Epigenetic analysis of cfDNA from human

The following sections contain tutorials on how to analyze the cfDNA data from

Teo, Y. V. et al. Cell-free DNA as a biomarker of aging. Aging Cell, e12890, doi:10.1111/acel.12890 (2018).

###Tutorial###

**Dependencies**

The following versions are tested for the pipeline:
samtools/0.1.18
bwa/0.7.15
bedtools/2.25.0
bowtie2/2.3.0
```
module load samtools/0.1.18  
module load bwa/0.7.15  
module load R  
module load bedtools/2.25.0  
module load bowtie2/2.3.0
```

Alignment of cfDNA fastq file to hs37d5.fa using BWA-MEM (set the $outdir and $sample accordingly)

```
bwa mem -t 16 -M $index ${sample}_1.fq.gz ${sample}_2.fq.gz > ${outdir}/${sample}_bwa_hs37d5.sam
```
Sort, remove duplicates and index
```
samtools view -bS ${outdir}/${sample}_bwa_hs37d5.sam|samtools sort -o ${outdir}/${sample}_bwa_hs37d5_sorted.bam  
samtools index ${outdir}/${sample}_bwa_hs37d5_sorted.bam ${outdir}/${sample}_bwa_hs37d5_sorted.bam.bai  

java -jar picard.jar MarkDuplicates I=${outdir}/${sample}_bwa_hs37d5_sorted.bam \  
REMOVE_DUPLICATES=true ASSUME_SORTED=coordinate METRICS_FILE=${outdir}/${sample}_metrics.txt \  
OUTPUT=${outdir}/${sample}_bwa_hs37d5_sorted_removedup.bam  


samtools index ${outdir}/${sample}_bwa_hs37d5_sorted_removedup.bam \  
${outdir}/${sample}_bwa_hs37d5_sorted_removedup.bam.bai  
```

Filter bam file
Remove hs37d5 decoy, MT, X, Y and human virus reads for danpos analysis
Remove secondary alignment and retain only proper reads
```
samtools view -b ${outdir}/${sample}_bwa_hs37d5_sorted_removedup.bam 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 GL000191.1 GL000192.1 GL000193.1 GL000194.1 GL000195.1 GL000196.1 GL000197.1 GL000198.1 GL000199.1 GL000200.1 GL000201.1 GL000202.1 GL000203.1 GL000204.1 GL000205.1 GL000206.1 GL000207.1 GL000208.1 GL000209.1 GL000210.1 GL000211.1 GL000212.1 GL000213.1 GL000214.1 GL000215.1 GL000216.1 GL000217.1 GL000218.1 GL000219.1 GL000220.1 GL000221.1 GL000222.1 GL000223.1 GL000224.1 GL000225.1 GL000226.1 GL000227.1 GL000228.1 GL000229.1 GL000230.1 GL000231.1 GL000232.1 GL000233.1 GL000234.1 GL000235.1 GL000236.1 GL000237.1 GL000238.1 GL000239.1 GL000240.1 GL000241.1 GL000242.1 GL000243.1 GL000244.1 GL000245.1 GL000246.1 GL000247.1 GL000248.1 GL000249.1 -F 256 -f2 -q 1|samtools sort -o ${outdir}/${sample}_bwa_hs37d5_danpos_sorted_proper_rmdup.bam  
samtools index ${outdir}/${sample}_bwa_hs37d5_danpos_sorted_proper_rmdup.bam

```

Remove reads that were soft-clipped > 75bp
```
python filter_cigar.py ${outdir}/${sample}_bwa_hs37d5_danpos_sorted_proper_rmdup.bam ${outdir}/${sample}_bwa_hs37d5_danpos_sorted_proper_rmdup_cigar.bam
```

Determine insert size (Figure 1) 
```
file=${outdir}/${sample}_bwa_hs37d5_danpos_sorted_proper_rmdup_cigar.bam  
java -jar picard.jar CollectInsertSizeMetrics INPUT=${file} \  
OUTPUT=${file%_bwa_hs37d5_danpos_sorted_proper_rmdup_cigar.bam}_insertsize.txt HISTOGRAM_FILE=${file%_bwa_hs37d5_danpos_sorted_proper_rmdup_cigar.bam}_insertsizeHIST.pdf

```
Identify cfDNA signals for every 10bp of the genome using danpos-2.2.2 **(Figure 1)**
```
file=${outdir}/${sample}_bwa_hs37d5_danpos_sorted_proper_rmdup_cigar.bam  
python danpos.py dpos $file \  
-m 1 -o danposresult10bp
```
 
Identify nucleosome signals at TSS, TTS and CTCF (Adjust sample output from dpos accordingly) **(Figure 4)**
ENCFF833FTF_availablechr.bed contains CTCF binding sites of GM12878 (obtained from GEO: GSM733752)
```
python danpos.py profile danposresult10bp/pooled/${sample}.Fnor.smooth.wig --genefile_paths ./Data/TranscriptAnno-GRCh37.75_availablechr_TSS_TTS.txt \  
--genomic_sites TSS --name TSS --heatmap 0  

python danpos.py profile danposresult10bp/pooled/${sample}.Fnor.smooth.wig --genefile_paths ./Data/TranscriptAnno-GRCh37.75_availablechr_TSS_TTS.txt \  
--genomic_sites TTS --name TTS --heatmap 0  

python danpos.py profile danposresult10bp/pooled/${sample}.Fnor.smooth.wig --genefile_paths ./Data/ENCFF833FTF_availablechr.bed \  
--genomic_sites center,region --name CTCF --heatmap 0 
```

Downsample reads to look at 100kb regions between different samples (based on the paper, we downsampled to 46301323 reads)
```
java -jar picard.jar SortSam I=${file} O=${outdir}/${sample}s_sorted_proper_rmdup_tmp.bam SORT_ORDER=queryname  

samtools view ${outdir}/${sample}_sorted_proper_rmdup_tmp.bam|cut -f1 |sort -u|shuf|head -n 46301323 > ${outdir}/${sample}_tmp  

java -jar picard.jar FilterSamReads \  
INPUT=${outdir}/${sample}_sorted_proper_rmdup_tmp.bam FILTER=includeReadList READ_LIST_FILE=${outdir}/${sample}_tmp OUTPUT=${outdir}/${sample}_downsampled_46301323.bam  
samtools sort ${outdir}/${sample}_downsampled_46301323.bam -o ${outdir}/${sample}_downsampled_46301323_sorted.bam

```
Intersect the coverage of cfDNA with the 100kb subcompartments obtained from GM12878 HiC dataset (GSE63525_GM12878_subcompartments.100kb.sorted.bed is from https://github.com/shendurelab/cfDNA/tree/master/peak_density)
```
bedtools intersect -a ${outdir}/${sample}_downsampled_46301323_sorted.bam -b <(sed 's/chr//' ./Data/GSE63525_GM12878_subcompartments.100kb.sorted.bed) -wb -bed \  
|sed 's/\/1//'|sed 's/\/2//'|cut -f4,13,14,15,16|sort -u|cut -f2,3,4,5|sort|uniq -c|tr '[:blank:]' \\t |sed 's/^[ \t]*//;s/[ \t]*$//' > ${outdir}/${sample}_full_compartment.bed
```

Plots from **Figure 2 and 3** are obtained with this example script: CompiledRCodeFig2and3.R

   
The python scripts used in **Figure 5** is modified using the original script from (https://github.com/nerettilab/RepEnrich2)to suit the analysis of cfDNA dataset, specifically focused on L1HS and AluY repeat elements.
hg19_repeatmasker.txt and RepEnrich2_setup_hg19/ can be prepared and setup by following https://github.com/nerettilab/RepEnrich2
```
bowtie2 -q -p 20 -x ${hs37d5Index} -1 ${sample}_1.fq.gz -2 ${sample}_2.fq.gz -S ${sample}.sam  
samtools view -bS ${sample}.sam > ${sample}.bam  
rm /users/yteo/data/yteo/BGI_cfDNA_cleanrun/TE/${sample}.sam  
python ModifiedRepEnrich2/RepEnrich2_subset_YT.py ${sample}.bam 30 ${file} --pairedend TRUE  
python ModifiedRepEnrich2/RepEnrich2_YT.py hg19_repeatmasker.txt ${sample} ${sample} RepEnrich2_setup_hg19/ ${file}_multimap_R1.fastq --fastqfile2 ${sample}_multimap_R2.fastq ${sample}_unique.bam --cpus 20 --pairedend TRUE
```
You'll find three output files for each repeat element  
\*_ALUY_multi_frac.txt  
\*_ALUY_consensuscoord_count.txt  
\*_ALUY_consensuscoord_num.txt  
\*_L1HS_multi_frac.txt  
\*_L1HS_consensuscoord_count.txt  
\*_L1HS_consensuscoord_num.txt  

Subsequently, run RFig5.R  
 
The tissues-of-origin (**Figure 6**) was identified using the method in  
Snyder, M. W., Kircher, M., Hill, A. J., Daza, R. M., & Shendure, J. (2016). Cell‐free DNA comprises an in vivo nucleosome footprint that informs its tissues‐of‐origin. Cell, 164(1–2), 57–68. https://doi.org/10. 1016/j.cell.2015.11.050  
github link: https://github.com/shendurelab/cfDNA  

We replaced the gene expression of tissues from GTEx instead of cancer cell lines. The sample tissues were filtered as explained in Teo et al. 2018  


