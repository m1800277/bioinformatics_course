#!/bin/bash

trimmomatic PE \
-threads 4 \
-phred33 \
$1 $2 \
-baseout ~/ngs_course/dnaseq_pipeline/data/trimmed_fastq/trimmed_data \
ILLUMINACLIP:/home/ubuntu/anaconda3/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic-0.39-2/adapters/NexteraPE-PE.fa:2:30:10 \
TRAILING:25 MINLEN:50
# Trimmomatic used to remove sequencer specific reads below a predetermined threshold, making it easier to align the sequence to the transcriptome since poor quality bases are removed
# Options used in ths command:
# -threads determines the number of processors Trimmomatic runs, in this case 4
# -phred refers to the quality of the reads,either being -phred33 or -phred 64, in this case -phred33
# $1 and $2 is used to indicate that the two inputs for the trimmomatic command are the two files put when entering the bash command to start the pipeline, i.e., the NGS0001.R1.fastq.gz and NGS0001.R2.fastq.gz
# -baseout is the location of the output, i.e., the trimmed reads. In this case, they are output into the trimmed_fastq subdirectory
# ILLUMINACLIP is used to remove the sequences from the read that are from adapters or other illumina specific sequences
# TRAILING and LEADINGare two options used to remove bases off of the end and begining of the read respecively, only if they are below the threshold quality specified, in this case 25 and 50 respectively

fastqc -t 4 /home/ubuntu/ngs_course/dnaseq_pipeline/data/trimmed_fastq/trimmed_data_1P \
/home/ubuntu/ngs_course/dnaseq_pipeline/data/trimmed_fastq/trimmed_data_2P
# fastqc is used to perform quality control on the trimmed reads, allowng for a quick assessment on whether further analysis of the read is worth doing
# Options used in this command:
# -t determines the amount of threads to be processed simultaneously, in this case 4
# fastqc is being performed on the 2 trimmed paired reads, with their full directory being provided

mkdir ~/ngs_course/dnaseq_pipeline/data/aligned_data
# making a new directory for the aligned data files, not necessary but is useful for data management

bwa mem -t 4 -v 1 -R '@RG\tID:HWI-D0011.50.H7AP8ADXX.1.NGS0001\tSM:NGS0001\tPL:ILLUMINA\tLB:nextera-ngs0001-blood\tDT:2017-02-23\tPU:HWI-D00119' -I 250,50 ~/ngs_course/dnaseq_pipeline/data/reference/hg19.fa.gz ~/ngs_course/dnaseq_pipeline/data/trimmed_fastq/trimmed_data_1P ~/ngs_course/dnaseq_pipeline/data/trimmed_fastq/trimmed_data_2P > ~/ngs_course/dnaseq_pipeline/data/aligned_data/NGS0001.sam
# bwa is used for read alignments, and bwa mem is prefered since it performs better
# Options used:
# -t used to specify the number of threads, in this case 4
# -v used for the verbosity level, in this case 1 for error
# -R used for the read group header information, with the following headers:
# tID is the read group identifier, used to tag read group information and must be unique
# tSM is the name of the sample sequence, where samples wit the same name are treated as being the same sample
# tPL is the platform used to produce the reads, in this case ILLUMINA
# tLB is the identity of the DNA preparation library, used to determine whether marked duplicates are present
# tDT is the date
# tPU is the platform unit, with information on the flowcell barcode, lane, and a lbrary specific identifier
# -I is used to specify the mean and standard deviation, in this case 250 and 50 respectively
# The reference genome used is the hg19.fa.gz file, which has already been indexed as part of the project setup phase
# The input files are the quality trimmed paired reads
# The output is a .sam file, which is a Sequence Alignment Map

rm ~/ngs_course/dnaseq_pipeline/data/trimmed_fastq/*data*
# Not a necessary step, but is done to conserve limited allocated memory. These files will no longer be required to finish the pipeline

samtools view -h -b -S ~/ngs_course/dnaseq_pipeline/data/aligned_data/NGS0001.sam > ~/ngs_course/dnaseq_pipeline/data/aligned_data/NGS0001.bam
# Converting the .sam file to a .bam file, the binary version of a .sam file
# Options used:
# -h used to include a header
# -b used to indicate that it is converting to a .bam file
# -S not necessary, used just in case for compatibility reasons
# .sam file is input, .bam file is output

rm ~/ngs_course/dnaseq_pipeline/data/aligned_data/NGS0001.sam
# Not necessary, but removing the .sam file to save memory, as the .sam file can be quite large and it is not needed anymore

samtools sort ~/ngs_course/dnaseq_pipeline/data/aligned_data/NGS0001.bam > ~/ngs_course/dnaseq_pipeline/data/aligned_data/NGS0001_sorted.bam
# sort is used to align the data by the leftmost coordinates
# .bam file is input, aligned .bam file is output

rm ~/ngs_course/dnaseq_pipeline/data/aligned_data/NGS0001.bam
# not necessary, but removng unaligned .bam file to save memory, as it is not useful anymore

samtools index ~/ngs_course/dnaseq_pipeline/data/aligned_data/NGS0001_sorted.bam
# index used to index and sort the .bam file
# aligned .bam is input, .bai files are output

picard MarkDuplicates -I ~/ngs_course/dnaseq_pipeline/data/aligned_data/NGS0001_sorted.bam -O ~/ngs_course/dnaseq_pipeline/data/aligned_data/NGS0001_sorted_marked.bam -M ~/ngs_course/dnaseq_pipeline/data/aligned_data/marked_dup_metrics.txt
# Marking duplicated reads from the sorted .bam file using picard, by examining aligned records from a library to the.bam file, locating any duplicate reads
# Options used:
# -I to indicate input file, in this case the sorted .bam file
# -O to indicate the output file, in this case a sorted marked .bam file
# -M to indicate the .txt file where marked duplicates are listed

rm ~/ngs_course/dnaseq_pipeline/data/aligned_data/NGS0001_sorted.bam
# removing the obsolete sorted .bam file in order to conserve memory

samtools index ~/ngs_course/dnaseq_pipeline/data/aligned_data/NGS0001_sorted_marked.bam
# index used to index and sort the sorted marked .bam file

samtools view -F 1796 -q 20 -o ~/ngs_course/dnaseq_pipeline/data/aligned_data/NGS0001_sorted_filtered.bam ~/ngs_course/dnaseq_pipeline/data/aligned_data/NGS0001_sorted_marked.bam
# Filtering out the duplicate reads using samtools view
# Options used:
# -F used to indicate that any alignments flagged will not be output. The 1796 used to flag 4 criteria; unmapped reads, non-primary alignments, failed quality check reads, and duplicate reads
# -q used to set the minimum MAPQ score to 20
# -o used to indicate the output file, in this case a sorted filtered .bam file
# Input file is the sorted marked .bam file

rm ~/ngs_course/dnaseq_pipeline/data/aligned_data/NGS0001_sorted_marked.bam
# removing the sorted marked file in order to save memory

samtools index ~/ngs_course/dnaseq_pipeline/data/aligned_data/NGS0001_sorted_filtered.bam
# Index used to index and sort the filtered .bam file

samtools flagstat ~/ngs_course/dnaseq_pipeline/data/aligned_data/NGS0001_sorted_filtered.bam
# flagstat is used to count the number of read alignments for eah of the 13 different flag types
# input is the filtered .bam file

samtools idxstats ~/ngs_course/dnaseq_pipeline/data/aligned_data/NGS0001_sorted_filtered.bam
# idxstats is used to report the summary statistics of the allignment, by printing the contents of the index file
# Input is the filtered .bam file

bedtools coverage -a ~/ngs_course/dnaseq_pipeline/data/aligned_data/NGS0001_sorted_filtered.bam \
-b ~/ngs_course/dnaseq_pipeline/data/annotation.bed
# bedtools coverage is used to generate depth of coverage statistics
# Options used:
# -a is the filtered .bam file as an input
# -b is the .bed file which the .bm file is to be compared with

picard CollectInsertSizeMetrics -I ~/ngs_course/dnaseq_pipeline/data/aligned_data/NGS0001_sorted_filtered.bam -O ~/ngs_course/dnaseq_pipeline/data/aligned_data/insert_size_metrics.txt -H ~/ngs_course/dnaseq_pipeline/data/aligned_data/insert_size_histogram.pdf -M 0.5
# CollectInsetSizeMetrics is used to generate the standard alignment statistics related to the insert size
# Options used:
# -I used to indicate the input file, in this case the filtered .bam file
# -O used to indicate the output file, in this case a .txt file containing the metrics
# -H is used to indicate the output histogram, a .pdf file
# -M is used to indicate the minimum percentage of data displayed in the histogram, in ths case 50%

zcat ~/ngs_course/dnaseq_pipeline/data/reference/hg19.fa.gz > ~/ngs_course/dnaseq_pipeline/data/reference/hg19.fa
# Uncompressing the reference data, it was initially compressed to save memory

samtools faidx ~/ngs_course/dnaseq_pipeline/data/reference/hg19.fa
# Indexing the uncompressed reference

freebayes --bam ~/ngs_course/dnaseq_pipeline/data/aligned_data/NGS0001_sorted_filtered.bam --fasta-reference ~/ngs_course/dnaseq_pipeline/data/reference/hg19.fa --vcf ~/ngs_course/dnaseq_pipeline/results/NGS0001.vcf
# freebayes used to call the variants of the data by comparing it to the reference data
# Options used:
# --bam used to set the .bam file where the variants will be compared
# --fasta-reference sets the freference data file,, in this case the uncompressed hg19.fa file
# --vcf sets the output file in a *.vcf file format, which is a Variant Call File

bgzip ~/ngs_course/dnaseq_pipeline/results/NGS0001.vcf
# Compressing the .vcf file, preparing it for use further down the pipeline

tabix -p vcf ~/ngs_course/dnaseq_pipeline/results/NGS0001.vcf.gz
# tabix used to index the compressed .vcf file
# Opeions used:
# -p used to format the indexing, in this case formatting it to index a .vcf file
# input file is the compressed .vcf file

vcffilter -f "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" \
~/ngs_course/dnaseq_pipeline/results/NGS0001.vcf.gz > ~/ngs_course/dnaseq_pipeline/results/NGS0001_filtered.vcf
# vcffilter is used to filter out the bad calls of the .vcf file, using several freebayers information fields generated when creating the .vcf file
# Options used are:
# -f is used to specify the filter which is applied to the .vcf file, removing bad copies
# QUAL > 1 states that low quality calls with a score less than 1 are removed
# AO > 10 states that observations should be at lease 10 log units; SAF argument means the observation is the number of observations on the forward strand, in this case at least 10; the SAR argument means the observation is the number of observations on the reverse strand, in this case more than nothing; RPR are the reads placed right and RPL are the reads placed left, in this case it is set to greater than one for both
# input file is the compressed .vcf file

bedtools intersect -header -wa -a ~/ngs_course/dnaseq_pipeline/results/NGS0001_filtered.vcf -b ~/ngs_course/dnaseq_pipeline/data/annotation.bed \
> ~/ngs_course/dnaseq_pipeline/results/NGS0001_filtered_annotation.vcf
# bedtools intersect is used to filter the filtered .vcf file to only display regions of data observed in the .bed file as well
# Options used:
# -header used to print out the header from the file
# -wa used to write the initial entry from file a for every overlap
# -a determines the file to be compared, in this case the filtered .vcf file
# -b determines the file a is compared to, in this case the .bed file
# output is a .vcf file

bgzip ~/ngs_course/dnaseq_pipeline/results/NGS0001_filtered_annotation.vcf
# compressing the file output from the bedtools intersect operation, preparing it for use further down the pipeline

tabix -p vcf ~/ngs_course/dnaseq_pipeline/results/NGS0001_filtered_annotation.vcf.gz
# tabix is used to index the compressed .vcf file
# Options used:
# -p used to format the indexing, in ths case fomatting it to index a .vcf file

~/annovar/convert2annovar.pl -format vcf4 ~/ngs_course/dnaseq_pipeline/results/NGS0001_filtered_annotated.vcf.gz > ~/ngs_course/dnaseq_pipeline/results/NGS0001_filtered_annotated.avinput
# annovar is an annotating software which annotates every variant in the .vcf file, in this case the filtered .vcf file containing only the reads that are also present in the .bed file
# convert2annovar.pl used to convert the .vcf file to a .pl file, compatible with annovar for annotation
# Options used:
# -format used to indicate that the input file is a .vcf file
# output file is a .avinput format file

~/annovar/table_annovar.pl ~/ngs_course/dnaseq_pipeline/results/NGS0001_filtered_annotated.avinput humandb/ -buildver hg19 \
-out ~/ngs_course/dnaseq_pipeline/results/NGS0001_filtered_annotated -remove \
-protocol refGene,ensGene,clinvar_20180603,exac03,dbnsfp31a_interpro -operation g,g,f,f,f -otherinfo -nastring . -csvout
# table_annovar generates a table containing the annotated variants which were determined to be useful with all the previous filtering
# Options used:
# -buildver used to indicare the build version of the reference genome, in this case hg19
# -out determines the output file
# -protocol designates the databases annovar uses in its annotations. The following databases were used; refGene, ensGene, clinvar_20180603, exac03, dbnsfp31a_interpro

