#!/bin/bash

#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N HW_WGS3_1GB_workbook
#PBS -l walltime=06:00:00
#PBS -l select=1:ncpus=24:mem=20GB 
#PBS -q defaultQ
#PBS -o WGS3_1GB_workbook.txt

# Submit job
## qsub -P RDS-FSC-Heartworm_MLR-RW WGS3_1GB_workbook.sh

### Reminder: Make sure you have your seqs and ref seqs ready in the folder


# FastQC
cd /project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_3/1GB/fastqc

module load fastqc

INPUTDIR="/project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_3/1GB/seq"
NCPU=2
OUTDIR="/project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_3/1GB/fastqc"

fastqc -t $NCPU -o $OUTDIR $INPUTDIR/*.fq

## This worked. Downloaded webpage outputs onto my Z drive and viewed this in normal web browser. Compared these outputs to the Galaxy output (of the same sample) and they looked identical.


# Trimmomatic
cd /project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_3/1GB/trim

module load trimmomatic

module load parallel

parallel --colsep "\t" 'java -jar /usr/local/trimmomatic/0.38/trimmomatic-0.38.jar PE /project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_3/1GB/seq/{1}.fq /project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_3/1GB/seq/{2}.fq /project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_3/1GB/trim/{1}_trimpaired.fq /project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_3/1GB/trim/{1}_trimunpaired.fq /project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_3/1GB/trim/{2}_trimpaired.fq /project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_3/1GB/trim/{2}_trimunpaired.fq AVGQUAL:30 MINLEN:150' :::: /project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_3/1GB/samples.txt


# BWA Index
## The first step of using BWA is to make an index of the reference genome in fasta format.
cd /project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_3/1GB/map
module load bwa
bwa index -p index_nDi.2.2 -a bwtsw /project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_3/1GB/refseq/nDi.2.2.fna


# Map with BWA-MEM to D. immitis reference genome
parallel --colsep "\t" 'bwa mem index_nDi.2.2 /project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_3/1GB/trim/{1}_trimpaired.fq /project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_3/1GB/trim/{2}_trimpaired.fq -t 2 -v 1 > {3}_bwa_mem_alignments.bam' :::: /project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_3/1GB/samples.txt


# Sortsam
cd /project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_3/1GB/sort
module load picard/2.18.23
parallel --colsep "\t" 'java -jar /usr/local/picard/2.18.23/picard.jar SortSam \
      I=/project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_3/1GB/map/{3}_bwa_mem_alignments.bam \
      O={3}_sorted.bam \
      SORT_ORDER=coordinate QUIET=true VERBOSITY=ERROR VALIDATION_STRINGENCY=LENIENT' :::: /project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_3/1GB/samples.txt


# Samtools flagstat
cd /project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_3/1GB/flagstat
module load samtools/1.9
parallel --colsep "\t" 'samtools flagstat /project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_3/1GB/sort/{3}_sorted.bam' :::: /project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_3/1GB/samples.txt

# The above code runs it parallel, but it didn't label which stats results belonged to which sample. So, run each sample individually (below).


# JS6342_DSFP220004094-1a_HLGLMDSX3_L1
samtools flagstat /project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_3/1GB/sort/JS6342_DSFP220004094-1a_HLGLMDSX3_L1_sorted.bam

# JS6342_DSFP220004094-1a_HLGNNDSX3_L1
samtools flagstat /project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_3/1GB/sort/JS6342_DSFP220004094-1a_HLGNNDSX3_L1_sorted.bam

# JS6343_DSFP220004095-1a_HJJFMDSX3_L2
samtools flagstat /project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_3/1GB/sort/JS6343_DSFP220004095-1a_HJJFMDSX3_L2_sorted.bam

# JS6344_DSFP220004096-1a_HJJFMDSX3_L2
samtools flagstat /project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_3/1GB/sort/JS6344_DSFP220004096-1a_HJJFMDSX3_L2_sorted.bam

# JS6345_DSFP220004097-1a_HJJFMDSX3_L2
samtools flagstat /project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_3/1GB/sort/JS6345_DSFP220004097-1a_HJJFMDSX3_L2_sorted.bam

# JS6346_DSFP220004098-1a_HJJFMDSX3_L3
samtools flagstat /project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_3/1GB/sort/JS6346_DSFP220004098-1a_HJJFMDSX3_L3_sorted.bam

# JS6347_DSFP220004099-1a_HJJFMDSX3_L2
samtools flagstat /project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_3/1GB/sort/JS6347_DSFP220004099-1a_HJJFMDSX3_L2_sorted.bam

# JS6348_DSFP220004100-1a_HLGNNDSX3_L2
samtools flagstat /project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_3/1GB/sort/JS6348_DSFP220004100-1a_HLGNNDSX3_L2_sorted.bam

# JS6349_DSFP220004101-1a_HJJFMDSX3_L2
samtools flagstat /project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_3/1GB/sort/JS6349_DSFP220004101-1a_HJJFMDSX3_L2_sorted.bam

# JS6350_DSFP220004102-1a_HJJFMDSX3_L2
samtools flagstat /project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_3/1GB/sort/JS6350_DSFP220004102-1a_HJJFMDSX3_L2_sorted.bam

# JS6351_DSFP220004103-1a_HJJFMDSX3_L2
samtools flagstat /project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_3/1GB/sort/JS6351_DSFP220004103-1a_HJJFMDSX3_L2_sorted.bam

# JS6352_DSFP220004104-1a_HJJFMDSX3_L2
samtools flagstat /project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_3/1GB/sort/JS6352_DSFP220004104-1a_HJJFMDSX3_L2_sorted.bam

# JS6353_DSFP220004105-1a_HJJFMDSX3_L2
samtools flagstat /project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_3/1GB/sort/JS6353_DSFP220004105-1a_HJJFMDSX3_L2_sorted.bam

# JS6354_DSFP220004106-1a_HJJFMDSX3_L2
samtools flagstat /project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_3/1GB/sort/JS6354_DSFP220004106-1a_HJJFMDSX3_L2_sorted.bam

# JS6355_DSFP220004107-1a_HJJFMDSX3_L2
samtools flagstat /project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_3/1GB/sort/JS6355_DSFP220004107-1a_HJJFMDSX3_L2_sorted.bam

# JS6356_DSFP220004108-1a_HJJFMDSX3_L2
samtools flagstat /project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_3/1GB/sort/JS6356_DSFP220004108-1a_HJJFMDSX3_L2_sorted.bam

# JS6357_DSFP220004109-1a_HJJFMDSX3_L2
samtools flagstat /project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_3/1GB/sort/JS6357_DSFP220004109-1a_HJJFMDSX3_L2_sorted.bam

# JS6358_DSFP220004110-1a_HJJFMDSX3_L2
samtools flagstat /project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_3/1GB/sort/JS6358_DSFP220004110-1a_HJJFMDSX3_L2_sorted.bam

# JS6359_DSFP220004111-1a_HLGNNDSX3_L1
samtools flagstat /project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_3/1GB/sort/JS6359_DSFP220004111-1a_HLGNNDSX3_L1_sorted.bam

# JS6360_DSFP220004112-1a_HLGNNDSX3_L1
samtools flagstat /project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_3/1GB/sort/JS6360_DSFP220004112-1a_HLGNNDSX3_L1_sorted.bam

# Normalize fasta
cd /project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_3/1GB/refseq
java -jar /usr/local/picard/2.18.23/picard.jar NormalizeFasta \
      I=nDi.2.2.fna \
      O=normalized_nDi.2.2.fasta

# Make index
samtools faidx normalized_nDi.2.2.fasta

# Generate pileup
cd /project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_3/1GB/pileup

parallel --colsep "\t" 'samtools mpileup -f /project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_3/1GB/refseq/normalized_nDi.2.2.fasta /project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_3/1GB/sort/{3}_sorted.bam > {3}_pileup.bcf' :::: /project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_3/1GB/samples.txt






