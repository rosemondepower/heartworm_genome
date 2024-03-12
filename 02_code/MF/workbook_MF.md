# Dirofilaria immitis WGS Lab Book - Microfilaria

In this project, I have performed whole-genome amplification and whole-genome sequencing on individual microfilaria. The goal is to map my fastq sequences and determine the % of heartworm and dog DNA and thus demonstrate how this method is not appropriate for population genetics analysis. I have already performed this on Galaxy, now I just want to repeat it on Artemis.

## Merge fastq files for the same sample

```bash
#!/bin/bash

# PBS directives 
#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N fastq_merge
#PBS -l select=1:ncpus=16:mem=60GB
#PBS -l walltime=12:00:00
#PBS -m abe
#PBS -q defaultQ
#PBS -o fastq_merge.txt
#PBS -M rosemonde.power@sydney.edu.au

cd /scratch/RDS-FSC-Heartworm_MLR-RW/MF/fastq
zcat JS6089_FDSW210306864-1r_H5YLMDSX2_L3_1.fq.gz JS6089_FDSW210306864-1r_HWKN2DSXY_L4_1.fq.gz | gzip > JS6089_merged_1.fq.gz
zcat JS6089_FDSW210306864-1r_H5YLMDSX2_L3_2.fq.gz JS6089_FDSW210306864-1r_HWKN2DSXY_L4_2.fq.gz | gzip > JS6089_merged_2.fq.gz

zcat JS6090_FDSW210306865-1r_H5YLMDSX2_L3_1.fq.gz JS6090_FDSW210306865-1r_HWKN2DSXY_L4_1.fq.gz | gzip > JS6090_merged_1.fq.gz
zcat JS6090_FDSW210306865-1r_H5YLMDSX2_L3_2.fq.gz JS6090_FDSW210306865-1r_HWKN2DSXY_L4_2.fq.gz | gzip > JS6090_merged_2.fq.gz

zcat JS6091_FDSW210306866-1r_H5YLMDSX2_L3_1.fq.gz JS6091_FDSW210306866-1r_HWKN2DSXY_L4_1.fq.gz | gzip > JS6091_merged_1.fq.gz
zcat JS6091_FDSW210306866-1r_H5YLMDSX2_L3_2.fq.gz JS6091_FDSW210306866-1r_HWKN2DSXY_L4_2.fq.gz | gzip > JS6091_merged_2.fq.gz

zcat JS6092_FDSW210306867-1r_H5YLMDSX2_L3_1.fq.gz JS6092_FDSW210306867-1r_HWKN2DSXY_L4_1.fq.gz | gzip > JS6092_merged_1.fq.gz
zcat JS6092_FDSW210306867-1r_H5YLMDSX2_L3_2.fq.gz JS6092_FDSW210306867-1r_HWKN2DSXY_L4_2.fq.gz | gzip > JS6092_merged_2.fq.gz

zcat JS6093_FDSW210306868-1r_H5YLMDSX2_L3_1.fq.gz JS6093_FDSW210306868-1r_HWKN2DSXY_L4_1.fq.gz | gzip > JS6093_merged_1.fq.gz
zcat JS6093_FDSW210306868-1r_H5YLMDSX2_L3_2.fq.gz JS6093_FDSW210306868-1r_HWKN2DSXY_L4_2.fq.gz | gzip > JS6093_merged_2.fq.gz

zcat JS6094_FDSW210306869-1r_H5TWCDSX2_L3_1.fq.gz JS6094_FDSW210306869-1r_H5YLMDSX2_L3_1.fq.gz JS6094_FDSW210306869-1r_HW7W2DSXY_L1_1.fq.gz | gzip > JS6094_merged_1.fq.gz
zcat JS6094_FDSW210306869-1r_H5TWCDSX2_L3_2.fq.gz JS6094_FDSW210306869-1r_H5YLMDSX2_L3_2.fq.gz JS6094_FDSW210306869-1r_HW7W2DSXY_L1_2.fq.gz | gzip > JS6094_merged_2.fq.gz

zcat JS6095_FDSW210306870-1r_H5YLMDSX2_L3_1.fq.gz JS6095_FDSW210306870-1r_HWKN2DSXY_L4_1.fq.gz | gzip > JS6095_merged_1.fq.gz
zcat JS6095_FDSW210306870-1r_H5YLMDSX2_L3_2.fq.gz JS6095_FDSW210306870-1r_HWKN2DSXY_L4_2.fq.gz | gzip > JS6095_merged_2.fq.gz

zcat JS6096_FDSW210306871-1r_H5YLMDSX2_L3_1.fq.gz JS6096_FDSW210306871-1r_HWKN2DSXY_L4_1.fq.gz | gzip > JS6096_merged_1.fq.gz
zcat JS6096_FDSW210306871-1r_H5YLMDSX2_L3_2.fq.gz JS6096_FDSW210306871-1r_HWKN2DSXY_L4_2.fq.gz | gzip > JS6096_merged_2.fq.gz

zcat JS6097_FDSW210306872-1r_H5YLMDSX2_L3_1.fq.gz JS6097_FDSW210306872-1r_HWKN2DSXY_L4_1.fq.gz | gzip > JS6097_merged_1.fq.gz
zcat JS6097_FDSW210306872-1r_H5YLMDSX2_L3_2.fq.gz JS6097_FDSW210306872-1r_HWKN2DSXY_L4_2.fq.gz | gzip > JS6097_merged_2.fq.gz
```

### Check that files merged correctly

I can check that I have the same number of reads before & after merging. I can do this by counting the number of lines. Collect info into Excel sheet.

**Forward reads**

```bash
module load parallel/20160222

# Count Forward reads
cd /scratch/RDS-FSC-Heartworm_MLR-RW/MF/fastq/raw

# Raw data files
for f in *_1.fq.gz; do echo $f;zcat $f|wc -l ; done > /scratch/RDS-FSC-Heartworm_MLR-RW/MF/analysis/count/count_raw_1.txt

cd /scratch/RDS-FSC-Heartworm_MLR-RW/MF/analysis/count

# Every 2nd line moved to new column
sed 'N;s/\n/ /g' count_raw_1.txt | column -t > raw_1.txt

# Get list of unique sample name. Make sample list file.
awk '{print $1}' OFS="\t" raw_1.txt | cut -c1-6 | uniq > samples_1.txt

# Finds all files for each individual sample. Saves to new file for each sample.
parallel --colsep "\t" 'grep {1} raw_1.txt > {1}_raw_1.txt' :::: samples_1.txt
## Academic tradition requires you to cite works you base your article on.
##When using programs that use GNU Parallel to process data for publication
##please cite:
##  O. Tange (2011): GNU Parallel - The Command-Line Power Tool,
##  ;login: The USENIX Magazine, February 2011:42-47.


# prints sample ID and total at the bottom. Extract last line.
for f in JS*_raw_1.txt; do awk '{sum+=$2;print $1" "$2} END {print "'$f'", sum}' $f | tail -1 > total_$f; done

# combine files
cat total_JS*_raw_1.txt > total_raw_1.txt

# remove files I don't need anymore
rm *JS*.txt


# Merged data files
cd /scratch/RDS-FSC-Heartworm_MLR-RW/MF/fastq/merged

# Forward reads
for f in *_1.fq.gz; do echo $f;zcat $f|wc -l ; done > /scratch/RDS-FSC-Heartworm_MLR-RW/MF/analysis/count/count_merged_1.txt

cd /scratch/RDS-FSC-Heartworm_MLR-RW/MF/analysis/count

# Every 2nd line moved to new column
sed 'N;s/\n/ /g' count_merged_1.txt | column -t > merged_1.txt


# Join the raw & merged stats for FORWARD reads
paste total_raw_1.txt merged_1.txt | column -s $'\t' -t > total_both_1.txt

# Make txt file into csv file
mv total_both_1.txt total_both_1.csv
```


**Reverse reads**

```bash
# Count Reverse reads
cd /scratch/RDS-FSC-Heartworm_MLR-RW/MF/fastq/raw

# Raw data files
for f in *_2.fq.gz; do echo $f;zcat $f|wc -l ; done > /scratch/RDS-FSC-Heartworm_MLR-RW/MF/analysis/count/count_raw_2.txt

cd /scratch/RDS-FSC-Heartworm_MLR-RW/MF/analysis/count

# Every 2nd line moved to new column
sed 'N;s/\n/ /g' count_raw_2.txt | column -t > raw_2.txt

# Get list of unique sample name. Make sample list file.
awk '{print $1}' OFS="\t" raw_2.txt | cut -c1-6 | uniq > samples_2.txt

# Finds all files for each individual sample. Saves to new file for each sample.
parallel --colsep "\t" 'grep {1} raw_2.txt > {1}_raw_2.txt' :::: samples_2.txt

# prints sample ID and total at the bottom. Extract last line.
for f in JS*_raw_2.txt; do awk '{sum+=$2;print $1" "$2} END {print "'$f'", sum}' $f | tail -1 > total_$f; done

# combine files
cat total_JS*_raw_2.txt > total_raw_2.txt

# remove files I don't need anymore
rm *JS*.txt


# Merged data files
cd /scratch/RDS-FSC-Heartworm_MLR-RW/MF/fastq/merged

# Forward reads
for f in *_2.fq.gz; do echo $f;zcat $f|wc -l ; done > /scratch/RDS-FSC-Heartworm_MLR-RW/MF/analysis/count/count_merged_2.txt

cd /scratch/RDS-FSC-Heartworm_MLR-RW/MF/analysis/count

# Every 2nd line moved to new column
sed 'N;s/\n/ /g' count_merged_2.txt | column -t > merged_2.txt


# Join the raw & merged stats for FORWARD reads
paste total_raw_2.txt merged_2.txt | column -s $'\t' -t > total_both_2.txt

# Make txt file into csv file
mv total_both_2.txt total_both_2.csv
```
Now we have 2 excel files: 1. Raw vs merged FORWARD reads and 2. Raw vs merged REVERSE reads. Compare the total numbers and ensure that they match up so we didn't lose any data in the merging process.

After inspecting the tables, everything matches up. We have the same number of lines in the raw and merged fastq files. We can continue with the analysis using the merged files.


## FastQC


We want to get some stats on the raw data.

```bash
#!/bin/bash

# PBS directives 
#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N fastqc
#PBS -l select=1:ncpus=1:mem=10GB
#PBS -l walltime=00:30:00
#PBS -m e
#PBS -q defaultQ
#PBS -o fastqc.txt
#PBS -M rosemonde.power@sydney.edu.au
#PBS -J 1-18

# Submit job
## qsub ../fastqc.sh

# Load modules
module load fastqc/0.11.8

# FastQC
cd /scratch/RDS-FSC-Heartworm_MLR-RW/MF/analysis

NCPU=2
OUTDIR=/scratch/RDS-FSC-Heartworm_MLR-RW/MF/analysis/fastqc

config=/scratch/RDS-FSC-Heartworm_MLR-RW/MF/analysis/info.txt
sample=$(awk -v taskID=$PBS_ARRAY_INDEX '$1==taskID {print $2}' $config) 

echo "sample is: $sample"

fastqc -t $NCPU -o ${OUTDIR} /scratch/RDS-FSC-Heartworm_MLR-RW/MF/fastq/merged/${sample}.fq.gz
```


## MultiQC 

```bash
#!/bin/bash

# PBS directives 
#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N multiqc
#PBS -l select=1:ncpus=1:mem=25GB
#PBS -l walltime=03:00:00
#PBS -m e
#PBS -q defaultQ
#PBS -o multiqc.txt
#PBS -M rosemonde.power@sydney.edu.au

# Submit job
## qsub ../multiqc.sh
cd /project/RDS-FSC-Heartworm_MLR-RW/MultiQC

module load git/2.25.0
module load python/3.9.15

# Run MultiQC to combine all of the FastQC reports for the raw fastq files
multiqc /scratch/RDS-FSC-Heartworm_MLR-RW/MF/analysis/fastqc/raw -o /scratch/RDS-FSC-Heartworm_MLR-RW/MF/analysis/fastqc/raw
```


## Trimming

### Trimmomatic

```bash
#!/bin/bash

# PBS directives 
#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N trimmomatic
#PBS -l select=1:ncpus=2:mem=30GB
#PBS -l walltime=04:00:00
#PBS -m e
#PBS -q defaultQ
#PBS -o trimmomatic.txt
#PBS -J 1-9

# qsub ../trimmomatic.pbs

config=/scratch/RDS-FSC-Heartworm_MLR-RW/MF/analysis/trimmomatic/info.txt
sample=$(awk -v taskID=$PBS_ARRAY_INDEX '$1==taskID {print $2}' $config) 
NCPU=2

echo "sample is: $sample"

cd /scratch/RDS-FSC-Heartworm_MLR-RW/MF/analysis/trimmomatic

# Load modules
module load trimmomatic/0.39

# Run Trimmomatic
java -jar /usr/local/trimmomatic/0.39/trimmomatic-0.39.jar PE \
-threads $NCPU \
-phred33 \
/scratch/RDS-FSC-Heartworm_MLR-RW/MF/fastq/merged/${sample}_merged_1.fq.gz \
/scratch/RDS-FSC-Heartworm_MLR-RW/MF/fastq/merged/${sample}_merged_2.fq.gz \
${sample}_1_trimpaired.fq.gz ${sample}_1_trimunpaired.fq.gz \
${sample}_2_trimpaired.fq.gz ${sample}_2_trimunpaired.fq.gz \
SLIDINGWINDOW:10:20 MINLEN:50

# SLIDINGWINDOW:10:20 means it will scan the read with a 10-base wide sliding window, cutting when the average quality per base drops below 20.
```
All new Illumina uses phred33. Only need to worry about phred 64 if you've got pretty old data...


## FastQC & Multi-QC AFTER TRIMMING

Check to see how the data looks after trimming.

```bash
#!/bin/bash

# PBS directives 
#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N fastQC_trimmed
#PBS -l select=1:ncpus=1:mem=20GB
#PBS -l walltime=02:00:00
#PBS -m e
#PBS -q defaultQ
#PBS -o fastQC_trimmed.txt
#PBS -J 1-18

# Submit job
## qsub ../fastqc_trimmed.pbs

# Load modules
module load fastqc/0.11.8

# FastQC
cd /scratch/RDS-FSC-Heartworm_MLR-RW/MF/analysis
NCPU=1
OUTDIR=/scratch/RDS-FSC-Heartworm_MLR-RW/MF/analysis/fastqc/trimmed

config=/scratch/RDS-FSC-Heartworm_MLR-RW/MF/analysis/trimmomatic/info.txt
sample=$(awk -v taskID=$PBS_ARRAY_INDEX '$1==taskID {print $2}' $config) 

fastqc -t $NCPU -o ${OUTDIR} /scratch/RDS-FSC-Heartworm_MLR-RW/MF/analysis/trimmomatic/${sample}_trimpaired.fq.gz
```

```bash
#!/bin/bash

# PBS directives 
#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N multiqc_trimmed
#PBS -l select=1:ncpus=1:mem=10GB
#PBS -l walltime=03:00:00
#PBS -m e
#PBS -q defaultQ
#PBS -o multiqc_trimmed.txt

# Submit job
# qsub ../multiqc_trimmed.sh

# Run MultiQC to combine all of the FastQC reports for the trimmed fastq files

cd /project/RDS-FSC-Heartworm_MLR-RW/MultiQC

# Load modules
module load git/2.25.0
module load python/3.9.15

multiqc /scratch/RDS-FSC-Heartworm_MLR-RW/MF/analysis/fastqc/trimmed -o /scratch/RDS-FSC-Heartworm_MLR-RW/MF/analysis/fastqc/trimmed
```


### Map trimmed reads to combined D. immitis & Wol & dog genome


```bash
#!/bin/bash

# PBS directives 
#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N mapping
#PBS -l select=1:ncpus=16:mem=50GB
#PBS -l walltime=02:00:00
#PBS -m abe
#PBS -q defaultQ
#PBS -o mapping.txt
#PBS -J 1-9

# qsub ../mapping.pbs

config=/scratch/RDS-FSC-Heartworm_MLR-RW/MF/analysis/mapping/info.txt
sample=$(awk -v taskID=$PBS_ARRAY_INDEX '$1==taskID {print $2}' $config) 
NCPU=16

# Set working directory
cd /scratch/RDS-FSC-Heartworm_MLR-RW/MF/analysis/mapping

# Load modules
module load bwa/0.7.17

# map the reads, with a separate mapping job for each sample
bwa mem \
-R '@RG\tID:${SAMPLE_NAME}\tSM:${SAMPLE_NAME}\tPL:illumina' \
/project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_ALL/batch1/analysis/mapping/reference_di_wol_dog.fa \
-t $NCPUS \
../trimmomatic/${sample}_1_trimpaired.fq.gz \
../trimmomatic/${sample}_2_trimpaired.fq.gz \
> /scratch/RDS-FSC-Heartworm_MLR-RW/MF/analysis/mapping/${sample}.tmp.sam

# Set the num_threads param to directly scale with the number of cpus using the PBS environment variable "${NCPUS}).
```

Convert to bam & sort the mapped reads:

```bash
#!/bin/bash

# PBS directives 
#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N mapping_sort
#PBS -l select=1:ncpus=1:mem=100GB
#PBS -l walltime=05:00:00
#PBS -m e
#PBS -q defaultQ
#PBS -o mapping_sort.txt
#PBS -J 1-9

# qsub ../mapping_sort.pbs

config=/scratch/RDS-FSC-Heartworm_MLR-RW/MF/analysis/mapping/info.txt
sample=$(awk -v taskID=$PBS_ARRAY_INDEX '$1==taskID {print $2}' $config) 
NCPU=1

echo "sample is: $sample"

# Set working directory
cd /scratch/RDS-FSC-Heartworm_MLR-RW/MF/analysis/mapping

# Load modules
module load samtools/1.17

# Mapping stats
samtools flagstat ${sample}.tmp.sam > flagstat1/${sample}_flagstat1.txt
	
# convert the sam to bam format
samtools view -q 15 -b -o ${sample}.tmp.bam ${sample}.tmp.sam

# sort the mapped reads in the bam file
samtools sort ${sample}.tmp.bam -o ${sample}.sorted.bam
 
# index the sorted bam
samtools index ${sample}.sorted.bam

# Mapping stats after filtering
samtools flagstat ${sample}.sorted.bam > flagstat2/${sample}_flagstat2.txt
```


Combine flagstat files for all samples so it's easier to read.

```bash
#!/bin/bash

# PBS directives 
#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N multiqc_flagstat1
#PBS -l select=1:ncpus=1:mem=1GB
#PBS -l walltime=00:10:00
#PBS -m e
#PBS -q defaultQ
#PBS -o multiqc_flagstat1.txt

# Submit job
# qsub ../multiqc_flagstat1.sh

cd /project/RDS-FSC-Heartworm_MLR-RW/MultiQC

# Load modules
module load git/2.25.0
module load python/3.9.15

multiqc /scratch/RDS-FSC-Heartworm_MLR-RW/MF/analysis/mapping/flagstat1/*_flagstat1.txt -o /scratch/RDS-FSC-Heartworm_MLR-RW/MF/analysis/mapping/flagstat1
```

```bash
#!/bin/bash

# PBS directives 
#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N multiqc_flagstat2
#PBS -l select=1:ncpus=1:mem=1GB
#PBS -l walltime=00:10:00
#PBS -m e
#PBS -q defaultQ
#PBS -o multiqc_flagstat2.txt

# Submit job
# qsub ../multiqc_flagstat2.sh

cd /project/RDS-FSC-Heartworm_MLR-RW/MultiQC

# Load modules
module load git/2.25.0
module load python/3.9.15

multiqc /scratch/RDS-FSC-Heartworm_MLR-RW/MF/analysis/mapping/flagstat2/*_flagstat2.txt -o /scratch/RDS-FSC-Heartworm_MLR-RW/MF/analysis/mapping/flagstat2
```



## Extract reads that mapped to the *D. immitis* genome

If I mapped to the *D. immitis* and dog genomes separately, there could be reads that mapped to both genomes. To avoid this, I mapped to the combined D. immitis/dog genome. I can now extract the reads that mapped to only the *D. immitis* genome and use this for downstream analyses.

### D. immitis without Wolbachia:

```bash
#!/bin/bash

# PBS directives 
#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N mapping_extract_di
#PBS -l select=1:ncpus=1:mem=20GB
#PBS -l walltime=06:00:00
#PBS -m e
#PBS -q defaultQ
#PBS -o mapping_extract_di.txt
#PBS -J 1-9

# qsub ../mapping_extract_di.pbs

# Set working directory
cd /scratch/RDS-FSC-Heartworm_MLR-RW/MF/analysis/mapping

# Load modules
module load samtools/1.17

config=/scratch/RDS-FSC-Heartworm_MLR-RW/MF/analysis/mapping/info.txt
sample=$(awk -v taskID=$PBS_ARRAY_INDEX '$1==taskID {print $2}' $config) 
NCPU=1

# Extract reads that only mapped to D. immitis.
samtools view -b -h -L /project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_ALL/batch1/analysis/mapping/dimmitis_WSI_2.2_noWb.bed /scratch/RDS-FSC-Heartworm_MLR-RW/MF/analysis/mapping/${sample}.sorted.bam > /scratch/RDS-FSC-Heartworm_MLR-RW/MF/analysis/mapping/${sample}_extract_di.bam
# Should still be in sorted form
# -b flag makes sure the output is bam
# -h flag includes the header in SAM output

samtools view /scratch/RDS-FSC-Heartworm_MLR-RW/MF/analysis/mapping/${sample}_extract_di.bam | head

# I do not have to sort the bam file again, it should still be sorted.

## QC
# How many D. immitis reads were extracted?

samtools flagstat /scratch/RDS-FSC-Heartworm_MLR-RW/MF/analysis/mapping/${sample}_extract_di.bam > /scratch/RDS-FSC-Heartworm_MLR-RW/MF/analysis/mapping/extract_di_flagstat/${sample}_extract_di_flagstat.txt
```


### Wolbachia

```bash
#!/bin/bash

# PBS directives 
#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N mapping_extract_Wb
#PBS -l select=1:ncpus=1:mem=20GB
#PBS -l walltime=06:00:00
#PBS -m e
#PBS -q defaultQ
#PBS -o mapping_extract_Wb.txt
#PBS -J 1-9

# qsub ../mapping_extract_Wb.pbs

# Set working directory
cd /scratch/RDS-FSC-Heartworm_MLR-RW/MF/analysis/mapping

# Load modules
module load samtools/1.17

config=/scratch/RDS-FSC-Heartworm_MLR-RW/MF/analysis/mapping/info.txt
sample=$(awk -v taskID=$PBS_ARRAY_INDEX '$1==taskID {print $2}' $config) 
NCPU=1

# Extract reads that only mapped to D. immitis.
samtools view -b -h -L /project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_ALL/batch1/analysis/mapping/dimmitis_WSI_2.2_Wb.bed /scratch/RDS-FSC-Heartworm_MLR-RW/MF/analysis/mapping/${sample}.sorted.bam > /scratch/RDS-FSC-Heartworm_MLR-RW/MF/analysis/mapping/${sample}_extract_Wb.bam
# Should still be in sorted form
# -b flag makes sure the output is bam
# -h flag includes the header in SAM output

samtools view /scratch/RDS-FSC-Heartworm_MLR-RW/MF/analysis/mapping/${sample}_extract_Wb.bam | head

# I do not have to sort the bam file again, it should still be sorted.

## QC
# How many D. immitis reads were extracted?

samtools flagstat /scratch/RDS-FSC-Heartworm_MLR-RW/MF/analysis/mapping/${sample}_extract_Wb.bam > /scratch/RDS-FSC-Heartworm_MLR-RW/MF/analysis/mapping/extract_Wb_flagstat/${sample}_extract_Wb_flagstat.txt
```


### Dog

I can now extract the reads that mapped to only the dog genome to see how much contamination there is.

Get bed file for dog genome:

```bash
#!/bin/bash

# PBS directives 
#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N mapping_bed_dog
#PBS -l select=1:ncpus=1:mem=10GB
#PBS -l walltime=01:00:30
#PBS -m e
#PBS -q defaultQ
#PBS -o mapping_bed_dog.txt

# qsub ../mapping_bed_dog.sh

# Set working directory
cd /project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_ALL/batch1/analysis/mapping

# Load modules
module load samtools/1.17

# Index the reference file (from Steve's paper) using samtools faidx
samtools faidx /project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_ALL/batch1/reference/GCA_014441545.1_ROS_Cfam_1.0_genomic.fna

# Get the scaffolds/positions.
head /project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_ALL/batch1/reference/GCA_014441545.1_ROS_Cfam_1.0_genomic.fna
# Column 1 is the chromosome/scaffold, column 2 is how long it is, then there's some other info.

# Get chromosome, then start and end positions
awk '{print $1, "1", $2}' OFS="\t" /project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_ALL/batch1/reference/GCA_014441545.1_ROS_Cfam_1.0_genomic.fna.fai | head

# Save this info as a bed file
awk '{print $1, "1", $2}' OFS="\t" /project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_ALL/batch1/reference/GCA_014441545.1_ROS_Cfam_1.0_genomic.fna.fai > /project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_ALL/batch1/reference/GCA_014441545.1_ROS_Cfam_1.0_genomic.bed
# Now we have a nice bed file that has info telling us where things are
```

Now extract dog reads:

```bash
#!/bin/bash

# PBS directives 
#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N mapping_extract_dog
#PBS -l select=1:ncpus=1:mem=20GB
#PBS -l walltime=06:00:00
#PBS -m e
#PBS -q defaultQ
#PBS -o mapping_extract_dog.txt
#PBS -J 1-9

# qsub ../mapping_extract_dog.pbs

# Set working directory
cd /scratch/RDS-FSC-Heartworm_MLR-RW/MF/analysis/mapping

# Load modules
module load samtools/1.17

config=/scratch/RDS-FSC-Heartworm_MLR-RW/MF/analysis/mapping/info.txt
sample=$(awk -v taskID=$PBS_ARRAY_INDEX '$1==taskID {print $2}' $config) 
NCPU=1

# Extract reads that only mapped to D. immitis.
samtools view -b -h -L /project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_ALL/batch1/analysis/mapping/GCA_014441545.1_ROS_Cfam_1.0_genomic.bed /scratch/RDS-FSC-Heartworm_MLR-RW/MF/analysis/mapping/${sample}.sorted.bam > /scratch/RDS-FSC-Heartworm_MLR-RW/MF/analysis/mapping/${sample}_extract_dog.bam
# Should still be in sorted form
# -b flag makes sure the output is bam
# -h flag includes the header in SAM output

samtools view /scratch/RDS-FSC-Heartworm_MLR-RW/MF/analysis/mapping/${sample}_extract_dog.bam | head

# I do not have to sort the bam file again, it should still be sorted.

## QC
# How many D. immitis reads were extracted?

samtools flagstat /scratch/RDS-FSC-Heartworm_MLR-RW/MF/analysis/mapping/${sample}_extract_dog.bam > /scratch/RDS-FSC-Heartworm_MLR-RW/MF/analysis/mapping/extract_dog_flagstat/${sample}_extract_dog_flagstat.txt
```

Combine flagstat files for all samples so it's easier to read.

```bash
#!/bin/bash

# PBS directives 
#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N multiqc_extract_flagstat
#PBS -l select=1:ncpus=1:mem=1GB
#PBS -l walltime=00:30:00
#PBS -m e
#PBS -q defaultQ
#PBS -o multiqc_extract_flagstat.txt

# Submit job
# qsub ../multiqc_extract_flagstat.sh

cd /project/RDS-FSC-Heartworm_MLR-RW/MultiQC

# Load modules
module load git/2.25.0
module load python/3.9.15

multiqc /scratch/RDS-FSC-Heartworm_MLR-RW/MF/analysis/mapping/extract_di_flagstat/*_extract_di_flagstat.txt -o /scratch/RDS-FSC-Heartworm_MLR-RW/MF/analysis/mapping/extract_di_flagstat

multiqc /scratch/RDS-FSC-Heartworm_MLR-RW/MF/analysis/mapping/extract_Wb_flagstat/*_extract_Wb_flagstat.txt -o /scratch/RDS-FSC-Heartworm_MLR-RW/MF/analysis/mapping/extract_Wb_flagstat

multiqc /scratch/RDS-FSC-Heartworm_MLR-RW/MF/analysis/mapping/extract_dog_flagstat/*_extract_dog_flagstat.txt -o /scratch/RDS-FSC-Heartworm_MLR-RW/MF/analysis/mapping/extract_dog_flagstat
```