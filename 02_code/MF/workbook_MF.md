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

# Extract reads that only mapped to Wolbachia
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


## Coverage

Adopted the code from Javier's paper - acknowledge this in methods.

```bash
#!/bin/bash

# PBS directives 
#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N coverage
#PBS -l select=3:ncpus=1:mem=50GB
#PBS -l walltime=24:00:00
#PBS -m abe
#PBS -q defaultQ
#PBS -o coverage.txt

# qsub ../coverage.sh

WORKING_DIR=/scratch/RDS-FSC-Heartworm_MLR-RW/MF/analysis/mapping
cd ${WORKING_DIR}

WINDOW='100000'

module load bamtools/2.5.1
module load bedtools/2.31.0
module load samtools/1.17

for i in /scratch/RDS-FSC-Heartworm_MLR-RW/MF/analysis/mapping/*sorted.bam; do

bamtools header -in ${i} | grep "^@SQ" | awk -F'[:\t]' '{printf $3"\t"1"\t"$5"\n"}' OFS="\t" > ${i%.bam}.chr.bed
bamtools header -in ${i} | grep "^@SQ" | awk -F'[:\t]' '{printf $3"\t"$5"\n"}' OFS="\t" > ${i%.bam}.chr.genome

bedtools makewindows -g ${i%.bam}.chr.genome -w ${WINDOW} > ${i%.bam}.${WINDOW}_window.bed

samtools bedcov -Q 20 ${i%.bam}.chr.bed ${i} | awk -F'\t' '{printf $1"\t"$2"\t"$3"\t"$4"\t"$4/($3-$2)"\n"}' OFS="\t" > ${i%.bam}.chr.cov
samtools bedcov -Q 20 ${i%.bam}.${WINDOW}_window.bed ${i} | awk -F'\t' '{printf $1"\t"$2"\t"$3"\t"$4"\t"$4/($3-$2)"\n"}' OFS="\t" > ${i%.bam}.${WINDOW}_window.cov

rm ${i%.bam}.chr.bed ${i%.bam}.${WINDOW}_window.bed ${i%.bam}.chr.genome;

done

for i in *.chr.cov; do 

printf "${i}\n" > ${i}.tmp | awk '{print $5}' OFS="\t" ${i} >> ${i}.tmp;

done

paste *.tmp > coverage_stats.summary
```
```bash
# moved all relevant coverage files into coverage folder
cd /scratch/RDS-FSC-Heartworm_MLR-RW/MF/analysis/mapping/coverage
# remove temporary files
rm *.tmp
```


### Generate quantitative stats on coverage for supplementary tables etc
Extract mtDNA, Wb and nuclear (mean & stddev) data

For nuclear, we will select only the defined Chr (chrX and chr1 to chr4)

```bash
#!/bin/bash

# PBS directives 
#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N coverage_stats
#PBS -l select=1:ncpus=1:mem=4GB
#PBS -l walltime=00:10:00
#PBS -m e
#PBS -q defaultQ
#PBS -o coverage_stats.txt

# qsub ../coverage_stats.sh

cd /scratch/RDS-FSC-Heartworm_MLR-RW/MF/analysis/mapping/coverage
# Load modules
module load datamash/1.7

# extract mtDNA and nuclear (mean & stddev) data
for i in *sorted.chr.cov; do
	name=${i%.sorted.chr.cov};
	nuc=$(grep -v "scaffold\|Wb\|Mt\|CM025\|JAA" ${i%.sorted.chr.cov}.sorted.100000_window.cov | datamash mean 5 sstdev 5 );
	mtDNA=$(grep "chrMtDNA" ${i} | cut -f5 );
	Wb=$(grep 'chrWb' ${i} | cut -f5 ); 
	echo -e "${name}\t${nuc}\t${mtDNA}\t${Wb}";
done > 'mito_wolb_cov.stats'
```
Transferred all the relevant files into the R_analysis folder on my computer for further analysis in R. Now we'll generate some plots and stats.

Coverage in R:

```R
# HW WGS Coverage

# load libraries
library(tidyverse)
library(ggpubr)
library(ggsci)
library(stringr)

setwd("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/HW_WGS_1/Artemis/coverage")

#first, I have to read the nuclear cov stat and estimate the mean and sd
#then, to add it to 'mito_wolb_cov.stats

nuc_mito_wb_cov <- read.table('mito_wolb_cov.stats', header = F) %>% as_tibble()

colnames(nuc_mito_wb_cov) <- c('ID', 'nuc_cov', 'sd_nuc_cov', 'mito_cov', 'wb_cov')

write_csv(nuc_mito_wb_cov, 'nuc_mit_wb_cov.csv')

# nuclear, mitochondrial and Wb DNA coverage ratio

n_m <- ggplot(nuc_mito_wb_cov, aes(x=ID, y=mito_cov/nuc_cov)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(title = "Nuc. to mito. genome coverage ratio", y = "Coverage Ratio")
n_m

n_wb <- ggplot(nuc_mito_wb_cov, aes(x=ID, y=wb_cov/nuc_cov)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(title = "Nuc. to Wolb. genome coverage ratio", y = "Coverage Ratio")
n_wb

m_wb <- ggplot(nuc_mito_wb_cov, aes(x=ID, y=mito_cov/wb_cov)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(title = "Mito. to Wolb. genome coverage ratio", y = "Coverage Ratio")
m_wb

ggarrange(n_m, n_wb, m_wb, ncol = 3)
ggsave("cov_ratios.png", height=6, width=20)


# list file names
file_names.window <- list.files(path = "C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/HW_WGS_1/Artemis/coverage",pattern = ".sorted.100000_window.cov")

# load data using file names, and make a formatted data frame
setwd("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/HW_WGS_1/Artemis/coverage")

data <- purrr::map_df(file_names.window, function(x) {
  data <- read.delim(x, header = F, sep="\t")
  data$V1 <- str_replace(data$V1, 'dirofilaria_immitis_', '')
  data <- tibble::rowid_to_column(data, "NUM")
  cbind(sample_name = gsub(".sorted.100000_window.cov","",x), data)
})
colnames(data) <- c("ID", "NUM", "CHR", "START", "END", 
                    "RAW_COVERAGE", "PROPORTION_COVERAGE")



# D. immitis coverage

# remove scaffolds, mitochondrial and wolbachia genome
data_nuc <- dplyr::filter(data, !grepl("scaffold|MtDNA|Wb|JAAUVH|CM025",CHR))
# Also remove chromosomes called JAAUVH010000344 and CM025130.1 etc. - they are part of the dog genome. Removing them here is totally fine.

# data$SEX <- str_detect(data$SCAF,"Trichuris_trichiura_1_")


# plot the general cov for each sample
ggplot(data_nuc, aes(NUM, PROPORTION_COVERAGE/(median(PROPORTION_COVERAGE)), group = ID, col = CHR)) +
  geom_point(size=0.5) +
  labs( x = "Genome position" , y = "Relative coverage per 100kb window") +
  theme_bw() + theme(strip.text.x = element_text(size = 6)) +
  facet_wrap(~ID, scales = "free_y")

ggsave("ALL_genomewide_prop_median_coverage_allsamples.png", height=11.25, width=15)


# this shows the proportion of coverage relative to the median coverage of all samples. Let's also just look at the proportion coverage by itself (not in relation to anything else).
ggplot(data_nuc, aes(NUM, PROPORTION_COVERAGE, group = ID, col = CHR)) +
  geom_point(size=0.5) +
  labs( x = "Genome position" , y = "Raw coverage per 100kb window") +
  theme_bw() + theme(strip.text.x = element_text(size = 6)) +
  facet_wrap(~ID, scales = "free_y")


# include line for the mean
mean_coverage_per_sample <- aggregate(PROPORTION_COVERAGE ~ ID, data = data_nuc, FUN = mean)

data_nuc <- merge(data_nuc, mean_coverage_per_sample, by = "ID", suffixes = c("", "_mean"))

ggplot(data_nuc, aes(NUM, PROPORTION_COVERAGE, group = ID, col = CHR)) +
  geom_point(size=0.5) +
  geom_hline(aes(yintercept = PROPORTION_COVERAGE_mean), linetype = "dashed", color = "black") +
  labs( x = "Genome position" , y = "Raw coverage per 100kb window") +
  theme_bw() + theme(strip.text.x = element_text(size = 6)) +
  facet_wrap(~ID, scales = "free_y")

ggsave("ALL_genomewide_prop_coverage_allsamples.png", height=11.25, width=15)





# Let's see only the chrX to explore the sex of the sample
#Plotting with the chr1 helps to see differences
data_nuc %>%
  filter(., grepl("chrX|chr1",CHR)) %>%
  ggplot(aes(NUM, PROPORTION_COVERAGE/(median(PROPORTION_COVERAGE)), group = ID, col = CHR)) +
  geom_point(size=0.2) +
  labs( x = "Genome position" , y = "Relative coverage per 100kb window") +
  theme_bw() + theme(strip.text.x = element_text(size = 6)) +
  facet_wrap(~ID, scales = "free_y")

ggsave("chrXtochr1_genomewide_coverage_allsamples.png", height=11.25, width=15)



###### Wolbachia coverage

# remove D. immitis and dog genome
data_wb <- dplyr::filter(data, !grepl("scaffold|MtDNA|chrX|chr1|chr2|chr3|chr4|JAAUVH|CM025",CHR))

# plot cov for each sample
ggplot(data_wb, aes(NUM, PROPORTION_COVERAGE/(median(PROPORTION_COVERAGE)), group = ID, col = CHR)) +
  geom_point(size=0.5) +
  labs( x = "Genome position" , y = "Relative coverage per 100kb window") +
  theme_bw() + theme(strip.text.x = element_text(size = 6)) +
  facet_wrap(~ID, scales = "free_y")

ggsave("wb_prop_median_coverage_allsamples.png", height=11.25, width=15)

# plot cov for proportionsnwithout relation to the median
ggplot(data_wb, aes(NUM, PROPORTION_COVERAGE), group = ID, col = CHR) +
  geom_point(size=0.5) +
  labs( x = "Genome position" , y = "Relative coverage per 100kb window") +
  theme_bw() + theme(strip.text.x = element_text(size = 6)) +
  facet_wrap(~ID, scales = "free_y")

ggsave("wb_prop_coverage_allsamples.png", height=11.25, width=15)


###### mtDNA

# remove D. immitis and dog genome
data_mtDNA <- dplyr::filter(data, !grepl("scaffold|Wb|chrX|chr1|chr2|chr3|chr4|JAAUVH|CM025",CHR))

# plot cov for each sample
ggplot(data_mtDNA, aes(NUM, PROPORTION_COVERAGE/(median(PROPORTION_COVERAGE)), group = ID, col = CHR)) +
  geom_point(size=0.5) +
  labs( x = "Genome position" , y = "Relative coverage per 100kb window") +
  theme_bw() + theme(strip.text.x = element_text(size = 6)) +
  facet_wrap(~ID, scales = "free_y")

ggsave("mtDNA_prop_median_coverage_allsamples.png", height=11.25, width=15)

# plot cov for proportionsnwithout relation to the median
ggplot(data_mtDNA, aes(NUM, PROPORTION_COVERAGE), group = ID, col = CHR) +
  geom_point(size=0.5) +
  labs( x = "Genome position" , y = "Relative coverage per 100kb window") +
  theme_bw() + theme(strip.text.x = element_text(size = 6)) +
  facet_wrap(~ID, scales = "free_y")

ggsave("mtDNA_prop_coverage_allsamples.png", height=11.25, width=15)
## the window sizes are just too large to plot the mitochondrial and Wolbachia coverage. 
```
