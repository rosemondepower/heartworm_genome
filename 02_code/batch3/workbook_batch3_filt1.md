# Dirofilaria immitis WGS Lab Book - Batch3

### Rose Power USYD 2023

## Public data

**D. repens data**

There were 9 SRR runs for the same sample bioproject(SAMN11159215). I can probably merge these together if it's all the same sample?

```bash
#!/bin/bash

# PBS directives 
#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N D.repens
#PBS -l select=1:ncpus=4:mem=30GB
#PBS -l walltime=10:00:00
#PBS -m abe
#PBS -q defaultQ
#PBS -o D.repens.txt
#PBS -J 1-9

config=/scratch/RDS-FSC-Heartworm_MLR-RW/batch3/fastq/D.repens/info.txt
sample=$(awk -v taskID=$PBS_ARRAY_INDEX '$1==taskID {print $1}' $config) 

cd /scratch/RDS-FSC-Heartworm_MLR-RW/batch3/fastq/D.repens/

module load sratoolkit/3.0.3

fastq-dump --gzip ${sample}
```

```bash
zcat SRR*.fastq.gz | gzip > D.repens_merged.fastq.gz

# Check that the number of lines matches
zcat $f | wc -l
```
A bunch of files downloaded for D. repens. Merge them all together?

```bash
#!/bin/bash

# PBS directives 
#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N D.repens_merge
#PBS -l select=1:ncpus=4:mem=30GB
#PBS -l walltime=05:00:00
#PBS -m abe
#PBS -q defaultQ
#PBS -o D.repens_merge.txt
#PBS -M rosemonde.power@sydney.edu.au

cd /scratch/RDS-FSC-Heartworm_MLR-RW/batch3/fastq/D.repens
zcat SRR* | gzip > Drepens_merged.fastq.gz
```

Check if it has the same number of lines.

```bash
module load parallel/20160222

# Count Forward reads
cd /scratch/RDS-FSC-Heartworm_MLR-RW/batch3/fastq/D.repens

# Raw data files
for f in SRR*.fastq.gz; do echo $f;zcat $f|wc -l ; done > /scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/count/count_raw_Drepens.txt

cd /scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/count

# Every 2nd line moved to new column
sed 'N;s/\n/ /g' count_raw_Drepens.txt | column -t > raw_Drepens.txt

# Get list of unique sample name. Make sample list file.
awk '{print $1}' OFS="\t" raw_Drepens.txt | cut -c1-10 | uniq > samples_Drepens.txt

# Finds all files for each individual sample. Saves to new file for each sample.
parallel --colsep "\t" 'grep {1} raw_Drepens.txt > {1}_raw_Drepens.txt' :::: samples_Drepens.txt

# prints sample ID and total at the bottom. Extract last line.
for f in SRR*_raw_Drepens.txt; do awk '{sum+=$2;print $1" "$2} END {print "'$f'", sum}' $f | tail -1 > total_$f; done

# combine files
cat total_SRR*_raw_Drepens.txt > total_raw_Drepens.txt

# remove files I don't need anymore
rm SRR*.txt
rm *SRR*.txt


# Merged data files
cd /scratch/RDS-FSC-Heartworm_MLR-RW/batch3/fastq/D.repens

# Forward reads
zcat Drepens_merged.fastq.gz |wc -l > /scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/count/count_merged_Drepens.txt
# And it's a match
```


## Merging files for the same sample

Running it as a job just kept going on forever...? Do this on the terminal instead.
Even on the terminal it kept going forever. Maybe specify the exact files to cat so it's not on an infinite loop?

```bash
# Ran on the terminal
cd /scratch/RDS-FSC-Heartworm_MLR-RW/batch3/fastq/novogene/JS6648
zcat JS6648_DKDN230043077-1A_H7MFFDSX7_L1_1.fq.gz JS6648_DKDN230043077-1A_HC2NMDSX7_L2_1.fq.gz | gzip > JS6648_merged_1.fq.gz #done
zcat JS6648_DKDN230043077-1A_H7MFFDSX7_L1_2.fq.gz JS6648_DKDN230043077-1A_HC2NMDSX7_L2_2.fq.gz | gzip > JS6648_merged_2.fq.gz #done


#!/bin/bash

# PBS directives 
#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N JS6656_merge
#PBS -l select=1:ncpus=4:mem=60GB
#PBS -l walltime=02:00:00
#PBS -m abe
#PBS -q defaultQ
#PBS -o JS6656_merge.txt
#PBS -M rosemonde.power@sydney.edu.au

cd /scratch/RDS-FSC-Heartworm_MLR-RW/batch3/fastq/novogene/JS6656
zcat JS6656_DKDN230043085-1A_H7MFFDSX7_L1_1.fq.gz JS6656_DKDN230043085-1A_HC2NMDSX7_L2_1.fq.gz | gzip > JS6656_merged_1.fq.gz
zcat JS6656_DKDN230043085-1A_H7MFFDSX7_L1_2.fq.gz JS6656_DKDN230043085-1A_HC2NMDSX7_L2_2.fq.gz | gzip > JS6656_merged_2.fq.gz
#done

#!/bin/bash

# PBS directives 
#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N JS6657_merge
#PBS -l select=1:ncpus=4:mem=60GB
#PBS -l walltime=02:00:00
#PBS -m abe
#PBS -q defaultQ
#PBS -o JS6657_merge.txt
#PBS -M rosemonde.power@sydney.edu.au

cd /scratch/RDS-FSC-Heartworm_MLR-RW/batch3/fastq/novogene/JS6657
zcat JS6657_DKDN230043086-1A_H7MFFDSX7_L1_1.fq.gz JS6657_DKDN230043086-1A_HC2NMDSX7_L2_1.fq.gz | gzip > JS6657_merged_1.fq.gz
zcat JS6657_DKDN230043086-1A_H7MFFDSX7_L1_2.fq.gz JS6657_DKDN230043086-1A_HC2NMDSX7_L2_2.fq.gz | gzip > JS6657_merged_2.fq.gz
#done


#!/bin/bash

# PBS directives 
#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N JS6659_merge
#PBS -l select=1:ncpus=4:mem=60GB
#PBS -l walltime=02:00:00
#PBS -m abe
#PBS -q defaultQ
#PBS -o JS6659_merge.txt
#PBS -M rosemonde.power@sydney.edu.au

cd /scratch/RDS-FSC-Heartworm_MLR-RW/batch3/fastq/novogene/JS6659
zcat JS6659_DKDN230043088-1A_H7MFFDSX7_L2_1.fq.gz JS6659_DKDN230043088-1A_H7WVCDSX7_L1_1.fq.gz | gzip > JS6659_merged_1.fq.gz
zcat JS6659_DKDN230043088-1A_H7MFFDSX7_L2_2.fq.gz JS6659_DKDN230043088-1A_H7WVCDSX7_L1_2.fq.gz | gzip > JS6659_merged_2.fq.gz
#done
```


### Check that files merged correctly

I can check that I have the same number of reads before & after merging. I can do this by counting the number of lines. Collect info into Excel sheet.

**Forward reads**

```bash
module load parallel/20160222

# Count Forward reads
cd /scratch/RDS-FSC-Heartworm_MLR-RW/batch3/fastq/raw

# Raw data files
for f in *_1.fq.gz; do echo $f;zcat $f|wc -l ; done > /scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/count/count_raw_1.txt

cd /scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/count

# Every 2nd line moved to new column
sed 'N;s/\n/ /g' count_raw_1.txt | column -t > raw_1.txt

# Get list of unique sample name. Make sample list file.
awk '{print $1}' OFS="\t" raw_1.txt | cut -c1-6 | uniq > samples_1.txt

# Finds all files for each individual sample. Saves to new file for each sample.
parallel --colsep "\t" 'grep {1} raw_1.txt > {1}_raw_1.txt' :::: samples_1.txt

# prints sample ID and total at the bottom. Extract last line.
for f in JS*_raw_1.txt; do awk '{sum+=$2;print $1" "$2} END {print "'$f'", sum}' $f | tail -1 > total_$f; done

# combine files
cat total_JS*_raw_1.txt > total_raw_1.txt

# remove files I don't need anymore
rm JS*.txt
rm *JS*.txt


# Merged data files
cd /scratch/RDS-FSC-Heartworm_MLR-RW/batch3/fastq/merged

# Forward reads
for f in *_1.fq.gz; do echo $f;zcat $f|wc -l ; done > /scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/count/count_merged_1.txt

cd /scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/count

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
cd /scratch/RDS-FSC-Heartworm_MLR-RW/batch3/fastq/raw

# Raw data files
for f in *_2.fq.gz; do echo $f;zcat $f|wc -l ; done > /scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/count/count_raw_2.txt

cd /scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/count

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
rm JS*.txt
rm *JS*.txt


# Merged data files
cd /scratch/RDS-FSC-Heartworm_MLR-RW/batch3/fastq/merged

# Forward reads
for f in *_2.fq.gz; do echo $f;zcat $f|wc -l ; done > /scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/count/count_merged_2.txt

cd /scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/count

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
#PBS -J 1-38

# Submit job
## qsub ../fastqc.pbs

# Load modules
module load fastqc/0.11.8

# FastQC
cd /scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis

NCPU=2
OUTDIR=/scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/fastqc

config=/scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/info.txt
sample=$(awk -v taskID=$PBS_ARRAY_INDEX '$1==taskID {print $2}' $config) 

echo "sample is: $sample"

fastqc -t $NCPU -o ${OUTDIR} /scratch/RDS-FSC-Heartworm_MLR-RW/batch3/fastq/${sample}.fq.gz
```


## MultiQC - need to do this still

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
## qsub ../multiqc.pbs
cd /project/RDS-FSC-Heartworm_MLR-RW/MultiQC

module load git
module load python/3.9.15

# Run MultiQC to combine all of the FastQC reports for the raw fastq files
multiqc /scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/fastqc/raw -o /scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/fastqc/raw
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
#PBS -J 1-19

# qsub ../trimmomatic.pbs

config=/scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/trimmomatic/info.txt
sample=$(awk -v taskID=$PBS_ARRAY_INDEX '$1==taskID {print $2}' $config) 
NCPU=2

echo "sample is: $sample"

cd /scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/trimmomatic

# Load modules
module load trimmomatic/0.38

# Run Trimmomatic
java -jar /usr/local/trimmomatic/0.38/trimmomatic-0.38.jar PE \
-threads $NCPU \
-phred33 \
/scratch/RDS-FSC-Heartworm_MLR-RW/batch3/fastq/${sample}_1.fq.gz \
/scratch/RDS-FSC-Heartworm_MLR-RW/batch3/fastq/${sample}_2.fq.gz \
${sample}_1_trimpaired.fq.gz ${sample}_1_trimunpaired.fq.gz \
${sample}_2_trimpaired.fq.gz ${sample}_2_trimunpaired.fq.gz \
SLIDINGWINDOW:10:20 MINLEN:50

# SLIDINGWINDOW:10:20 means it will scan the read with a 10-base wide sliding window, cutting when the average quality per base drops below 20.

# Instead of SLIDINGWINDOW, in my previous practice code I used 'AVGQUAL:30 MINLEN:150'.
```
Assume all files are phred33 quality encoded? The SRR files have '????' in the quality scores so I'll have to check this somehow..
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
#PBS -J 1-38

# Submit job
## qsub ../fastqc_trimmed.pbs

# Load modules
module load fastqc/0.11.8

# FastQC
cd /scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis
NCPU=1
OUTDIR=/scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/fastqc/trimmed

config=/scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/info.txt
sample=$(awk -v taskID=$PBS_ARRAY_INDEX '$1==taskID {print $2}' $config) 

fastqc -t $NCPU -o ${OUTDIR} /scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/trimmomatic/${sample}_trimpaired.fq.gz
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
# qsub ../multiqc_trimmed.pbs

# Run MultiQC to combine all of the FastQC reports for the trimmed fastq files

cd /project/RDS-FSC-Heartworm_MLR-RW/MultiQC

# Load modules
module load git/2.25.0
module load python/3.9.15

multiqc /scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/fastqc/trimmed -o /scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/fastqc/trimmed
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
#PBS -J 1-19

# qsub ../mapping.pbs

config=/scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping/info.txt
sample=$(awk -v taskID=$PBS_ARRAY_INDEX '$1==taskID {print $2}' $config) 
NCPU=16

# Set working directory
cd /scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping

# Load modules
module load bwa/0.7.17

# map the reads, with a separate mapping job for each sample
bwa mem /project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_ALL/batch1/analysis/mapping/reference_di_wol_dog.fa \
-t $NCPUS \
../trimmomatic/${sample}_1_trimpaired.fq.gz \
../trimmomatic/${sample}_2_trimpaired.fq.gz \
> /scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping/${sample}.tmp.sam

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
#PBS -J 1-19

# qsub ../mapping_sort.pbs

config=/scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping/info.txt
sample=$(awk -v taskID=$PBS_ARRAY_INDEX '$1==taskID {print $2}' $config) 
NCPU=1

echo "sample is: $sample"

# Set working directory
cd /scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping

# Load modules
module load samtools/1.9

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
# qsub ../multiqc_flagstat1.pbs

cd /project/RDS-FSC-Heartworm_MLR-RW/MultiQC

# Load modules
module load git/2.25.0
module load python/3.9.15

multiqc /scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping/flagstat1/*_flagstat1.txt -o /scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping/flagstat1
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
# qsub ../multiqc_flagstat2.pbs

cd /project/RDS-FSC-Heartworm_MLR-RW/MultiQC

# Load modules
module load git/2.25.0
module load python/3.9.15

multiqc /scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping/flagstat2/*_flagstat2.txt -o /scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping/flagstat2
```



## Extract reads that mapped to the *D. immitis* genome

If I mapped to the *D. immitis* and dog genomes separately, there could be reads that mapped to both genomes. To avoid this, I mapped to the combined D. immitis/dog genome. I can now extract the reads that mapped to only the *D. immitis* genome and use this for downstream analyses.


```bash
#!/bin/bash

# PBS directives 
#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N mapping_extract
#PBS -l select=1:ncpus=1:mem=20GB
#PBS -l walltime=06:00:00
#PBS -m e
#PBS -q defaultQ
#PBS -o mapping_extract.txt
#PBS -J 1-19

# qsub ../mapping_extract.pbs

# Set working directory
cd /scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping

# Load modules
module load samtools/1.9

config=/scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping/info.txt
sample=$(awk -v taskID=$PBS_ARRAY_INDEX '$1==taskID {print $2}' $config) 
NCPU=1

# Extract reads that only mapped to D. immitis.
samtools view -b -h -L /project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_ALL/batch1/analysis/mapping/dimmitis_WSI_2.2.bed /scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping/${sample}.sorted.bam > /scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping/${sample}_extract.bam
# Should still be in sorted form
# -b flag makes sure the output is bam
# -h flag includes the header in SAM output

samtools view /scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping/${sample}_extract.bam | head

# I do not have to sort the bam file again, it should still be sorted.

## QC
# How many D. immitis reads were extracted?

samtools flagstat /scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping/${sample}_extract.bam > /scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping/extract_flagstat/${sample}_extract_flagstat.txt
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
# qsub ../multiqc_extract_flagstat.pbs

cd /project/RDS-FSC-Heartworm_MLR-RW/MultiQC

# Load modules
module load git/2.25.0
module load python/3.9.15

multiqc /scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping/extract_flagstat/*_extract_flagstat.txt -o /scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping/extract_flagstat
```



## Index the extracted bam files

```bash
#!/bin/bash

# PBS directives 
#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N mapping_extract_index
#PBS -l select=1:ncpus=1:mem=10GB
#PBS -l walltime=02:00:00
#PBS -m e
#PBS -q defaultQ
#PBS -o mapping_extract_index.txt
#PBS -J 1-19

#qsub ../mapping_extract_index.pbs

cd /scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping

# Load modules
module load samtools/1.9

config=/scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping/info.txt
sample=$(awk -v taskID=$PBS_ARRAY_INDEX '$1==taskID {print $2}' $config) 
NCPU=1

samtools index ${sample}_extract.bam
```


## Adding read groups to bam files

I didn't add read groups during the mapping step, but luckily I can use samtools addreplacerg to add them in after mapping.

I could've done this step at the mapping stage this time but I wanted to keep things consistent with the previous batches.

```bash
#!/bin/bash

# PBS directives 
#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N read_groups
#PBS -l select=1:ncpus=1:mem=20GB
#PBS -l walltime=02:00:00
#PBS -m e
#PBS -q defaultQ
#PBS -o read_groups.txt
#PBS -J 1-19

# qsub ../read_groups.pbs

WORKING_DIR=//scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping
cd "${WORKING_DIR}"

config=/scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping/info.txt
SAMPLE_NAME=$(awk -v taskID=$PBS_ARRAY_INDEX '$1==taskID {print $2}' $config)

# Load modules
module load samtools/1.9

samtools addreplacerg -r "@RG\tRG:${SAMPLE_NAME}\tID:${SAMPLE_NAME}\tSM:${SAMPLE_NAME}" -o /scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping/bams/${SAMPLE_NAME}_rg.bam  /scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping/${SAMPLE_NAME}_extract.bam
```


## Index the read group bam files

```bash
#!/bin/bash

# PBS directives 
#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N rg_index
#PBS -l select=1:ncpus=2:mem=10GB
#PBS -l walltime=01:00:00
#PBS -m e
#PBS -q defaultQ
#PBS -o rg_index.txt
#PBS -J 1-19

# qsub ../rg_index.pbs


WORKING_DIR=/scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping/bams
cd "${WORKING_DIR}"

config=/scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping/info.txt

sample=$(awk -v taskID=$PBS_ARRAY_INDEX '$1==taskID {print $2}' $config) 
bam=/scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping/bams/${sample}_rg.bam

echo "sample is: $sample"
echo "bam is: $bam"

# Load modules
module load samtools/1.9

# index bam files
samtools index ${bam}
```


## Calling variants in a faster way using arrays

Now run on all samples:

```bash
#!/bin/bash

# PBS directives 
#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N variant_calling
#PBS -l select=1:ncpus=2:mem=20GB
#PBS -l walltime=250:00:00
#PBS -m e
#PBS -q defaultQ
#PBS -o variant_calling.txt
#PBS -M rosemonde.power@sydney.edu.au
#PBS -J 1-19


# qsub ../variant_calling.pbs


WORKING_DIR=/scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping/bams
cd "${WORKING_DIR}"

config=/scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping/info.txt
ref=/project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_ALL/batch1/analysis/mapping/reference_di_wol_dog.fa

sample=$(awk -v taskID=$PBS_ARRAY_INDEX '$1==taskID {print $2}' $config) 
bam=/scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping/bams/${sample}_rg.bam
vcf=/scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping/vcf/${sample}.g.vcf.gz

# Load modules
#module load gatk/4.2.1.0
module load gatk/4.1.4.1
module load samtools/1.9


# make gvcf per sample
gatk --java-options "-Xmx8g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
HaplotypeCaller \
-R ${ref} \
-I ${bam} \
-O ${vcf} \
-ERC GVCF
```



