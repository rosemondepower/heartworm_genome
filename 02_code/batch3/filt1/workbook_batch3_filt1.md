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

### Extra 3:

```bash
#!/bin/bash

# PBS directives 
#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N extra3_fastqc
#PBS -l select=1:ncpus=1:mem=10GB
#PBS -l walltime=00:30:00
#PBS -m e
#PBS -q defaultQ
#PBS -o extra3_fastqc.txt
#PBS -M rosemonde.power@sydney.edu.au
#PBS -J 1-6

# Submit job
## qsub ../extra3_fastqc.pbs

# Load modules
module load fastqc/0.11.8

# FastQC
cd /scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis

NCPU=2
OUTDIR=/scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/fastqc/raw

config=/scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/extra3_info.txt
sample=$(awk -v taskID=$PBS_ARRAY_INDEX '$1==taskID {print $2}' $config) 

echo "sample is: $sample"

fastqc -t $NCPU -o ${OUTDIR} /scratch/RDS-FSC-Heartworm_MLR-RW/batch3/fastq/${sample}.fq.gz
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


### Extra 3

```bash
#!/bin/bash

# PBS directives 
#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N extra3_trimmomatic
#PBS -l select=1:ncpus=2:mem=30GB
#PBS -l walltime=04:00:00
#PBS -m e
#PBS -q defaultQ
#PBS -o extra3_trimmomatic.txt
#PBS -J 1-3

# qsub ../extra3_trimmomatic.pbs

config=/scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/trimmomatic/extra3_info.txt
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

### Extra 3

```bash
#!/bin/bash

# PBS directives 
#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N extra3_fastQC_trimmed
#PBS -l select=1:ncpus=1:mem=20GB
#PBS -l walltime=02:00:00
#PBS -m e
#PBS -q defaultQ
#PBS -o extra3_fastQC_trimmed.txt
#PBS -J 1-6

# Submit job
## qsub ../extra3_fastqc_trimmed.pbs

# Load modules
module load fastqc/0.11.8

# FastQC
cd /scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis
NCPU=1
OUTDIR=/scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/fastqc/trimmed

config=/scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/extra3_info.txt
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

### Extra 3

```bash
#!/bin/bash

# PBS directives 
#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N extra3_mapping
#PBS -l select=1:ncpus=16:mem=50GB
#PBS -l walltime=02:00:00
#PBS -m abe
#PBS -q defaultQ
#PBS -o extra3_mapping.txt
#PBS -J 1-3

# qsub ../extra3_mapping.pbs

config=/scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping/extra3_info.txt
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

### Extra 3

```bash
#!/bin/bash

# PBS directives 
#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N extra3_mapping_sort
#PBS -l select=1:ncpus=1:mem=100GB
#PBS -l walltime=05:00:00
#PBS -m e
#PBS -q defaultQ
#PBS -o extra3_mapping_sort.txt
#PBS -J 1-3

# qsub ../extra3_mapping_sort.pbs

config=/scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping/extra3_info.txt
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

### Extra 3

```bash
#!/bin/bash

# PBS directives 
#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N extra3_mapping_extract
#PBS -l select=1:ncpus=1:mem=20GB
#PBS -l walltime=06:00:00
#PBS -m e
#PBS -q defaultQ
#PBS -o extra3_mapping_extract.txt
#PBS -J 1-3

# qsub ../extra3_mapping_extract.pbs

# Set working directory
cd /scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping

# Load modules
module load samtools/1.9

config=/scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping/extra3_info.txt
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

### Extra 3

```bash
#!/bin/bash

# PBS directives 
#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N extra3_mapping_extract_index
#PBS -l select=1:ncpus=1:mem=10GB
#PBS -l walltime=02:00:00
#PBS -m e
#PBS -q defaultQ
#PBS -o extra3_mapping_extract_index.txt
#PBS -J 1-3

#qsub ../extra3_mapping_extract_index.pbs

cd /scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping

# Load modules
module load samtools/1.9

config=/scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping/extra3_info.txt
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

### Extra 3

```bash
#!/bin/bash

# PBS directives 
#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N extra3_read_groups
#PBS -l select=1:ncpus=1:mem=20GB
#PBS -l walltime=02:00:00
#PBS -m e
#PBS -q defaultQ
#PBS -o extra3_read_groups.txt
#PBS -J 1-3

# qsub ../extra3_read_groups.pbs

WORKING_DIR=//scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping
cd "${WORKING_DIR}"

config=/scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping/extra3_info.txt
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

### Extra 3

```bash
#!/bin/bash

# PBS directives 
#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N extra3_rg_index
#PBS -l select=1:ncpus=2:mem=10GB
#PBS -l walltime=01:00:00
#PBS -m e
#PBS -q defaultQ
#PBS -o extra3_rg_index.txt
#PBS -J 1-3

# qsub ../extra3_rg_index.pbs


WORKING_DIR=/scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping/bams
cd "${WORKING_DIR}"

config=/scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping/extra3_info.txt

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

### Extra 3

```bash
#!/bin/bash

# PBS directives 
#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N extra3_variant_calling
#PBS -l select=1:ncpus=2:mem=20GB
#PBS -l walltime=250:00:00
#PBS -m e
#PBS -q defaultQ
#PBS -o extra3_variant_calling.txt
#PBS -M rosemonde.power@sydney.edu.au
#PBS -J 1-3


# qsub ../extra3_variant_calling.pbs


WORKING_DIR=/scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping/bams
cd "${WORKING_DIR}"

config=/scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping/extra3_info.txt
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



## Joint call variants

Merge the GVCF files we generated with HaplotypeCaller into a single GVCF file.

```bash
#!/bin/bash

# PBS directives 
#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N jointcall_variants.pbs
#PBS -l select=1:ncpus=1:mem=50GB
#PBS -l walltime=20:00:00
#PBS -m abe
#PBS -q defaultQ
#PBS -o jointcall_variants.txt
#PBS -M rosemonde.power@sydney.edu.au

# qsub ../jointcall_variants.pbs

# Perform joint genotyping on one or more samples pre-called with HaplotypeCaller

WORKING_DIR=/scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping/vcf
cd "${WORKING_DIR}"

cohort=Dirofilaria_immitis_Dec2023
ref=/project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_ALL/batch1/analysis/mapping/reference_di_wol_dog.fa

gvcf=/scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping/vcf/${cohort}.g.vcf.gz
vcf=/scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping/vcf/${cohort}.vcf.gz

# Load modules
module load gatk/4.2.1.0

# collect all sample g.vcfs (from all batches) into a list, to make input for CombineGVCFs
ls /scratch/RDS-FSC-Heartworm_MLR-RW/batch1/analysis/vcf/*.g.vcf.gz /scratch/RDS-FSC-Heartworm_MLR-RW/batch2/analysis/mapping/vcf/*.g.vcf.gz /scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping/vcf/*.g.vcf.gz > /scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping/vcf/gvcf.list

args=$(while read line; do
  echo "-V ${line}"
done < /scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping/vcf/gvcf.list)

echo $args

# create cohort gvcf
gatk --java-options "-Xmx28g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
        CombineGVCFs \
        -R ${ref} \
        ${args} \
        -O ${gvcf}


# Genotype cohort vcf
gatk --java-options "-Xmx28g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
        GenotypeGVCFs \
        -R ${ref} \
        -V ${gvcf} \
        -O ${vcf}
```


## SNPs QC

### Querying SNP and INDEL QC profiles to determine thresholds for filters

Adopted from Javier's paper.

```bash
#!/bin/bash

# PBS directives 
#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N snps_qc
#PBS -l select=1:ncpus=1:mem=10GB
#PBS -l walltime=00:20:00
#PBS -m abe
#PBS -q defaultQ
#PBS -o snps_qc.txt
#PBS -M rosemonde.power@sydney.edu.au

# qsub ../snps_qc.pbs


# load gatk
module load gatk/4.1.4.1

WORKING_DIR=/scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping/filter

# set reference, vcf, and mitochondrial and Wb contig
REFERENCE=/project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_ALL/batch1/analysis/mapping/dimmitis_WSI_2.2.fa
VCF=/scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping/filter/Dirofilaria_immitis_Dec2023.vcf.gz
MIT_CONTIG=dirofilaria_immitis_chrMtDNA
WB_CONTIG=dirofilaria_immitis_chrWb


# Make .dict file for reference sequence
## gatk CreateSequenceDictionary -R ${REFERENCE} -O /project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_ALL/data/analysis/mapping/dimmitis_WSI_2.2.dict # already did this previously

cd ${WORKING_DIR}

# select nuclear SNPs
gatk SelectVariants \
--reference ${REFERENCE} \
--variant ${VCF} \
--select-type-to-include SNP \
--exclude-intervals ${MIT_CONTIG} \
--exclude-intervals ${WB_CONTIG} \
--output ${VCF%.vcf.gz}.nuclearSNPs.vcf

# select nuclear INDELs
gatk SelectVariants \
--reference ${REFERENCE} \
--variant ${VCF} \
--select-type-to-include INDEL \
--exclude-intervals ${MIT_CONTIG} \
--exclude-intervals ${WB_CONTIG} \
--output ${VCF%.vcf.gz}.nuclearINDELs.vcf

# select mitochondrial SNPs
gatk SelectVariants \
--reference ${REFERENCE} \
--variant ${VCF} \
--select-type-to-include SNP \
--intervals ${MIT_CONTIG} \
--output ${VCF%.vcf.gz}.mitoSNPs.vcf

# select mitochondrial INDELs
gatk SelectVariants \
--reference ${REFERENCE} \
--variant ${VCF} \
--select-type-to-include INDEL \
--intervals ${MIT_CONTIG} \
--output ${VCF%.vcf.gz}.mitoINDELs.vcf

# select WB SNPs
gatk SelectVariants \
--reference ${REFERENCE} \
--variant ${VCF} \
--select-type-to-include SNP \
--intervals ${WB_CONTIG} \
--output ${VCF%.vcf.gz}.WbSNPs.vcf

# select WB INDELs
gatk SelectVariants \
--reference ${REFERENCE} \
--variant ${VCF} \
--select-type-to-include INDEL \
--intervals ${WB_CONTIG} \
--output ${VCF%.vcf.gz}.WbINDELs.vcf

# make a table of nuclear SNP data
gatk VariantsToTable \
--reference ${REFERENCE} \
--variant Dirofilaria_immitis_Dec2023.nuclearSNPs.vcf \
--fields CHROM --fields POS --fields QUAL --fields QD --fields DP --fields MQ --fields MQRankSum --fields FS --fields ReadPosRankSum --fields SOR \
--output GVCFall_nuclearSNPs.table

# make a table of nuclear INDEL data
gatk VariantsToTable \
--reference ${REFERENCE} \
--variant Dirofilaria_immitis_Dec2023.nuclearINDELs.vcf \
--fields CHROM --fields POS --fields QUAL --fields QD --fields DP --fields MQ --fields MQRankSum --fields FS --fields ReadPosRankSum --fields SOR \
--output GVCFall_nuclearINDELs.table

# make a table of mito SNP data
gatk VariantsToTable \
--reference ${REFERENCE} \
--variant Dirofilaria_immitis_Dec2023.mitoSNPs.vcf \
--fields CHROM --fields POS --fields QUAL --fields QD --fields DP --fields MQ --fields MQRankSum --fields FS --fields ReadPosRankSum --fields SOR \
--output GVCFall_mitoSNPs.table

# make a table of mito INDEL data
gatk VariantsToTable \
--reference ${REFERENCE} \
--variant Dirofilaria_immitis_Dec2023.mitoINDELs.vcf \
--fields CHROM --fields POS --fields QUAL --fields QD --fields DP --fields MQ --fields MQRankSum --fields FS --fields ReadPosRankSum --fields SOR \
--output GVCFall_mitoINDELs.table

# make a table of Wb SNP data
gatk VariantsToTable \
--reference ${REFERENCE} \
--variant Dirofilaria_immitis_Dec2023.WbSNPs.vcf \
--fields CHROM --fields POS --fields QUAL --fields QD --fields DP --fields MQ --fields MQRankSum --fields FS --fields ReadPosRankSum --fields SOR \
--output GVCFall_WbSNPs.table

# make a table of Wb INDEL data
gatk VariantsToTable \
--reference ${REFERENCE} \
--variant Dirofilaria_immitis_Dec2023.WbINDELs.vcf \
--fields CHROM --fields POS --fields QUAL --fields QD --fields DP --fields MQ --fields MQRankSum --fields FS --fields ReadPosRankSum --fields SOR \
--output GVCFall_WbINDELs.table
```




### Make some density plots of the data and get quantiles in R

```R
# Make some density plots of the data and get quantiles in R

# load libraries
library(patchwork)
require(data.table)
library(tidyverse)
library(gridExtra)

setwd("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/batch3/filter1/snps_qc")

set.seed(123)

# Import data
## Nuclear SNPs & indels
VCF_nuclear_snps <- fread('GVCFall_nuclearSNPs.table', header = TRUE, fill=TRUE, na.strings=c("","NA"), sep = "\t")
VCF_nuclear_snps <- sample_frac(VCF_nuclear_snps, 0.2) # select fraction of rows
VCF_nuclear_indels <- fread('GVCFall_nuclearINDELs.table', header = TRUE, fill=TRUE, na.strings=c("","NA"), sep = "\t")
VCF_nuclear_indels <- sample_frac(VCF_nuclear_indels, 0.2)
dim(VCF_nuclear_snps) # returns the dimensions of the data frame
dim(VCF_nuclear_indels)
VCF_nuclear <- rbind(VCF_nuclear_snps, VCF_nuclear_indels) # joins multiple rows to form a single batch
VCF_nuclear$Variant <- factor(c(rep("SNPs", dim(VCF_nuclear_snps)[1]), rep("Indels", dim(VCF_nuclear_indels)[1]))) # make a new column saying whether it's a SNP or indel

## Mitochondrial SNPs & indels
VCF_mito_snps <- fread('GVCFall_mitoSNPs.table', header = TRUE, fill=TRUE, na.strings=c("","NA"), sep = "\t")
VCF_mito_indels <- fread('GVCFall_mitoINDELs.table', header = TRUE, fill=TRUE, na.strings=c("","NA"), sep = "\t")
dim(VCF_mito_snps)
dim(VCF_mito_indels)
VCF_mito <- rbind(VCF_mito_snps, VCF_mito_indels)
VCF_mito$Variant <- factor(c(rep("SNPs", dim(VCF_mito_snps)[1]), rep("Indels", dim(VCF_mito_indels)[1])))

## Wolbachia SNPs & indels
VCF_wb_snps <- fread('GVCFall_WbSNPs.table', header = TRUE, fill=TRUE, na.strings=c("","NA"), sep = "\t")
VCF_wb_indels <- fread('GVCFall_WbINDELs.table', header = TRUE, fill=TRUE, na.strings=c("","NA"), sep = "\t")
dim(VCF_wb_snps)
dim(VCF_wb_indels)
VCF_wb <- rbind(VCF_wb_snps, VCF_wb_indels)
VCF_wb$Variant <- factor(c(rep("SNPs", dim(VCF_wb_snps)[1]), rep("Indels", dim(VCF_wb_indels)[1])))

# Set colours for the plots
snps <- '#A9E2E4'
indels <- '#F4CCCA'

# make function which makes plots
fun_variant_summaries <- function(data, title){
  # gatk hardfilter: SNP & INDEL QUAL < 0
  QUAL_quant <- quantile(data$QUAL, c(.01,.99), na.rm=T)
  
  QUAL <-
    ggplot(data, aes(x=log10(QUAL), fill=Variant)) +
    geom_density(alpha=.3) +
    geom_vline(xintercept=0, size=0.7, col="red") +
    geom_vline(xintercept=c(log10(QUAL_quant[2]), log10(QUAL_quant[3])), size=0.7, col="blue") +
    #xlim(0,10000) +
    theme_bw() +
    labs(title=paste0(title,": QUAL"))
  
  
  # DP doesnt have a hardfilter
  DP_quant <- quantile(data$DP, c(.01,.99), na.rm=T)
  
  DP <-
    ggplot(data, aes(x=log10(DP), fill=Variant)) +
    geom_density(alpha=0.3) +
    geom_vline(xintercept=log10(DP_quant), col="blue") +
    theme_bw() +
    labs(title=paste0(title,": DP"))
  
  # gatk hardfilter: SNP & INDEL QD < 2
  QD_quant <- quantile(data$QD, c(.01,.99), na.rm=T)
  
  QD <-
    ggplot(data, aes(x=QD, fill=Variant)) +
    geom_density(alpha=.3) +
    geom_vline(xintercept=2, size=0.7, col="red") +
    geom_vline(xintercept=QD_quant, size=0.7, col="blue") +
    theme_bw() +
    labs(title=paste0(title,": QD"))
  
  # gatk hardfilter: SNP FS > 60, INDEL FS > 200
  FS_quant <- quantile(data$FS, c(.01,.99), na.rm=T)
  
  FS <-
    ggplot(data, aes(x=log10(FS), fill=Variant)) +
    geom_density(alpha=.3) +
    geom_vline(xintercept=c(log10(60), log10(200)), size=0.7, col="red") +
    geom_vline(xintercept=log10(FS_quant), size=0.7, col="blue") +
    #xlim(0,250) +
    theme_bw() +
    labs(title=paste0(title,": FS"))
  
  # gatk hardfilter: SNP & INDEL MQ < 30
  MQ_quant <- quantile(data$MQ, c(.01,.99), na.rm=T)
  
  MQ <-
    ggplot(data, aes(x=MQ, fill=Variant)) + geom_density(alpha=.3) +
    geom_vline(xintercept=40, size=0.7, col="red") +
    geom_vline(xintercept=MQ_quant, size=0.7, col="blue") +
    theme_bw() +
    labs(title=paste0(title,": MQ"))
  
  # gatk hardfilter: SNP MQRankSum < -20
  MQRankSum_quant <- quantile(data$MQRankSum, c(.01,.99), na.rm=T)
  
  MQRankSum <-
    ggplot(data, aes(x=log10(MQRankSum), fill=Variant)) + geom_density(alpha=.3) +
    geom_vline(xintercept=log10(-20), size=0.7, col="red") +
    geom_vline(xintercept=log10(MQRankSum_quant), size=0.7, col="blue") +
    theme_bw() +
    labs(title=paste0(title,": MQRankSum"))
  
  
  # gatk hardfilter: SNP SOR < 4 , INDEL SOR > 10
  SOR_quant <- quantile(data$SOR, c(.01, .99), na.rm=T)
  
  SOR <-
    ggplot(data, aes(x=SOR, fill=Variant)) +
    geom_density(alpha=.3) +
    geom_vline(xintercept=c(4, 10), size=1, colour = c(snps,indels)) +
    geom_vline(xintercept=SOR_quant, size=0.7, col="blue") +
    theme_bw() +
    labs(title=paste0(title,": SOR"))
  
  # gatk hardfilter: SNP ReadPosRankSum <-10 , INDEL ReadPosRankSum < -20
  ReadPosRankSum_quant <- quantile(data$ReadPosRankSum, c(.01,.99), na.rm=T)
  
  ReadPosRankSum <-
    ggplot(data, aes(x=ReadPosRankSum, fill=Variant)) +
    geom_density(alpha=.3) +
    geom_vline(xintercept=c(-10,10,-20,20), size=1, colour = c(snps,snps,indels,indels)) +
    xlim(-10, 10) +
    geom_vline(xintercept=ReadPosRankSum_quant, size=0.7, col="blue") +
    theme_bw() +
    labs(title=paste0(title,": ReadPosRankSum"))
  
  
  plot <- QUAL + DP + QD + FS + MQ + MQRankSum + SOR + ReadPosRankSum + plot_layout(ncol=2)
  
  print(plot)
  
  ggsave(paste0("plot_",title,"_variant_summaries.png"), height=20, width=15, type="cairo")
  
  
  # generate a table of quantiles for each variant feature
  QUAL_quant <- data %>% group_by(Variant) %>% summarise(quants = list(quantile(QUAL, probs = c(0.01,0.05,0.95,0.99),na.rm=T))) %>% unnest_wider(quants)
  QUAL_quant$name <- "QUAL"
  DP_quant <- data %>% group_by(Variant) %>% summarise(quants = list(quantile(DP, probs = c(0.01,0.05,0.95,0.99),na.rm=T))) %>% unnest_wider(quants)
  DP_quant$name <- "DP"
  QD_quant <- data %>% group_by(Variant) %>% summarise(quants = list(quantile(QD, probs = c(0.01,0.05,0.95,0.99),na.rm=T))) %>% unnest_wider(quants)
  QD_quant$name <- "QD"
  FS_quant <- data %>% group_by(Variant) %>% summarise(quants = list(quantile(FS, probs = c(0.01,0.05,0.95,0.99),na.rm=T))) %>% unnest_wider(quants)
  FS_quant$name <- "FS"
  MQ_quant <- data %>% group_by(Variant) %>% summarise(quants = list(quantile(MQ, probs = c(0.01,0.05,0.95,0.99),na.rm=T))) %>% unnest_wider(quants)
  MQ_quant$name <- "MQ"
  MQRankSum_quant <- data %>% group_by(Variant) %>% summarise(quants = list(quantile(MQRankSum, probs = c(0.01,0.05,0.95,0.99),na.rm=T))) %>% unnest_wider(quants)
  MQRankSum_quant$name <- "MQRankSum"
  SOR_quant <- data %>% group_by(Variant) %>% summarise(quants = list(quantile(SOR, probs = c(0.01,0.05,0.95,0.99),na.rm=T))) %>% unnest_wider(quants)
  SOR_quant$name <- "SOR"
  ReadPosRankSum_quant <- data %>% group_by(Variant) %>% summarise(quants = list(quantile(ReadPosRankSum, probs = c(0.01,0.05,0.95,0.99),na.rm=T))) %>% unnest_wider(quants)
  ReadPosRankSum_quant$name <- "ReadPosRankSum"
  
  quantiles <- bind_rows(QUAL_quant,DP_quant, QD_quant, FS_quant, MQ_quant, MQRankSum_quant, SOR_quant, ReadPosRankSum_quant)
  quantiles$name <- c("QUAL_Indels","QUAL_SNPs","DP_indels","DP_SNPs", "QD_indels","QD_SNPs", "FS_indels","FS_SNPs", "MQ_indels","MQ_SNPs", "MQRankSum_indels","MQRankSum_SNPs", "SOR_indels","SOR_SNPs","ReadPosRankSum_indels","ReadPosRankSum_SNPs")
  
  png(paste0("table_",title,"_variant_quantiles.png"), width=1000,height=500,bg = "white")
  print(quantiles)
  grid.table(quantiles)
  dev.off()
  
}

# run nuclear variants
fun_variant_summaries(VCF_nuclear,"nuclear")
# run mitochondrial variants
fun_variant_summaries(VCF_mito,"mitochondrial")
# run wolbachia variants
fun_variant_summaries(VCF_wb,"wolbachia")






# Plot for mtDNA
##The mitochondrial genome has fewer SNPs so we might want to plot the individual points instead of density. This will give us a better idea of how to filter this data.

  # gatk hardfilter: SNP & INDEL QUAL < 0
  QUAL <-
    ggplot(VCF_mito) +
    geom_point( aes(x=log10(QUAL), y=Variant, colour=Variant)) +
    #xlim(0,10000) +
    theme_bw() +
    labs(title=paste0("mitochondrial",": QUAL"))
  QUAL
  
  # DP doesnt have a hardfilter
  DP <-
    ggplot(VCF_mito) +
    geom_point( aes(x=log10(DP), y=Variant, colour=Variant)) +
    theme_bw() +
    labs(title=paste0("mitochondrial",": DP"))
  DP
  
  # gatk hardfilter: SNP & INDEL QD < 2
  QD <-
    ggplot(VCF_mito) +
    geom_point( aes(x=QD, y=Variant, colour=Variant)) +
    theme_bw() +
    labs(title=paste0("mitochondrial",": QD"))
  QD
  # Removed 1 rows containing missing values (`geom_point()`). 
  
  # gatk hardfilter: SNP FS > 60, INDEL FS > 200
  FS <-
    ggplot(VCF_mito) +
    geom_point( aes(x=log10(FS), y=Variant, colour=Variant)) +
    #xlim(0,250) +
    theme_bw() +
    labs(title=paste0("mitochondrial",": FS"))
  FS
  
  # gatk hardfilter: SNP & INDEL MQ < 30
  MQ <-
    ggplot(VCF_mito) +
    geom_point( aes(x=MQ, y=Variant, colour=Variant)) +
    theme_bw() +
    labs(title=paste0("mitochondrial",": MQ"))
  MQ
  
  # gatk hardfilter: SNP MQRankSum < -20
  MQRankSum <-
    ggplot(VCF_mito) +
    geom_point( aes(x=log10(MQRankSum), y=Variant, colour=Variant)) +
    theme_bw() +
    labs(title=paste0("mitochondrial",": MQRankSum"))
  MQRankSum
  # Removed 100 rows containing missing values (`geom_point()`).  
  
  
  # gatk hardfilter: SNP SOR < 4 , INDEL SOR > 10
  SOR <-
    ggplot(VCF_mito) +
    geom_point( aes(x=SOR, y=Variant, colour=Variant)) +
    theme_bw() +
    labs(title=paste0("mitochondrial",": SOR"))
  SOR
  
  # gatk hardfilter: SNP ReadPosRankSum <-10 , INDEL ReadPosRankSum < -20
  ReadPosRankSum <-
    ggplot(VCF_mito) +
    geom_point( aes(x=ReadPosRankSum, y=Variant, colour=Variant)) +
    theme_bw() +
    labs(title=paste0("mitochondrial",": ReadPosRankSum"))
  ReadPosRankSum
  # Removed 8 rows containing missing values (`geom_point()`). 
  
  
  plot <- QUAL + DP + QD + FS + MQ + MQRankSum + SOR + ReadPosRankSum + plot_layout(ncol=2)
  
  print(plot)
  
  ggsave(paste0("plot_mitochondrial_indiv_variant_summaries.png"), height=20, width=15, type="cairo")
  
  
  
  
  # Repeat the same thing as above but use geom_jitter instead of geom_point.
  # gatk hardfilter: SNP & INDEL QUAL < 0
  QUAL <-
    ggplot(VCF_mito) +
    geom_jitter( aes(x=log10(QUAL), y=Variant, colour=Variant)) +
    #xlim(0,10000) +
    theme_bw() +
    labs(title=paste0("mitochondrial",": QUAL"))
  QUAL
  
  # DP doesnt have a hardfilter
  DP <-
    ggplot(VCF_mito) +
    geom_jitter( aes(x=log10(DP), y=Variant, colour=Variant)) +
    theme_bw() +
    labs(title=paste0("mitochondrial",": DP"))
  DP
  
  # gatk hardfilter: SNP & INDEL QD < 2
  QD <-
    ggplot(VCF_mito) +
    geom_jitter( aes(x=QD, y=Variant, colour=Variant)) +
    theme_bw() +
    labs(title=paste0("mitochondrial",": QD"))
  QD
  # Removed 1 rows containing missing values (`geom_point()`). 
  
  # gatk hardfilter: SNP FS > 60, INDEL FS > 200
  FS <-
    ggplot(VCF_mito) +
    geom_jitter( aes(x=log10(FS), y=Variant, colour=Variant)) +
    #xlim(0,250) +
    theme_bw() +
    labs(title=paste0("mitochondrial",": FS"))
  FS
  
  # gatk hardfilter: SNP & INDEL MQ < 30
  MQ <-
    ggplot(VCF_mito) +
    geom_jitter( aes(x=MQ, y=Variant, colour=Variant)) +
    theme_bw() +
    labs(title=paste0("mitochondrial",": MQ"))
  MQ
  
  # gatk hardfilter: SNP MQRankSum < -20
  MQRankSum <-
    ggplot(VCF_mito) +
    geom_jitter( aes(x=log10(MQRankSum), y=Variant, colour=Variant)) +
    theme_bw() +
    labs(title=paste0("mitochondrial",": MQRankSum"))
  MQRankSum
  # Removed 100 rows containing missing values (`geom_point()`). 
  
  
  # gatk hardfilter: SNP SOR < 4 , INDEL SOR > 10
  SOR <-
    ggplot(VCF_mito) +
    geom_jitter( aes(x=SOR, y=Variant, colour=Variant)) +
    theme_bw() +
    labs(title=paste0("mitochondrial",": SOR"))
  SOR
  
  # gatk hardfilter: SNP ReadPosRankSum <-10 , INDEL ReadPosRankSum < -20
  ReadPosRankSum <-
    ggplot(VCF_mito) +
    geom_jitter( aes(x=ReadPosRankSum, y=Variant, colour=Variant)) +
    theme_bw() +
    labs(title=paste0("mitochondrial",": ReadPosRankSum"))
  ReadPosRankSum
  # Removed 8 rows containing missing values (`geom_point()`). 
  
  
  plot <- QUAL + DP + QD + FS + MQ + MQRankSum + SOR + ReadPosRankSum + plot_layout(ncol=2)
  
  print(plot)
  
  ggsave(paste0("plot_mitochondrial_indiv_variant_summaries_jitter.png"), height=20, width=15, type="cairo")
```


### Attending to the quantiles, thresholds for specific parameters are established

```bash
#!/bin/bash

# PBS directives 
#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N snps_qc2
#PBS -l select=1:ncpus=1:mem=10GB
#PBS -l walltime=00:20:00
#PBS -m e
#PBS -q defaultQ
#PBS -o snps_qc2.txt

# qsub ../snps_qc2.pbs


# load gatk
module load gatk/4.1.4.1

WORKING_DIR=/scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping/filter

# set reference, vcf
REFERENCE=/project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_ALL/batch1/analysis/mapping/dimmitis_WSI_2.2.fa

cd ${WORKING_DIR}

#Nuclear
gatk VariantFiltration \
--reference ${REFERENCE} \
--variant Dirofilaria_immitis_Dec2023.nuclearSNPs.vcf \
--filter-expression 'QUAL < 32 || DP < 178 || DP > 15337 || MQ < 27.92 || SOR > 8.396 || QD < 0.310 || FS > 98.425 || MQRankSum < -7.714 || ReadPosRankSum < -9.931 || ReadPosRankSum > 3.590' \
--filter-name "SNP_filtered" \
--output Dirofilaria_immitis_Dec2023.nuclearSNPs.filtered.vcf

gatk VariantFiltration \
--reference ${REFERENCE} \
--variant Dirofilaria_immitis_Dec2023.nuclearINDELs.vcf \
--filter-expression 'QUAL < 38 || DP < 629 || DP > 12230 || MQ < 36.28 || SOR > 8.121 || QD < 0.410 || FS > 98.740 || MQRankSum < -4.361 || ReadPosRankSum < -10.940 || ReadPosRankSum > 3.160' \
--filter-name "INDEL_filtered" \
--output Dirofilaria_immitis_Dec2023.nuclearINDELs.filtered.vcf

#Mitochondrial
gatk VariantFiltration \
--reference ${REFERENCE} \
--variant Dirofilaria_immitis_Dec2023.mitoSNPs.vcf \
--filter-expression ' QUAL < 72 || DP < 112843 || DP > 455045 || MQ < 44.84 || SOR > 12.735 || QD < 0.386 || FS > 161.202 || MQRankSum < -11.087 || ReadPosRankSum < -29.759 || ReadPosRankSum > 3.735' \
--filter-name "SNP_filtered" \
--output Dirofilaria_immitis_Dec2023.mitoSNPs.filtered.vcf

gatk VariantFiltration \
--reference ${REFERENCE} \
--variant Dirofilaria_immitis_Dec2023.mitoINDELs.vcf \
--filter-expression 'QUAL < 39 || DP < 138404 || DP > 426512 || MQ < 45.26 || SOR > 7.321 || QD < 0.025 || FS > 53.312 || MQRankSum < -13.171 || ReadPosRankSum < -18.135 || ReadPosRankSum > 5.310' \
--filter-name "INDEL_filtered" \
--output Dirofilaria_immitis_Dec2023.mitoINDELs.filtered.vcf

#Wolbachia
gatk VariantFiltration \
--reference ${REFERENCE} \
--variant Dirofilaria_immitis_Dec2023.WbSNPs.vcf \
--filter-expression ' QUAL < 31 || DP < 55226 || DP > 88657 || MQ < 40.00 || SOR > 8.254 || QD < 0.320 || FS > 116.62 || MQRankSum < -6.558 || ReadPosRankSum < -11.518 || ReadPosRankSum > 4.478' \
--filter-name "SNP_filtered" \
--output Dirofilaria_immitis_Dec2023.WbSNPs.filtered.vcf

gatk VariantFiltration \
--reference ${REFERENCE} \
--variant Dirofilaria_immitis_Dec2023.WbINDELs.vcf \
--filter-expression 'QUAL < 37 || DP < 63884 || DP > 92402 || MQ < 42.18 || SOR > 8.254 || QD < 0.360 || FS > 135.318 || MQRankSum < -4.425 || ReadPosRankSum < -11.870 || ReadPosRankSum > 2.840' \
--filter-name "INDEL_filtered" \
--output Dirofilaria_immitis_Dec2023.WbINDELs.filtered.vcf

# once done, count the filtered sites - funny use of "|" allows direct markdown table format
echo -e "| Filtered_VCF | Variants_PASS | Variants_FILTERED |\n| -- | -- | -- | " > filter.stats

for i in *filtered.vcf; do
     name=${i}; pass=$( grep -E 'PASS' ${i} | wc -l ); filter=$( grep -E 'filter' ${i} | wc -l );
     echo -e "| ${name} | ${pass} | ${filter} |" >> filter.stats
done
```
This is the summary of the filtered variants ('filter.stats'):

| Filtered_VCF | Variants_PASS | Variants_FILTERED |
| -- | -- | -- | 
| Dirofilaria_immitis_Dec2023.mitoINDELs.filtered.vcf | 43 | 16 |
| Dirofilaria_immitis_Dec2023.mitoSNPs.filtered.vcf | 90 | 24 |
| Dirofilaria_immitis_Dec2023.nuclearINDELs.filtered.vcf | 554174 | 69143 |
| Dirofilaria_immitis_Dec2023.nuclearSNPs.filtered.vcf | 546844 | 71421 |
| Dirofilaria_immitis_Dec2023.WbINDELs.filtered.vcf | 5671 | 716 |
| Dirofilaria_immitis_Dec2023.WbSNPs.filtered.vcf | 4340 | 509 |


Filtering looks ok. Didn't lose too many SNPs.

Usually we would merge the SNP and INDEL files together to make a joined VCF. However, I want to disregard the indels moving forward and only focus on the SNPs. Some downstream tools don't like having indels in there.


### Filter genotypes based on x3 depth per genotype

```bash
#!/bin/bash

# PBS directives 
#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N snps_qc3
#PBS -l select=1:ncpus=1:mem=10GB
#PBS -l walltime=00:30:00
#PBS -m e
#PBS -q defaultQ
#PBS -o snps_qc3.txt

# qsub ../snps_qc3.pbs

# load gatk
module load gatk/4.1.4.1

WORKING_DIR=/scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping/filter

# set reference, vcf
REFERENCE=/project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_ALL/batch1/analysis/mapping/dimmitis_WSI_2.2.fa
VCF=/scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping/filter/Dirofilaria_immitis_Dec2023.vcf.gz

cd ${WORKING_DIR}

#Nuclear 
gatk VariantFiltration \
--reference ${REFERENCE} \
--variant ${VCF%.vcf.gz}.nuclearSNPs.filtered.vcf \
--genotype-filter-expression ' DP < 3 '  \
--genotype-filter-name "DP_lt3" \
--output ${VCF%.vcf.gz}.nuclearSNPs.DPfiltered.vcf

gatk SelectVariants \
--reference ${REFERENCE} \
--variant ${VCF%.vcf.gz}.nuclearSNPs.DPfiltered.vcf \
--set-filtered-gt-to-nocall \
--output ${VCF%.vcf.gz}.nuclearSNPs.DPfilterNoCall.vcf

#Mito
gatk VariantFiltration \
--reference ${REFERENCE} \
--variant ${VCF%.vcf.gz}.mitoSNPs.filtered.vcf \
--genotype-filter-expression ' DP < 3 '  \
--genotype-filter-name "DP_lt3" \
--output ${VCF%.vcf.gz}.mitoSNPs.DPfiltered.vcf

gatk SelectVariants \
--reference ${REFERENCE} \
--variant ${VCF%.vcf.gz}.mitoSNPs.DPfiltered.vcf \
--set-filtered-gt-to-nocall \
--output ${VCF%.vcf.gz}.mitoSNPs.DPfilterNoCall.vcf

#wolbachia
gatk VariantFiltration \
--reference ${REFERENCE} \
--variant ${VCF%.vcf.gz}.WbSNPs.filtered.vcf \
--genotype-filter-expression ' DP < 3 '  \
--genotype-filter-name "DP_lt3" \
--output ${VCF%.vcf.gz}.WbSNPs.DPfiltered.vcf

gatk SelectVariants \
--reference ${REFERENCE} \
--variant ${VCF%.vcf.gz}.WbSNPs.DPfiltered.vcf \
--set-filtered-gt-to-nocall \
--output ${VCF%.vcf.gz}.WbSNPs.DPfilterNoCall.vcf
```


### Now we apply a set of standard filters for population genomics

```bash
#!/bin/bash

# PBS directives 
#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N snps_qc4
#PBS -l select=1:ncpus=1:mem=4GB
#PBS -l walltime=00:10:00
#PBS -m e
#PBS -q defaultQ
#PBS -o snps_qc4.txt

# qsub ../snps_qc4.pbs

# load gatk
module load gatk/4.1.4.1
module load vcftools/0.1.14

WORKING_DIR=/scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping/filter

# set vcf
VCF=/scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping/filter/Dirofilaria_immitis_Dec2023.vcf.gz


cd ${WORKING_DIR}

#Nuclear variants
vcftools \
--vcf ${VCF%.vcf.gz}.nuclearSNPs.DPfilterNoCall.vcf \
--remove-filtered-geno-all \
--remove-filtered-all \
--min-alleles 2 \
--max-alleles 2 \
--hwe 1e-06 \
--maf 0.02 \
--recode \
--recode-INFO-all \
--out ${VCF%.vcf.gz}.nuclear_SNPs.final
#After filtering, kept 83 out of 83 Individuals
#Outputting VCF file...
#After filtering, kept 256535 out of a possible 618259 Sites

# Changed --maf from 0.05 to 0.02 to follow Javier's code

#--- nuclear SNPs
vcftools --vcf Dirofilaria_immitis_Dec2023.nuclear_SNPs.final.recode.vcf --remove-indels
#After filtering, kept 83 out of 83 Individuals
#After filtering, kept 256535 out of a possible 256535 Sites

#--- nuclear  INDELs
vcftools --vcf Dirofilaria_immitis_Dec2023.nuclear_SNPs.final.recode.vcf --keep-only-indels
#After filtering, kept 83 out of 83 Individuals
#After filtering, kept 0 out of a possible 256535 Sites




#Mitochondrial variants
vcftools \
--vcf ${VCF%.vcf.gz}.mitoSNPs.DPfilterNoCall.vcf \
--remove-filtered-geno-all \
--remove-filtered-all \
--min-alleles 2 \
--max-alleles 2 \
--maf 0.02 \
--recode \
--recode-INFO-all \
--out ${VCF%.vcf.gz}.mito_SNPs.final
#After filtering, kept 83 out of 83 Individuals
#Outputting VCF file...
#After filtering, kept 46 out of a possible 108 Sites


#--- mito SNPs
vcftools --vcf Dirofilaria_immitis_Dec2023.mito_SNPs.final.recode.vcf --remove-indels
#After filtering, kept 83 out of 83 Individuals
#After filtering, kept 46 out of a possible 46 Sites

#--- mito INDELs
vcftools --vcf Dirofilaria_immitis_Dec2023.mito_SNPs.final.recode.vcf --keep-only-indels
#After filtering, kept 83 out of 83 Individuals
#After filtering, kept 0 out of a possible 46 Sites



#wolbachia variants
vcftools \
--vcf ${VCF%.vcf.gz}.WbSNPs.DPfilterNoCall.vcf \
--remove-filtered-geno-all \
--remove-filtered-all \
--min-alleles 2 \
--max-alleles 2 \
--maf 0.02 \
--recode \
--recode-INFO-all \
--out ${VCF%.vcf.gz}.Wb_SNPs.final
#After filtering, kept 83 out of 83 Individuals
#Outputting VCF file...
#After filtering, kept 467 out of a possible 4843 Sites


#--- Wb SNPs
vcftools --vcf Dirofilaria_immitis_Dec2023.Wb_SNPs.final.recode.vcf --remove-indels
#After filtering, kept 83 out of 83 Individuals
#After filtering, kept 467 out of a possible 467 Sites

#--- Wb INDELs
vcftools --vcf Dirofilaria_immitis_Dec2023.Wb_SNPs.final.recode.vcf --keep-only-indels
#After filtering, kept 83 out of 83 Individuals
#After filtering, kept 0 out of a possible 467 Sites
```

### Now, we are filtering by missingness

```bash
#!/bin/bash

# PBS directives 
#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N snps_qc5
#PBS -l select=1:ncpus=1:mem=10GB
#PBS -l walltime=00:10:00
#PBS -m e
#PBS -q defaultQ
#PBS -o snps_qc5.txt

# qsub ../snps_qc5.pbs

# load modules
module load gatk/4.1.4.1
module load vcftools/0.1.14

WORKING_DIR=/scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping/filter

# set vcf
VCF=/scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping/filter/Dirofilaria_immitis_Dec2023.vcf.gz

cd ${WORKING_DIR}

#determine missingness per individual
vcftools --vcf ${VCF%.vcf.gz}.nuclear_SNPs.final.recode.vcf --out nuclear --missing-indv
vcftools --vcf ${VCF%.vcf.gz}.mito_SNPs.final.recode.vcf --out mito --missing-indv
vcftools --vcf ${VCF%.vcf.gz}.Wb_SNPs.final.recode.vcf --out wb --missing-indv
```

### Check the missingess in R

```R
 # Check missingness
  data_nuclear <- read.delim("nuclear.imiss", header=T)
  data_mito <- read.delim("mito.imiss", header=T)
  data_wb <- read.delim("wb.imiss", header=T)
  
  #creating the function - per sample
  fun_plot_missingness <- function(data,title) {
    
    plot <- ggplot(data, aes(INDV, 1-F_MISS)) +
      geom_boxplot(color = 'brown') +
      geom_point(size = 1, color = 'brown4') +
      theme_bw() +
      labs(x="Sample ID", y="Proportion of total variants present (1-missingness)")+
      ggtitle(title) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    print(plot)
    ggsave(paste0("plot_missingness_figure_",title,".png"))
  }
  
  # plotting for each dataset
  fun_plot_missingness(data_nuclear, "nuclear_variants")

  fun_plot_missingness(data_mito,"mitochondrial_variants")

  fun_plot_missingness(data_wb, "wb_variants")
```

In Javier's code, he generated a different sample list for each database and evaluated the max missingness. For now, I will just keep all the samples except Pavia ERR034942-3 (these are paired end, I could've combined them earlier but I already have ERR034941 which is a mate pair library of the same sample). Also now have the extra 3 (JS6665, JS6670, JS6675).

```bash
# For nuclear (n=81) - nuclear_samplelist.keep
JS6277
JS6278
JS6279
JS6280
JS6281
JS6342
JS6343
JS6344
JS6345
JS6346
JS6347
JS6349
JS6350
JS6351
JS6352
JS6353
JS6354
JS6355
JS6356
JS6357
JS6358
JS6359
JS6360
JS6368
JS6369
JS6370
SRR13154013
SRR13154014
SRR13154015
SRR13154016
SRR13154017
ERR034940
ERR034941
JS6597_DKDN230025568-1A_HHMLHDSX7_L1
JS6598_DKDN230025569-1A_HHMLHDSX7_L1
JS6599_DKDN230025570-1A_HHMLHDSX7_L1
JS6600_DKDN230025571-1A_HHMLHDSX7_L1
JS6601_DKDN230025572-1A_HHMLHDSX7_L1
JS6602_DKDN230025573-1A_HHMLHDSX7_L1
JS6603_DKDN230025574-1A_HHMLHDSX7_L1
JS6604_DKDN230025575-1A_HHMLHDSX7_L1
JS6605_DKDN230025576-1A_HHMLHDSX7_L1
JS6606_DKDN230025577-1A_HHMLHDSX7_L1
JS6607_DKDN230025578-1A_HHMLHDSX7_L4
JS6608_DKDN230025579-1A_HHMLHDSX7_L4
JS6609_DKDN230025580-1A_HHMLHDSX7_L4
JS6610_DKDN230025581-1A_HHMLHDSX7_L1
JS6611_DKDN230025582-1A_HHMLHDSX7_L1
JS6612_DKDN230025583-1A_HHMLHDSX7_L1
JS6613_DKDN230025584-1A_HHMLHDSX7_L1
JS6614_DKDN230025585-1A_HHMLHDSX7_L1
JS6615_DKDN230025586-1A_HHMLHDSX7_L1
JS6616_DKDN230025587-1A_HHMLHDSX7_L1
JS6617_DKDN230025588-1A_HHMLHDSX7_L1
SRR10533236
SRR10533237
SRR10533238
SRR10533239
SRR10533240
JS6648_merged
JS6649_DKDN230043078-1A_HC2NMDSX7_L2
JS6650_DKDN230043079-1A_HC2NMDSX7_L2
JS6651_DKDN230043080-1A_HC2NMDSX7_L2
JS6652_DKDN230043081-1A_HC2NMDSX7_L2
JS6653_DKDN230043082-1A_HC2NMDSX7_L2
JS6654_DKDN230043083-1A_HC2NMDSX7_L2
JS6655_DKDN230043084-1A_HC2NMDSX7_L2
JS6656_merged
JS6657_merged
JS6658_DKDN230043087-1A_HC2NMDSX7_L2
JS6659_merged
JS6660_DKDN230043089-1A_H7WVCDSX7_L1
JS6664_DKDN230043093-1A_H7WVCDSX7_L1
JS6666_DKDN230043095-1A_H7WVCDSX7_L1
JS6667_DKDN230043096-1A_H7WVCDSX7_L4
JS6668_DKDN230043097-1A_H7WVCDSX7_L1
JS6671_DKDN230043100-1A_H7WVCDSX7_L1
JS6678_DKDN230043107-1A_H7WVCDSX7_L1
JS6665_DKDN230043094-1A_HFKCLDSX7_L1
JS6670_DKDN230043099-1A_HFKCLDSX7_L1
JS6675_DKDN230043104-1A_HFKCLDSX7_L1


# For mithochondiral (n=81) - mito_samplelist.keep
#Same as above

# For wb (n=81) - wb_samplelist.keep
#same as above
```


### Let's check different thresholds for each dataset

```bash
#!/bin/bash

# PBS directives 
#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N snps_qc6
#PBS -l select=1:ncpus=1:mem=10GB
#PBS -l walltime=00:20:00
#PBS -m e
#PBS -q defaultQ
#PBS -o snps_qc6.txt

# qsub ../snps_qc6.pbs

# load gatk
module load gatk/4.1.4.1
module load vcftools/0.1.14

WORKING_DIR=/scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping/filter

# set vcf
VCF=/scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping/filter/Dirofilaria_immitis_Dec2023.vcf.gz

cd ${WORKING_DIR}

# For nuclear variants
for i in 0.7 0.8 0.9 1; do
     vcftools --vcf ${VCF%.vcf.gz}.nuclear_SNPs.final.recode.vcf --keep nuclear_samplelist.keep --max-missing ${i} ;
done

# max-missing = 0.7
#After filtering, kept 81 out of 83 Individuals
#After filtering, kept 251542 out of a possible 256535 Sites

# max-missing = 0.8
#After filtering, kept 81 out of 83 Individuals
#After filtering, kept 249005 out of a possible 256535 Sites

# max-missing = 0.9
# After filtering, kept 81 out of 83 Individuals
#After filtering, kept 240259 out of a possible 256535 Sites

# max-missing = 1
# After filtering, kept 81 out of 83 Individuals
#After filtering, kept 10479 out of a possible 256535 Sites


# For mito variants
for i in 0.7 0.8 0.9 1; do
     vcftools --vcf ${VCF%.vcf.gz}.mito_SNPs.final.recode.vcf --keep mito_samplelist.keep --max-missing ${i} ;
done

# max-missing = 0.7
# After filtering, kept 81 out of 83 Individuals
#After filtering, kept 46 out of a possible 46 Sites

# max-missing = 0.8
# After filtering, kept 81 out of 83 Individuals
#After filtering, kept 46 out of a possible 46 Sites

# max-missing = 0.9
# After filtering, kept 81 out of 83 Individuals
#After filtering, kept 46 out of a possible 46 Sites

# max-missing = 1
# After filtering, kept 81 out of 83 Individuals
#After filtering, kept 45 out of a possible 46 Sites


# For Wb variants
for i in 0.7 0.8 0.9 1; do
     vcftools --vcf ${VCF%.vcf.gz}.Wb_SNPs.final.recode.vcf --keep wb_samplelist.keep --max-missing ${i} ;
done

# max-missing = 0.7
# After filtering, kept 81 out of 83 Individuals
#After filtering, kept 467 out of a possible 467 Sites

# max-missing = 0.8
# After filtering, kept 81 out of 83 Individuals
#After filtering, kept 467 out of a possible 467 Sites

# max-missing = 0.9
# After filtering, kept 81 out of 83 Individuals
#After filtering, kept 415 out of a possible 467 Sites

# max-missing = 1
# After filtering, kept 81 out of 83 Individuals
#After filtering, kept 1 out of a possible 467 Sites
```

Selecting a max missingness of 0.9 for nuclear, 1 for mito and 0.9 for Wb is sensible.

```bash
#!/bin/bash

# PBS directives 
#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N snps_qc7
#PBS -l select=1:ncpus=1:mem=10GB
#PBS -l walltime=00:10:00
#PBS -m e
#PBS -q defaultQ
#PBS -o snps_qc7.txt

# qsub ../snps_qc7.pbs

# load gatk
module load gatk/4.1.4.1
module load vcftools/0.1.14

WORKING_DIR=/scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping/filter

# set vcf
VCF=/scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping/filter/Dirofilaria_immitis_Dec2023.vcf.gz

cd ${WORKING_DIR}

# For nuclear
vcftools --vcf ${VCF%.vcf.gz}.nuclear_SNPs.final.recode.vcf \
     --keep nuclear_samplelist.keep \
     --max-missing 0.9 \
     --recode --recode-INFO-all \
     --out FINAL_SETS/nuclear_samples3x_missing0.9

# For mito
vcftools --vcf ${VCF%.vcf.gz}.mito_SNPs.final.recode.vcf \
     --keep mito_samplelist.keep \
     --max-missing 1 \
     --recode --recode-INFO-all \
     --out FINAL_SETS/mito_samples3x_missing1

# For wb
vcftools --vcf ${VCF%.vcf.gz}.Wb_SNPs.final.recode.vcf \
     --keep wb_samplelist.keep \
     --max-missing 0.9 \
     --recode --recode-INFO-all \
     --out FINAL_SETS/wb_samples3x_missing0.9
```


### Also, we will select only the variants in the chr 1 to chr4, avoiding the chrX and the scaffolds

```bash
#!/bin/bash

# PBS directives 
#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N snps_qc8
#PBS -l select=1:ncpus=1:mem=4GB
#PBS -l walltime=00:10:00
#PBS -m e
#PBS -q defaultQ
#PBS -o snps_qc8.txt

# qsub ../snps_qc8.pbs

# load gatk
module load vcftools/0.1.14

cd /scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping/filter

vcftools --vcf FINAL_SETS/nuclear_samples3x_missing0.9.recode.vcf \
--chr dirofilaria_immitis_chr1 \
--chr dirofilaria_immitis_chr2 \
--chr dirofilaria_immitis_chr3 \
--chr dirofilaria_immitis_chr4 \
--recode --out FINAL_SETS/nuclear_samples3x_missing0.9.chr1to4
# After filtering, kept 81 out of 81 Individuals
#Outputting VCF file...
#After filtering, kept 175415 out of a possible 240259 Sites


vcftools --vcf FINAL_SETS/nuclear_samples3x_missing0.9.chr1to4.recode.vcf --remove-indels
# After filtering, kept 81 out of 81 Individuals
#After filtering, kept 175415 out of a possible 175415 Sites


vcftools --vcf FINAL_SETS/nuclear_samples3x_missing0.9.chr1to4.recode.vcf --keep-only-indels
# After filtering, kept 81 out of 81 Individuals
#After filtering, kept 0 out of a possible 175415 Sites- this makes sense because I removed the indels earlier and only focused on the SNPs.
```

Rename samples to more meaningful names:

```bash
#!/bin/bash

# PBS directives 
#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N vcf_rename
#PBS -l select=1:ncpus=2:mem=20GB
#PBS -l walltime=00:20:00
#PBS -m abe
#PBS -q defaultQ
#PBS -o vcf_rename.txt
#PBS -M rosemonde.power@sydney.edu.au

# qsub ../vcf_rename.pbs

# load modules
module load bcftools/1.11

cd /scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping/filter/FINAL_SETS

bcftools reheader -s rename.txt nuclear_samples3x_missing0.9.chr1to4.recode.vcf -o nuclear_samples3x_missing0.9.chr1to4.recode.RENAMED.vcf
```

rename.txt:
ERR034940	ERR034940
ERR034941	ERR034941
JS6277	Dog2.1
JS6278	Dog2.2
JS6279	Dog2.3
JS6280	Dog2.4
JS6281	JS6281
JS6342	JS6342
JS6343	JS6343
JS6344	JS6344
JS6345	JS6345
JS6346	Fox1
JS6347	Fox2
JS6349	JS6349
JS6350	JS6350
JS6351	JS6351
JS6352	JS6352
JS6353	JS6353
JS6354	JS6354
JS6355	JS6355
JS6356	JS6356
JS6357	JS6357
JS6358	JS6358
JS6359	JS6359
JS6360	JS6360
JS6368	JS6368
JS6369	JS6369
JS6370	JS6370
JS6597_DKDN230025568-1A_HHMLHDSX7_L1	JS6597
JS6598_DKDN230025569-1A_HHMLHDSX7_L1	JS6598
JS6599_DKDN230025570-1A_HHMLHDSX7_L1	JS6599
JS6600_DKDN230025571-1A_HHMLHDSX7_L1	JS6600
JS6601_DKDN230025572-1A_HHMLHDSX7_L1	Dog3.1
JS6602_DKDN230025573-1A_HHMLHDSX7_L1	Dog3.2
JS6603_DKDN230025574-1A_HHMLHDSX7_L1	Dog3.3
JS6604_DKDN230025575-1A_HHMLHDSX7_L1	JS6604
JS6605_DKDN230025576-1A_HHMLHDSX7_L1	JS6605
JS6606_DKDN230025577-1A_HHMLHDSX7_L1	JS6606
JS6607_DKDN230025578-1A_HHMLHDSX7_L4	JS6344-REP
JS6608_DKDN230025579-1A_HHMLHDSX7_L4	Dog2.3-REP
JS6609_DKDN230025580-1A_HHMLHDSX7_L4	Dog1.5-REP
JS6610_DKDN230025581-1A_HHMLHDSX7_L1	JS6610
JS6611_DKDN230025582-1A_HHMLHDSX7_L1	JS6611
JS6612_DKDN230025583-1A_HHMLHDSX7_L1	JS6612
JS6613_DKDN230025584-1A_HHMLHDSX7_L1	JS6613
JS6614_DKDN230025585-1A_HHMLHDSX7_L1	JS6614
JS6615_DKDN230025586-1A_HHMLHDSX7_L1	JS6615
JS6616_DKDN230025587-1A_HHMLHDSX7_L1	JS6616
JS6617_DKDN230025588-1A_HHMLHDSX7_L1	JS6617
SRR10533236	SRR10533236
SRR10533237	SRR10533237
SRR10533238	SRR10533238
SRR10533239	SRR10533239
SRR10533240	SRR10533240
SRR13154013	Dog1.1
SRR13154014	Dog1.2
SRR13154015	Dog1.3
SRR13154016	Dog1.4
SRR13154017	Dog1.5
JS6648_merged	JS6648
JS6649_DKDN230043078-1A_HC2NMDSX7_L2	JS6649
JS6650_DKDN230043079-1A_HC2NMDSX7_L2	JS6650
JS6651_DKDN230043080-1A_HC2NMDSX7_L2	JS6651
JS6652_DKDN230043081-1A_HC2NMDSX7_L2	JS6652
JS6653_DKDN230043082-1A_HC2NMDSX7_L2	JS6653
JS6654_DKDN230043083-1A_HC2NMDSX7_L2	JS6654
JS6655_DKDN230043084-1A_HC2NMDSX7_L2	JS6655
JS6656_merged	JS6656
JS6657_merged	JS6657
JS6658_DKDN230043087-1A_HC2NMDSX7_L2	JS6658
JS6659_merged	JS6659
JS6660_DKDN230043089-1A_H7WVCDSX7_L1	JS6660
JS6664_DKDN230043093-1A_H7WVCDSX7_L1	JS6664
JS6666_DKDN230043095-1A_H7WVCDSX7_L1	JS6666
JS6667_DKDN230043096-1A_H7WVCDSX7_L4	JS6667
JS6668_DKDN230043097-1A_H7WVCDSX7_L1	JS6668
JS6671_DKDN230043100-1A_H7WVCDSX7_L1	JS6671
JS6678_DKDN230043107-1A_H7WVCDSX7_L1	JS6678
JS6665_DKDN230043094-1A_HFKCLDSX7_L1	JS6665
JS6670_DKDN230043099-1A_HFKCLDSX7_L1	JS6670
JS6675_DKDN230043104-1A_HFKCLDSX7_L1	JS6675



Check that it worked:

```bash
# List all sample names in the original VCF
bcftools query -l /scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping/filter/FINAL_SETS/nuclear_samples3x_missing0.9.chr1to4.recode.RENAMED.vcf
```

Yep all the sample names are changed.


## Select variants for chr1, chr2, chr3 and chr4 (SEPARATELY)

```bash
#!/bin/bash

# PBS directives 
#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N snps_chr1-4_separate
#PBS -l select=1:ncpus=1:mem=2GB
#PBS -l walltime=00:20:00
#PBS -m e
#PBS -q defaultQ
#PBS -o snps_chr1-4_separate.txt

# qsub ../snps_chr1-4_separate.pbs

# load gatk
module load vcftools/0.1.14

cd /scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping/filter/FINAL_SETS

# chr1
vcftools --vcf nuclear_samples3x_missing0.9.chr1to4.recode.RENAMED.vcf \
--chr dirofilaria_immitis_chr1 \
--recode --out nuclear_samples3x_missing0.9.chr1
#After filtering, kept 81 out of 81 Individuals
#Outputting VCF file...
#After filtering, kept 46816 out of a possible 175415 Sites

# chr2
vcftools --vcf nuclear_samples3x_missing0.9.chr1to4.recode.RENAMED.vcf \
--chr dirofilaria_immitis_chr2 \
--recode --out nuclear_samples3x_missing0.9.chr2
# After filtering, kept 81 out of 81 Individuals
#Outputting VCF file...
#After filtering, kept 38423 out of a possible 175415 Sites

# chr3
vcftools --vcf nuclear_samples3x_missing0.9.chr1to4.recode.RENAMED.vcf \
--chr dirofilaria_immitis_chr3 \
--recode --out nuclear_samples3x_missing0.9.chr3
# After filtering, kept 81 out of 81 Individuals
#Outputting VCF file...
#After filtering, kept 46302 out of a possible 175415 Sites

# chr4
vcftools --vcf nuclear_samples3x_missing0.9.chr1to4.recode.RENAMED.vcf \
--chr dirofilaria_immitis_chr4 \
--recode --out nuclear_samples3x_missing0.9.chr4
# After filtering, kept 81 out of 81 Individuals
#Outputting VCF file...
#After filtering, kept 43874 out of a possible 175415 Sites
```