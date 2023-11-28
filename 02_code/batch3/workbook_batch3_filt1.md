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



## FastQC & Multi-QC on the merged files

```bash
#!/bin/bash

# PBS directives 
#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N fastqc_merged
#PBS -l select=2:ncpus=2:mem=30GB
#PBS -l walltime=05:00:00
#PBS -m abe
#PBS -q defaultQ
#PBS -o fastqc_merged.txt
#PBS -M rosemonde.power@sydney.edu.au

# Submit job
## qsub fastqc_merged.pbs

# Load modules
module load fastqc/0.11.8

# FastQC
cd /scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/fastqc

INPUTDIR="/scratch/RDS-FSC-Heartworm_MLR-RW/batch3/fastq"
NCPU=24
OUTDIR="/scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/fastqc"

fastqc -t $NCPU -o $OUTDIR $INPUTDIR/*.fq.gz
fastqc -t $NCPU -o $OUTDIR $INPUTDIR/*.fastq.gz
```


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









```bash
#!/bin/bash

# PBS directives 
#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N fastqc_daisy
#PBS -l select=1:ncpus=24:mem=40GB
#PBS -l walltime=01:00:00
#PBS -m e
#PBS -q defaultQ
#PBS -o fastqc_daisy.txt

# Submit job
## qsub fastqc_daisy.pbs

# Load modules
module load fastqc/0.11.8

# FastQC
cd /project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_ALL/data/analysis/fastqc/merged

INPUTDIR="/project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_ALL/data/fastq/merged"
NCPU=24
OUTDIR="/project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_ALL/data/analysis/fastqc/merged"

fastqc -t $NCPU -o $OUTDIR $INPUTDIR/*.fastq.gz
```

```bash
#!/bin/bash

# PBS directives 
#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N multiqc_merged
#PBS -l select=1:ncpus=1:mem=2GB
#PBS -l walltime=00:03:00
#PBS -m e
#PBS -q defaultQ
#PBS -o multiqc_merged.txt

# Submit job
# qsub ../multiqc_merged.pbs

# Run MultiQC to combine all of the FastQC reports for the fastq files (fastq files of samples with same library were merged)

cd /project/RDS-FSC-Heartworm_MLR-RW/MultiQC

# Load modules
module load git/2.25.0
module load python/3.9.15

multiqc /project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_ALL/data/analysis/fastqc/merged -o /project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_ALL/data/analysis/fastqc/merged
```

![](output/multiqc_merged/fastqc_sequence_counts_plot.png)

Some samples have much fewer reads. This makes sense because some samples were only sequenced for 1GB due to low quality.


![](output/multiqc_merged/fastqc_per_base_sequence_quality_plot.png)
![](output/multiqc_merged/fastqc_per_sequence_quality_scores_plot.png)

Quality looks pretty good.


![](output/multiqc_merged/fastqc_per_sequence_gc_content_plot.png)

JS6278 (this was one of Wilson's HWs) and JS6348 (this was the D. roemeri sample) failed the GC content check. 