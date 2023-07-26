# Dirofilaria immitis WGS Lab Book - Extra Data

### Rose Power USYD 2023

In this workbook, I will obtain the additional public data from Javier's paper and analyse this data to get the GVCF files. Then I will do joint variant calling using my original 31 samples and this extra data.


## Get public data mentioned in Javier's paper

```bash
#!/bin/bash

# PBS directives 
#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N shin_godel
#PBS -l select=1:ncpus=16:mem=50GB
#PBS -l walltime=03:00:00
#PBS -m e
#PBS -q defaultQ
#PBS -o shin_godel.txt

# Submit job
## qsub ../shin_godel.pbs

cd /scratch/RDS-FSC-Heartworm_MLR-RW/mapping/steve_data/fastq

# (Shi 2020) data
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/039/SRR10533239/SRR10533239_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/036/SRR10533236/SRR10533236_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/038/SRR10533238/SRR10533238_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/040/SRR10533240/SRR10533240_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/040/SRR10533240/SRR10533240_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/039/SRR10533239/SRR10533239_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/037/SRR10533237/SRR10533237_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/036/SRR10533236/SRR10533236_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/037/SRR10533237/SRR10533237_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/038/SRR10533238/SRR10533238_2.fastq.gz

# (Godel 2012)
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR034/ERR034941/ERR034941_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR034/ERR034943/ERR034943_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR034/ERR034940/ERR034940_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR034/ERR034940/ERR034940_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR034/ERR034942/ERR034942_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR034/ERR034942/ERR034942_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR034/ERR034941/ERR034941_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR034/ERR034943/ERR034943_1.fastq.gz
```

## FastQC & Multi-QC

We want to get some stats on the raw data.

```bash
#!/bin/bash

# PBS directives 
#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N fastqc
#PBS -l select=1:ncpus=1:mem=10GB
#PBS -l walltime=00:20:00
#PBS -m e
#PBS -q defaultQ
#PBS -o fastqc.txt
#PBS -M rosemonde.power@sydney.edu.au
#PBS -J 1-18

# Submit job
## qsub ../fastqc.pbs

# Load modules
module load fastqc/0.11.8

# FastQC
cd /scratch/RDS-FSC-Heartworm_MLR-RW/mapping/steve_data/analysis

NCPU=2
OUTDIR=/scratch/RDS-FSC-Heartworm_MLR-RW/mapping/steve_data/analysis/fastqc/raw

config=/scratch/RDS-FSC-Heartworm_MLR-RW/mapping/steve_data/analysis/fastqc/raw/info.txt
sample=$(awk -v taskID=$PBS_ARRAY_INDEX '$1==taskID {print $2}' $config) 

echo "sample is: $sample"

fastqc -t $NCPU -o ${OUTDIR} /scratch/RDS-FSC-Heartworm_MLR-RW/mapping/steve_data/fastq/${sample}.fastq.gz
```

## MultiQC

```bash
#!/bin/bash

# PBS directives 
#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N multiqc
#PBS -l select=1:ncpus=1:mem=4GB
#PBS -l walltime=00:10:00
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
multiqc /scratch/RDS-FSC-Heartworm_MLR-RW/mapping/steve_data/analysis/fastqc/raw -o /scratch/RDS-FSC-Heartworm_MLR-RW/mapping/steve_data/analysis/fastqc/raw
```


## Trimming

### Trimmomatic

```bash
#!/bin/bash

# PBS directives 
#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N trimmomatic
#PBS -l select=1:ncpus=2:mem=20GB
#PBS -l walltime=00:45:00
#PBS -m e
#PBS -q defaultQ
#PBS -o trimmomatic.txt
#PBS -J 1-9

# qsub ../trimmomatic.pbs

config=/scratch/RDS-FSC-Heartworm_MLR-RW/mapping/steve_data/analysis/info.txt
sample=$(awk -v taskID=$PBS_ARRAY_INDEX '$1==taskID {print $2}' $config) 
NCPU=2

echo "sample is: $sample"

cd /scratch/RDS-FSC-Heartworm_MLR-RW/mapping/steve_data/analysis/trimmomatic

# Load modules
module load trimmomatic/0.38

# Run Trimmomatic
java -jar /usr/local/trimmomatic/0.38/trimmomatic-0.38.jar PE \
-threads $NCPU \
-phred33 \
/scratch/RDS-FSC-Heartworm_MLR-RW/mapping/steve_data/fastq/${sample}_1.fastq.gz \
/scratch/RDS-FSC-Heartworm_MLR-RW/mapping/steve_data/fastq/${sample}_2.fastq.gz \
${sample}_1_trimpaired.fq.gz ${sample}_1_trimunpaired.fq.gz \
${sample}_2_trimpaired.fq.gz ${sample}_2_trimunpaired.fq.gz \
SLIDINGWINDOW:10:20 MINLEN:50

# SLIDINGWINDOW:10:20 means it will scan the read with a 10-base wide sliding window, cutting when the average quality per base drops below 20.

# Instead of SLIDINGWINDOW, in my previous practice code I used 'AVGQUAL:30 MINLEN:150'.
```
Assume all files are phred33 quality encoded? The SRR files have '????' in the quality scores so I'll have to check this somehow..

