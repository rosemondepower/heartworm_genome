# Using Nimbus supercomputer

## Log into Nimbus

```bash
ssh -i `name_of_your_key`.pem `login_name`@###.###.###.###
```

## Transfer files into Nimbus

```bash
# from local terminal, type:
scp -i `name_of_your_key`.pem Z:/PRJ-Heartworm_MLR/project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_ALL/data/analysis/mapping/local_file `ubuntu`@###.###.###.###:`instance_file_path`

# transferring whole directory
scp -r -i heartworm.pem Z:/PRJ-Heartworm_MLR/scratch/RDS-FSC-Heartworm_MLR-RW/mapping/bams/bams_rg ubuntu@###.###.###.###8:data/
```

But it took 2+ hours to transfer 1 bam file into my Nimbus directory...


## Using containers

```bash
# Export the paths for your data directory(ies) to Singularity, so that they can be readable by the container.
export SINGULARITY_BINDPATH=/data
echo 'export SINGULARITY_BINDPATH=/data' >> ~/.bashrc

shpc show quay.io/biocontainers/gatk
module avail gatk
module load quay.io/biocontainers/gatk4/4.4.0.0--py36hdfd78af_0/module
gatk --version
module list


gatk --java-options "-Xmx4g" HaplotypeCaller \
-R reference_di_wol_dog.fa \
-I JS6277_rg.bam \
-O JS6277.g.vcf.gz \
-ERC GVCF


# redirect logs
gatk --java-options "-Xmx4g" HaplotypeCaller \
-R reference_di_wol_dog.fa \
-I JS6277_rg.bam \
-O JS6277.g.vcf.gz \
-ERC GVCF > output.log 2> err.log

sample=
for i in ${sample}
do
    gatk --java-options "-Xmx4g" HaplotypeCaller \
    -R reference_di_wol_dog.fa \
    -I ${i}_rg.bam \
    -O ${i}.g.vcf.gz \
    -ERC GVCF
done


```


- learn how to get output error and text files (> and 2>)
- learn how to run things in the background
- learn how to run multiple samples in a loop



```bash
#!/bin/bash

# PBS directives 
#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N variant_calling
#PBS -l select=1:ncpus=24:mem=80GB
#PBS -l walltime=72:00:00
#PBS -m e
#PBS -q defaultQ
#PBS -o variant_calling.txt
#PBS -J 1-2

# qsub ../variant_calling.pbs


WORKING_DIR=/scratch/RDS-FSC-Heartworm_MLR-RW/mapping/bams
cd "${WORKING_DIR}"

config=/project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_ALL/data/analysis/info.txt
ref=/project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_ALL/data/analysis/mapping/reference_di_wol_dog.fa

sample=$(awk -v taskID=$PBS_ARRAY_INDEX '$1==taskID {print $2}' $config) 
bam=/scratch/RDS-FSC-Heartworm_MLR-RW/mapping/bams/${sample}_rg.bam
vcf=/scratch/RDS-FSC-Heartworm_MLR-RW/mapping/vcf/${sample}.g.vcf.gz

echo "sample is: $sample"
echo "bam is: $bam"

# Load modules
#module load gatk/4.2.1.0
module load gatk/4.1.4.1
module load samtools/1.9

# index bam files
samtools index ${bam}

# make gvcf per sample
gatk --java-options "-Xmx8g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
HaplotypeCaller \
-R ${ref} \
-I ${bam} \
-O ${vcf} \
-ERC GVCF
#--native-pair-hmm-threads ${NCPUS} # can maybe try using this to multithread and speed things up
```

