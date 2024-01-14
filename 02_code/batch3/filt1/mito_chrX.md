## SNPs in Mitochondria and ChrX

## Select chrX

```bash
#!/bin/bash

# PBS directives 
#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N snps_chrX
#PBS -l select=1:ncpus=1:mem=4GB
#PBS -l walltime=00:10:00
#PBS -m e
#PBS -q defaultQ
#PBS -o snps_chrX.txt

# qsub ../snps_chrX.pbs

# load gatk
module load vcftools/0.1.14
module load bcftools/1.11

cd /scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping/filter

vcftools --vcf FINAL_SETS/nuclear_samples3x_missing0.9.recode.vcf \
--chr dirofilaria_immitis_chrX \
--recode --out FINAL_SETS/nuclear_samples3x_missing0.9.chrX
#After filtering, kept 81 out of 81 Individuals
#Outputting VCF file...
#After filtering, kept 64158 out of a possible 240259 Sites

vcftools --vcf FINAL_SETS/nuclear_samples3x_missing0.9.chrX.recode.vcf --remove-indels
#After filtering, kept 81 out of 81 Individuals
#After filtering, kept 64158 out of a possible 64158 Sites

vcftools --vcf FINAL_SETS/nuclear_samples3x_missing0.9.chrX.recode.vcf --keep-only-indels
#After filtering, kept 81 out of 81 Individuals
#After filtering, kept 0 out of a possible 64158 Sites

#Rename samples to more meaningful names:

cd FINAL_SETS

bcftools reheader -s rename.txt nuclear_samples3x_missing0.9.chrX.recode.vcf -o nuclear_samples3x_missing0.9.chrX.recode.RENAMED.vcf
```

Check that it worked:

```bash
# List all sample names in the original VCF
bcftools query -l /scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping/filter/FINAL_SETS/nuclear_samples3x_missing0.9.chrX.recode.RENAMED.vcf
```

Yep all the sample names are changed.



## Rename mito

```bash
#!/bin/bash

# PBS directives 
#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N vcf_rename_mt
#PBS -l select=1:ncpus=2:mem=20GB
#PBS -l walltime=00:20:00
#PBS -m abe
#PBS -q defaultQ
#PBS -o vcf_rename_mt.txt
#PBS -M rosemonde.power@sydney.edu.au

# qsub ../vcf_rename_mt.pbs

# load modules
module load bcftools/1.11

cd /scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping/filter/FINAL_SETS

bcftools reheader -s rename.txt mito_samples3x_missing1.recode.vcf -o mito_samples3x_missing1.recode.RENAMED.vcf
```

Check that it worked:

```bash
# List all sample names in the original VCF
bcftools query -l /scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping/filter/FINAL_SETS/mito_samples3x_missing1.recode.RENAMED.vcf
```

Yep all the sample names are changed.

## PCA

