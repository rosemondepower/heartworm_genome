
# Convert vcf file to eigenstrat format

ADMIXTOOLS requires eigenstrat input files. To convert a vcf file to eigenstrat, we will use a conversion script written by Joana. That script just requires the name of the vcf file (NB - the script can use uncompressed vcfs too). 

```bash
module load bsub.py/0.42.1
module load compilers/gcc/12.2
module load vcftools/0.1.16-c4
module load perl-velvetoptimiser/2.2.6

cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/BATCH3/02_VARIANTS/ADMIXTOOLS

bash convertVCFtoEigenstrat.sh nuclear_samples3x_missing0.9.chr1to4.recode.RENAMED.vcf --renameScaff

# convertVCFtoEigenstrat.sh downloaded from: https://github.com/joanam/scripts/blob/master/convertVCFtoEigenstrat.sh

```