# SNP concordance between replicate samples

Replicate samples are:
- AUS_BNE_AD_003 & AUS_BNE_AD_003_R
- AUS_SYD_AD_005 & AUS_SYD_AD_005_R
- AUS_SYD_AD_008 & AUS_SYD_AD_008_R
- PAN_PUE_AD_004 & PAN_PUE_AD_004_R

## Get VCF files ready

Get 1 VCF file containing the first replicate of each sample. Then get another VCF file containing the second replicates.

```bash
module load PaM/environment
module load bsub.py/0.42.1
module load vcftools/0.1.16-c4
module load bcftools/1.17--h3cc50cf_1

cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/05_ANALYSIS
mkdir SNP_CONCORDANCE
cd SNP_CONCORDANCE

# replicate 1
vcftools --vcf /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/04_VARIANTS/FILTER1/NO_OUTGROUPS/FINAL_SETS/nuclear_samples3x_missing0.9.chr1to4.recode.vcf \
--indv AUS_BNE_AD_003 \
--indv AUS_SYD_AD_005 \
--indv AUS_SYD_AD_008 \
--indv PAN_PUE_AD_004 \
--recode --out nuclear_samples3x_missing0.9.chr1to4.REP1
# After filtering, kept 4 out of 128 Individuals
# After filtering, kept 196550 out of a possible 196550 Sites

# check sample names
grep "#CHROM" nuclear_samples3x_missing0.9.chr1to4.REP1.recode.vcf

# replicate 2
vcftools --vcf /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/04_VARIANTS/FILTER1/NO_OUTGROUPS/FINAL_SETS/nuclear_samples3x_missing0.9.chr1to4.recode.vcf \
--indv AUS_BNE_AD_003_R \
--indv AUS_SYD_AD_005_R \
--indv AUS_SYD_AD_008_R \
--indv PAN_PUE_AD_004_R \
--recode --out nuclear_samples3x_missing0.9.chr1to4.REP2
# After filtering, kept 4 out of 128 Individuals
# After filtering, kept 196550 out of a possible 196550 Sites

# check sample names
grep "#CHROM" nuclear_samples3x_missing0.9.chr1to4.REP2.recode.vcf



# calculate discordance between replicate samples
vcftools --vcf nuclear_samples3x_missing0.9.chr1to4.REP1.recode.vcf \
--diff nuclear_samples3x_missing0.9.chr1to4.REP2.recode.vcf \
--diff-indv-map names.txt \
--diff-indv-discordance
# This compares the original input file (REP1) to another specified VCF (REP2). 
# --diff-indv-map allows me to specify an indivdual's name in file 1 and file 2
# --diff-indv-discordance calculates discordance on a per-individual basis.
# names.txt:
'
AUS_BNE_AD_003_R	AUS_BNE_AD_003
AUS_SYD_AD_005_R	AUS_SYD_AD_005
AUS_SYD_AD_008_R	AUS_SYD_AD_008
PAN_PUE_AD_004_R	PAN_PUE_AD_004
'

# output:
INDV	N_COMMON_CALLED	N_DISCORD	DISCORDANCE
AUS_BNE_AD_003	196354	295	0.00150239 = < 0.2% discordance
AUS_SYD_AD_005	196058	977	0.00498322 = < 0.5% discordance
AUS_SYD_AD_008	196194	352	0.00179414 = < 0.2% discordance
PAN_PUE_AD_004	196434	324	0.00164941 = < 0.2% discordance

```