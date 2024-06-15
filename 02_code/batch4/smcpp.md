# SMC++ for estimating size history of populations

## Batch 4 data

```bash
module load bsub.py/0.42.1
module load smcpp/1.15.3-c1
module load common-apps/htslib/1.9.229
module load bcftools/1.14--h88f3f91_0
module load vcftools/0.1.16-c4

cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/04_VARIANTS/FILTER2/NO_OUTGROUPS/SMCPP

# prep VCF file
bgzip -c /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/04_VARIANTS/FILTER2/NO_OUTGROUPS/FINAL_SETS/nuclear_samples3x_missing0.9.chr1to4.recode.vcf > nuclear_samples3x_missing0.9.chr1to4.recode.vcf.gz
tabix -p vcf nuclear_samples3x_missing0.9.chr1to4.recode.vcf.gz
VCF=/lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/04_VARIANTS/FILTER2/NO_OUTGROUPS/SMCPP/nuclear_samples3x_missing0.9.chr1to4.recode.vcf.gz

## There may be missing genotypes designated as '.', whereas smc++ expects it to be './.' - so see if there are any lines like this:
zcat ${VCF} | \
  bcftools view | \
  grep -v '^#' | \
  awk '{
    found = 0;
    for (i=10; i<=90; i++) {
      split($i, gt, ":");
      if (gt[1] == "." && length(gt[1]) == 1 && found == 0) {
        print "Found \".\" in line", NR, "in column", i;
        found = 1;
      }
    }
  }'
# There are.

# view the actual lines
zcat ${VCF} | \
  bcftools view | \
  grep -v '^#' | \
  awk '{
    for (i=10; i<=90; i++) {
      split($i, gt, ":");
      if (gt[1] == "." && length(gt[1]) == 1) {
        print "Found \".\" in line", NR ":", $0;
        break;
      }
    }
  }'

## Replace them as './.'
zcat ${VCF} | \
  awk '
  BEGIN {OFS="\t"}
  /^#/ {print; next}
  {
    for (i=10; i<=NF; i++) {
      split($i, gt, ":");
      if (gt[1] == "." && length(gt[1]) == 1) {
        gt[1] = "./."
      }
      $i = gt[1]
      for (j=2; j<=length(gt); j++) {
        $i = $i ":" gt[j]
      }
    }
    print
  }' | bgzip > smcpp.vcf.gz

# check if it worked
zcat smcpp.vcf.gz | \
  bcftools view | \
  grep -v '^#' | \
  awk '{
    found = 0;
    for (i=10; i<=90; i++) {
      split($i, gt, ":");
      if (gt[1] == "." && length(gt[1]) == 1 && found == 0) {
        print "Found \".\" in line", NR, "in column", i;
        found = 1;
      }
    }
  }'

zcat smcpp.vcf.gz | \
  bcftools view | \
  grep -v '^#' | \
  awk '{
    for (i=10; i<=90; i++) {
      split($i, gt, ":");
      if (gt[1] == "." && length(gt[1]) == 1) {
        print "Found \".\" in line", NR ":", $0;
        break;
      }
    }
  }'
# Yep, none anymore! Continue with smc++ analysis.
tabix -p vcf smcpp.vcf.gz
```

# Plot per population 

bsub.py --queue long 10 run_smcpp "run_smcpp.sh"

```bash
# Load modules
module load smcpp/1.15.3-c1
module load common-apps/htslib/1.9.229
module load vcftools/0.1.16-c4

cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/04_VARIANTS/FILTER2/NO_OUTGROUPS/SMCPP
mkdir DATA

# Get sample names for each population
ASIA=$(awk -F'_' '$1 == "MYS" || $1 == "THA" {print}' nuclear_samplelist.keep | paste -sd ',')
AUS=$(awk -F'_' '$1 == "AUS" {print}' nuclear_samplelist.keep | paste -sd ',')
CAM=$(awk -F'_' '$1 == "PAN" || $1 == "CRI" {print}' nuclear_samplelist.keep | paste -sd ',')
EUR=$(awk -F'_' '$1 == "GRC" || $1 == "ITA" || $1 == "ROU" {print}' nuclear_samplelist.keep | paste -sd ',')
USA=$(awk -F'_' '$1 == "USA" {print}' nuclear_samplelist.keep | paste -sd ',')
declare -A populations
populations=(
  [ASIA]=$ASIA
  [AUS]=$AUS
  [CAM]=$CAM
  [EUR]=$EUR
  [USA]=$USA
)
# nuclear_samplelist.keep - the samples outputted in the final vcf. Removed repeat samples.

# ASIA
vcftools --gzvcf smcpp.vcf.gz \
--indv MYS_SEL_AD_001 \
--indv THA_BKK_AD_001 \
--indv THA_BKK_AD_002 \
--indv THA_BKK_AD_003 \
--indv THA_BKK_AD_004 \
--indv THA_BKK_AD_005 \
--indv THA_BKK_AD_006 \
--max-missing 1 --recode --out ASIA
bgzip -f ASIA.recode.vcf
tabix ASIA.recode.vcf.gz

# Convert VCF to the SMC++ input format with vcf2smc
for chr in {1..4}; do
  smc++ vcf2smc ASIA.recode.vcf.gz DATA/ASIA.chr${chr}.smc.gz dirofilaria_immitis_chr${chr} ASIA:${ASIA};
done
# each call to vcf2smc processes a single contig. VCFs containing multiple contigs should be processed via multiple separate runs.

# Fit the model using Estimate
## # mutation rate (C.elegans) = 2.7e-9 (https://doi.org/10.1073/pnas.0904895106)
smc++ estimate --timepoints 1 1000000 -o ASIA/ 2.7e-9 DATA/ASIA.*.smc.gz

# plot
## Use generation time of D. immitis ~ 5 yrs
smc++ plot -g 5 -c SMCPP_ASIA.png ASIA/model.final.json


# AUS
vcftools --gzvcf smcpp.vcf.gz \
--indv AUS_BNE_AD_001 \
--indv AUS_BNE_AD_002 \
--indv AUS_BNE_AD_003 \
--indv AUS_BNE_AD_004 \
--indv AUS_BNE_AD_006 \
--indv AUS_BNE_AD_008 \
--indv AUS_BNE_AD_009 \
--indv AUS_CNS_AD_001 \
--indv AUS_CNS_AD_002 \
--indv AUS_LHR_AD_001 \
--indv AUS_ROK_AD_001 \
--indv AUS_SYD_AD_001 \
--indv AUS_SYD_AD_002 \
--indv AUS_SYD_AD_003 \
--indv AUS_SYD_AD_004 \
--indv AUS_SYD_AD_005 \
--indv AUS_SYD_AD_006 \
--indv AUS_SYD_AD_007 \
--indv AUS_SYD_AD_008 \
--indv AUS_SYD_AD_009 \
--indv AUS_SYD_AD_010 \
--indv AUS_SYD_AD_011 \
--indv AUS_SYD_AD_012 \
--indv AUS_SYD_AD_013 \
--indv AUS_SYD_AD_014 \
--indv AUS_SYD_AD_015 \
--indv AUS_TVS_AD_001 \
--indv AUS_TVS_AD_002 \
--indv AUS_TVS_AD_003 \
--indv AUS_TVS_AD_004 \
--indv AUS_TVS_AD_005 \
--indv AUS_TVS_AD_006 \
--indv AUS_TVS_AD_007 \
--indv AUS_TVS_AD_008 \
--indv AUS_TVS_AD_010 \
--indv AUS_TVS_AD_011 \
--indv AUS_TVS_AD_012 \
--indv AUS_TVS_AD_013 \
--indv AUS_TVS_AD_014 \
--indv AUS_TVS_AD_015 \
--indv AUS_TVS_AD_016 \
--indv AUS_TVS_AD_017 \
--indv AUS_TVS_AD_018 \
--indv AUS_TVS_AD_019 \
--max-missing 1 --recode --out AUS
bgzip -f AUS.recode.vcf
tabix AUS.recode.vcf.gz

# Convert VCF to the SMC++ input format with vcf2smc
for chr in {1..4}; do
  smc++ vcf2smc AUS.recode.vcf.gz DATA/AUS.chr${chr}.smc.gz dirofilaria_immitis_chr${chr} AUS:${AUS};
done
# each call to vcf2smc processes a single contig. VCFs containing multiple contigs should be processed via multiple separate runs.

# Fit the model using Estimate
## # mutation rate (C.elegans) = 2.7e-9 (https://doi.org/10.1073/pnas.0904895106)
smc++ estimate --timepoints 1 1000000 -o AUS/ 2.7e-9 DATA/AUS.*.smc.gz

# plot
## Use generation time of D. immitis ~ 5 yrs
smc++ plot -g 5 -c SMCPP_AUS.png AUS/model.final.json



# CAM
vcftools --gzvcf smcpp.vcf.gz \
--indv CRI_SJO_AD_001 \
--indv PAN_BOC_AD_001 \
--indv PAN_BOC_AD_002 \
--indv PAN_BOC_AD_003 \
--indv PAN_BOC_AD_004 \
--indv PAN_BOC_AD_005 \
--indv PAN_PUE_AD_001 \
--indv PAN_PUE_AD_002 \
--indv PAN_PUE_AD_003 \
--indv PAN_PUE_AD_004 \
--indv PAN_PUE_AD_005 \
--indv PAN_PUE_AD_006 \
--indv PAN_SLO_AD_001 \
--indv PAN_SLO_AD_002 \
--indv PAN_SLO_AD_003 \
--max-missing 1 --recode --out CAM
bgzip -f CAM.recode.vcf
tabix CAM.recode.vcf.gz

# Convert VCF to the SMC++ input format with vcf2smc
for chr in {1..4}; do
  smc++ vcf2smc CAM.recode.vcf.gz DATA/CAM.chr${chr}.smc.gz dirofilaria_immitis_chr${chr} CAM:${CAM};
done
# each call to vcf2smc processes a single contig. VCFs containing multiple contigs should be processed via multiple separate runs.

# Fit the model using Estimate
## # mutation rate (C.elegans) = 2.7e-9 (https://doi.org/10.1073/pnas.0904895106)
smc++ estimate --timepoints 1 1000000 -o CAM/ 2.7e-9 DATA/CAM.*.smc.gz

# plot
## Use generation time of D. immitis ~ 5 yrs
smc++ plot -g 5 -c SMCPP_CAM.png CAM/model.final.json




# EUR
vcftools --gzvcf smcpp.vcf.gz \
--indv CRI_SJO_AD_001 \
--indv GRC_XAN_AD_001 \
--indv GRC_XAN_AD_002 \
--indv GRC_XAN_AD_003 \
--indv GRC_XAN_AD_004 \
--indv GRC_XAN_AD_005 \
--indv GRC_XAN_AD_006 \
--indv GRC_XAN_AD_007 \
--indv GRC_XAN_AD_008 \
--indv GRC_XAN_AD_009 \
--indv GRC_XAN_AD_010 \
--indv GRC_XAN_AD_011 \
--indv GRC_XAN_AD_012 \
--indv ITA_NEA_AD_001 \
--indv ITA_NEA_AD_002 \
--indv ITA_NEA_AD_003 \
--indv ITA_PAV_AD_001 \
--indv ROU_BUC_AD_001 \
--indv ROU_COM_AD_001 \
--indv ROU_GIU_AD_001 \
--max-missing 1 --recode --out EUR
bgzip -f EUR.recode.vcf
tabix EUR.recode.vcf.gz

# Convert VCF to the SMC++ input format with vcf2smc
for chr in {1..4}; do
  smc++ vcf2smc EUR.recode.vcf.gz DATA/EUR.chr${chr}.smc.gz dirofilaria_immitis_chr${chr} EUR:${EUR};
done
# each call to vcf2smc processes a single contig. VCFs containing multiple contigs should be processed via multiple separate runs.

# Fit the model using Estimate
## # mutation rate (C.elegans) = 2.7e-9 (https://doi.org/10.1073/pnas.0904895106)
smc++ estimate --timepoints 1 1000000 -o EUR/ 2.7e-9 DATA/EUR.*.smc.gz

# plot
## Use generation time of D. immitis ~ 5 yrs
smc++ plot -g 5 -c SMCPP_EUR.png EUR/model.final.json





# USA
vcftools --gzvcf smcpp.vcf.gz \
--indv USA_FLO_AD_001 \
--indv USA_FLO_AD_002 \
--indv USA_FLO_AD_003 \
--indv USA_FLO_AD_004 \
--indv USA_FLO_AD_005 \
--indv USA_FLO_AD_006 \
--indv USA_FLO_AD_007 \
--indv USA_FLO_AD_008 \
--indv USA_FLO_AD_009 \
--indv USA_FLO_AD_010 \
--indv USA_FLO_AD_011 \
--indv USA_FLO_AD_012 \
--indv USA_FLO_AD_013 \
--indv USA_FLO_AD_014 \
--indv USA_FLO_AD_015 \
--indv USA_FLO_AD_016 \
--indv USA_GEO_AD_001 \
--indv USA_ILL_AD_001 \
--indv USA_ILL_AD_002 \
--indv USA_LOU_AD_001 \
--indv USA_LOU_AD_002 \
--indv USA_MIS_AD_001 \
--indv USA_MIS_AD_002 \
--indv USA_TEX_AD_001 \
--indv USA_TEX_AD_002 \
--indv USA_TEX_AD_003 \
--indv USA_TEX_AD_004 \
--indv USA_TEX_AD_005 \
--indv USA_TEX_AD_006 \
--indv USA_TEX_AD_007 \
--indv USA_TEX_AD_008 \
--indv USA_TEX_AD_009 \
--indv USA_TEX_AD_010 \
--indv USA_TEX_AD_011 \
--indv USA_TEX_AD_012 \
--indv USA_TEX_AD_013 \
--indv USA_TEX_AD_014 \
--indv USA_TEX_AD_015 \
--max-missing 1 --recode --out USA
bgzip -f USA.recode.vcf
tabix USA.recode.vcf.gz

# Convert VCF to the SMC++ input format with vcf2smc
for chr in {1..4}; do
  smc++ vcf2smc USA.recode.vcf.gz DATA/USA.chr${chr}.smc.gz dirofilaria_immitis_chr${chr} USA:${USA};
done
# each call to vcf2smc processes a single contig. VCFs containing multiple contigs should be processed via multiple separate runs.

# Fit the model using Estimate
## # mutation rate (C.elegans) = 2.7e-9 (https://doi.org/10.1073/pnas.0904895106)
smc++ estimate --timepoints 1 1000000 -o USA/ 2.7e-9 DATA/USA.*.smc.gz

# plot
## Use generation time of D. immitis ~ 5 yrs
smc++ plot -g 5 -c SMCPP_USA.png USA/model.final.json
```



########################################################################################
```bash
# Clean split models

The split command fits two-population clean split models (assumes no ongoing gene flow between the two populations after they diverged)

```bash
# Create datasets containing the joint frequency spectrum for both populations. Run split.
## ASIA & AUS
for chr in {1..4}; do
mkdir split/ASIA_AUS/chr${chr}
smc++ vcf2smc ${VCF_MOD} split/ASIA_AUS/chr${chr}/chr${chr}_ASIA_AUS.smc.gz dirofilaria_immitis_chr${chr} ASIA:${ASIA} AUS:${AUS}
smc++ vcf2smc ${VCF_MOD} split/ASIA_AUS/chr${chr}/chr${chr}_AUS_ASIA.smc.gz dirofilaria_immitis_chr${chr} AUS:${AUS} ASIA:${ASIA}
smc++ split --timepoints 0 100000 -o split/ASIA_AUS/chr${chr} chr${chr}_ASIA_model/model.final.json chr${chr}_AUS_model.json/model.final.json split/ASIA_AUS/chr${chr}/*.smc.gz
smc++ plot split/ASIA_AUS/chr${chr}/chr${chr}_ASIA_AUS.pdf split/ASIA_AUS/chr${chr}/model.final.json;
done

## ASIA & CAM
for chr in {1..4}; do
mkdir split/ASIA_CAM/chr${chr}
smc++ vcf2smc ${VCF_MOD} split/ASIA_CAM/chr${chr}/chr${chr}_ASIA_CAM.smc.gz dirofilaria_immitis_chr${chr} ASIA:${ASIA} CAM:${CAM}
smc++ vcf2smc ${VCF_MOD} split/ASIA_CAM/chr${chr}/chr${chr}_CAM_ASIA.smc.gz dirofilaria_immitis_chr${chr} CAM:${CAM} ASIA:${ASIA}
smc++ split -o split/ASIA_CAM/chr${chr} chr${chr}_ASIA_model.json/model.final.json chr${chr}_CAM_model.json/model.final.json split/ASIA_CAM/chr${chr}/*.smc.gz
smc++ plot split/ASIA_CAM/chr${chr}/chr${chr}_ASIA_CAM.pdf split/ASIA_CAM/chr${chr}/model.final.json;
done

## ASIA & EUR

## ASIA & USA

## AUS & CAM

## AUS & EUR

## AUS & USA

## CAM & EUR

## CAM & USA

## EUR & USA





```