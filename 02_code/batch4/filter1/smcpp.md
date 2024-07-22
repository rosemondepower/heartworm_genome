# SMC++ for estimating size history of populations

## Batch 4 data

## Prep VCF file for smcpp

```bash
module load bsub.py/0.42.1
module load smcpp/1.15.3-c1
module load common-apps/htslib/1.9.229
module load bcftools/1.14--h88f3f91_0
module load vcftools/0.1.16-c4

cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/05_ANALYSIS/SMCPP

# prep VCF file
ln -s /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/04_VARIANTS/FILTER1/NO_OUTGROUPS/FINAL_SETS/nuclear_samples3x_missing0.9.chr1to4.recode.vcf
bgzip -c nuclear_samples3x_missing0.9.chr1to4.recode.vcf > nuclear_samples3x_missing0.9.chr1to4.recode.vcf.gz
tabix -p vcf nuclear_samples3x_missing0.9.chr1to4.recode.vcf.gz
VCF=/lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/05_ANALYSIS/SMCPP/nuclear_samples3x_missing0.9.chr1to4.recode.vcf.gz

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

## Plot per population 

### timepoint = 1 to 1 million years, generation time = 5 years

bsub.py --queue long 10 run_smcpp_t1m_g5 "run_smcpp_t1m_g5.sh"

```bash
# Load modules
module load smcpp/1.15.3-c1
module load common-apps/htslib/1.9.229
module load vcftools/0.1.16-c4

cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/05_ANALYSIS/SMCPP
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
--indv AUS_SYD_AD_017 \
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
Put this output into a folder called 't1m/g5'. Now try out a few other parameters.

## 21/7/24 up to here



### timepoint = 1 to 1 million years, generation time = 1 year

```bash
module load bsub.py/0.42.1
cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/04_VARIANTS/FILTER2/NO_OUTGROUPS/SMCPP/t1m/g1
bsub.py --queue long 10 run_smcpp_t1m_g1 "run_smcpp_t1m_g1.sh"
```

```bash
# Load modules
module load smcpp/1.15.3-c1

cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/04_VARIANTS/FILTER2/NO_OUTGROUPS/SMCPP/t1m/g1


# ASIA
# plot
## Use generation time of D. immitis ~ 1 yrs
smc++ plot -g 1 -c SMCPP_ASIA_t1m_g1.png ../g5/ASIA/model.final.json

# AUS
# plot
## Use generation time of D. immitis ~ 1 yrs
smc++ plot -g 1 -c SMCPP_AUS_t1m_g1.png ../g5/AUS/model.final.json

# CAM
# plot
## Use generation time of D. immitis ~ 1 yrs
smc++ plot -g 1 -c SMCPP_CAM_t1m_g1.png ../g5/CAM/model.final.json

# EUR
# plot
## Use generation time of D. immitis ~ 1 yrs
smc++ plot -g 1 -c SMCPP_EUR_t1m_g1.png ../g5/EUR/model.final.json

# USA
# plot
## Use generation time of D. immitis ~ 1 yrs
smc++ plot -g 1 -c SMCPP_USA_t1m_g1.png ../g5/USA/model.final.json
```

### timepoint = 1 to 1 million years, generation time = 2.5 years

```bash
cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/04_VARIANTS/FILTER2/NO_OUTGROUPS/SMCPP/t1m/g2.5
bsub.py --queue long 10 run_smcpp_t1m_g2.5 "run_smcpp_t1m_g2.5.sh"
```

```bash
# Load modules
module load smcpp/1.15.3-c1

cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/04_VARIANTS/FILTER2/NO_OUTGROUPS/SMCPP/t1m/g2.5

# ASIA
# plot
## Use generation time of D. immitis ~ 2.5 yrs
smc++ plot -g 2.5 -c SMCPP_ASIA_t1m_g2.5.png ../g5/ASIA/model.final.json

# AUS
# plot
## Use generation time of D. immitis ~ 2.5 yrs
smc++ plot -g 2.5 -c SMCPP_AUS_t1m_g2.5.png ../g5/AUS/model.final.json

# CAM
# plot
## Use generation time of D. immitis ~ 2.5 yrs
smc++ plot -g 2.5 -c SMCPP_CAM_t1m_g2.5.png ../g5/CAM/model.final.json

# EUR
# plot
## Use generation time of D. immitis ~ 2.5 yrs
smc++ plot -g 2.5 -c SMCPP_EUR_t1m_g2.5.png ../g5/EUR/model.final.json

# USA
# plot
## Use generation time of D. immitis ~ 2.5 yrs
smc++ plot -g 2.5 -c SMCPP_USA_t1m_g2.5.png ../g5/USA/model.final.json
```

### timepoint = 1 to 1.5 million years, generation time = 1, 2.5, 5 years

```bash
cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/04_VARIANTS/FILTER2/NO_OUTGROUPS/SMCPP/t1.5m
bsub.py --queue long 10 run_smcpp_t1.5m "run_smcpp_t1.5m.sh"
```

```bash
# Load modules
module load smcpp/1.15.3-c1

cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/04_VARIANTS/FILTER2/NO_OUTGROUPS/SMCPP/t1.5m
mkdir g1
mkdir g2.5
mkdir g5

# nuclear_samplelist.keep - the samples outputted in the final vcf. Removed repeat samples.

# ASIA
# Fit the model using Estimate
## # mutation rate (C.elegans) = 2.7e-9 (https://doi.org/10.1073/pnas.0904895106)
smc++ estimate --timepoints 1 1500000 -o ASIA/ 2.7e-9 ../t1m/g5/DATA/ASIA.*.smc.gz
# plot
## Use generation time of D. immitis ~ 1, 2.5 & 5 yrs
smc++ plot -g 1 -c g1/SMCPP_ASIA_t1.5m_g1.png ASIA/model.final.json
smc++ plot -g 2.5 -c g2.5/SMCPP_ASIA_t1.5m_g2.5.png ASIA/model.final.json
smc++ plot -g 5 -c g5/SMCPP_ASIA_t1.5m_g5.png ASIA/model.final.json

# AUS
# Fit the model using Estimate
## # mutation rate (C.elegans) = 2.7e-9 (https://doi.org/10.1073/pnas.0904895106)
smc++ estimate --timepoints 1 1500000 -o AUS/ 2.7e-9 ../t1m/g5/DATA/AUS.*.smc.gz
# plot
## Use generation time of D. immitis ~ 5 yrs
smc++ plot -g 1 -c g1/SMCPP_AUS_t1.5m_g1.png AUS/model.final.json
smc++ plot -g 2.5 -c g2.5/SMCPP_AUS_t1.5m_g2.5.png AUS/model.final.json
smc++ plot -g 5 -c g5/SMCPP_AUS_t1.5m_g5.png AUS/model.final.json

# CAM
# Fit the model using Estimate
## # mutation rate (C.elegans) = 2.7e-9 (https://doi.org/10.1073/pnas.0904895106)
smc++ estimate --timepoints 1 1500000 -o CAM/ 2.7e-9 ../t1m/g5/DATA/CAM.*.smc.gz
# plot
## Use generation time of D. immitis ~ 5 yrs
smc++ plot -g 1 -c g1/SMCPP_CAM_t1.5m_g1.png CAM/model.final.json
smc++ plot -g 2.5 -c g2.5/SMCPP_CAM_t1.5m_g2.5.png CAM/model.final.json
smc++ plot -g 5 -c g5/SMCPP_CAM_t1.5m_g5.png CAM/model.final.json

# EUR
# Fit the model using Estimate
## # mutation rate (C.elegans) = 2.7e-9 (https://doi.org/10.1073/pnas.0904895106)
smc++ estimate --timepoints 1 1500000 -o EUR/ 2.7e-9 ../t1m/g5/DATA/EUR.*.smc.gz
# plot
## Use generation time of D. immitis ~ 5 yrs
smc++ plot -g 1 -c g1/SMCPP_EUR_t1.5m_g1.png EUR/model.final.json
smc++ plot -g 2.5 -c g2.5/SMCPP_EUR_t1.5m_g2.5.png EUR/model.final.json
smc++ plot -g 5 -c g5/SMCPP_EUR_t1.5m_g5.png EUR/model.final.json

# USA
# Fit the model using Estimate
## # mutation rate (C.elegans) = 2.7e-9 (https://doi.org/10.1073/pnas.0904895106)
smc++ estimate --timepoints 1 1500000 -o USA/ 2.7e-9 ../t1m/g5/DATA/USA.*.smc.gz
# plot
## Use generation time of D. immitis ~ 5 yrs
smc++ plot -g 1 -c g1/SMCPP_USA_t1.5m_g1.png USA/model.final.json
smc++ plot -g 2.5 -c g2.5/SMCPP_USA_t1.5m_g2.5.png USA/model.final.json
smc++ plot -g 5 -c g5/SMCPP_USA_t1.5m_g5.png USA/model.final.json
```


### timepoint = 1 to 1 million years, generation time = 2.5 years, -c 1kbp

The -c parameter will treat runs of homozygosity longer than -c bp as missing.

cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/04_VARIANTS/FILTER2/NO_OUTGROUPS/SMCPP/c1kbp
bsub.py --queue long 10 run_smcpp_c1kbp "run_smcpp_c1kbp.sh"

```bash
# Load modules
module load smcpp/1.15.3-c1

cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/04_VARIANTS/FILTER2/NO_OUTGROUPS/SMCPP/c1kbp
mkdir DATA

# Get sample names for each population
ASIA=$(awk -F'_' '$1 == "MYS" || $1 == "THA" {print}' ../nuclear_samplelist.keep | paste -sd ',')
AUS=$(awk -F'_' '$1 == "AUS" {print}' ../nuclear_samplelist.keep | paste -sd ',')
CAM=$(awk -F'_' '$1 == "PAN" || $1 == "CRI" {print}' ../nuclear_samplelist.keep | paste -sd ',')
EUR=$(awk -F'_' '$1 == "GRC" || $1 == "ITA" || $1 == "ROU" {print}' ../nuclear_samplelist.keep | paste -sd ',')
USA=$(awk -F'_' '$1 == "USA" {print}' ../nuclear_samplelist.keep | paste -sd ',')
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
# Convert VCF to the SMC++ input format with vcf2smc
for chr in {1..4}; do
  smc++ vcf2smc -c 1000 ../ASIA.recode.vcf.gz DATA/ASIA.chr${chr}.smc.gz dirofilaria_immitis_chr${chr} ASIA:${ASIA};
done
# each call to vcf2smc processes a single contig. VCFs containing multiple contigs should be processed via multiple separate runs.

# Fit the model using Estimate
## # mutation rate (C.elegans) = 2.7e-9 (https://doi.org/10.1073/pnas.0904895106)
smc++ estimate --timepoints 1 1000000 -o ASIA/ 2.7e-9 DATA/ASIA.*.smc.gz

# plot
## Use generation time of D. immitis ~ 2.5 yrs
smc++ plot -g 2.5 -c SMCPP_ASIA_c1kbp.png ASIA/model.final.json


# AUS
# Convert VCF to the SMC++ input format with vcf2smc
for chr in {1..4}; do
  smc++ vcf2smc -c 1000 ../AUS.recode.vcf.gz DATA/AUS.chr${chr}.smc.gz dirofilaria_immitis_chr${chr} AUS:${AUS};
done
# each call to vcf2smc processes a single contig. VCFs containing multiple contigs should be processed via multiple separate runs.

# Fit the model using Estimate
## # mutation rate (C.elegans) = 2.7e-9 (https://doi.org/10.1073/pnas.0904895106)
smc++ estimate --timepoints 1 1000000 -o AUS/ 2.7e-9 DATA/AUS.*.smc.gz

# plot
## Use generation time of D. immitis ~ 2.5 yrs
smc++ plot -g 2.5 -c SMCPP_AUS_c1kbp.png AUS/model.final.json



# CAM
# Convert VCF to the SMC++ input format with vcf2smc
for chr in {1..4}; do
  smc++ vcf2smc -c 1000 ../CAM.recode.vcf.gz DATA/CAM.chr${chr}.smc.gz dirofilaria_immitis_chr${chr} CAM:${CAM};
done
# each call to vcf2smc processes a single contig. VCFs containing multiple contigs should be processed via multiple separate runs.

# Fit the model using Estimate
## # mutation rate (C.elegans) = 2.7e-9 (https://doi.org/10.1073/pnas.0904895106)
smc++ estimate --timepoints 1 1000000 -o CAM/ 2.7e-9 DATA/CAM.*.smc.gz

# plot
## Use generation time of D. immitis ~ 2.5 yrs
smc++ plot -g 2.5 -c SMCPP_CAM_c1kbp.png CAM/model.final.json




# EUR
# Convert VCF to the SMC++ input format with vcf2smc
for chr in {1..4}; do
  smc++ vcf2smc -c 1000 ../EUR.recode.vcf.gz DATA/EUR.chr${chr}.smc.gz dirofilaria_immitis_chr${chr} EUR:${EUR};
done
# each call to vcf2smc processes a single contig. VCFs containing multiple contigs should be processed via multiple separate runs.

# Fit the model using Estimate
## # mutation rate (C.elegans) = 2.7e-9 (https://doi.org/10.1073/pnas.0904895106)
smc++ estimate --timepoints 1 1000000 -o EUR/ 2.7e-9 DATA/EUR.*.smc.gz

# plot
## Use generation time of D. immitis ~ 2.5 yrs
smc++ plot -g 2.5 -c SMCPP_EUR_c1kbp.png EUR/model.final.json





# USA
# Convert VCF to the SMC++ input format with vcf2smc
for chr in {1..4}; do
  smc++ vcf2smc -c 1000 ../USA.recode.vcf.gz DATA/USA.chr${chr}.smc.gz dirofilaria_immitis_chr${chr} USA:${USA};
done
# each call to vcf2smc processes a single contig. VCFs containing multiple contigs should be processed via multiple separate runs.

# Fit the model using Estimate
## # mutation rate (C.elegans) = 2.7e-9 (https://doi.org/10.1073/pnas.0904895106)
smc++ estimate --timepoints 1 1000000 -o USA/ 2.7e-9 DATA/USA.*.smc.gz

# plot
## Use generation time of D. immitis ~ 2.5 yrs
smc++ plot -g 2.5 -c SMCPP_USA_c1kbp.png USA/model.final.json
```


## Split

The split command fits two-population clean split models (assumes no ongoing gene flow between the two populations after they diverged).

## Timepoint 1 to 1,000,000, generation time ~2.5 years

cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/04_VARIANTS/FILTER2/NO_OUTGROUPS/SMCPP/split/g2.5/run1
bsub.py --queue long 10 run_smcpp_split_g2.5_run1 "run_smcpp_split_g2.5_run1.sh"

```bash
# Load modules
module load smcpp/1.15.3-c1
module load common-apps/htslib/1.9.229
module load vcftools/0.1.16-c4

cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/04_VARIANTS/FILTER2/NO_OUTGROUPS/SMCPP/split/g2.5/run1

# Create necessary directories
mkdir DATA
mkdir DATA/ASIA_AUS DATA/ASIA_CAM DATA/ASIA_EUR DATA/ASIA_USA
mkdir DATA/AUS_CAM DATA/AUS_EUR DATA/AUS_USA
mkdir DATA/CAM_EUR DATA/CAM_USA DATA/EUR_USA

# Get sample names for each population
ASIA=$(awk -F'_' '$1 == "MYS" || $1 == "THA" {print $0}' ../../../nuclear_samplelist.keep | paste -sd ',')
AUS=$(awk -F'_' '$1 == "AUS" {print $0}' ../../../nuclear_samplelist.keep | paste -sd ',')
CAM=$(awk -F'_' '$1 == "PAN" || $1 == "CRI" {print $0}' ../../../nuclear_samplelist.keep | paste -sd ',')
EUR=$(awk -F'_' '$1 == "GRC" || $1 == "ITA" || $1 == "ROU" {print $0}' ../../../nuclear_samplelist.keep | paste -sd ',')
USA=$(awk -F'_' '$1 == "USA" {print $0}' ../../../nuclear_samplelist.keep | paste -sd ',')

declare -A populations
populations=(
  [ASIA]="$ASIA"
  [AUS]="$AUS"
  [CAM]="$CAM"
  [EUR]="$EUR"
  [USA]="$USA"
)

# Function to create dataset and run smc++ commands
process_pair() {
  local pop1=$1
  local pop2=$2
  local out_prefix="../../${pop1}_${pop2}"
  local dir_prefix="DATA/${pop1}_${pop2}"
  
  # Generate VCF for the pair
  vcftools --gzvcf ../../../smcpp.vcf.gz \
    $(echo ${populations[$pop1]} | tr ',' '\n' | awk '{print "--indv "$0}') \
    $(echo ${populations[$pop2]} | tr ',' '\n' | awk '{print "--indv "$0}') \
    --max-missing 1 --recode --out $out_prefix

  bgzip -f ${out_prefix}.recode.vcf
  tabix ${out_prefix}.recode.vcf.gz

  for chr in {1..4}; do
    smc++ vcf2smc ${out_prefix}.recode.vcf.gz ${dir_prefix}/${pop1}_${pop2}.chr${chr}.smc.gz dirofilaria_immitis_chr${chr} ${pop1}:${populations[$pop1]} ${pop2}:${populations[$pop2]}
    smc++ vcf2smc ${out_prefix}.recode.vcf.gz ${dir_prefix}/${pop2}_${pop1}.chr${chr}.smc.gz dirofilaria_immitis_chr${chr} ${pop2}:${populations[$pop2]} ${pop1}:${populations[$pop1]}
  done

  smc++ split --timepoints 1 1000000 -o ${pop1}_${pop2} 2.7e-9 ../../../t1m/g5/${pop1}/model.final.json ../../../t1m/g5/${pop2}/model.final.json ${dir_prefix}/*.smc.gz
  smc++ plot -g 2.5 -c SMCPP_${pop1}_${pop2}_g2.5.png ${pop1}_${pop2}/model.final.json
}

# Process each population pair
process_pair ASIA AUS
process_pair ASIA CAM
process_pair ASIA EUR
process_pair ASIA USA
process_pair AUS CAM
process_pair AUS EUR
process_pair AUS USA
process_pair CAM EUR
process_pair CAM USA
process_pair EUR USA
```


cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/04_VARIANTS/FILTER2/NO_OUTGROUPS/SMCPP/split/g2.5/run2
bsub.py --done "run_smcpp_split_g2.5_run1" --queue long 10 run_smcpp_split_g2.5_run2 "run_smcpp_split_g2.5_run2.sh"

```bash
# Load modules
module load smcpp/1.15.3-c1
module load common-apps/htslib/1.9.229
module load vcftools/0.1.16-c4

cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/04_VARIANTS/FILTER2/NO_OUTGROUPS/SMCPP/split/g2.5/run2

# Create necessary directories
mkdir DATA
mkdir DATA/AUS_ASIA DATA/CAM_ASIA DATA/EUR_ASIA DATA/USA_ASIA
mkdir DATA/CAM_AUS DATA/EUR_AUS DATA/USA_AUS
mkdir DATA/EUR_CAM DATA/USA_CAM DATA/USA_EUR

# Get sample names for each population
ASIA=$(awk -F'_' '$1 == "MYS" || $1 == "THA" {print $0}' ../../../nuclear_samplelist.keep | paste -sd ',')
AUS=$(awk -F'_' '$1 == "AUS" {print $0}' ../../../nuclear_samplelist.keep | paste -sd ',')
CAM=$(awk -F'_' '$1 == "PAN" || $1 == "CRI" {print $0}' ../../../nuclear_samplelist.keep | paste -sd ',')
EUR=$(awk -F'_' '$1 == "GRC" || $1 == "ITA" || $1 == "ROU" {print $0}' ../../../nuclear_samplelist.keep | paste -sd ',')
USA=$(awk -F'_' '$1 == "USA" {print $0}' ../../../nuclear_samplelist.keep | paste -sd ',')

declare -A populations
populations=(
  [ASIA]="$ASIA"
  [AUS]="$AUS"
  [CAM]="$CAM"
  [EUR]="$EUR"
  [USA]="$USA"
)

# Function to create dataset and run smc++ commands
process_pair() {
  local pop1=$1
  local pop2=$2
  local out_prefix="../../${pop1}_${pop2}"
  local dir_prefix="DATA/${pop1}_${pop2}"
  
  # Generate VCF for the pair
  vcftools --gzvcf ../../../smcpp.vcf.gz \
    $(echo ${populations[$pop1]} | tr ',' '\n' | awk '{print "--indv "$0}') \
    $(echo ${populations[$pop2]} | tr ',' '\n' | awk '{print "--indv "$0}') \
    --max-missing 1 --recode --out $out_prefix

  bgzip -f ${out_prefix}.recode.vcf
  tabix ${out_prefix}.recode.vcf.gz

  for chr in {1..4}; do
    smc++ vcf2smc ${out_prefix}.recode.vcf.gz ${dir_prefix}/${pop1}_${pop2}.chr${chr}.smc.gz dirofilaria_immitis_chr${chr} ${pop1}:${populations[$pop1]} ${pop2}:${populations[$pop2]}
    smc++ vcf2smc ${out_prefix}.recode.vcf.gz ${dir_prefix}/${pop2}_${pop1}.chr${chr}.smc.gz dirofilaria_immitis_chr${chr} ${pop2}:${populations[$pop2]} ${pop1}:${populations[$pop1]}
  done

  smc++ split --timepoints 1 1000000 -o ${pop1}_${pop2} 2.7e-9 ../../../t1m/g5/${pop1}/model.final.json ../../../t1m/g5/${pop2}/model.final.json ${dir_prefix}/*.smc.gz
  smc++ plot -g 2.5 -c SMCPP_${pop1}_${pop2}_g2.5.png ${pop1}_${pop2}/model.final.json
}

# Process each population pair
process_pair AUS ASIA
process_pair CAM ASIA
process_pair EUR ASIA
process_pair USA ASIA
process_pair CAM AUS
process_pair EUR AUS
process_pair USA AUS
process_pair EUR CAM
process_pair USA CAM
process_pair USA EUR
```