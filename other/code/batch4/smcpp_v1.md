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

### timepoint = 1 to 1mill, generation time t = 2, 4 & 6 years

bsub.py --queue long 20 run_smcpp_t1m "run_smcpp_t1m.sh"

```bash
#!/bin/bash

# Load modules
module load smcpp/1.15.3-c1
module load common-apps/htslib/1.9.229
module load vcftools/0.1.16-c4

cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/05_ANALYSIS/SMCPP
mkdir DATA

# Get sample names for each population
ASIA=$(awk -F'_' '$1 == "MYS" || $1 == "THA" {print}' nuclear_samplelist.keep | paste -sd ',')
AUS=$(awk -F'_' '$1 == "AUS" {print}' nuclear_samplelist.keep | paste -sd ',')
CENAM=$(awk -F'_' '$1 == "PAN" || $1 == "CRI" {print}' nuclear_samplelist.keep | paste -sd ',')
EUR=$(awk -F'_' '$1 == "GRC" || $1 == "ITA" || $1 == "ROU" {print}' nuclear_samplelist.keep | paste -sd ',')
USA=$(awk -F'_' '$1 == "USA" {print}' nuclear_samplelist.keep | paste -sd ',')
declare -A populations
populations=(
  [ASIA]=$ASIA
  [AUS]=$AUS
  [CENAM]=$CENAM
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
# mutation rate (C.elegans) = 2.7e-9 (https://doi.org/10.1073/pnas.0904895106)
mkdir t1m
cd t1m
mkdir g2 g4 g6

smc++ estimate --timepoints 1 1000000 -o ASIA/ 2.7e-9 ../DATA/ASIA.*.smc.gz

# plot
## Use generation time of 2 years
smc++ plot -g 2 -c g2/SMCPP_ASIA_t1m_g2.pdf ASIA/model.final.json

## Use generation time of 4 years
smc++ plot -g 4 -c g4/SMCPP_ASIA_t1m_g4.pdf ASIA/model.final.json

## Use generation time of 6 years
smc++ plot -g 6 -c g6/SMCPP_ASIA_t1m_g6.pdf ASIA/model.final.json


# AUS
cd ..
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
# mutation rate (C.elegans) = 2.7e-9 (https://doi.org/10.1073/pnas.0904895106)
cd t1m
smc++ estimate --timepoints 1 1000000 -o AUS/ 2.7e-9 ../DATA/AUS.*.smc.gz

# plot
## Use generation time of 2 years
smc++ plot -g 2 -c g2/SMCPP_AUS_t1m_g2.pdf AUS/model.final.json

## Use generation time of 4 years
smc++ plot -g 4 -c g4/SMCPP_AUS_t1m_g4.pdf AUS/model.final.json

## Use generation time of 6 years
smc++ plot -g 6 -c g6/SMCPP_AUS_t1m_g6.pdf AUS/model.final.json



# CENAM
cd ..
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
--max-missing 1 --recode --out CENAM
bgzip -f CENAM.recode.vcf
tabix CENAM.recode.vcf.gz

# Convert VCF to the SMC++ input format with vcf2smc
for chr in {1..4}; do
  smc++ vcf2smc CENAM.recode.vcf.gz DATA/CENAM.chr${chr}.smc.gz dirofilaria_immitis_chr${chr} CENAM:${CENAM};
done
# each call to vcf2smc processes a single contig. VCFs containing multiple contigs should be processed via multiple separate runs.

# Fit the model using Estimate
# mutation rate (C.elegans) = 2.7e-9 (https://doi.org/10.1073/pnas.0904895106)
cd t1m
smc++ estimate --timepoints 1 1000000 -o CENAM/ 2.7e-9 ../DATA/CENAM.*.smc.gz

# plot
## Use generation time of 2 years
smc++ plot -g 2 -c g2/SMCPP_CENAM_t1m_g2.pdf CENAM/model.final.json

## Use generation time of 4 years
smc++ plot -g 4 -c g4/SMCPP_CENAM_t1m_g4.pdf CENAM/model.final.json

## Use generation time of 6 years
smc++ plot -g 6 -c g6/SMCPP_CENAM_t1m_g6.pdf CENAM/model.final.json




# EUR
cd ..
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
# mutation rate (C.elegans) = 2.7e-9 (https://doi.org/10.1073/pnas.0904895106)
cd t1m
smc++ estimate --timepoints 1 1000000 -o EUR/ 2.7e-9 ../DATA/EUR.*.smc.gz

# plot
## Use generation time of 2 years
smc++ plot -g 2 -c g2/SMCPP_EUR_t1m_g2.pdf EUR/model.final.json

## Use generation time of 4 years
smc++ plot -g 4 -c g4/SMCPP_EUR_t1m_g4.pdf EUR/model.final.json

## Use generation time of 6 years
smc++ plot -g 6 -c g6/SMCPP_EUR_t1m_g6.pdf EUR/model.final.json





# USA
cd ..
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
# mutation rate (C.elegans) = 2.7e-9 (https://doi.org/10.1073/pnas.0904895106)
cd t1m
smc++ estimate --timepoints 1 1000000 -o USA/ 2.7e-9 ../DATA/USA.*.smc.gz

# plot
## Use generation time of 2 years
smc++ plot -g 2 -c g2/SMCPP_USA_t1m_g2.pdf USA/model.final.json

## Use generation time of 4 years
smc++ plot -g 4 -c g4/SMCPP_USA_t1m_g4.pdf USA/model.final.json

## Use generation time of 6 years
smc++ plot -g 6 -c g6/SMCPP_USA_t1m_g6.pdf USA/model.final.json
```

Successfully completed.

```R
# smc++ - plotting all populations together
# timepoint = 1 to 1mill

# load libraries
library(tidyverse)
library(patchwork)
library(ggsci)
library(RColorBrewer)

####################################################################
# Generation = 2 yrs
####################################################################

setwd("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/batch4/FILTER1/NO_OUTGROUPS/smcpp/t1m/g2")

# load data
ASIA_data <- read.delim("SMCPP_ASIA_t1m_g2.csv", header = T , sep = ",")
ASIA_data$ID <- "ASIA"
AUS_data <- read.delim("SMCPP_AUS_t1m_g2.csv", header = T , sep = ",")
AUS_data$ID <- "AUS"
CENAM_data <- read.delim("SMCPP_CENAM_t1m_g2.csv", header = T , sep = ",")
CENAM_data$ID <- "CENAM"
EUR_data <- read.delim("SMCPP_EUR_t1m_g2.csv", header = T , sep = ",")
EUR_data$ID <- "EUR"
USA_data <- read.delim("SMCPP_USA_t1m_g2.csv", header = T , sep = ",")
USA_data$ID <- "USA"

data <- bind_rows(ASIA_data, AUS_data, CENAM_data, EUR_data, USA_data)

pop_colours <- c("ASIA" = "hotpink", "AUS" = "cornflowerblue", "CENAM" = "purple", "EUR" = "forestgreen", "USA" = "tomato2")

ggplot(data,aes(x,y,col=ID)) +
  geom_line(size=1) +
  labs(title = "Generation time = 2 years", x = "Years before present", y = "Effective population size (Ne)", col="Population") +
  theme_bw() +
  scale_x_log10(labels = prettyNum) +
  ylim(0,5e6) +
  scale_colour_manual(values = pop_colours)
# In scale_x_log10(labels = prettyNum) :
# log-10 transformation introduced infinite values.

ggsave("plot_smcpp_t1m_g2.png")
ggsave("plot_smcpp_t1m_g2.pdf", height = 4, width = 5, useDingbats = FALSE)


ggplot(data,aes(x,y,col=ID)) +
  geom_rect(aes(xmin=25000,ymin=0,xmax=35000,ymax=1.0E5), fill="grey80", col=NA) +
  geom_line(size=1) +
  labs(title = "Generation time = 2 years", x = "Years before present", y = "Effective population size (Ne)", col="Population") +
  theme_bw() +
  scale_x_log10(labels = prettyNum) +
  ylim(0, 1e5) +
  scale_colour_manual(values = pop_colours)
# 1: In scale_x_log10(labels = prettyNum) :
# log-10 transformation introduced infinite values.
# 2: Removed 19 rows containing missing values or values outside the scale range (`geom_line()`). 

ggsave("plot_smcpp_t1m_g2_cut.png")
ggsave("plot_smcpp_t1m_g2_cut.pdf", height = 4, width = 5, useDingbats = FALSE)




####################################################################
# Generation = 4 yrs
####################################################################

setwd("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/batch4/FILTER1/NO_OUTGROUPS/smcpp/t1m/g4")

# load data
ASIA_data <- read.delim("SMCPP_ASIA_t1m_g4.csv", header = T , sep = ",")
ASIA_data$ID <- "ASIA"
AUS_data <- read.delim("SMCPP_AUS_t1m_g4.csv", header = T , sep = ",")
AUS_data$ID <- "AUS"
CENAM_data <- read.delim("SMCPP_CENAM_t1m_g4.csv", header = T , sep = ",")
CENAM_data$ID <- "CENAM"
EUR_data <- read.delim("SMCPP_EUR_t1m_g4.csv", header = T , sep = ",")
EUR_data$ID <- "EUR"
USA_data <- read.delim("SMCPP_USA_t1m_g4.csv", header = T , sep = ",")
USA_data$ID <- "USA"

data <- bind_rows(ASIA_data, AUS_data, CENAM_data, EUR_data, USA_data)

ggplot(data,aes(x,y,col=ID)) +
  geom_line(size=1) +
  labs(title = "Generation time = 4 years", x = "Years before present", y = "Effective population size (Ne)", col="Population") +
  theme_bw() +
  scale_x_log10(labels = prettyNum) +
  ylim(0,5e6) +
  scale_colour_manual(values = pop_colours)
# In scale_x_log10(labels = prettyNum) :
# log-10 transformation introduced infinite values.

ggsave("plot_smcpp_t1m_g4.png")
ggsave("plot_smcpp_t1m_g4.pdf", height = 4, width = 5, useDingbats = FALSE)


ggplot(data,aes(x,y,col=ID)) +
  geom_rect(aes(xmin=25000,ymin=0,xmax=35000,ymax=1.0E5), fill="grey80", col=NA) +
  geom_line(size=1) +
  labs(title = "Generation time = 4 years", x = "Years before present", y = "Effective population size (Ne)", col="Population") +
  theme_bw() +
 scale_x_log10(labels = prettyNum) +
  ylim(0, 1e5) +
  scale_colour_manual(values = pop_colours)
# 1: In scale_x_log10(labels = prettyNum) :
# log-10 transformation introduced infinite values.
# 2: Removed 19 rows containing missing values or values outside the scale range (`geom_line()`). 

ggsave("plot_smcpp_t1m_g4_cut.png")
ggsave("plot_smcpp_t1m_g4_cut.pdf", height = 4, width = 5, useDingbats = FALSE)



####################################################################
# Generation = 6 yrs
####################################################################

setwd("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/batch4/FILTER1/NO_OUTGROUPS/smcpp/t1m/g6")

# load data
ASIA_data <- read.delim("SMCPP_ASIA_t1m_g6.csv", header = T , sep = ",")
ASIA_data$ID <- "ASIA"
AUS_data <- read.delim("SMCPP_AUS_t1m_g6.csv", header = T , sep = ",")
AUS_data$ID <- "AUS"
CENAM_data <- read.delim("SMCPP_CENAM_t1m_g6.csv", header = T , sep = ",")
CENAM_data$ID <- "CENAM"
EUR_data <- read.delim("SMCPP_EUR_t1m_g6.csv", header = T , sep = ",")
EUR_data$ID <- "EUR"
USA_data <- read.delim("SMCPP_USA_t1m_g6.csv", header = T , sep = ",")
USA_data$ID <- "USA"

data <- bind_rows(ASIA_data, AUS_data, CENAM_data, EUR_data, USA_data)

ggplot(data,aes(x,y,col=ID)) +
  geom_line(size=1) +
  labs(title = "Generation time = 6 years", x = "Years before present", y = "Effective population size (Ne)", col="Population") +
  theme_bw() +
  scale_x_log10(labels = prettyNum) +
  ylim(0,5e6) +
  scale_colour_manual(values = pop_colours)
# In scale_x_log10(labels = prettyNum) :
# log-10 transformation introduced infinite values.

ggsave("plot_smcpp_t1m_g6.png")
ggsave("plot_smcpp_t1m_g6.pdf", height = 4, width = 5, useDingbats = FALSE)


ggplot(data,aes(x,y,col=ID)) +
  geom_rect(aes(xmin=25000,ymin=0,xmax=35000,ymax=1.0E5), fill="grey80", col=NA) +
  geom_line(size=1) +
  labs(title = "Generation time = 6 years", x = "Years before present", y = "Effective population size (Ne)", col="Population") +
  theme_bw() +
  scale_x_log10(labels = prettyNum) +
  ylim(0, 1e5) +
  scale_colour_manual(values = pop_colours)
# 1: In scale_x_log10(labels = prettyNum) :
# log-10 transformation introduced infinite values.
# 2: Removed 19 rows containing missing values or values outside the scale range (`geom_line()`). 

ggsave("plot_smcpp_t1m_g6_cut.png")
ggsave("plot_smcpp_t1m_g6_cut.pdf", height = 4, width = 5, useDingbats = FALSE)
```



### timepoint = 1 to 5 mill, generation time = 1, 2, 4 & 6 years

```bash
bsub.py --queue long 20 run_smcpp_t5m "run_smcpp_t5m.sh"
```

```bash
#!/bin/bash

# Load modules
module load smcpp/1.15.3-c1

cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/05_ANALYSIS/SMCPP

# Get sample names for each population
ASIA=$(awk -F'_' '$1 == "MYS" || $1 == "THA" {print}' nuclear_samplelist.keep | paste -sd ',')
AUS=$(awk -F'_' '$1 == "AUS" {print}' nuclear_samplelist.keep | paste -sd ',')
CENAM=$(awk -F'_' '$1 == "PAN" || $1 == "CRI" {print}' nuclear_samplelist.keep | paste -sd ',')
EUR=$(awk -F'_' '$1 == "GRC" || $1 == "ITA" || $1 == "ROU" {print}' nuclear_samplelist.keep | paste -sd ',')
USA=$(awk -F'_' '$1 == "USA" {print}' nuclear_samplelist.keep | paste -sd ',')
declare -A populations
populations=(
  [ASIA]=$ASIA
  [AUS]=$AUS
  [CENAM]=$CENAM
  [EUR]=$EUR
  [USA]=$USA
)

# set up directories
mkdir t5m
cd t5m
mkdir g2 g4 g6
mkdir g1

# ASIA
# Fit the model using Estimate
## # mutation rate (C.elegans) = 2.7e-9 (https://doi.org/10.1073/pnas.0904895106)
smc++ estimate --timepoints 1 5000000 -o ASIA/ 2.7e-9 ../DATA/ASIA.*.smc.gz
# plot
## Use generation time of D. immitis 2, 4 & 6 yrs
smc++ plot -g 2 -c g2/SMCPP_ASIA_t5m_g2.pdf ASIA/model.final.json
smc++ plot -g 4 -c g4/SMCPP_ASIA_t5m_g4.pdf ASIA/model.final.json
smc++ plot -g 6 -c g6/SMCPP_ASIA_t5m_g6.pdf ASIA/model.final.json
# do gen=1 year
smc++ plot -g 1 -c g1/SMCPP_ASIA_t5m_g1.pdf ASIA/model.final.json

# AUS
# Fit the model using Estimate
## # mutation rate (C.elegans) = 2.7e-9 (https://doi.org/10.1073/pnas.0904895106)
smc++ estimate --timepoints 1 5000000 -o AUS/ 2.7e-9 ../DATA/AUS.*.smc.gz
# plot
## Use generation time of D. immitis 2, 4 & 6 yrs
smc++ plot -g 2 -c g2/SMCPP_AUS_t5m_g2.pdf AUS/model.final.json
smc++ plot -g 4 -c g4/SMCPP_AUS_t5m_g4.pdf AUS/model.final.json
smc++ plot -g 6 -c g6/SMCPP_AUS_t5m_g6.pdf AUS/model.final.json
# do gen=1 year
smc++ plot -g 1 -c g1/SMCPP_AUS_t5m_g1.pdf AUS/model.final.json

# CENAM
# Fit the model using Estimate
## # mutation rate (C.elegans) = 2.7e-9 (https://doi.org/10.1073/pnas.0904895106)
smc++ estimate --timepoints 1 5000000 -o CENAM/ 2.7e-9 ../DATA/CENAM.*.smc.gz
# plot
## Use generation time of D. immitis 2, 4 & 6 yrs
smc++ plot -g 2 -c g2/SMCPP_CENAM_t5m_g2.pdf CENAM/model.final.json
smc++ plot -g 4 -c g4/SMCPP_CENAM_t5m_g4.pdf CENAM/model.final.json
smc++ plot -g 6 -c g6/SMCPP_CENAM_t5m_g6.pdf CENAM/model.final.json
# do gen=1 year
smc++ plot -g 1 -c g1/SMCPP_CENAM_t5m_g1.pdf CENAM/model.final.json

# EUR
# Fit the model using Estimate
## # mutation rate (C.elegans) = 2.7e-9 (https://doi.org/10.1073/pnas.0904895106)
smc++ estimate --timepoints 1 5000000 -o EUR/ 2.7e-9 ../DATA/EUR.*.smc.gz
# plot
## Use generation time of D. immitis 2, 4 & 6 yrs
smc++ plot -g 2 -c g2/SMCPP_EUR_t5m_g2.pdf EUR/model.final.json
smc++ plot -g 4 -c g4/SMCPP_EUR_t5m_g4.pdf EUR/model.final.json
smc++ plot -g 6 -c g6/SMCPP_EUR_t5m_g6.pdf EUR/model.final.json
# do gen=1 year
smc++ plot -g 1 -c g1/SMCPP_EUR_t5m_g1.pdf EUR/model.final.json

# USA
# Fit the model using Estimate
## # mutation rate (C.elegans) = 2.7e-9 (https://doi.org/10.1073/pnas.0904895106)
smc++ estimate --timepoints 1 5000000 -o USA/ 2.7e-9 ../DATA/USA.*.smc.gz
# plot
## Use generation time of D. immitis 2, 4 & 6 yrs
smc++ plot -g 2 -c g2/SMCPP_USA_t5m_g2.pdf USA/model.final.json
smc++ plot -g 4 -c g4/SMCPP_USA_t5m_g4.pdf USA/model.final.json
smc++ plot -g 6 -c g6/SMCPP_USA_t5m_g6.pdf USA/model.final.json
# do gen=1 year
smc++ plot -g 1 -c g1/SMCPP_USA_t5m_g1.pdf USA/model.final.json
```

```R
# smc++ - plotting all populations together
# timepoint = 1 to 5 mill

# load libraries
library(tidyverse)
library(patchwork)
library(ggsci)
library(RColorBrewer)


####################################################################
# Generation = 1 yr
####################################################################

setwd("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/batch4/FILTER1/NO_OUTGROUPS/smcpp/t5m/g1")

# load data
ASIA_data <- read.delim("SMCPP_ASIA_t5m_g1.csv", header = T , sep = ",")
ASIA_data$ID <- "ASIA"
AUS_data <- read.delim("SMCPP_AUS_t5m_g1.csv", header = T , sep = ",")
AUS_data$ID <- "AUS"
CENAM_data <- read.delim("SMCPP_CENAM_t5m_g1.csv", header = T , sep = ",")
CENAM_data$ID <- "CENAM"
EUR_data <- read.delim("SMCPP_EUR_t5m_g1.csv", header = T , sep = ",")
EUR_data$ID <- "EUR"
USA_data <- read.delim("SMCPP_USA_t5m_g1.csv", header = T , sep = ",")
USA_data$ID <- "USA"

data <- bind_rows(ASIA_data, AUS_data, CENAM_data, EUR_data, USA_data)

pop_colours <- c("ASIA" = "hotpink", "AUS" = "cornflowerblue", "CENAM" = "purple", "EUR" = "forestgreen", "USA" = "tomato2")


ggplot(data,aes(x,y,col=ID)) +
  geom_rect(aes(xmin=25000,ymin=0,xmax=35000,ymax=200000), fill="grey80", col=NA) +
  geom_line(size=1) +
  labs(title = "Generation time = 1 year", x = "Years before present", y = "Effective population size (Ne)", col="Population") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_x_log10(labels = prettyNum) +
  ylim(0, 200000) +
  scale_colour_manual(values = pop_colours)
# 1: In scale_x_log10(labels = prettyNum) :
# log-10 transformation introduced infinite values.
# 2: Removed 19 rows containing missing values or values outside the scale range (`geom_line()`). 

ggsave("plot_smcpp_t5m_g1.png")
ggsave("plot_smcpp_t5m_g1.pdf", height = 4, width = 5, useDingbats = FALSE)


####################################################################
# Generation = 2 yrs
####################################################################

setwd("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/batch4/FILTER1/NO_OUTGROUPS/smcpp/t5m/g2")

# load data
ASIA_data <- read.delim("SMCPP_ASIA_t5m_g2.csv", header = T , sep = ",")
ASIA_data$ID <- "ASIA"
AUS_data <- read.delim("SMCPP_AUS_t5m_g2.csv", header = T , sep = ",")
AUS_data$ID <- "AUS"
CENAM_data <- read.delim("SMCPP_CENAM_t5m_g2.csv", header = T , sep = ",")
CENAM_data$ID <- "CENAM"
EUR_data <- read.delim("SMCPP_EUR_t5m_g2.csv", header = T , sep = ",")
EUR_data$ID <- "EUR"
USA_data <- read.delim("SMCPP_USA_t5m_g2.csv", header = T , sep = ",")
USA_data$ID <- "USA"

data <- bind_rows(ASIA_data, AUS_data, CENAM_data, EUR_data, USA_data)


ggplot(data,aes(x,y,col=ID)) +
  geom_rect(aes(xmin=25000,ymin=0,xmax=35000,ymax=300000), fill="grey80", col=NA) +
  geom_line(size=1) +
  labs(title = "Generation time = 2 years", x = "Years before present", y = "Effective population size (Ne)", col="Population") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_x_log10(labels = prettyNum) +
  ylim(0, 300000) +
  scale_colour_manual(values = pop_colours)
# 1: In scale_x_log10(labels = prettyNum) :
# log-10 transformation introduced infinite values.
# 2: Removed 19 rows containing missing values or values outside the scale range (`geom_line()`). 

ggsave("plot_smcpp_t5m_g2.png")
ggsave("plot_smcpp_t5m_g2.pdf", height = 4, width = 5, useDingbats = FALSE)


ggplot(data,aes(x,y,col=ID)) +
  geom_rect(aes(xmin=25000,ymin=0,xmax=35000,ymax=100000), fill="grey80", col=NA) +
  geom_line(size=1) +
  labs(title = "Generation time = 2 years", x = "Years before present", y = "Effective population size (Ne)", col="Population") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_x_log10(labels = prettyNum) +
  ylim(0, 100000) +
  scale_colour_manual(values = pop_colours)

ggsave("plot_smcpp_t5m_g2_cut.png")
ggsave("plot_smcpp_t5m_g2_cut.pdf", height = 4, width = 5, useDingbats = FALSE)


####################################################################
# Generation = 4 yrs
####################################################################

setwd("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/batch4/FILTER1/NO_OUTGROUPS/smcpp/t5m/g4")

# load data
ASIA_data <- read.delim("SMCPP_ASIA_t5m_g4.csv", header = T , sep = ",")
ASIA_data$ID <- "ASIA"
AUS_data <- read.delim("SMCPP_AUS_t5m_g4.csv", header = T , sep = ",")
AUS_data$ID <- "AUS"
CENAM_data <- read.delim("SMCPP_CENAM_t5m_g4.csv", header = T , sep = ",")
CENAM_data$ID <- "CENAM"
EUR_data <- read.delim("SMCPP_EUR_t5m_g4.csv", header = T , sep = ",")
EUR_data$ID <- "EUR"
USA_data <- read.delim("SMCPP_USA_t5m_g4.csv", header = T , sep = ",")
USA_data$ID <- "USA"

data <- bind_rows(ASIA_data, AUS_data, CENAM_data, EUR_data, USA_data)


ggplot(data,aes(x,y,col=ID)) +
  geom_rect(aes(xmin=25000,ymin=0,xmax=35000,ymax=300000), fill="grey80", col=NA) +
  geom_line(size=1) +
  labs(title = "Generation time = 4 years", x = "Years before present", y = "Effective population size (Ne)", col="Population") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
 scale_x_log10(labels = prettyNum) +
  ylim(0, 300000) +
  scale_colour_manual(values = pop_colours)
# 1: In scale_x_log10(labels = prettyNum) :
# log-10 transformation introduced infinite values.
# 2: Removed 19 rows containing missing values or values outside the scale range (`geom_line()`). 

ggsave("plot_smcpp_t5m_g4.png")
ggsave("plot_smcpp_t5m_g4.pdf", height = 4, width = 5, useDingbats = FALSE)


ggplot(data,aes(x,y,col=ID)) +
  geom_rect(aes(xmin=25000,ymin=0,xmax=35000,ymax=100000), fill="grey80", col=NA) +
  geom_line(size=1) +
  labs(title = "Generation time = 4 years", x = "Years before present", y = "Effective population size (Ne)", col="Population") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_x_log10(labels = prettyNum) +
  ylim(0, 100000) +
  scale_colour_manual(values = pop_colours)

ggsave("plot_smcpp_t5m_g4_cut.png")
ggsave("plot_smcpp_t5m_g4_cut.pdf", height = 4, width = 5, useDingbats = FALSE)



####################################################################
# Generation = 6 yrs
####################################################################

setwd("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/batch4/FILTER1/NO_OUTGROUPS/smcpp/t5m/g6")

# load data
ASIA_data <- read.delim("SMCPP_ASIA_t5m_g6.csv", header = T , sep = ",")
ASIA_data$ID <- "ASIA"
AUS_data <- read.delim("SMCPP_AUS_t5m_g6.csv", header = T , sep = ",")
AUS_data$ID <- "AUS"
CENAM_data <- read.delim("SMCPP_CENAM_t5m_g6.csv", header = T , sep = ",")
CENAM_data$ID <- "CENAM"
EUR_data <- read.delim("SMCPP_EUR_t5m_g6.csv", header = T , sep = ",")
EUR_data$ID <- "EUR"
USA_data <- read.delim("SMCPP_USA_t5m_g6.csv", header = T , sep = ",")
USA_data$ID <- "USA"

data <- bind_rows(ASIA_data, AUS_data, CENAM_data, EUR_data, USA_data)

ggplot(data,aes(x,y,col=ID)) +
  geom_rect(aes(xmin=25000,ymin=0,xmax=35000,ymax=300000), fill="grey80", col=NA) +
  geom_line(size=1) +
  labs(title = "Generation time = 6 years", x = "Years before present", y = "Effective population size (Ne)", col="Population") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_x_log10(labels = prettyNum) +
  ylim(0,300000) +
  scale_colour_manual(values = pop_colours)
# In scale_x_log10(labels = prettyNum) :
# log-10 transformation introduced infinite values.

ggsave("plot_smcpp_t5m_g6.png")
ggsave("plot_smcpp_t5m_g6.pdf", height = 4, width = 5, useDingbats = FALSE)


ggplot(data,aes(x,y,col=ID)) +
  geom_rect(aes(xmin=25000,ymin=0,xmax=35000,ymax=100000), fill="grey80", col=NA) +
  geom_line(size=1) +
  labs(title = "Generation time = 6 years", x = "Years before present", y = "Effective population size (Ne)", col="Population") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_x_log10(labels = prettyNum) +
  ylim(0, 100000) +
  scale_colour_manual(values = pop_colours)
# 1: In scale_x_log10(labels = prettyNum) :
# log-10 transformation introduced infinite values.
# 2: Removed 19 rows containing missing values or values outside the scale range (`geom_line()`). 

ggsave("plot_smcpp_t5m_g6_cut.png")
ggsave("plot_smcpp_t5m_g6_cut.pdf", height = 4, width = 5, useDingbats = FALSE)
```




### timepoint = 1 to 1 mill, generation time = 4 years, -c 1kbp

The -c parameter will treat runs of homozygosity longer than -c bp as missing.

```bash
bsub.py --queue long 20 run_smcpp_c1kbp "run_smcpp_c1kbp.sh"
```

```bash
#!/bin/bash

# Load modules
module load smcpp/1.15.3-c1

cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/05_ANALYSIS/SMCPP/c1kbp
mkdir DATA

# Get sample names for each population
ASIA=$(awk -F'_' '$1 == "MYS" || $1 == "THA" {print}' ../nuclear_samplelist.keep | paste -sd ',')
AUS=$(awk -F'_' '$1 == "AUS" {print}' ../nuclear_samplelist.keep | paste -sd ',')
CENAM=$(awk -F'_' '$1 == "PAN" || $1 == "CRI" {print}' ../nuclear_samplelist.keep | paste -sd ',')
EUR=$(awk -F'_' '$1 == "GRC" || $1 == "ITA" || $1 == "ROU" {print}' ../nuclear_samplelist.keep | paste -sd ',')
USA=$(awk -F'_' '$1 == "USA" {print}' ../nuclear_samplelist.keep | paste -sd ',')
declare -A populations
populations=(
  [ASIA]=$ASIA
  [AUS]=$AUS
  [CENAM]=$CENAM
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
# mutation rate (C.elegans) = 2.7e-9 (https://doi.org/10.1073/pnas.0904895106)
smc++ estimate --timepoints 1 1000000 -o ASIA/ 2.7e-9 DATA/ASIA.*.smc.gz

# plot
## Use generation time of D. immitis ~ 4yrs
smc++ plot -g 4 -c SMCPP_ASIA_c1kbp.pdf ASIA/model.final.json


# AUS
# Convert VCF to the SMC++ input format with vcf2smc
for chr in {1..4}; do
  smc++ vcf2smc -c 1000 ../AUS.recode.vcf.gz DATA/AUS.chr${chr}.smc.gz dirofilaria_immitis_chr${chr} AUS:${AUS};
done
# each call to vcf2smc processes a single contig. VCFs containing multiple contigs should be processed via multiple separate runs.

# Fit the model using Estimate
# mutation rate (C.elegans) = 2.7e-9 (https://doi.org/10.1073/pnas.0904895106)
smc++ estimate --timepoints 1 1000000 -o AUS/ 2.7e-9 DATA/AUS.*.smc.gz

# plot
## Use generation time of D. immitis ~ 4 yrs
smc++ plot -g 4 -c SMCPP_AUS_c1kbp.pdf AUS/model.final.json



# CENAM
# Convert VCF to the SMC++ input format with vcf2smc
for chr in {1..4}; do
  smc++ vcf2smc -c 1000 ../CENAM.recode.vcf.gz DATA/CENAM.chr${chr}.smc.gz dirofilaria_immitis_chr${chr} CENAM:${CENAM};
done
# each call to vcf2smc processes a single contig. VCFs containing multiple contigs should be processed via multiple separate runs.

# Fit the model using Estimate
# mutation rate (C.elegans) = 2.7e-9 (https://doi.org/10.1073/pnas.0904895106)
smc++ estimate --timepoints 1 1000000 -o CENAM/ 2.7e-9 DATA/CENAM.*.smc.gz

# plot
## Use generation time of D. immitis ~ 4 yrs
smc++ plot -g 4 -c SMCPP_CENAM_c1kbp.pdf CENAM/model.final.json




# EUR
# Convert VCF to the SMC++ input format with vcf2smc
for chr in {1..4}; do
  smc++ vcf2smc -c 1000 ../EUR.recode.vcf.gz DATA/EUR.chr${chr}.smc.gz dirofilaria_immitis_chr${chr} EUR:${EUR};
done
# each call to vcf2smc processes a single contig. VCFs containing multiple contigs should be processed via multiple separate runs.

# Fit the model using Estimate
# mutation rate (C.elegans) = 2.7e-9 (https://doi.org/10.1073/pnas.0904895106)
smc++ estimate --timepoints 1 1000000 -o EUR/ 2.7e-9 DATA/EUR.*.smc.gz

# plot
## Use generation time of D. immitis ~ 4 yrs
smc++ plot -g 4 -c SMCPP_EUR_c1kbp.pdf EUR/model.final.json





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
## Use generation time of D. immitis ~ 4 yrs
smc++ plot -g 4 -c SMCPP_USA_c1kbp.pdf USA/model.final.json
```

```R
# smc++ - plotting all populations together
# timepoint = 1 to 1 mill, g=4 years, -c 1kbp

# load libraries
library(tidyverse)
library(patchwork)
library(ggsci)
library(RColorBrewer)


####################################################################
# Generation = 4 years
####################################################################

setwd("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/batch4/FILTER1/NO_OUTGROUPS/smcpp/c1kbp")

# load data
ASIA_data <- read.delim("SMCPP_ASIA_c1kbp.csv", header = T , sep = ",")
ASIA_data$ID <- "ASIA"
AUS_data <- read.delim("SMCPP_AUS_c1kbp.csv", header = T , sep = ",")
AUS_data$ID <- "AUS"
CENAM_data <- read.delim("SMCPP_CENAM_c1kbp.csv", header = T , sep = ",")
CENAM_data$ID <- "CENAM"
EUR_data <- read.delim("SMCPP_EUR_c1kbp.csv", header = T , sep = ",")
EUR_data$ID <- "EUR"
USA_data <- read.delim("SMCPP_USA_c1kbp.csv", header = T , sep = ",")
USA_data$ID <- "USA"

data <- bind_rows(ASIA_data, AUS_data, CENAM_data, EUR_data, USA_data)

pop_colours <- c("ASIA" = "hotpink", "AUS" = "cornflowerblue", "CENAM" = "purple", "EUR" = "forestgreen", "USA" = "tomato2")


ggplot(data,aes(x,y,col=ID)) +
  geom_rect(aes(xmin=25000,ymin=0,xmax=35000,ymax=5e6), fill="grey80", col=NA) +
  geom_line(size=1) +
  labs(title = "Generation time = 4 years, Runs of homozygosity > 1kbp = missing", x = "Years before present", y = "Effective population size (Ne)", col="Population") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_x_log10(labels = prettyNum) +
  ylim(0, 5e6) +
  scale_colour_manual(values = pop_colours)
# 1: In scale_x_log10(labels = prettyNum) :
# log-10 transformation introduced infinite values.

ggsave("plot_smcpp_c1kbp.png")
ggsave("plot_smcpp_c1kbp.pdf", height = 4, width = 5, useDingbats = FALSE)



ggplot(data,aes(x,y,col=ID)) +
  geom_rect(aes(xmin=25000,ymin=0,xmax=35000,ymax=3e5), fill="grey80", col=NA) +
  geom_line(size=1) +
  labs(title = "Generation time = 4 years, Runs of homozygosity > 1kbp = missing", x = "Years before present", y = "Effective population size (Ne)", col="Population") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_x_log10(labels = prettyNum) +
  ylim(0, 3e5) +
  scale_colour_manual(values = pop_colours)

ggsave("plot_smcpp_c1kbp_cut.png")
ggsave("plot_smcpp_c1kbp_cut.pdf", height = 4, width = 5, useDingbats = FALSE)
```





## Split

The split command fits two-population clean split models (assumes no ongoing gene flow between the two populations after they diverged).

## Timepoint 1 to 5mill generation time ~4 years

```bash
cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/05_ANALYSIS/SMCPP
bsub.py --queue long --threads 4 20 run_smcpp_split_t5m_g4 "run_smcpp_split_t5m_g4.sh"
```

```bash
#!/bin/bash

# Load modules
module load smcpp/1.15.3-c1
module load common-apps/htslib/1.9.229
module load vcftools/0.1.16-c4

cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/05_ANALYSIS/SMCPP
mkdir SPLIT
cd SPLIT

# Create necessary directories
mkdir DATA
mkdir DATA/ASIA_AUS DATA/ASIA_CENAM DATA/ASIA_EUR DATA/ASIA_USA
mkdir DATA/AUS_CENAM DATA/AUS_EUR DATA/AUS_USA
mkdir DATA/CENAM_EUR DATA/CENAM_USA DATA/EUR_USA

# Get sample names for each population
ASIA=$(awk -F'_' '$1 == "MYS" || $1 == "THA" {print $0}' ../nuclear_samplelist.keep | paste -sd ',')
AUS=$(awk -F'_' '$1 == "AUS" {print $0}' ../nuclear_samplelist.keep | paste -sd ',')
CENAM=$(awk -F'_' '$1 == "PAN" || $1 == "CRI" {print $0}' ../nuclear_samplelist.keep | paste -sd ',')
EUR=$(awk -F'_' '$1 == "GRC" || $1 == "ITA" || $1 == "ROU" {print $0}' ../nuclear_samplelist.keep | paste -sd ',')
USA=$(awk -F'_' '$1 == "USA" {print $0}' ../nuclear_samplelist.keep | paste -sd ',')

declare -A populations
populations=(
  [ASIA]="$ASIA"
  [AUS]="$AUS"
  [CENAM]="$CENAM"
  [EUR]="$EUR"
  [USA]="$USA"
)

# Function to create dataset and run smc++ commands
process_pair() {
  local pop1=$1
  local pop2=$2
  local out_prefix="./${pop1}_${pop2}"
  local dir_prefix="DATA/${pop1}_${pop2}"
  
  # Generate VCF for the pair
  vcftools --gzvcf ../smcpp.vcf.gz \
    $(echo ${populations[$pop1]} | tr ',' '\n' | awk '{print "--indv "$0}') \
    $(echo ${populations[$pop2]} | tr ',' '\n' | awk '{print "--indv "$0}') \
    --max-missing 1 --recode --out $out_prefix

  bgzip -f ${out_prefix}.recode.vcf
  tabix ${out_prefix}.recode.vcf.gz

  # create datasets containing the joint frequency spectrum for both pops
  for chr in {1..4}; do
    smc++ vcf2smc ${out_prefix}.recode.vcf.gz ${dir_prefix}/${pop1}_${pop2}.chr${chr}.smc.gz dirofilaria_immitis_chr${chr} ${pop1}:${populations[$pop1]} ${pop2}:${populations[$pop2]}
    smc++ vcf2smc ${out_prefix}.recode.vcf.gz ${dir_prefix}/${pop2}_${pop1}.chr${chr}.smc.gz dirofilaria_immitis_chr${chr} ${pop2}:${populations[$pop2]} ${pop1}:${populations[$pop1]}
  done

  # copy smc++ data for individual pops
  cp ../DATA/${pop1}*.smc.gz ${dir_prefix}
  cp ../DATA/${pop2}*.smc.gz ${dir_prefix}

  # run split
  smc++ split --timepoints 1 5000000 -o ${pop1}_${pop2} ../t5m/${pop1}/model.final.json ../t5m/${pop2}/model.final.json ${dir_prefix}/*.smc.gz
  smc++ plot -g 4 -c SMCPP_${pop1}_${pop2}_t5m_g4.pdf ${pop1}_${pop2}/model.final.json
}

# Process each population pair
process_pair ASIA AUS
process_pair ASIA CENAM
process_pair ASIA EUR
process_pair ASIA USA
process_pair AUS CENAM
process_pair AUS EUR
process_pair AUS USA
process_pair CENAM EUR
process_pair CENAM USA
process_pair EUR USA
```

```R
# smc++ - plotting all populations together
# timepoint = 1 to 5 million
## include population splits

# load libraries
library(tidyverse)
library(patchwork)
library(ggsci)
library(RColorBrewer)


####################################################################
# Generation = 4 yr
####################################################################

setwd("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/batch4/FILTER1/NO_OUTGROUPS/smcpp/t5m/g4")

# load data
ASIA_data <- read.delim("SMCPP_ASIA_t5m_g4.csv", header = T , sep = ",")
ASIA_data$ID <- "ASIA"
AUS_data <- read.delim("SMCPP_AUS_t5m_g4.csv", header = T , sep = ",")
AUS_data$ID <- "AUS"
CENAM_data <- read.delim("SMCPP_CENAM_t5m_g4.csv", header = T , sep = ",")
CENAM_data$ID <- "CENAM"
EUR_data <- read.delim("SMCPP_EUR_t5m_g4.csv", header = T , sep = ",")
EUR_data$ID <- "EUR"
USA_data <- read.delim("SMCPP_USA_t5m_g4.csv", header = T , sep = ",")
USA_data$ID <- "USA"

data <- bind_rows(ASIA_data, AUS_data, CENAM_data, EUR_data, USA_data)

pop_colours <- c("ASIA" = "hotpink", "AUS" = "cornflowerblue", "CENAM" = "purple", "EUR" = "forestgreen", "USA" = "tomato2")

splits <- read.csv("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/batch4/FILTER1/NO_OUTGROUPS/smcpp/SPLIT/t5m/g4/split_dates_g4.csv")


ggplot(data,aes(x,y,col=ID)) +
  geom_rect(aes(xmin=25000,ymin=0,xmax=35000,ymax=150000), fill="grey80", col=NA) +
  geom_line(size=2) +
  geom_vline(data = splits, aes(xintercept = end_x), linetype = "dotted", color = "black", size = 0.6) +
  labs(x = "Years before present", y = "Effective population size (Ne)", col="Population") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_x_log10(labels = prettyNum) +
  ylim(0, 155000) +
  scale_colour_manual(values = pop_colours) +
  theme(text = element_text(size = 20),
        legend.position = "right",
        legend.justification = c(0, 1),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(100, 10, 10, 20)
        )
# 1: In scale_x_log10(labels = prettyNum) :
# log-10 transformation introduced infinite values.

ggsave("plot_smcpp_t5m_g4_splits.tif", height = 8, width = 16, dpi = 600)
ggsave("plot_smcpp_t5m_g4_splits.pdf", height = 8, width = 16, useDingbats = FALSE, dpi = 600)




# Plotting the splits values to double check 
# Ensure the splits dataframe contains the necessary columns
splits <- splits %>%
  mutate(label = paste0("x = ", prettyNum(end_x, big.mark = ",")))

# Plotting with text labels for vertical lines
ggplot(data, aes(x, y, col = ID)) +
  geom_rect(aes(xmin = 25000, ymin = 0, xmax = 35000, ymax = 150000), fill = "grey80", col = NA) +
  geom_line(size = 2) +
  geom_vline(data = splits, aes(xintercept = end_x), linetype = "dotted", color = "black", size = 0.6) +
  geom_text_repel(data = splits, aes(x = end_x, y = 155000, label = label), 
            angle = 90, vjust = -0.5, hjust = 1, size = 3, inherit.aes = FALSE) +
  labs(x = "Years before present", y = "Effective population size (Ne)", col = "Population") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_x_log10(labels = prettyNum) +
  ylim(0, 155000) +
  scale_colour_manual(values = pop_colours)
```


## Timepoint 1y to 5mill generation time ~1 year

```bash
cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/05_ANALYSIS/SMCPP
bsub.py --queue long 20 run_smcpp_split_t5m_g1 "run_smcpp_split_t5m_g1.sh"
```

```bash
#!/bin/bash

# Load modules
module load smcpp/1.15.3-c1
module load common-apps/htslib/1.9.229
module load vcftools/0.1.16-c4

cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/05_ANALYSIS/SMCPP/SPLIT

# Get sample names for each population
ASIA=$(awk -F'_' '$1 == "MYS" || $1 == "THA" {print $0}' ../nuclear_samplelist.keep | paste -sd ',')
AUS=$(awk -F'_' '$1 == "AUS" {print $0}' ../nuclear_samplelist.keep | paste -sd ',')
CENAM=$(awk -F'_' '$1 == "PAN" || $1 == "CRI" {print $0}' ../nuclear_samplelist.keep | paste -sd ',')
EUR=$(awk -F'_' '$1 == "GRC" || $1 == "ITA" || $1 == "ROU" {print $0}' ../nuclear_samplelist.keep | paste -sd ',')
USA=$(awk -F'_' '$1 == "USA" {print $0}' ../nuclear_samplelist.keep | paste -sd ',')

declare -A populations
populations=(
  [ASIA]="$ASIA"
  [AUS]="$AUS"
  [CENAM]="$CENAM"
  [EUR]="$EUR"
  [USA]="$USA"
)

# Function to create dataset and run smc++ commands
process_pair() {
  local pop1=$1
  local pop2=$2

  smc++ plot -g 1 -c g1/SMCPP_${pop1}_${pop2}_t5m_g1.pdf ${pop1}_${pop2}/model.final.json
}

# Process each population pair
process_pair ASIA AUS
process_pair ASIA CENAM
process_pair ASIA EUR
process_pair ASIA USA
process_pair AUS CENAM
process_pair AUS EUR
process_pair AUS USA
process_pair CENAM EUR
process_pair CENAM USA
process_pair EUR USA
```


```R
# smc++ - plotting all populations together
# timepoint = 1 to 5 million
## include population splits

# load libraries
library(tidyverse)
library(patchwork)
library(ggsci)
library(RColorBrewer)
library(ggrepel)


####################################################################
# Generation = 4 yr
####################################################################

setwd("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/batch4/FILTER1/NO_OUTGROUPS/smcpp/t5m/g1")

# load data
ASIA_data <- read.delim("SMCPP_ASIA_t5m_g1.csv", header = T , sep = ",")
ASIA_data$ID <- "ASIA"
AUS_data <- read.delim("SMCPP_AUS_t5m_g1.csv", header = T , sep = ",")
AUS_data$ID <- "AUS"
CENAM_data <- read.delim("SMCPP_CENAM_t5m_g1.csv", header = T , sep = ",")
CENAM_data$ID <- "CENAM"
EUR_data <- read.delim("SMCPP_EUR_t5m_g1.csv", header = T , sep = ",")
EUR_data$ID <- "EUR"
USA_data <- read.delim("SMCPP_USA_t5m_g1.csv", header = T , sep = ",")
USA_data$ID <- "USA"

data <- bind_rows(ASIA_data, AUS_data, CENAM_data, EUR_data, USA_data)

pop_colours <- c("ASIA" = "hotpink", "AUS" = "cornflowerblue", "CENAM" = "purple", "EUR" = "forestgreen", "USA" = "tomato2")

splits <- read.csv("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/batch4/FILTER1/NO_OUTGROUPS/smcpp/SPLIT/t5m/g1/split_dates_g1.csv")


ggplot(data,aes(x,y,col=ID)) +
  geom_rect(aes(xmin=25000,ymin=0,xmax=35000,ymax=150000), fill="grey80", col=NA) +
  geom_line(size=2) +
  geom_vline(data = splits, aes(xintercept = end_x), linetype = "dotted", color = "black", size = 0.6) +
  labs(x = "Years before present", y = "Effective population size (Ne)", col="Population") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_x_log10(labels = prettyNum) +
  ylim(0, 155000) +
  scale_colour_manual(values = pop_colours) +
  theme(text = element_text(size = 20),
        legend.position = "right",
        legend.justification = c(0, 1),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(100, 10, 10, 20)
        )
# 1: In scale_x_log10(labels = prettyNum) :
# log-10 transformation introduced infinite values.

ggsave("plot_smcpp_t5m_g1_splits.tif", height = 8, width = 16, dpi = 600)
ggsave("plot_smcpp_t5m_g1_splits.pdf", height = 8, width = 16, useDingbats = FALSE, dpi = 600)




# Plotting the splits values to double check 
# Ensure the splits dataframe contains the necessary columns
splits <- splits %>%
  mutate(label = paste0("x = ", prettyNum(end_x, big.mark = ",")))

# Plotting with text labels for vertical lines
ggplot(data, aes(x, y, col = ID)) +
  geom_rect(aes(xmin = 25000, ymin = 0, xmax = 35000, ymax = 150000), fill = "grey80", col = NA) +
  geom_line(size = 2) +
  geom_vline(data = splits, aes(xintercept = end_x), linetype = "dotted", color = "black", size = 0.6) +
  geom_text_repel(data = splits, aes(x = end_x, y = 155000, label = label), 
            angle = 90, vjust = -0.5, hjust = 1, size = 3, inherit.aes = FALSE) +
  labs(x = "Years before present", y = "Effective population size (Ne)", col = "Population") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_x_log10(labels = prettyNum) +
  ylim(0, 155000) +
  scale_colour_manual(values = pop_colours)


```
