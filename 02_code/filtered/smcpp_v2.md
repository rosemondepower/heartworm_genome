# SMC++ for estimating size history of populations - FINAL

Run:
- smc++ with split by region (dog hosts only)
- smc++ by host

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



## Run smc++ by region again but for DOG HOSTS ONLY

### timepoint = 1 to 5 mill, generation time = 1, 2, 4 & 6 years

```bash
# load modules
module load PaM/environment
module load bsub.py/0.42.1
module load smcpp/1.15.3-c1
# module load common-apps/htslib/1.9.229 # not available
module load htslib-1.19/perl-5.38.0
module load vcftools/0.1.16-c4

cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/05_ANALYSIS/SMCPP/DOG_ONLY
chmod a+x run_smcpp_t5m_DOG_ONLY.sh
bsub.py --queue long 20 run_smcpp_t5m_DOG_ONLY "run_smcpp_t5m_DOG_ONLY.sh"
```

```bash
#!/bin/bash

# Load modules
module load PaM/environment
module load smcpp/1.15.3-c1
# module load common-apps/htslib/1.9.229 # not available
module load htslib-1.19/perl-5.38.0
module load vcftools/0.1.16-c4

cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/05_ANALYSIS/SMCPP/DOG_ONLY

# Get sample names for each population
ASIA=$(awk -F'_' '$1 == "MYS" || $1 == "THA" {print}' DOG_ONLY_samplelist.keep | paste -sd ',')
AUS=$(awk -F'_' '$1 == "AUS" {print}' DOG_ONLY_samplelist.keep | paste -sd ',')
CENAM=$(awk -F'_' '$1 == "PAN" || $1 == "CRI" {print}' DOG_ONLY_samplelist.keep | paste -sd ',')
EUR=$(awk -F'_' '$1 == "GRC" || $1 == "ITA" || $1 == "ROU" {print}' DOG_ONLY_samplelist.keep | paste -sd ',')
USA=$(awk -F'_' '$1 == "USA" {print}' DOG_ONLY_samplelist.keep | paste -sd ',')
declare -A populations
populations=(
  [ASIA]=$ASIA
  [AUS]=$AUS
  [CENAM]=$CENAM
  [EUR]=$EUR
  [USA]=$USA
)
# DOG_ONLY_samplelist.keep - is from when I was running pixy - it's basically all the samples in the final VCF for dog hosts only and minus repeat samples.

mkdir DATA

# ASIA
vcftools --gzvcf ../smcpp.vcf.gz \
--indv MYS_SEL_AD_001 \
--indv THA_BKK_AD_001 \
--indv THA_BKK_AD_002 \
--indv THA_BKK_AD_003 \
--indv THA_BKK_AD_004 \
--indv THA_BKK_AD_005 \
--indv THA_BKK_AD_006 \
--max-missing 1 --recode --out ASIA_DOG_ONLY
bgzip -f ASIA_DOG_ONLY.recode.vcf
tabix ASIA_DOG_ONLY.recode.vcf.gz

# Convert VCF to the SMC++ input format with vcf2smc
for chr in {1..4}; do
  smc++ vcf2smc ASIA_DOG_ONLY.recode.vcf.gz DATA/ASIA_DOG_ONLY.chr${chr}.smc.gz dirofilaria_immitis_chr${chr} ASIA:${ASIA};
done
# each call to vcf2smc processes a single contig. VCFs containing multiple contigs should be processed via multiple separate runs.

# Fit the model using Estimate
## # mutation rate (C.elegans) = 2.7e-9 (https://doi.org/10.1073/pnas.0904895106)
mkdir t5m
cd t5m
mkdir g1 g2 g4 g6

smc++ estimate --timepoints 1 5000000 -o ASIA_DOG_ONLY/ 2.7e-9 ../DATA/ASIA_DOG_ONLY.*.smc.gz

# plot
## Use generation time of D. immitis 1, 2, 4 & 6 yrs
smc++ plot -g 1 -c g1/SMCPP_ASIA_DOG_ONLY_t5m_g1.pdf ASIA_DOG_ONLY/model.final.json
smc++ plot -g 2 -c g2/SMCPP_ASIA_DOG_ONLY_t5m_g2.pdf ASIA_DOG_ONLY/model.final.json
smc++ plot -g 4 -c g4/SMCPP_ASIA_DOG_ONLY_t5m_g4.pdf ASIA_DOG_ONLY/model.final.json
smc++ plot -g 6 -c g6/SMCPP_ASIA_DOG_ONLY_t5m_g6.pdf ASIA_DOG_ONLY/model.final.json


# AUS
cd ..
vcftools --gzvcf ../smcpp.vcf.gz \
--indv AUS_BNE_AD_001 \
--indv AUS_BNE_AD_002 \
--indv AUS_BNE_AD_003 \
--indv AUS_BNE_AD_004 \
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
--max-missing 1 --recode --out AUS_DOG_ONLY
bgzip -f AUS_DOG_ONLY.recode.vcf
tabix AUS_DOG_ONLY.recode.vcf.gz

# Convert VCF to the SMC++ input format with vcf2smc
for chr in {1..4}; do
  smc++ vcf2smc AUS_DOG_ONLY.recode.vcf.gz DATA/AUS_DOG_ONLY.chr${chr}.smc.gz dirofilaria_immitis_chr${chr} AUS:${AUS};
done
# each call to vcf2smc processes a single contig. VCFs containing multiple contigs should be processed via multiple separate runs.


# Fit the model using Estimate
## # mutation rate (C.elegans) = 2.7e-9 (https://doi.org/10.1073/pnas.0904895106)
cd t5m
smc++ estimate --timepoints 1 5000000 -o AUS_DOG_ONLY/ 2.7e-9 ../DATA/AUS_DOG_ONLY.*.smc.gz

# plot
## Use generation time of D. immitis 1, 2, 4 & 6 yrs
smc++ plot -g 1 -c g1/SMCPP_AUS_DOG_ONLY_t5m_g1.pdf AUS_DOG_ONLY/model.final.json
smc++ plot -g 2 -c g2/SMCPP_AUS_DOG_ONLY_t5m_g2.pdf AUS_DOG_ONLY/model.final.json
smc++ plot -g 4 -c g4/SMCPP_AUS_DOG_ONLY_t5m_g4.pdf AUS_DOG_ONLY/model.final.json
smc++ plot -g 6 -c g6/SMCPP_AUS_DOG_ONLY_t5m_g6.pdf AUS_DOG_ONLY/model.final.json




# CENAM
cd ..
vcftools --gzvcf ../smcpp.vcf.gz \
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
--max-missing 1 --recode --out CENAM_DOG_ONLY
bgzip -f CENAM_DOG_ONLY.recode.vcf
tabix CENAM_DOG_ONLY.recode.vcf.gz

# Convert VCF to the SMC++ input format with vcf2smc
for chr in {1..4}; do
  smc++ vcf2smc CENAM_DOG_ONLY.recode.vcf.gz DATA/CENAM_DOG_ONLY.chr${chr}.smc.gz dirofilaria_immitis_chr${chr} CENAM:${CENAM};
done
# each call to vcf2smc processes a single contig. VCFs containing multiple contigs should be processed via multiple separate runs.


# Fit the model using Estimate
## # mutation rate (C.elegans) = 2.7e-9 (https://doi.org/10.1073/pnas.0904895106)
cd t5m
smc++ estimate --timepoints 1 5000000 -o CENAM_DOG_ONLY/ 2.7e-9 ../DATA/CENAM_DOG_ONLY.*.smc.gz

# plot
## Use generation time of D. immitis 1, 2, 4 & 6 yrs
smc++ plot -g 1 -c g1/SMCPP_CENAM_DOG_ONLY_t5m_g1.pdf CENAM_DOG_ONLY/model.final.json
smc++ plot -g 2 -c g2/SMCPP_CENAM_DOG_ONLY_t5m_g2.pdf CENAM_DOG_ONLY/model.final.json
smc++ plot -g 4 -c g4/SMCPP_CENAM_DOG_ONLY_t5m_g4.pdf CENAM_DOG_ONLY/model.final.json
smc++ plot -g 6 -c g6/SMCPP_CENAM_DOG_ONLY_t5m_g6.pdf CENAM_DOG_ONLY/model.final.json




# EUR
cd ..
vcftools --gzvcf ../smcpp.vcf.gz \
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
--indv ITA_PAV_AD_001 \
--max-missing 1 --recode --out EUR_DOG_ONLY
bgzip -f EUR_DOG_ONLY.recode.vcf
tabix EUR_DOG_ONLY.recode.vcf.gz

# Convert VCF to the SMC++ input format with vcf2smc
for chr in {1..4}; do
  smc++ vcf2smc EUR_DOG_ONLY.recode.vcf.gz DATA/EUR_DOG_ONLY.chr${chr}.smc.gz dirofilaria_immitis_chr${chr} EUR:${EUR};
done
# each call to vcf2smc processes a single contig. VCFs containing multiple contigs should be processed via multiple separate runs.

# Fit the model using Estimate
## # mutation rate (C.elegans) = 2.7e-9 (https://doi.org/10.1073/pnas.0904895106)
cd t5m
smc++ estimate --timepoints 1 5000000 -o EUR_DOG_ONLY/ 2.7e-9 ../DATA/EUR_DOG_ONLY.*.smc.gz

# plot
## Use generation time of D. immitis 1, 2, 4 & 6 yrs
smc++ plot -g 1 -c g1/SMCPP_EUR_DOG_ONLY_t5m_g1.pdf EUR_DOG_ONLY/model.final.json
smc++ plot -g 2 -c g2/SMCPP_EUR_DOG_ONLY_t5m_g2.pdf EUR_DOG_ONLY/model.final.json
smc++ plot -g 4 -c g4/SMCPP_EUR_DOG_ONLY_t5m_g4.pdf EUR_DOG_ONLY/model.final.json
smc++ plot -g 6 -c g6/SMCPP_EUR_DOG_ONLY_t5m_g6.pdf EUR_DOG_ONLY/model.final.json


# USA
cd ..
vcftools --gzvcf ../smcpp.vcf.gz \
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
--indv USA_TEX_AD_002 \
--indv USA_TEX_AD_003 \
--indv USA_TEX_AD_004 \
--indv USA_TEX_AD_005 \
--indv USA_TEX_AD_006 \
--indv USA_TEX_AD_007 \
--indv USA_TEX_AD_008 \
--indv USA_TEX_AD_009 \
--indv USA_TEX_AD_010 \
--indv USA_TEX_AD_013 \
--indv USA_TEX_AD_014 \
--max-missing 1 --recode --out USA_DOG_ONLY
bgzip -f USA_DOG_ONLY.recode.vcf
tabix USA_DOG_ONLY.recode.vcf.gz

# Convert VCF to the SMC++ input format with vcf2smc
for chr in {1..4}; do
  smc++ vcf2smc USA_DOG_ONLY.recode.vcf.gz DATA/USA_DOG_ONLY.chr${chr}.smc.gz dirofilaria_immitis_chr${chr} USA:${USA};
done
# each call to vcf2smc processes a single contig. VCFs containing multiple contigs should be processed via multiple separate runs.

# Fit the model using Estimate
## # mutation rate (C.elegans) = 2.7e-9 (https://doi.org/10.1073/pnas.0904895106)
cd t5m
smc++ estimate --timepoints 1 5000000 -o USA_DOG_ONLY/ 2.7e-9 ../DATA/USA_DOG_ONLY.*.smc.gz

# plot
## Use generation time of D. immitis 1, 2, 4 & 6 yrs
smc++ plot -g 1 -c g1/SMCPP_USA_DOG_ONLY_t5m_g1.pdf USA_DOG_ONLY/model.final.json
smc++ plot -g 2 -c g2/SMCPP_USA_DOG_ONLY_t5m_g2.pdf USA_DOG_ONLY/model.final.json
smc++ plot -g 4 -c g4/SMCPP_USA_DOG_ONLY_t5m_g4.pdf USA_DOG_ONLY/model.final.json
smc++ plot -g 6 -c g6/SMCPP_USA_DOG_ONLY_t5m_g6.pdf USA_DOG_ONLY/model.final.json
```
Successfully completed

### Stacked plot for all populations, comparing the various generation times

```R
# smc++ - plotting all populations together
# timepoint = 1 to 5 million years

# load libraries
library(tidyverse)
library(patchwork)
library(ggsci)
library(RColorBrewer)

setwd("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/batch4/FILTER1/NO_OUTGROUPS/smcpp/DOG_ONLY/t5m")

# Directories for each generation time
gen_times <- list(
  "1" = "C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/batch4/FILTER1/NO_OUTGROUPS/smcpp/DOG_ONLY/t5m/g1",
  "2" = "C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/batch4/FILTER1/NO_OUTGROUPS/smcpp/DOG_ONLY/t5m/g2",
  "4" = "C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/batch4/FILTER1/NO_OUTGROUPS/smcpp/DOG_ONLY/t5m/g4",
  "6" = "C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/batch4/FILTER1/NO_OUTGROUPS/smcpp/DOG_ONLY/t5m/g6"
)


data_list <- list()

# get data
for (gen in names(gen_times)) {
  # Set working directory
  setwd(gen_times[[gen]])
  
  # Set filenames
  file_paths <- list(
    "ASIA" = paste0("SMCPP_ASIA_DOG_ONLY_t5m_g", gen, ".csv"),
    "AUS" = paste0("SMCPP_AUS_DOG_ONLY_t5m_g", gen, ".csv"),
    "CENAM" = paste0("SMCPP_CENAM_DOG_ONLY_t5m_g", gen, ".csv"),
    "EUR" = paste0("SMCPP_EUR_DOG_ONLY_t5m_g", gen, ".csv"),
    "USA" = paste0("SMCPP_USA_DOG_ONLY_t5m_g", gen, ".csv")
  )

  data_list[[gen]] <- list()
  
  # Load each file
  for (region in names(file_paths)) {
    file_path <- file_paths[[region]]
      region_data <- read.delim(file_path, header = T, sep = ",")
      region_data$ID <- region
      region_data$Generation <- as.numeric(gen)
      data_list[[gen]][[region]] <- region_data
    }
  }


# Combine data
final_data <- do.call(rbind, lapply(data_list, function(x) do.call(rbind, x)))

pop_colours <- c("ASIA" = "hotpink", "AUS" = "cornflowerblue", "CENAM" = "purple", "EUR" = "forestgreen", "USA" = "tomato2")

# Plot all generation times together for comparison. Do stacked plot.
plot <- ggplot(final_data, aes(x = x, y = y, col = ID)) +
  geom_rect(aes(xmin=25000,ymin=0,xmax=35000,ymax=125000), fill="grey80", col=NA) + # dog domestication
  geom_rect(aes(xmin=116000,ymin=0,xmax=130000,ymax=125000), fill="grey80", col=NA) + # last interglacial
  geom_rect(aes(xmin=58000,ymin=0,xmax=72000,ymax=125000), fill="grey80", col=NA) + # MIS 4
  geom_line(size = 1) +
  guides(color = guide_legend(override.aes = list(size = 10, linewidth=3))) +
  facet_grid(Generation~., scales = "free_y", labeller = labeller(Generation = c("1" = "1 year", "2" = "2 years", "4" = "4 years", "6" = "6 years"))) +
  labs(x = "Years before present", y = expression("Effective population size" ~ (N[e])), col = "Population") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 18),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18)) +
  scale_x_log10(labels = prettyNum) +
  ylim(0,125000) +
  scale_colour_manual(values = pop_colours)

# Save the plot
ggsave("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/batch4/FILTER1/NO_OUTGROUPS/smcpp/DOG_ONLY/t5m/plot_smcpp_t5m_DOG_ONLY.tif", plot, dpi = 300, height = 10, width = 8)
ggsave("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/batch4/FILTER1/NO_OUTGROUPS/smcpp/DOG_ONLY/t5m/plot_smcpp_t5m_DOG_ONLY.png", plot, dpi = 300, height = 10, width = 8)
```



### Split: Timepoint 1 to 5mill generation time, g = 1 and 4 years

The split command fits two-population clean split models (assumes no ongoing gene flow between the two populations after they diverged).


```bash
# load modules
module load PaM/environment
module load bsub.py/0.42.1
module load smcpp/1.15.3-c1
# module load common-apps/htslib/1.9.229 # not available
module load htslib-1.19/perl-5.38.0
module load vcftools/0.1.16-c4

cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/05_ANALYSIS/SMCPP/DOG_ONLY
chmod a+x run_smcpp_split_t5m_g4_DOG_ONLY.sh
bsub.py --queue long 20 run_smcpp_split_t5m_g4_DOG_ONLY "run_smcpp_split_t5m_g4_DOG_ONLY.sh"
```

```bash
#!/bin/bash

cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/05_ANALYSIS/SMCPP/DOG_ONLY
mkdir SPLIT
cd SPLIT

# Create necessary directories
mkdir DATA
mkdir g1 g4
mkdir DATA/ASIA_AUS_DOG_ONLY DATA/ASIA_CENAM_DOG_ONLY DATA/ASIA_EUR_DOG_ONLY DATA/ASIA_USA_DOG_ONLY
mkdir DATA/AUS_CENAM_DOG_ONLY DATA/AUS_EUR_DOG_ONLY DATA/AUS_USA_DOG_ONLY
mkdir DATA/CENAM_EUR_DOG_ONLY DATA/CENAM_USA_DOG_ONLY DATA/EUR_USA_DOG_ONLY

# Get sample names for each population
ASIA=$(awk -F'_' '$1 == "MYS" || $1 == "THA" {print $0}' ../DOG_ONLY_samplelist.keep | paste -sd ',')
AUS=$(awk -F'_' '$1 == "AUS" {print $0}' ../DOG_ONLY_samplelist.keep | paste -sd ',')
CENAM=$(awk -F'_' '$1 == "PAN" || $1 == "CRI" {print $0}' ../DOG_ONLY_samplelist.keep | paste -sd ',')
EUR=$(awk -F'_' '$1 == "GRC" || $1 == "ITA" || $1 == "ROU" {print $0}' ../DOG_ONLY_samplelist.keep | paste -sd ',')
USA=$(awk -F'_' '$1 == "USA" {print $0}' ../DOG_ONLY_samplelist.keep | paste -sd ',')

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
  local out_prefix="./${pop1}_${pop2}_DOG_ONLY"
  local dir_prefix="DATA/${pop1}_${pop2}_DOG_ONLY"
  
  # Generate VCF for the pair
  vcftools --gzvcf ../../smcpp.vcf.gz \
    $(echo ${populations[$pop1]} | tr ',' '\n' | awk '{print "--indv "$0}') \
    $(echo ${populations[$pop2]} | tr ',' '\n' | awk '{print "--indv "$0}') \
    --max-missing 1 --recode --out $out_prefix

  bgzip -f ${out_prefix}.recode.vcf
  tabix ${out_prefix}.recode.vcf.gz

  # create datasets containing the joint frequency spectrum for both pops
  for chr in {1..4}; do
    smc++ vcf2smc ${out_prefix}.recode.vcf.gz ${dir_prefix}/${pop1}_${pop2}_DOG_ONLY.chr${chr}.smc.gz dirofilaria_immitis_chr${chr} ${pop1}:${populations[$pop1]} ${pop2}:${populations[$pop2]}
    smc++ vcf2smc ${out_prefix}.recode.vcf.gz ${dir_prefix}/${pop2}_${pop1}_DOG_ONLY.chr${chr}.smc.gz dirofilaria_immitis_chr${chr} ${pop2}:${populations[$pop2]} ${pop1}:${populations[$pop1]}
  done

  # copy smc++ data for individual pops
  cp ../DATA/${pop1}*.smc.gz ${dir_prefix}
  cp ../DATA/${pop2}*.smc.gz ${dir_prefix}

  # run split
  smc++ split --timepoints 1 5000000 -o ${pop1}_${pop2}_DOG_ONLY ../t5m/${pop1}_DOG_ONLY/model.final.json ../t5m/${pop2}_DOG_ONLY/model.final.json ${dir_prefix}/*.smc.gz
  smc++ plot -g 1 -c g1/SMCPP_${pop1}_${pop2}_t5m_g1_DOG_ONLY.pdf ${pop1}_${pop2}_DOG_ONLY/model.final.json
  smc++ plot -g 4 -c g4/SMCPP_${pop1}_${pop2}_t5m_g4_DOG_ONLY.pdf ${pop1}_${pop2}_DOG_ONLY/model.final.json
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
Successfully completed


### Plot the split times for g=4 years

```R
# smc++ Split - combine all split dates into a single csv file
# timepoint = 1 to 5 million years, g = 4
## DOG SAMPLES ONLY

# load libraries
library(tidyverse)
library(patchwork)
library(ggsci)
library(RColorBrewer)

setwd("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/batch4/FILTER1/NO_OUTGROUPS/smcpp/DOG_ONLY/t5m/g4/SPLIT")

# Read all datasets
ASIA_AUS_data <- read.delim("SMCPP_ASIA_AUS_t5m_g4_DOG_ONLY.csv", header = TRUE, sep = ",")
ASIA_CENAM_data <- read.delim("SMCPP_ASIA_CENAM_t5m_g4_DOG_ONLY.csv", header = TRUE, sep = ",")
ASIA_EUR_data <- read.delim("SMCPP_ASIA_EUR_t5m_g4_DOG_ONLY.csv", header = TRUE, sep = ",")
ASIA_USA_data <- read.delim("SMCPP_ASIA_USA_t5m_g4_DOG_ONLY.csv", header = TRUE, sep = ",")
AUS_CENAM_data <- read.delim("SMCPP_AUS_CENAM_t5m_g4_DOG_ONLY.csv", header = TRUE, sep = ",")
AUS_EUR_data <- read.delim("SMCPP_AUS_EUR_t5m_g4_DOG_ONLY.csv", header = TRUE, sep = ",")
AUS_USA_data <- read.delim("SMCPP_AUS_USA_t5m_g4_DOG_ONLY.csv", header = TRUE, sep = ",")
CENAM_EUR_data <- read.delim("SMCPP_CENAM_EUR_t5m_g4_DOG_ONLY.csv", header = TRUE, sep = ",")
CENAM_USA_data <- read.delim("SMCPP_CENAM_USA_t5m_g4_DOG_ONLY.csv", header = TRUE, sep = ",")
EUR_USA_data <- read.delim("SMCPP_EUR_USA_t5m_g4_DOG_ONLY.csv", header = TRUE, sep = ",")

# Add pair identifiers and subset IDs
datasets <- list(
  ASIA_AUS = ASIA_AUS_data,
  ASIA_CENAM = ASIA_CENAM_data,
  ASIA_EUR = ASIA_EUR_data,
  ASIA_USA = ASIA_USA_data,
  AUS_CENAM = AUS_CENAM_data,
  AUS_EUR = AUS_EUR_data,
  AUS_USA = AUS_USA_data,
  CENAM_EUR = CENAM_EUR_data,
  CENAM_USA = CENAM_USA_data,
  EUR_USA = EUR_USA_data
)

for (name in names(datasets)) {
  datasets[[name]] <- datasets[[name]] %>%
    mutate(Pair = name, subset_id = name)
}

# Combine all datasets
all_data <- bind_rows(datasets)

get_second_pop_end_point <- function(data) {
  switch_point <- which(data$label != data$label[1])[1]
  second_pop_data <- data[switch_point:nrow(data), ]
  end_x <- max(second_pop_data$x, na.rm = TRUE)
  return(end_x)
}

# Apply the function to each dataset and store the results
end_points <- data.frame(
  subset_id = names(datasets),
  end_x = sapply(datasets, get_second_pop_end_point)
)

write.csv(end_points, "split_dates_g4.csv", row.names=FALSE)
```

```R
# smc++ - plotting all populations together
# timepoint = 1 to 5 million years

# load libraries
library(tidyverse)
library(patchwork)
library(ggsci)
library(RColorBrewer)

setwd("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/batch4/FILTER1/NO_OUTGROUPS/smcpp/DOG_ONLY/t5m")

# Directories for each generation time
gen_times <- list(
  "1" = "C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/batch4/FILTER1/NO_OUTGROUPS/smcpp/DOG_ONLY/t5m/g1",
  "2" = "C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/batch4/FILTER1/NO_OUTGROUPS/smcpp/DOG_ONLY/t5m/g2",
  "4" = "C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/batch4/FILTER1/NO_OUTGROUPS/smcpp/DOG_ONLY/t5m/g4",
  "6" = "C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/batch4/FILTER1/NO_OUTGROUPS/smcpp/DOG_ONLY/t5m/g6"
)


data_list <- list()

# get data
for (gen in names(gen_times)) {
  # Set working directory
  setwd(gen_times[[gen]])
  
  # Set filenames
  file_paths <- list(
    "ASIA" = paste0("SMCPP_ASIA_DOG_ONLY_t5m_g", gen, ".csv"),
    "AUS" = paste0("SMCPP_AUS_DOG_ONLY_t5m_g", gen, ".csv"),
    "CENAM" = paste0("SMCPP_CENAM_DOG_ONLY_t5m_g", gen, ".csv"),
    "EUR" = paste0("SMCPP_EUR_DOG_ONLY_t5m_g", gen, ".csv"),
    "USA" = paste0("SMCPP_USA_DOG_ONLY_t5m_g", gen, ".csv")
  )

  data_list[[gen]] <- list()
  
  # Load each file
  for (region in names(file_paths)) {
    file_path <- file_paths[[region]]
      region_data <- read.delim(file_path, header = T, sep = ",")
      region_data$ID <- region
      region_data$Generation <- as.numeric(gen)
      data_list[[gen]][[region]] <- region_data
    }
  }


# Combine data
final_data <- do.call(rbind, lapply(data_list, function(x) do.call(rbind, x)))

pop_colours <- c("ASIA" = "hotpink", "AUS" = "cornflowerblue", "CENAM" = "purple", "EUR" = "forestgreen", "USA" = "tomato2")

# Plot all generation times together for comparison
plot <- ggplot(final_data, aes(x = x, y = y, col = ID)) +
  geom_rect(aes(xmin=25000,ymin=0,xmax=35000,ymax=125000), fill="grey80", col=NA) + # dog domestication
  geom_rect(aes(xmin=116000,ymin=0,xmax=130000,ymax=125000), fill="peachpuff1", col=NA) + # last interglacial
  geom_rect(aes(xmin=58000,ymin=0,xmax=72000,ymax=125000), fill="lightblue", col=NA) + # MIS 4
  geom_line(size = 1) +
  guides(color = guide_legend(override.aes = list(size = 8, linewidth=6))) +
  facet_grid(Generation~., scales = "free_y", labeller = labeller(Generation = c("1" = "1 year", "2" = "2 years", "4" = "4 years", "6" = "6 years"))) +
  labs(x = "Years before present", y = expression("Effective population size" ~ (N[e])), col = "Population") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 18),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18)) +
  scale_x_log10(labels = prettyNum) +
  ylim(0,125000) +
  scale_colour_manual(values = pop_colours)

# Save the plot
ggsave("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/batch4/FILTER1/NO_OUTGROUPS/smcpp/DOG_ONLY/t5m/plot_smcpp_t5m_DOG_ONLY.tif", plot, dpi = 300, height = 10, width = 8)
ggsave("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/batch4/FILTER1/NO_OUTGROUPS/smcpp/DOG_ONLY/t5m/plot_smcpp_t5m_DOG_ONLY.png", plot, dpi = 300, height = 10, width = 8)
```









## Run smc++ for all samples by HOST

### timepoint = 1 to 5 mill, generation time = 1, 2, 4 & 6 years

```bash
# load modules
module load PaM/environment
module load bsub.py/0.42.1
module load smcpp/1.15.3-c1
# module load common-apps/htslib/1.9.229 # not available
module load htslib-1.19/perl-5.38.0
module load vcftools/0.1.16-c4

cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/05_ANALYSIS/SMCPP/HOST
chmod a+x run_smcpp_t5m_HOST.sh
bsub.py --queue long 20 run_smcpp_t5m_HOST "run_smcpp_t5m_HOST.sh"
```

```bash
#!/bin/bash

# Load modules
module load PaM/environment
module load smcpp/1.15.3-c1
# module load common-apps/htslib/1.9.229 # not available
module load htslib-1.19/perl-5.38.0
module load vcftools/0.1.16-c4

cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/05_ANALYSIS/SMCPP/HOST

# Get sample names for each population
DOG=$(paste -sd ',' DOG_ONLY_samplelist.keep)
FOX="AUS_BNE_AD_006,ITA_NEA_AD_002,ITA_NEA_AD_003"
CAT="CRI_SJO_AD_001,USA_TEX_AD_011,USA_TEX_AD_012,USA_TEX_AD_015"
LEOPARD="ROU_BUC_AD_001"
WILDCAT="ROU_COM_AD_001"
JACKAL="ROU_GIU_AD_001"
FERRET="USA_TEX_AD_001"
declare -A populations
populations=(
  [DOG]=$DOG
  [FOX]=$FOX
  [CAT]=$CAT
  [LEOPARD]=$LEOPARD
  [WILDCAT]=$WILDCAT
  [JACKAL]=$JACKAL
  [FERRET]=$FERRET
)
# DOG_ONLY_samplelist.keep - is from when I was running pixy - it's basically all the samples in the final VCF for dog hosts only and minus repeat samples.

mkdir DATA

# DOG
vcftools --gzvcf ../smcpp.vcf.gz \
--indv AUS_BNE_AD_001 \
--indv AUS_BNE_AD_002 \
--indv AUS_BNE_AD_003 \
--indv AUS_BNE_AD_004 \
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
--indv ITA_PAV_AD_001 \
--indv MYS_SEL_AD_001 \
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
--indv THA_BKK_AD_001 \
--indv THA_BKK_AD_002 \
--indv THA_BKK_AD_003 \
--indv THA_BKK_AD_004 \
--indv THA_BKK_AD_005 \
--indv THA_BKK_AD_006 \
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
--indv USA_TEX_AD_002 \
--indv USA_TEX_AD_003 \
--indv USA_TEX_AD_004 \
--indv USA_TEX_AD_005 \
--indv USA_TEX_AD_006 \
--indv USA_TEX_AD_007 \
--indv USA_TEX_AD_008 \
--indv USA_TEX_AD_009 \
--indv USA_TEX_AD_010 \
--indv USA_TEX_AD_013 \
--indv USA_TEX_AD_014 \
--max-missing 1 --recode --out DOG
bgzip -f DOG.recode.vcf
tabix DOG.recode.vcf.gz

# Convert VCF to the SMC++ input format with vcf2smc
for chr in {1..4}; do
  smc++ vcf2smc DOG.recode.vcf.gz DATA/DOG.chr${chr}.smc.gz dirofilaria_immitis_chr${chr} DOG:${DOG};
done
# each call to vcf2smc processes a single contig. VCFs containing multiple contigs should be processed via multiple separate runs.

# Fit the model using Estimate
## # mutation rate (C.elegans) = 2.7e-9 (https://doi.org/10.1073/pnas.0904895106)
mkdir t5m
cd t5m
mkdir g1 g2 g4 g6

smc++ estimate --timepoints 1 5000000 -o DOG/ 2.7e-9 ../DATA/DOG.*.smc.gz

# plot
## Use generation time of D. immitis 1, 2, 4 & 6 yrs
smc++ plot -g 1 -c g1/SMCPP_DOG_t5m_g1.pdf DOG/model.final.json
smc++ plot -g 2 -c g2/SMCPP_DOG_t5m_g2.pdf DOG/model.final.json
smc++ plot -g 4 -c g4/SMCPP_DOG_t5m_g4.pdf DOG/model.final.json
smc++ plot -g 6 -c g6/SMCPP_DOG_t5m_g6.pdf DOG/model.final.json


# FOX
cd ..
vcftools --gzvcf ../smcpp.vcf.gz \
--indv AUS_BNE_AD_006 \
--indv ITA_NEA_AD_002 \
--indv ITA_NEA_AD_003 \
--max-missing 1 --recode --out FOX
bgzip -f FOX.recode.vcf
tabix FOX.recode.vcf.gz

# Convert VCF to the SMC++ input format with vcf2smc
for chr in {1..4}; do
  smc++ vcf2smc FOX.recode.vcf.gz DATA/FOX.chr${chr}.smc.gz dirofilaria_immitis_chr${chr} FOX:${FOX};
done
# each call to vcf2smc processes a single contig. VCFs containing multiple contigs should be processed via multiple separate runs.


# Fit the model using Estimate
## # mutation rate (C.elegans) = 2.7e-9 (https://doi.org/10.1073/pnas.0904895106)
cd t5m
smc++ estimate --timepoints 1 5000000 -o FOX/ 2.7e-9 ../DATA/FOX.*.smc.gz

# plot
## Use generation time of D. immitis 1, 2, 4 & 6 yrs
smc++ plot -g 1 -c g1/SMCPP_FOX_t5m_g1.pdf FOX/model.final.json
smc++ plot -g 2 -c g2/SMCPP_FOX_t5m_g2.pdf FOX/model.final.json
smc++ plot -g 4 -c g4/SMCPP_FOX_t5m_g4.pdf FOX/model.final.json
smc++ plot -g 6 -c g6/SMCPP_FOX_t5m_g6.pdf FOX/model.final.json




# CAT
cd ..
vcftools --gzvcf ../smcpp.vcf.gz \
--indv CRI_SJO_AD_001 \
--indv USA_TEX_AD_011 \
--indv USA_TEX_AD_012 \
--indv USA_TEX_AD_015 \
--max-missing 1 --recode --out CAT
bgzip -f CAT.recode.vcf
tabix CAT.recode.vcf.gz

# Convert VCF to the SMC++ input format with vcf2smc
for chr in {1..4}; do
  smc++ vcf2smc CAT.recode.vcf.gz DATA/CAT.chr${chr}.smc.gz dirofilaria_immitis_chr${chr} CAT:${CAT};
done
# each call to vcf2smc processes a single contig. VCFs containing multiple contigs should be processed via multiple separate runs.


# Fit the model using Estimate
## # mutation rate (C.elegans) = 2.7e-9 (https://doi.org/10.1073/pnas.0904895106)
cd t5m
smc++ estimate --timepoints 1 5000000 -o CAT/ 2.7e-9 ../DATA/CAT.*.smc.gz

# plot
## Use generation time of D. immitis 1, 2, 4 & 6 yrs
smc++ plot -g 1 -c g1/SMCPP_CAT_t5m_g1.pdf CAT/model.final.json
smc++ plot -g 2 -c g2/SMCPP_CAT_t5m_g2.pdf CAT/model.final.json
smc++ plot -g 4 -c g4/SMCPP_CAT_t5m_g4.pdf CAT/model.final.json
smc++ plot -g 6 -c g6/SMCPP_CAT_t5m_g6.pdf CAT/model.final.json




# LEOPARD
cd ..
vcftools --gzvcf ../smcpp.vcf.gz \
--indv ROU_BUC_AD_001 \
--max-missing 1 --recode --out LEOPARD
bgzip -f LEOPARD.recode.vcf
tabix LEOPARD.recode.vcf.gz

# Convert VCF to the SMC++ input format with vcf2smc
for chr in {1..4}; do
  smc++ vcf2smc LEOPARD.recode.vcf.gz DATA/LEOPARD.chr${chr}.smc.gz dirofilaria_immitis_chr${chr} LEOPARD:${LEOPARD};
done
# each call to vcf2smc processes a single contig. VCFs containing multiple contigs should be processed via multiple separate runs.

# Fit the model using Estimate
## # mutation rate (C.elegans) = 2.7e-9 (https://doi.org/10.1073/pnas.0904895106)
cd t5m
smc++ estimate --timepoints 1 5000000 -o LEOPARD/ 2.7e-9 ../DATA/LEOPARD.*.smc.gz

# plot
## Use generation time of D. immitis 1, 2, 4 & 6 yrs
smc++ plot -g 1 -c g1/SMCPP_LEOPARD_t5m_g1.pdf LEOPARD/model.final.json
smc++ plot -g 2 -c g2/SMCPP_LEOPARD_t5m_g2.pdf LEOPARD/model.final.json
smc++ plot -g 4 -c g4/SMCPP_LEOPARD_t5m_g4.pdf LEOPARD/model.final.json
smc++ plot -g 6 -c g6/SMCPP_LEOPARD_t5m_g6.pdf LEOPARD/model.final.json


# WILDCAT
cd ..
vcftools --gzvcf ../smcpp.vcf.gz \
--indv ROU_COM_AD_001 \
--max-missing 1 --recode --out WILDCAT
bgzip -f WILDCAT.recode.vcf
tabix WILDCAT.recode.vcf.gz

# Convert VCF to the SMC++ input format with vcf2smc
for chr in {1..4}; do
  smc++ vcf2smc WILDCAT.recode.vcf.gz DATA/WILDCAT.chr${chr}.smc.gz dirofilaria_immitis_chr${chr} WILDCAT:${WILDCAT};
done
# each call to vcf2smc processes a single contig. VCFs containing multiple contigs should be processed via multiple separate runs.

# Fit the model using Estimate
## # mutation rate (C.elegans) = 2.7e-9 (https://doi.org/10.1073/pnas.0904895106)
cd t5m
smc++ estimate --timepoints 1 5000000 -o WILDCAT/ 2.7e-9 ../DATA/WILDCAT.*.smc.gz

# plot
## Use generation time of D. immitis 1, 2, 4 & 6 yrs
smc++ plot -g 1 -c g1/SMCPP_WILDCAT_t5m_g1.pdf WILDCAT/model.final.json
smc++ plot -g 2 -c g2/SMCPP_WILDCAT_t5m_g2.pdf WILDCAT/model.final.json
smc++ plot -g 4 -c g4/SMCPP_WILDCAT_t5m_g4.pdf WILDCAT/model.final.json
smc++ plot -g 6 -c g6/SMCPP_WILDCAT_t5m_g6.pdf WILDCAT/model.final.json


# JACKAL
cd ..
vcftools --gzvcf ../smcpp.vcf.gz \
--indv ROU_GIU_AD_001 \
--max-missing 1 --recode --out JACKAL
bgzip -f JACKAL.recode.vcf
tabix JACKAL.recode.vcf.gz

# Convert VCF to the SMC++ input format with vcf2smc
for chr in {1..4}; do
  smc++ vcf2smc JACKAL.recode.vcf.gz DATA/JACKAL.chr${chr}.smc.gz dirofilaria_immitis_chr${chr} JACKAL:${JACKAL};
done
# each call to vcf2smc processes a single contig. VCFs containing multiple contigs should be processed via multiple separate runs.

# Fit the model using Estimate
## # mutation rate (C.elegans) = 2.7e-9 (https://doi.org/10.1073/pnas.0904895106)
cd t5m
smc++ estimate --timepoints 1 5000000 -o JACKAL/ 2.7e-9 ../DATA/JACKAL.*.smc.gz

# plot
## Use generation time of D. immitis 1, 2, 4 & 6 yrs
smc++ plot -g 1 -c g1/SMCPP_JACKAL_t5m_g1.pdf JACKAL/model.final.json
smc++ plot -g 2 -c g2/SMCPP_JACKAL_t5m_g2.pdf JACKAL/model.final.json
smc++ plot -g 4 -c g4/SMCPP_JACKAL_t5m_g4.pdf JACKAL/model.final.json
smc++ plot -g 6 -c g6/SMCPP_JACKAL_t5m_g6.pdf JACKAL/model.final.json


# FERRET
cd ..
vcftools --gzvcf ../smcpp.vcf.gz \
--indv USA_TEX_AD_001 \
--max-missing 1 --recode --out FERRET
bgzip -f FERRET.recode.vcf
tabix FERRET.recode.vcf.gz

# Convert VCF to the SMC++ input format with vcf2smc
for chr in {1..4}; do
  smc++ vcf2smc FERRET.recode.vcf.gz DATA/FERRET.chr${chr}.smc.gz dirofilaria_immitis_chr${chr} FERRET:${FERRET};
done
# each call to vcf2smc processes a single contig. VCFs containing multiple contigs should be processed via multiple separate runs.

# Fit the model using Estimate
## # mutation rate (C.elegans) = 2.7e-9 (https://doi.org/10.1073/pnas.0904895106)
cd t5m
smc++ estimate --timepoints 1 5000000 -o FERRET/ 2.7e-9 ../DATA/FERRET.*.smc.gz

# plot
## Use generation time of D. immitis 1, 2, 4 & 6 yrs
smc++ plot -g 1 -c g1/SMCPP_FERRET_t5m_g1.pdf FERRET/model.final.json
smc++ plot -g 2 -c g2/SMCPP_FERRET_t5m_g2.pdf FERRET/model.final.json
smc++ plot -g 4 -c g4/SMCPP_FERRET_t5m_g4.pdf FERRET/model.final.json
smc++ plot -g 6 -c g6/SMCPP_FERRET_t5m_g6.pdf FERRET/model.final.json
```
Successfully completed

### Stacked plot for all hosts, comparing the various generation times

```R
# smc++ - plotting all hosts together
# timepoint = 1 to 5 million years

# load libraries
library(tidyverse)
library(patchwork)
library(ggsci)
library(RColorBrewer)

setwd("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/batch4/FILTER1/NO_OUTGROUPS/smcpp/HOST/t5m")

# Directories for each generation time
gen_times <- list(
  "1" = "C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/batch4/FILTER1/NO_OUTGROUPS/smcpp/HOST/t5m/g1",
  "2" = "C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/batch4/FILTER1/NO_OUTGROUPS/smcpp/HOST/t5m/g2",
  "4" = "C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/batch4/FILTER1/NO_OUTGROUPS/smcpp/HOST/t5m/g4",
  "6" = "C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/batch4/FILTER1/NO_OUTGROUPS/smcpp/HOST/t5m/g6"
)


data_list <- list()

# get data
for (gen in names(gen_times)) {
  # Set working directory
  setwd(gen_times[[gen]])
  
  # Set filenames
  file_paths <- list(
    "DOG" = paste0("SMCPP_DOG_t5m_g", gen, ".csv"),
    "FOX" = paste0("SMCPP_FOX_t5m_g", gen, ".csv"),
    "CAT" = paste0("SMCPP_CAT_t5m_g", gen, ".csv"),
    "LEOPARD" = paste0("SMCPP_LEOPARD_t5m_g", gen, ".csv"),
    "WILDCAT" = paste0("SMCPP_WILDCAT_t5m_g", gen, ".csv"),
    "JACKAL" = paste0("SMCPP_JACKAL_t5m_g", gen, ".csv"),
    "FERRET" = paste0("SMCPP_FERRET_t5m_g", gen, ".csv")
  )

  data_list[[gen]] <- list()
  
  # Load each file
  for (region in names(file_paths)) {
    file_path <- file_paths[[region]]
      region_data <- read.delim(file_path, header = T, sep = ",")
      region_data$ID <- region
      region_data$Generation <- as.numeric(gen)
      data_list[[gen]][[region]] <- region_data
    }
  }


# Combine data
final_data <- do.call(rbind, lapply(data_list, function(x) do.call(rbind, x)))

# Plot all generation times together for comparison
plot <- ggplot(final_data, aes(x = x, y = y, col = ID)) +
  geom_rect(aes(xmin=25000,ymin=0,xmax=35000,ymax=500000), fill="grey80", col=NA) + # dog domestication
  geom_rect(aes(xmin=116000,ymin=0,xmax=130000,ymax=500000), fill="peachpuff1", col=NA) + # last interglacial
  geom_rect(aes(xmin=58000,ymin=0,xmax=72000,ymax=500000), fill="lightblue", col=NA) + # MIS 4
  geom_line(size = 1) +
  guides(color = guide_legend(override.aes = list(size = 8, linewidth=6))) +
  facet_grid(Generation~., scales = "free_y", labeller = labeller(Generation = c("1" = "1 year", "2" = "2 years", "4" = "4 years", "6" = "6 years"))) +
  labs(x = "Years before present", y = expression("Effective population size" ~ (N[e])), col = "Host") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 18),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18)) +
  scale_x_log10(labels = prettyNum) +
  #ylim(0,250000) +
  scale_color_jama()
plot

# Save the plot
ggsave("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/batch4/FILTER1/NO_OUTGROUPS/smcpp/HOST/t5m/plot_smcpp_t5m_HOST.tif", plot, dpi = 300, height = 10, width = 8)
ggsave("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/batch4/FILTER1/NO_OUTGROUPS/smcpp/HOST/t5m/plot_smcpp_t5m_HOST.png", plot, dpi = 300, height = 10, width = 8)
```