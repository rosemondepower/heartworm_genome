# LD - other code not used

Only using --geno-r2 seemed to work. When I tried using --hap-r2 for phased data, it said no SNPs were left after filtering. Illumina HiSeq produces unphased data by default, so it makes sense why --hap-r2 didn't work. To get phased data, additional steps are needed (usually long-read sequencing or specialised tools).

Try using various window sizes. Window size of 50 kbp would only compares sites within 50,000 bp of one another. NB: genome size of D. immitis is ~88,000,000. 

- 100 kbp window
- 50 kbp window
- 10 kbp window

Calculate LD for each group separately with different windows.

```bash
module load bsub.py/0.42.1
bsub.py 10 LD "LD.sh"
```

```bash
#!/bin/bash

# load modules
module load vcftools/0.1.16-c4

cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/05_ANALYSIS/LD

ln -s /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/04_VARIANTS/FILTER1/NO_OUTGROUPS/FINAL_SETS/nuclear_samples3x_missing0.9.chr1to4.recode.vcf

VCF=nuclear_samples3x_missing0.9.chr1to4.recode.vcf

# loop through each continent
# removed replicate samples

## 100 kb window
for i in ASIA AUS CENAM EUR USA; do
  vcftools --vcf ${VCF} \
    --geno-r2 \
    --ld-window-bp 100000 \
    --keep ${i}_samplelist.keep \
    --out LD_100000_${i} ;
done

## 50 kb window
for i in ASIA AUS CENAM EUR USA; do
  vcftools --vcf ${VCF} \
    --geno-r2 \
    --ld-window-bp 50000 \
    --keep ${i}_samplelist.keep \
    --out LD_50000_${i} ;
done


## 10 kb window
for i in ASIA AUS CENAM EUR USA; do
  vcftools --vcf ${VCF} \
    --geno-r2 \
    --ld-window-bp 10000 \
    --keep ${i}_samplelist.keep \
    --out LD_10000_${i} ;
done
```


## Plot LD

LD files are really big, run R on the farm.

Plot 50kb window:

```bash
module load cellgen/R/4.3.1
bsub.py 20 plot_LD_50000 "Rscript plot_LD_50000.R"
```

```R
version

setwd("/lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/05_ANALYSIS/LD")

library(ggplot2)
library(dplyr)



##########################################################
# plot 50 kbp window
##########################################################

# Function to process LD data
process_ld_data <- function(file_path) {
  # Read data
  ld_data <- read.table(file_path, header = TRUE)
  
  # Make distance column by calculating absolute distance between SNP1 and SNP2, converts this distance into kb
  ld_data <- ld_data %>% mutate(Distance = abs(POS2 - POS1) / 1000)
  
   # Separate into blocks
  ld_data_blocks <- ld_data %>%
    mutate(Bin = cut(Distance, breaks = seq(0, max(Distance) + 0.1, by = 0.1), include.lowest = TRUE)) %>%
    group_by(Bin) %>%
    summarize(Mean_R.2 = mean(R.2, na.rm = TRUE),
              Mean_Distance = mean(Distance))
  #make the x-axis labels the midpoint of each bin
  # each bin is 100bp or 0.1 kb

  # Add cohort identifier
  ld_data_blocks$Cohort <- gsub("LD_50000_|\\.geno\\.ld", "", basename(file_path))
  
  return(ld_data_blocks)
}

# File paths for the five cohorts
cohort_files <- c("LD_50000_USA.geno.ld", "LD_50000_EUR.geno.ld", "LD_50000_CENAM.geno.ld", "LD_50000_AUS.geno.ld", "LD_50000_ASIA.geno.ld")

# Process each cohort
cohort_data <- lapply(cohort_files, process_ld_data)

# Combine the results into a single data frame
all_data <- bind_rows(cohort_data)

all_data <- all_data %>%
  mutate(Location = sub("ld_50000_", "", Cohort))

# Define custom colors for each cohort
location_colours <- c("USA" = "firebrick3", "EUR" = "forestgreen", "CENAM" = "purple", "AUS" = "dodgerblue", "ASIA" = "hotpink")

# Plot
plot_50000 <- ggplot(all_data, aes(x = Mean_Distance, y = Mean_R.2, group = Location, color = Location)) +
  geom_point(alpha = 0) +  # Adjust alpha for point visibility
  geom_smooth(method = "loess", se = FALSE) +  # Optional: Add a smooth line for better visualization
  labs(x = "Distance (kb)", y = expression(R^2), title = "Linkage Disequilibrium", subtitle = "Window
  : 50,000 bp") +
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.key = element_rect(fill = "white")) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.05)) +
  scale_color_manual(values = location_colours)

plot_50000

# Save the plot
ggsave("plot_LD_50000.pdf", plot_50000, dpi = 300, height = 6, width = 8)
```

Plot 10 kb window:

```bash
bsub.py 20 plot_LD_10000 "Rscript plot_LD_10000.R"
```

```R
version

setwd("/lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/05_ANALYSIS/LD")

library(ggplot2)
library(dplyr)



##########################################################
# plot 10 kbp window
##########################################################

# Function to process LD data
process_ld_data <- function(file_path) {
  # Read data
  ld_data <- read.table(file_path, header = TRUE)
  
  # Make distance column by calculating absolute distance between SNP1 and SNP2, converts this distance into kb
  ld_data <- ld_data %>% mutate(Distance = abs(POS2 - POS1) / 1000)
  
   # Separate into blocks
  ld_data_blocks <- ld_data %>%
    mutate(Bin = cut(Distance, breaks = seq(0, max(Distance) + 0.1, by = 0.1), include.lowest = TRUE)) %>%
    group_by(Bin) %>%
    summarize(Mean_R.2 = mean(R.2, na.rm = TRUE),
              Mean_Distance = mean(Distance))
  #make the x-axis labels the midpoint of each bin
  # each bin is 100bp or 0.1 kb

  # Add cohort identifier
  ld_data_blocks$Cohort <- gsub("LD_10000_|\\.geno\\.ld", "", basename(file_path))
  
  return(ld_data_blocks)
}

# File paths for the five cohorts
cohort_files <- c("LD_10000_USA.geno.ld", "LD_10000_EUR.geno.ld", "LD_10000_CENAM.geno.ld", "LD_10000_AUS.geno.ld", "LD_10000_ASIA.geno.ld")

# Process each cohort
cohort_data <- lapply(cohort_files, process_ld_data)

# Combine the results into a single data frame
all_data <- bind_rows(cohort_data)

all_data <- all_data %>%
  mutate(Location = sub("ld_10000_", "", Cohort))

# Define custom colors for each cohort
location_colours <- c("USA" = "firebrick3", "EUR" = "forestgreen", "CENAM" = "purple", "AUS" = "dodgerblue", "ASIA" = "hotpink")

# Plot
plot_10000 <- ggplot(all_data, aes(x = Mean_Distance, y = Mean_R.2, group = Location, color = Location)) +
  geom_point(alpha = 0) +  # Adjust alpha for point visibility
  geom_smooth(method = "loess", se = FALSE) +  # Optional: Add a smooth line for better visualization
  labs(x = "Distance (kb)", y = expression(R^2), title = "Linkage Disequilibrium", subtitle = "Window
  : 10,000 bp") +
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.key = element_rect(fill = "white")) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.05)) +
  scale_color_manual(values = location_colours)

plot_10000

# Save the plot
ggsave("plot_LD_10000.pdf", plot_10000, dpi = 300, height = 6, width = 8)
```

Plot 100 kb window:

```bash
bsub.py 20 plot_LD_100000 "Rscript plot_LD_100000.R"
```

```R
version

setwd("/lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/05_ANALYSIS/LD")

library(ggplot2)
library(dplyr)



##########################################################
# plot 100 kbp window
##########################################################

# Function to process LD data
process_ld_data <- function(file_path) {
  # Read data
  ld_data <- read.table(file_path, header = TRUE)
  
  # Make distance column by calculating absolute distance between SNP1 and SNP2, converts this distance into kb
  ld_data <- ld_data %>% mutate(Distance = abs(POS2 - POS1) / 1000)
  
   # Separate into blocks
  ld_data_blocks <- ld_data %>%
    mutate(Bin = cut(Distance, breaks = seq(0, max(Distance) + 0.1, by = 0.1), include.lowest = TRUE)) %>%
    group_by(Bin) %>%
    summarize(Mean_R.2 = mean(R.2, na.rm = TRUE),
              Mean_Distance = mean(Distance))
  #make the x-axis labels the midpoint of each bin
  # each bin is 100bp or 0.1 kb

  # Add cohort identifier
  ld_data_blocks$Cohort <- gsub(".geno.ld", "", basename(file_path))
  
  return(ld_data)
}

# File paths for the five cohorts
cohort_files <- c("LD_100000_USA.geno.ld", "LD_100000_EUR.geno.ld", "LD_100000_CENAM.geno.ld", "LD_100000_AUS.geno.ld", "LD_100000_ASIA.geno.ld")

# Process each cohort
cohort_data <- lapply(cohort_files, process_ld_data)

# Combine the results into a single data frame
all_data <- bind_rows(cohort_data)

all_data <- all_data %>%
  mutate(Location = sub("ld_100000_", "", Cohort))

# Define custom colors for each cohort
location_colours <- c("USA" = "firebrick3", "EUR" = "forestgreen", "CENAM" = "purple", "AUS" = "dodgerblue", "ASIA" = "hotpink")

# Plot
plot_100000 <- ggplot(all_data, aes(x = Mean_Distance, y = Mean_R.2, group = Location, color = Location)) +
  geom_point(alpha = 0) +  # Adjust alpha for point visibility
  geom_smooth(method = "loess", se = FALSE) +  # Optional: Add a smooth line for better visualization
  labs(x = "Distance (kb)", y = expression(R^2), title = "Linkage Disequilibrium", subtitle = "Window
  : 100,000 bp") +
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.key = element_rect(fill = "white")) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.05)) +
  scale_color_manual(values = location_colours)

plot_100000

# Save the plot
ggsave("plot_LD_100000.pdf", plot_100000, dpi = 300, height = 6, width = 8)
```

Terminated - 100kb used too much memory



## SCRATCH

## LD using plink

### Get VCF per population

```bash
ln -s /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/05_ANALYSIS/SMCPP/AUS.recode.vcf.gz
```

### Get R2

Followed the steps from: https://speciationgenomics.github.io/ld_decay/

```bash
bsub.py 20 plink_LD.sh "plink_LD.sh"
```

```bash
#!/bin/bash

module load common-apps/htslib/1.9.229
module load vcftools/0.1.16-c4
module load plink/1.90b6.18--h516909a_0

VCF=/lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/05_ANALYSIS/LD/nuclear_samples3x_missing0.9.chr1to4.recode.vcf

vcftools --vcf ${VCF} \
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

VCF_AUS=AUS.recode.vcf.gz

plink --vcf ${VCF_AUS} --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--r2 gz --ld-window-kb 100 \
--ld-window-r2 0 \
--make-bed --out plink_LD
```
--r2 gz calculates r2 values for pairs of SNPs, output as compressed gz file
--ld-window 100 sets the number of SNPs in the sliding window
--ld-window-kb 1000 sets the size of the sliding window
--ld-window-r2 0 includes all pairs of SNPs in the output

### Calculating average LD across set distances

```bash
chmod a+x ld_decay_calc.py
bsub.py 20 ld_decay_calc "./ld_decay_calc.py -i plink_LD.ld.gz -o plink_LD"
```
Successfully completed.

Produced 2 files:
- pairwise SNPs arranged into distance classes
- average LD across 100 kb bins --> we'll take this and plot it in R to see the LD decay



### Decay plot

```R
## LD decay plotting script

rm(list = ls())

library(tidyverse)

# set path
my_bins <- "./plink_LD.ld_decay_bins"

# read in data
ld_bins <- read_tsv(my_bins)

# plot LD decay
ggplot(ld_bins, aes(distance, avg_R2)) + geom_line() +
  xlab("Distance (bp)") + ylab(expression(italic(r)^2))
```

