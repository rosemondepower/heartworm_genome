# Linkage disequilibrium

## Entire cohort

### get linkage disequilibrium statistics using vcftools

Only using --geno-r2 seemed to work. When I tried using --hap-r2 for phased data (I'm pretty sure mine are phased because the genotypes are like 0/0 etc?), it said no SNPs were left after filtering. But I was also encountering warnings re the PGT and PID format entries...

Actually - Illumina HiSeq produces unphased data by default, so it makes sense why --hap-r2 didn't work. To get phased data, additional steps are needed (usually long-read sequencing or specialised tools).





## Separate cohort into groups

Calculate LD for each group separately. Go with 10,000 bp window.

bsub.py 10 LD_vcftools_50000 "../LD_vcftools_50000.sh"

```bash
module load vcftools/0.1.16-c4

cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/BATCH3/02_VARIANTS/LD

vcftools --vcf ../FINAL_SETS/nuclear_samples3x_missing0.9.chr1to4.recode.RENAMED.vcf \
         --geno-r2 \
         --ld-window-bp 50000 \
         --keep USA_samples.txt \
         --out ld_50000_USA

vcftools --vcf ../FINAL_SETS/nuclear_samples3x_missing0.9.chr1to4.recode.RENAMED.vcf \
         --geno-r2 \
         --ld-window-bp 50000 \
         --keep CAM_samples.txt \
         --out ld_50000_CAM

vcftools --vcf ../FINAL_SETS/nuclear_samples3x_missing0.9.chr1to4.recode.RENAMED.vcf \
         --geno-r2 \
         --ld-window-bp 50000 \
         --keep EUR_samples.txt \
         --out ld_50000_EUR

vcftools --vcf ../FINAL_SETS/nuclear_samples3x_missing0.9.chr1to4.recode.RENAMED.vcf \
         --geno-r2 \
         --ld-window-bp 50000 \
         --keep ASIA_samples.txt \
         --out ld_50000_ASIA


vcftools --vcf ../FINAL_SETS/nuclear_samples3x_missing0.9.chr1to4.recode.RENAMED.vcf \
         --geno-r2 \
         --ld-window-bp 50000 \
         --keep AUS_samples.txt \
         --out ld_50000_AUS

# only compares sites within 50,000 bp of one another
```






















#scratch





## convert my vcf file to plink format

bsub.py 10 plink_convert "../plink_convert.sh"

```bash
# vcftools
module load vcftools/0.1.16-c4

cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/02_VARIANTS/LD

plink --vcf ../FINAL_SETS/nuclear_samples3x_missing0.9.chr1to4.recode.RENAMED.vcf --make-bed --allow-extra-chr --out plink_output 

# -allow-extra-chr because plink's default is the human genome. This parameter allows it to recognise other species' chromosome naming conventions.
```


## use plink to calculate parwise linkage disequilibrium

```bash
module load plink/1.90b6.18--h516909a_0
plink --bfile plink_output --r2 --ld-window-kb 1000 --ld-window 99999 --ld-window-r2 0.2 --out ld_results
# this calculated r2 (you can also use --ld for D'
# sliding window size of 1000 kb
# maximum physical distance of 99999 bases
# --ld-window-r2 sets the min r2 for SNP pairs to be included

```








Try using various window sizes. NB: genome size of D. immitis is ~88,000,000.

bsub.py 10 LD_vcftools_100000 "../LD_vcftools_100000.sh"

```bash
module load vcftools/0.1.16-c4

cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/BATCH3/02_VARIANTS/LD

vcftools --vcf ../FINAL_SETS/nuclear_samples3x_missing0.9.chr1to4.recode.RENAMED.vcf --geno-r2 --ld-window-bp 100000 --out ld_100000
# only compares sites within 100,000 bp of one another
```


bsub.py 10 LD_vcftools_50000 "../LD_vcftools_50000.sh"

```bash
module load vcftools/0.1.16-c4

cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/02_VARIANTS/LD

vcftools --vcf ../FINAL_SETS/nuclear_samples3x_missing0.9.chr1to4.recode.RENAMED.vcf --geno-r2 --ld-window-bp 50000 --out ld_50000
# only compares sites within 50,000 bp of one another
```

bsub.py 10 LD_vcftools_10000 "../LD_vcftools_10000.sh"

```bash
module load vcftools/0.1.16-c4

cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/02_VARIANTS/LD

vcftools --vcf ../FINAL_SETS/nuclear_samples3x_missing0.9.chr1to4.recode.RENAMED.vcf --geno-r2 --ld-window-bp 10000 --out ld_10000
# only compares sites within 10,000 bp of one another
```


### plotting LD in R

```R
library(ggplot2)

##########################################################
# No window
##########################################################

# Read in data
data_nowindow <- read.table("ld_nowindow.geno.ld", header = TRUE)

# plot
ggplot(data_nowindow, aes(x = abs(BP2 - BP1), y = R2)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(x = "Distance (bp)", y = "R^2", title = "Linkage Disequilibrium Decay Plot", subtitle = "No window")
  # abs(BP2 - BP1) gets the absolute physical distance between the two variants
  # geom_smooth fits a regression line to the points


##########################################################
# Window = 100,000
##########################################################

# Read in data
data_100000 <- read.table("ld_100000.geno.ld", header = TRUE)

# plot
ggplot(data_100000, aes(x = abs(BP2 - BP1), y = R2)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(x = "Distance (bp)", y = "R^2", title = "Linkage Disequilibrium Decay Plot", subtitle = "Window size = 100,000")
  # abs(BP2 - BP1) gets the absolute physical distance between the two variants
  # geom_smooth fits a regression line to the points



##########################################################
# Window = 50,000
##########################################################

# Read in data
data_50000 <- read.table("ld_50000.geno.ld", header = TRUE)

# plot
ggplot(data_50000, aes(x = abs(BP2 - BP1), y = R2)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(x = "Distance (bp)", y = "R^2", title = "Linkage Disequilibrium Decay Plot", subtitle = "Window size = 50,000")
  # abs(BP2 - BP1) gets the absolute physical distance between the two variants
  # geom_smooth fits a regression line to the points


##########################################################
# Window = 10,000
##########################################################

# Read in data
data_10000 <- read.table("ld_10000.geno.ld", header = TRUE)

# plot
ggplot(data_10000, aes(x = abs(BP2 - BP1), y = R2)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(x = "Distance (bp)", y = "R^2", title = "Linkage Disequilibrium Decay Plot", subtitle = "Window size = 10,000")
  # abs(BP2 - BP1) gets the absolute physical distance between the two variants
  # geom_smooth fits a regression line to the points
```