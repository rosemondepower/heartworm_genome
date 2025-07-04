# Dirofilaria immitis WGS Lab Book - Extra Data

# Filter #3

No filtering at all.



### Also, we will select only the variants in the chr 1 to chr4, avoiding the chrX and the scaffolds

```bash
#!/bin/bash

# PBS directives 
#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N snps_qc8
#PBS -l select=1:ncpus=1:mem=4GB
#PBS -l walltime=00:10:00
#PBS -m e
#PBS -q defaultQ
#PBS -o snps_qc8.txt

# qsub ../snps_qc8.pbs

# load gatk
module load vcftools/0.1.14

cd /scratch/RDS-FSC-Heartworm_MLR-RW/mapping/extra_data/analysis/mapping/filter/filter3

vcftools --vcf /scratch/RDS-FSC-Heartworm_MLR-RW/mapping/extra_data/analysis/mapping/filter/filter1/Dirofilaria_immitis_Sep2023.nuclearSNPs.vcf \
--chr dirofilaria_immitis_chr1 \
--chr dirofilaria_immitis_chr2 \
--chr dirofilaria_immitis_chr3 \
--chr dirofilaria_immitis_chr4 \
--recode --out nuclearSNPs.chr1to4
# After filtering, kept 61 out of 61 Individuals
#Outputting VCF file...
#After filtering, kept 361827 out of a possible 506734 Sites


vcftools --vcf nuclearSNPs.chr1to4.recode.vcf --remove-indels
# After filtering, kept 61 out of 61 Individuals
#After filtering, kept 361827 out of a possible 361827 Sites


vcftools --vcf nuclearSNPs.chr1to4.recode.vcf --keep-only-indels
# After filtering, kept 61 out of 61 Individuals
#After filtering, kept 0 out of a possible 361827 Sites - this makes sense because I removed the indels earlier and only focused on the SNPs.
```


Rename samples to more meaningful names:

```bash
#!/bin/bash

# PBS directives 
#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N vcf_rename
#PBS -l select=1:ncpus=2:mem=20GB
#PBS -l walltime=00:20:00
#PBS -m abe
#PBS -q defaultQ
#PBS -o vcf_rename.txt

# qsub ../vcf_rename.pbs

# load modules
module load bcftools/1.11

cd /scratch/RDS-FSC-Heartworm_MLR-RW/mapping/extra_data/analysis/mapping/filter/filter3

bcftools reheader -s rename.txt nuclearSNPs.chr1to4.recode.vcf -o nuclearSNPs.chr1to4.recode.RENAMED.vcf
```


Check that it worked:

```bash
# List all sample names in the original VCF
bcftools query -l /scratch/RDS-FSC-Heartworm_MLR-RW/mapping/extra_data/analysis/mapping/filter/filter2/FINAL_SETS/nuclear_samples3x_missing1.chr1to4.recode.RENAMED.vcf
```

Yep all the sample names are changed.


## PCA

```R
######################################################################

# PCA for filt3

######################################################################

library(tidyverse)
library(gdsfmt)
library(SNPRelate)
library(ggsci)
library(ggpubr)
library(reshape2)
library(viridis)
library(vcfR)
library(factoextra)
library(ggrepel)
library(ggtree)
library(poppr)
library(adegenet)
library(ape)
library(ggimage)


# Set working directory
setwd("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/extra_data/filter3/vcf/input")

#Preparing the data for plotting
# Set colours for different cities
scale_colour_pop <- function(...){
  ggplot2:::manual_scale(
    'colour', 
    values = setNames(
      c('cadetblue1',
        'steelblue1',
        'royalblue4',
        'lightcyan3',
        'darkblue',
        'orchid3',
        'springgreen4',
        'mediumpurple2',
        'orangered3',
        'indianred1',
        'coral2',
        'chocolate1',
        'darkorange3'), 
      c('Lockhart River Cooktown', 'Cairns', 
        'Townsville', 'Rockhampton',
        'Brisbane',
        'Sydney', 'Bangkok', 'Pavia', 'USA', 'Mississippi', 'Illinois', 'Louisiana', 'Georgia')), 
    ...
  )
}


#PCA on nuclear variants using genotypes
snpgdsClose(genofile) # you need this line to close any previous file open, otherwise it won't work if you want to re-run
vcf.in <- "C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/extra_data/filter3/vcf/input/nuclearSNPs.chr1to4.recode.RENAMED.vcf"
gds<-snpgdsVCF2GDS(vcf.in, "nucDNA.gds", method="biallelic.only")
genofile <- snpgdsOpen(gds)


pca <-snpgdsPCA(genofile, num.thread=2, autosome.only = F)
samples <- as.data.frame(pca$sample.id)
colnames(samples) <- "name"

# Metadata file that describes where the samples come from
metadata_file <- "C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/extra_data/location.csv"
metadata <- read.csv(metadata_file, header = TRUE)

data <- data.frame(sample.id = pca$sample.id,
                   EV1 = pca$eigenvect[,1],  
                   EV2 = pca$eigenvect[,2],
                   EV3 = pca$eigenvect[,3],
                   EV4 = pca$eigenvect[,4],
                   EV5 = pca$eigenvect[,5],
                   EV6 = pca$eigenvect[,6],
                   POPULATION = metadata$city,
                   SAMPLEID = metadata$sample_name,
                   stringsAsFactors = FALSE)

#Generating levels for population and super-population variables
data$POPULATION <- factor(data$POPULATION, 
                          levels = c('Lockhart River Cooktown', 'Cairns', 
                                     'Townsville', 'Rockhampton',
                                     'Brisbane',
                                     'Sydney', 'Bangkok', 'Pavia', 'USA', 'Mississippi', 'Illinois', 'Louisiana', 'Georgia'))



# Plot the eigenvectors
# Calculate the explained variance for non-NaN eigenvalues - some of the later eigenvectors are NaN.
## It is not unusual to have NaN (Not-a-Number) values among the eigenvalues in the context of Principal Component Analysis (PCA). NaN values in eigenvalues typically occur when there is missing or insufficient information in the data for certain principal components to be computed.

explained_variance <- 100 * pca$eigenval[!is.na(pca$eigenval)] / sum(pca$eigenval, na.rm = TRUE)

# Create a barplot
jpeg("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/extra_data/filter3/vcf/pca/barplot_eigenvalues.jpg", width = 800, height = 600)
barplot(explained_variance, col = "dark green", ylim = c(0, 20))
title(ylab = "Percent of variance explained") 
title(xlab = "Eigenvalues")
dev.off()





# Lets extract the variance associated with the top 4 PCs, so we can use them in our plots.
# Calculate the total variance (excluding NaN values)
eig.total <- sum(pca$eigenval, na.rm = TRUE)


PC1.variance <- formatC(head(pca$eigenval)[1]/eig.total * 100)
PC2.variance <- formatC(head(pca$eigenval)[2]/eig.total * 100)
PC3.variance <- formatC(head(pca$eigenval)[3]/eig.total * 100)
PC4.variance <- formatC(head(pca$eigenval)[4]/eig.total * 100)


PC1.variance
# 15.42
PC2.variance
# 11.19
PC3.variance
# 6.099
PC4.variance
# 5.787


# Make PCA plots

# PC1 vs PC2
PC1_PC2_plot <- ggplot(data, aes(x = EV1, y = EV2, color = POPULATION, label = SAMPLEID)) +
  theme_bw() +
  geom_point(alpha = 0.8, size = 2) +
  geom_text_repel(aes(label = SAMPLEID), box.padding = 0.05, point.padding = 0.01, segment.color = 'grey50', size = 3, hjust = 0, vjust = 0, max.overlaps = Inf, show.legend=FALSE) +
  labs(x = paste0("PC1 variance: ", round(pca$varprop[1] * 100, digits = 2), "%"),
       y = paste0("PC2 variance: ", round(pca$varprop[2] * 100, digits = 2), "%")) +
  scale_colour_pop() +
  theme(legend.position = "right") +  # Adjust legend position as needed
ggtitle("filter3: PC1 vs PC2")
  
PC1_PC2_plot

ggsave("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/extra_data/filter3/vcf/pca/PC1_PC2_plot.png", PC1_PC2_plot, height = 6, width = 8)

# PC1 vs PC3
PC1_PC3_plot <- ggplot(data, aes(x = EV1, y = EV3, color = POPULATION, label = SAMPLEID)) +
  theme_bw() +
  geom_point(alpha = 0.8, size = 2) +
  geom_text_repel(aes(label = SAMPLEID), box.padding = 0.05, point.padding = 0.01, segment.color = 'grey50', size = 3, hjust = 0, vjust = 0, max.overlaps = Inf, show.legend=FALSE) +
  labs(x = paste0("PC1 variance: ", round(pca$varprop[1] * 100, digits = 2), "%"),
       y = paste0("PC3 variance: ", round(pca$varprop[3] * 100, digits = 2), "%")) +
  scale_colour_pop() +
  theme(legend.position = "right") +  # Adjust legend position as needed
  ggtitle("filter3: PC1 vs PC3")

PC1_PC3_plot

ggsave("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/extra_data/filter3/vcf/pca/PC1_PC3_plot.png", PC1_PC3_plot, height = 6, width = 8)

# PC1 vs PC4
PC1_PC4_plot <- ggplot(data, aes(x = EV1, y = EV4, color = POPULATION, label = SAMPLEID)) +
  theme_bw() +
  geom_point(alpha = 0.8, size = 2) +
  geom_text_repel(aes(label = SAMPLEID), box.padding = 0.05, point.padding = 0.01, segment.color = 'grey50', size = 3, hjust = 0, vjust = 0, max.overlaps = Inf, show.legend=FALSE) +
  labs(x = paste0("PC1 variance: ", round(pca$varprop[1] * 100, digits = 2), "%"),
       y = paste0("PC4 variance: ", round(pca$varprop[4] * 100, digits = 2), "%")) +
  scale_colour_pop() +
  theme(legend.position = "right") +  # Adjust legend position as needed
  ggtitle("filter3: PC1 vs PC4")

PC1_PC4_plot

ggsave("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/extra_data/filter3/vcf/pca/PC1_PC4_plot.png", PC1_PC4_plot, height = 6, width = 8)

# PC1 vs PC5
PC1_PC5_plot <- ggplot(data, aes(x = EV1, y = EV5, color = POPULATION, label = SAMPLEID)) +
  theme_bw() +
  geom_point(alpha = 0.8, size = 2) +
  geom_text_repel(aes(label = SAMPLEID), box.padding = 0.05, point.padding = 0.01, segment.color = 'grey50', size = 3, hjust = 0, vjust = 0, max.overlaps = Inf, show.legend=FALSE) +
  labs(x = paste0("PC1 variance: ", round(pca$varprop[1] * 100, digits = 2), "%"),
       y = paste0("PC5 variance: ", round(pca$varprop[5] * 100, digits = 2), "%")) +
  scale_colour_pop() +
  theme(legend.position = "right") +  # Adjust legend position as needed
  ggtitle("filter3: PC1 vs PC5")

PC1_PC5_plot

ggsave("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/extra_data/filter3/vcf/pca/PC1_PC5_plot.png", PC1_PC5_plot, height = 6, width = 8)

# PC1 vs PC6
PC1_PC6_plot <- ggplot(data, aes(x = EV1, y = EV6, color = POPULATION, label = SAMPLEID)) +
  theme_bw() +
  geom_point(alpha = 0.8, size = 2) +
  geom_text_repel(aes(label = SAMPLEID), box.padding = 0.05, point.padding = 0.01, segment.color = 'grey50', size = 3, hjust = 0, vjust = 0, max.overlaps = Inf, show.legend=FALSE) +
  labs(x = paste0("PC1 variance: ", round(pca$varprop[1] * 100, digits = 2), "%"),
       y = paste0("PC6 variance: ", round(pca$varprop[6] * 100, digits = 2), "%")) +
  scale_colour_pop() +
  theme(legend.position = "right") +  # Adjust legend position as needed
  ggtitle("filter3: PC1 vs PC6")

PC1_PC6_plot

ggsave("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/extra_data/filter3/vcf/pca/PC1_PC6_plot.png", PC1_PC6_plot, height = 6, width = 8)

# PC2 vs PC3
PC2_PC3_plot <- ggplot(data, aes(x = EV2, y = EV3, color = POPULATION, label = SAMPLEID)) +
  theme_bw() +
  geom_point(alpha = 0.8, size = 2) +
  geom_text_repel(aes(label = SAMPLEID), box.padding = 0.05, point.padding = 0.01, segment.color = 'grey50', size = 3, hjust = 0, vjust = 0, max.overlaps = Inf, show.legend=FALSE) +
  labs(x = paste0("PC2 variance: ", round(pca$varprop[2] * 100, digits = 2), "%"),
       y = paste0("PC3 variance: ", round(pca$varprop[3] * 100, digits = 2), "%")) +
  scale_colour_pop() +
  theme(legend.position = "right") +  # Adjust legend position as needed
  ggtitle("filter3: PC2 vs PC3")

PC2_PC3_plot

ggsave("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/extra_data/filter3/vcf/pca/PC2_PC3_plot.png", PC2_PC3_plot, height = 6, width = 8)

# PC2 vs PC4
PC2_PC4_plot <- ggplot(data, aes(x = EV2, y = EV4, color = POPULATION, label = SAMPLEID)) +
  theme_bw() +
  geom_point(alpha = 0.8, size = 2) +
  geom_text_repel(aes(label = SAMPLEID), box.padding = 0.05, point.padding = 0.01, segment.color = 'grey50', size = 3, hjust = 0, vjust = 0, max.overlaps = Inf, show.legend=FALSE) +
  labs(x = paste0("PC2 variance: ", round(pca$varprop[2] * 100, digits = 2), "%"),
       y = paste0("PC4 variance: ", round(pca$varprop[4] * 100, digits = 2), "%")) +
  scale_colour_pop() +
  theme(legend.position = "right") +  # Adjust legend position as needed
  ggtitle("filter3: PC2 vs PC4")

PC2_PC4_plot

ggsave("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/extra_data/filter3/vcf/pca/PC2_PC4_plot.png", PC2_PC4_plot, height = 6, width = 8)

# PC3 vs PC4
PC3_PC4_plot <- ggplot(data, aes(x = EV3, y = EV4, color = POPULATION, label = SAMPLEID)) +
  theme_bw() +
  geom_point(alpha = 0.8, size = 2) +
  geom_text_repel(aes(label = SAMPLEID), box.padding = 0.05, point.padding = 0.01, segment.color = 'grey50', size = 3, hjust = 0, vjust = 0, max.overlaps = Inf, show.legend=FALSE) +
  labs(x = paste0("PC3 variance: ", round(pca$varprop[3] * 100, digits = 2), "%"),
       y = paste0("PC4 variance: ", round(pca$varprop[4] * 100, digits = 2), "%")) +
  scale_colour_pop() +
  theme(legend.position = "right") +  # Adjust legend position as needed
  ggtitle("filter3: PC3 vs PC4")

PC3_PC4_plot

ggsave("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/extra_data/filter3/vcf/pca/PC3_PC4_plot.png", PC3_PC4_plot, height = 6, width = 8)


######################################################################

# PCs on linear plane

######################################################################

# PC1
constant_y <- 0
pc1_line_data <- data.frame(PC1 = data$EV1, Y = constant_y)
PC1_line <- ggplot(pc1_line_data, aes(x = PC1, y = Y)) +
  geom_point(data = data, aes(x = EV1, y = rep(constant_y, nrow(data)), color = POPULATION), alpha = 0.8, size = 3) +
  geom_text_repel(data = data, aes(x = EV1, y = rep(constant_y, nrow(data)), label = SAMPLEID), box.padding = 0.05, point.padding = 0.01, segment.color = 'grey50', size = 2, hjust = 0, vjust = 0, max.overlaps = Inf, show.legend = FALSE) +
  theme_bw() +
  labs(x = "PC1", y = NULL) +
  scale_colour_pop() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "none") +
  theme(panel.grid = element_blank())

# Print the plot
print(PC1_line)
# Save plot
ggsave("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/extra_data/filter3/vcf/pca/PC1_line.png", height=2, width=8)

# PC2
constant_y <- 0
PC2_line_data <- data.frame(PC2 = data$EV2, Y = constant_y)
PC2_line <- ggplot(PC2_line_data, aes(x = PC2, y = Y)) +
  geom_point(data = data, aes(x = EV2, y = rep(constant_y, nrow(data)), color = POPULATION), alpha = 0.8, size = 3) +
  geom_text_repel(data = data, aes(x = EV2, y = rep(constant_y, nrow(data)), label = SAMPLEID), box.padding = 0.05, point.padding = 0.01, segment.color = 'grey50', size = 2, hjust = 0, vjust = 0, max.overlaps = Inf, show.legend = FALSE) +
  theme_bw() +
  labs(x = "PC2", y = NULL) +
  scale_colour_pop() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "none") +
  theme(panel.grid = element_blank())

# Print the plot
print(PC2_line)
# Save plot
ggsave("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/extra_data/filter3/vcf/pca/PC2_line.png", height=2, width=8)

# PC3
constant_y <- 0
PC3_line_data <- data.frame(PC3 = data$EV3, Y = constant_y)
PC3_line <- ggplot(PC3_line_data, aes(x = PC3, y = Y)) +
  geom_point(data = data, aes(x = EV3, y = rep(constant_y, nrow(data)), color = POPULATION), alpha = 0.8, size = 3) +
  geom_text_repel(data = data, aes(x = EV3, y = rep(constant_y, nrow(data)), label = SAMPLEID), box.padding = 0.05, point.padding = 0.01, segment.color = 'grey50', size = 2, hjust = 0, vjust = 0, max.overlaps = Inf, show.legend = FALSE) +
  theme_bw() +
  labs(x = "PC3", y = NULL) +
  scale_colour_pop() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "none") +
  theme(panel.grid = element_blank())

# Print the plot
print(PC3_line)
# Save plot
ggsave("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/extra_data/filter3/vcf/pca/PC3_line.png", height=2, width=8)

# PC4
constant_y <- 0
PC4_line_data <- data.frame(PC4 = data$EV4, Y = constant_y)
PC4_line <- ggplot(PC4_line_data, aes(x = PC4, y = Y)) +
  geom_point(data = data, aes(x = EV4, y = rep(constant_y, nrow(data)), color = POPULATION), alpha = 0.8, size = 3) +
  geom_text_repel(data = data, aes(x = EV4, y = rep(constant_y, nrow(data)), label = SAMPLEID), box.padding = 0.05, point.padding = 0.01, segment.color = 'grey50', size = 2, hjust = 0, vjust = 0, max.overlaps = Inf, show.legend = FALSE) +
  theme_bw() +
  labs(x = "PC4", y = NULL) +
  scale_colour_pop() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "none") +
  theme(panel.grid = element_blank())

# Print the plot
print(PC4_line)
# Save plot
ggsave("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/extra_data/filter3/vcf/pca/PC4_line.png", height=2, width=8)

# PC5
constant_y <- 0
PC5_line_data <- data.frame(PC5 = data$EV5, Y = constant_y)
PC5_line <- ggplot(PC5_line_data, aes(x = PC5, y = Y)) +
  geom_point(data = data, aes(x = EV5, y = rep(constant_y, nrow(data)), color = POPULATION), alpha = 0.8, size = 3) +
  geom_text_repel(data = data, aes(x = EV5, y = rep(constant_y, nrow(data)), label = SAMPLEID), box.padding = 0.05, point.padding = 0.01, segment.color = 'grey50', size = 2, hjust = 0, vjust = 0, max.overlaps = Inf, show.legend = FALSE) +
  theme_bw() +
  labs(x = "PC5", y = NULL) +
  scale_colour_pop() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "none") +
  theme(panel.grid = element_blank())

# Print the plot
print(PC5_line)
# Save plot
ggsave("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/extra_data/filter3/vcf/pca/PC5_line.png", height=2, width=8)

# PC6
constant_y <- 0
PC6_line_data <- data.frame(PC6 = data$EV6, Y = constant_y)
PC6_line <- ggplot(PC6_line_data, aes(x = PC6, y = Y)) +
  geom_point(data = data, aes(x = EV6, y = rep(constant_y, nrow(data)), color = POPULATION), alpha = 0.8, size = 3) +
  geom_text_repel(data = data, aes(x = EV6, y = rep(constant_y, nrow(data)), label = SAMPLEID), box.padding = 0.05, point.padding = 0.01, segment.color = 'grey50', size = 2, hjust = 0, vjust = 0, max.overlaps = Inf, show.legend = FALSE) +
  theme_bw() +
  labs(x = "PC6", y = NULL) +
  scale_colour_pop() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "none") +
  theme(panel.grid = element_blank())

# Print the plot
print(PC6_line)
# Save plot
ggsave("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/extra_data/filter3/vcf/pca/PC6_line.png", height=2, width=8)
```

## Select variants for chr1, chr2, chr3 and chr4 (SEPARATELY)

```bash
#!/bin/bash

# PBS directives 
#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N snps_chr1-4_separate
#PBS -l select=1:ncpus=1:mem=2GB
#PBS -l walltime=00:20:00
#PBS -m e
#PBS -q defaultQ
#PBS -o snps_chr1-4_separate.txt

# qsub ../snps_chr1-4_separate.pbs

# load gatk
module load vcftools/0.1.14

cd /scratch/RDS-FSC-Heartworm_MLR-RW/mapping/extra_data/analysis/mapping/filter/filter3

# chr1
vcftools --vcf nuclearSNPs.chr1to4.recode.RENAMED.vcf \
--chr dirofilaria_immitis_chr1 \
--recode --out nuclearSNPs.chr1
# After filtering, kept 61 out of 61 Individuals
#Outputting VCF file...
#After filtering, kept 97082 out of a possible 361827 Sites

# chr2
vcftools --vcf nuclearSNPs.chr1to4.recode.RENAMED.vcf \
--chr dirofilaria_immitis_chr2 \
--recode --out nuclearSNPs.chr2
# After filtering, kept 61 out of 61 Individuals
#Outputting VCF file...
#After filtering, kept 88518 out of a possible 361827 Sites

# chr3
vcftools --vcf nuclearSNPs.chr1to4.recode.RENAMED.vcf \
--chr dirofilaria_immitis_chr3 \
--recode --out nuclearSNPs.chr3
# After filtering, kept 61 out of 61 Individuals
#Outputting VCF file...
#After filtering, kept 91638 out of a possible 361827 Sites

# chr4
vcftools --vcf nuclearSNPs.chr1to4.recode.RENAMED.vcf \
--chr dirofilaria_immitis_chr4 \
--recode --out nuclearSNPs.chr4
# After filtering, kept 61 out of 61 Individuals
#Outputting VCF file...
#After filtering, kept 84589 out of a possible 361827 Sites
```

## PCA per chromosome

```R
######################################################################

# Per chromosome

######################################################################



# Chromosome 1

snpgdsClose(genofile_chr1)
vcf.in_chr1 <- "C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/extra_data/filter3/vcf/input/nuclearSNPs.chr1.recode.vcf"
gds_chr1<-snpgdsVCF2GDS(vcf.in_chr1, "nucDNA_chr1.gds", method="biallelic.only")
genofile_chr1 <- snpgdsOpen(gds_chr1)


pca_chr1 <-snpgdsPCA(genofile_chr1, num.thread=2, autosome.only = F)
samples_chr1 <- as.data.frame(pca_chr1$sample.id)
colnames(samples_chr1) <- "name"

data_chr1 <- data.frame(sample.id = pca_chr1$sample.id,
                        EV1_chr1 = pca_chr1$eigenvect[,1],  
                        EV2_chr1 = pca_chr1$eigenvect[,2],
                        EV3_chr1 = pca_chr1$eigenvect[,3],
                        EV4_chr1 = pca_chr1$eigenvect[,4],
                        EV5_chr1 = pca_chr1$eigenvect[,5],
                        EV6_chr1 = pca_chr1$eigenvect[,6],
                        POPULATION = metadata$city,
                        SAMPLEID = metadata$sample_name,
                        stringsAsFactors = FALSE)

#Generating levels for population and super-population variables
data_chr1$POPULATION <- factor(data_chr1$POPULATION, 
                               levels = c('Lockhart River Cooktown', 'Cairns', 
                                          'Townsville', 'Rockhampton',
                                          'Brisbane',
                                          'Sydney', 'Bangkok', 'Pavia', 'USA', 'Mississippi', 'Illinois', 'Louisiana', 'Georgia'))



# PC1 vs PC2
PC1_PC2_plot_chr1 <- ggplot(data_chr1, aes(x = EV1_chr1, y = EV2_chr1, color = POPULATION, label = SAMPLEID)) +
  theme_bw() +
  geom_point(alpha = 0.8, size = 2) +
  geom_text_repel(aes(label = SAMPLEID), box.padding = 0.05, point.padding = 0.01, segment.color = 'grey50', size = 3, hjust = 0, vjust = 0, max.overlaps = Inf, show.legend=FALSE) +
  labs(x = paste0("PC1 variance: ", round(pca_chr1$varprop[1] * 100, digits = 2), "%"),
       y = paste0("PC2 variance: ", round(pca_chr1$varprop[2] * 100, digits = 2), "%")) +
  scale_colour_pop() +
  theme(legend.position = "right") +  # Adjust legend position as needed
  ggtitle("filter3: PC1 vs PC2 - Chromosome 1")

PC1_PC2_plot_chr1

ggsave("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/extra_data/filter3/vcf/pca/PC1_PC2_plot_chr1.png", PC1_PC2_plot_chr1, height = 6, width = 8)



# Chromosome 2

snpgdsClose(genofile_chr2)
vcf.in_chr2 <- "C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/extra_data/filter3/vcf/input/nuclearSNPs.chr2.recode.vcf"
gds_chr2<-snpgdsVCF2GDS(vcf.in_chr2, "nucDNA_chr2.gds", method="biallelic.only")
genofile_chr2 <- snpgdsOpen(gds_chr2)


pca_chr2 <-snpgdsPCA(genofile_chr2, num.thread=2, autosome.only = F)
samples_chr2 <- as.data.frame(pca_chr2$sample.id)
colnames(samples_chr2) <- "name"

data_chr2 <- data.frame(sample.id = pca_chr2$sample.id,
                        EV1_chr2 = pca_chr2$eigenvect[,1],  
                        EV2_chr2 = pca_chr2$eigenvect[,2],
                        EV3_chr2 = pca_chr2$eigenvect[,3],
                        EV4_chr2 = pca_chr2$eigenvect[,4],
                        EV5_chr2 = pca_chr2$eigenvect[,5],
                        EV6_chr2 = pca_chr2$eigenvect[,6],
                        POPULATION = metadata$city,
                        SAMPLEID = metadata$sample_name,
                        stringsAsFactors = FALSE)

#Generating levels for population and super-population variables
data_chr2$POPULATION <- factor(data_chr2$POPULATION, 
                               levels = c('Lockhart River Cooktown', 'Cairns', 
                                          'Townsville', 'Rockhampton',
                                          'Brisbane',
                                          'Sydney', 'Bangkok', 'Pavia', 'USA', 'Mississippi', 'Illinois', 'Louisiana', 'Georgia'))



# PC1 vs PC2
PC1_PC2_plot_chr2 <- ggplot(data_chr2, aes(x = EV1_chr2, y = EV2_chr2, color = POPULATION, label = SAMPLEID)) +
  theme_bw() +
  geom_point(alpha = 0.8, size = 2) +
  geom_text_repel(aes(label = SAMPLEID), box.padding = 0.05, point.padding = 0.01, segment.color = 'grey50', size = 3, hjust = 0, vjust = 0, max.overlaps = Inf, show.legend=FALSE) +
  labs(x = paste0("PC1 variance: ", round(pca_chr2$varprop[1] * 100, digits = 2), "%"),
       y = paste0("PC2 variance: ", round(pca_chr2$varprop[2] * 100, digits = 2), "%")) +
  scale_colour_pop() +
  theme(legend.position = "right") +  # Adjust legend position as needed
  ggtitle("filter3: PC1 vs PC2 - Chromosome 2")

PC1_PC2_plot_chr2

ggsave("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/extra_data/filter3/vcf/pca/PC1_PC2_plot_chr2.png", PC1_PC2_plot_chr2, height = 6, width = 8)




# Chromosome 3

snpgdsClose(genofile_chr3)
vcf.in_chr3 <- "C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/extra_data/filter3/vcf/input/nuclearSNPs.chr3.recode.vcf"
gds_chr3<-snpgdsVCF2GDS(vcf.in_chr3, "nucDNA_chr3.gds", method="biallelic.only")
genofile_chr3 <- snpgdsOpen(gds_chr3)


pca_chr3 <-snpgdsPCA(genofile_chr3, num.thread=2, autosome.only = F)
samples_chr3 <- as.data.frame(pca_chr3$sample.id)
colnames(samples_chr3) <- "name"

data_chr3 <- data.frame(sample.id = pca_chr3$sample.id,
                        EV1_chr3 = pca_chr3$eigenvect[,1],  
                        EV2_chr3 = pca_chr3$eigenvect[,2],
                        EV3_chr3 = pca_chr3$eigenvect[,3],
                        EV4_chr3 = pca_chr3$eigenvect[,4],
                        EV5_chr3 = pca_chr3$eigenvect[,5],
                        EV6_chr3 = pca_chr3$eigenvect[,6],
                        POPULATION = metadata$city,
                        SAMPLEID = metadata$sample_name,
                        stringsAsFactors = FALSE)

#Generating levels for population and super-population variables
data_chr3$POPULATION <- factor(data_chr3$POPULATION, 
                               levels = c('Lockhart River Cooktown', 'Cairns', 
                                          'Townsville', 'Rockhampton',
                                          'Brisbane',
                                          'Sydney', 'Bangkok', 'Pavia', 'USA', 'Mississippi', 'Illinois', 'Louisiana', 'Georgia'))



# PC1 vs PC2
PC1_PC2_plot_chr3 <- ggplot(data_chr3, aes(x = EV1_chr3, y = EV2_chr3, color = POPULATION, label = SAMPLEID)) +
  theme_bw() +
  geom_point(alpha = 0.8, size = 2) +
  geom_text_repel(aes(label = SAMPLEID), box.padding = 0.05, point.padding = 0.01, segment.color = 'grey50', size = 3, hjust = 0, vjust = 0, max.overlaps = Inf, show.legend=FALSE) +
  labs(x = paste0("PC1 variance: ", round(pca_chr3$varprop[1] * 100, digits = 2), "%"),
       y = paste0("PC2 variance: ", round(pca_chr3$varprop[2] * 100, digits = 2), "%")) +
  scale_colour_pop() +
  theme(legend.position = "right") +  # Adjust legend position as needed
  ggtitle("filter3: PC1 vs PC2 - Chromosome 3")

PC1_PC2_plot_chr3

ggsave("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/extra_data/filter3/vcf/pca/PC1_PC2_plot_chr3.png", PC1_PC2_plot_chr3, height = 6, width = 8)


# Chromosome 4

snpgdsClose(genofile_chr4)

vcf.in_chr4 <- "C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/extra_data/filter3/vcf/input/nuclearSNPs.chr4.recode.vcf"
gds_chr4<-snpgdsVCF2GDS(vcf.in_chr4, "nucDNA_chr4.gds", method="biallelic.only")
genofile_chr4 <- snpgdsOpen(gds_chr4)


pca_chr4 <-snpgdsPCA(genofile_chr4, num.thread=2, autosome.only = F)
samples_chr4 <- as.data.frame(pca_chr4$sample.id)
colnames(samples_chr4) <- "name"

data_chr4 <- data.frame(sample.id = pca_chr4$sample.id,
                        EV1_chr4 = pca_chr4$eigenvect[,1],  
                        EV2_chr4 = pca_chr4$eigenvect[,2],
                        EV3_chr4 = pca_chr4$eigenvect[,3],
                        EV4_chr4 = pca_chr4$eigenvect[,4],
                        EV5_chr4 = pca_chr4$eigenvect[,5],
                        EV6_chr4 = pca_chr4$eigenvect[,6],
                        POPULATION = metadata$city,
                        SAMPLEID = metadata$sample_name,
                        stringsAsFactors = FALSE)

#Generating levels for population and super-population variables
data_chr4$POPULATION <- factor(data_chr4$POPULATION, 
                               levels = c('Lockhart River Cooktown', 'Cairns', 
                                          'Townsville', 'Rockhampton',
                                          'Brisbane',
                                          'Sydney', 'Bangkok', 'Pavia', 'USA', 'Mississippi', 'Illinois', 'Louisiana', 'Georgia'))



# PC1 vs PC2
PC1_PC2_plot_chr4 <- ggplot(data_chr4, aes(x = EV1_chr4, y = EV2_chr4, color = POPULATION, label = SAMPLEID)) +
  theme_bw() +
  geom_point(alpha = 0.8, size = 2) +
  geom_text_repel(aes(label = SAMPLEID), box.padding = 0.05, point.padding = 0.01, segment.color = 'grey50', size = 3, hjust = 0, vjust = 0, max.overlaps = Inf, show.legend=FALSE) +
  labs(x = paste0("PC1 variance: ", round(pca_chr4$varprop[1] * 100, digits = 2), "%"),
       y = paste0("PC2 variance: ", round(pca_chr4$varprop[2] * 100, digits = 2), "%")) +
  scale_colour_pop() +
  theme(legend.position = "right") +  # Adjust legend position as needed
  ggtitle("filter3: PC1 vs PC2 - Chromosome 4")

PC1_PC2_plot_chr4

ggsave("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/extra_data/filter3/vcf/pca/PC1_PC2_plot_chr4.png", PC1_PC2_plot_chr4, height = 6, width = 8)
```