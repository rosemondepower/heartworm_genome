# Wolbachia analysis

We have Wolbachia in our data. We can do PCA on our samples to see if Wolbachia follows the same pattern as the host (have they co-evolved together?)

## Select only the variants for Wb

```bash
#!/bin/bash

# PBS directives 
#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N wb
#PBS -l select=1:ncpus=1:mem=4GB
#PBS -l walltime=00:20:00
#PBS -m e
#PBS -q defaultQ
#PBS -o wb.txt

# qsub ../wb.pbs

# load gatk
module load vcftools/0.1.14

cd /scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping/filter

vcftools --vcf FINAL_SETS/wb_samples3x_missing0.9.recode.vcf --remove-indels
#After filtering, kept 78 out of 78 Individuals
#After filtering, kept 427 out of a possible 427 Sites


vcftools --vcf FINAL_SETS/wb_samples3x_missing0.9.recode.vcf --keep-only-indels
# After filtering, kept 78 out of 78 Individuals
#After filtering, kept 0 out of a possible 427 Sites- this makes sense because I removed the indels earlier and only focused on the SNPs.
```


## Rename samples to more meaningful names

```bash
#!/bin/bash

# PBS directives 
#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N wb_rename
#PBS -l select=1:ncpus=1:mem=10GB
#PBS -l walltime=00:20:00
#PBS -m abe
#PBS -q defaultQ
#PBS -o wb_rename.txt

# qsub ../wb_rename.pbs

# load modules
module load bcftools/1.11

cd /scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping/filter/FINAL_SETS

bcftools reheader -s rename.txt wb_samples3x_missing0.9.recode.vcf -o wb_samples3x_missing0.9.recode.RENAMED.vcf
```

Check that it worked:

```bash
# List all sample names in the original VCF
cd /scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping/filter/FINAL_SETS
bcftools query -l wb_samples3x_missing0.9.recode.RENAMED.vcf
```
Yep all the sample names are changed.




## PCA

```R
######################################################################

# PCA for filt1 - Wolbachia

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
library(RColorBrewer)
library(ggimage)


# Set working directory
setwd("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/extra_data/filter1/wolbachia/input")

#Preparing the data for plotting
# Set colours for different cities

red_palette <- brewer.pal(n = 9, name = "YlOrRd")
blue_palette <- brewer.pal(n = 9, name = "Blues")
green_palette <- brewer.pal(n = 9, name = "Greens")
purple_palette <- brewer.pal(n = 9, name = "Purples")
pink_palette <- brewer.pal(n=9, name = "PuRd")


scale_colour_pop <- function(...){
  ggplot2:::manual_scale(
    'colour', 
    values = setNames(
      c(blue_palette[9],
        blue_palette[8],
        blue_palette[7],
        blue_palette[6],
        blue_palette[5],
        blue_palette[4],
        red_palette[9],
        red_palette[8],
        red_palette[7],
        red_palette[6],
        red_palette[4],
        pink_palette[5],
        "purple4"), 
      c('Lockhart River', 'Cairns', 
        'Townsville', 'Rockhampton',
        'Brisbane',
        'Sydney', 'Illinois', 'Missouri', 'Georgia', 'Mississippi', 'Louisiana', 'Bangkok', 'Pavia')), 
    ...
  )
}








#PCA on nuclear variants using genotypes
snpgdsClose(genofile) # you need this line to close any previous file open, otherwise it won't work if you want to re-run
vcf.in <- "C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/extra_data/filter1/wolbachia/input/wb_samples3x_missing0.9.recode.RENAMED.vcf"
gds<-snpgdsVCF2GDS(vcf.in, "wbDNA.gds", method="biallelic.only")
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
                          levels = c('Illinois', 'Missouri', 'Georgia', 'Mississippi', 'Louisiana',
                                     'Pavia',
                                     'Bangkok',
                                     'Lockhart River', 'Cairns', 'Townsville', 'Rockhampton', 'Brisbane', 'Sydney'))



# Define base plot
base_plot <- ggplot() +
  theme_bw() +
  labs(x = "PC1 variance", y = "PC2 variance")

base_plot

# Set fixed limits for x and y axes
x_margin <- 0.1  # Adjust this value for the desired x-axis margin
y_margin <- 0.1  # Adjust this value for the desired y-axis margin
x_limits <- c(min(data$EV1) - x_margin, max(data$EV1) + x_margin)
y_limits <- c(min(data$EV2) - y_margin, max(data$EV2) + y_margin)


# Function for creating plots
create_plot <- function(data_subset) {
  p <- base_plot +
    geom_point(data = subset(data_subset), aes(EV1, EV2, col = POPULATION), alpha = 0.8, size = 3) +
    xlim(x_limits) + ylim(y_limits) +  # Set fixed limits for x and y axes
    scale_colour_pop() +
    labs(x = paste0("PC1 variance: ",round(pca$varprop[1]*100,digits=2),"%"),
         y = paste0("PC2 variance: ",round(pca$varprop[2]*100,digits=2),"%")) +
    theme(legend.key.width = unit(2, "cm"))
  return(p)
}

# List of data subsets (you can gradually add more samples)
subset_list <- list(
  subset(data, POPULATION %in% c("Illinois",  "Missouri", "Georgia", "Mississippi", "Louisiana")),
  subset(data, POPULATION %in% c("Illinois",  "Missouri", "Georgia", "Mississippi", "Louisiana", "Pavia")),
  subset(data, POPULATION %in% c("Illinois",  "Missouri", "Georgia", "Mississippi", "Louisiana", "Pavia", "Bangkok")),
  subset(data, POPULATION %in% c("Illinois",  "Missouri", "Georgia", "Mississippi", "Louisiana", "Pavia", "Bangkok", "Lockhart River", "Cairns", "Townsville", "Rockhampton", "Brisbane", "Sydney"))
)

# Set a fixed height and width for all plots
plot_height <- 6  # Adjust as needed
plot_width <- 8   # Adjust as needed

# Define a fixed legend width
legend_width <- 2  # Adjust as needed

# Define fixed x-axis limits
x_limits <- c(-0.4, 0.2)  # Adjust as needed


# Iterate through subsets and create/save plots
for (i in seq_along(subset_list)) {
  plot <- create_plot(subset_list[[i]])
  
  # Adjust the legend width in the plot
  plot <- plot + theme(legend.key.width = unit(legend_width, "cm"))
  
  # Set fixed x-axis limits
  plot <- plot + scale_x_continuous(limits = x_limits)
  
  ggsave(paste0("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/extra_data/filter1/wolbachia/pca/plot_", i, ".png"), plot, height = plot_height, width = plot_width)
}











# Label samples
# Function for creating plots
create_plot <- function(data_subset) {
  p <- base_plot +
    geom_point(data = subset(data_subset), aes(EV1, EV2, col = POPULATION), alpha = 0.8, size = 2) +
    geom_text_repel(data = subset(data_subset),
                    aes(EV1, EV2, label = POPULATION),
                    size = 3,
                    box.padding = 0.5, point.padding = 0.01, segment.color = 'grey50', size = 3, hjust = 0, vjust = 0, max.overlaps = Inf, show.legend=FALSE) +
    xlim(x_limits) + ylim(y_limits) +  # Set fixed limits for x and y axes
    scale_colour_pop() +
    labs(x = paste0("PC1 variance: ",round(pca$varprop[1]*100,digits=2),"%"),
         y = paste0("PC2 variance: ",round(pca$varprop[2]*100,digits=2),"%")) +
    theme(legend.key.width = unit(2, "cm"))
  return(p)
}

# List of data subsets (you can gradually add more samples)
subset_list <- list(
  subset(data, POPULATION %in% c("Illinois",  "Missouri", "Georgia", "Mississippi", "Louisiana")),
  subset(data, POPULATION %in% c("Illinois",  "Missouri", "Georgia", "Mississippi", "Louisiana", "Pavia")),
  subset(data, POPULATION %in% c("Illinois",  "Missouri", "Georgia", "Mississippi", "Louisiana", "Pavia", "Bangkok")),
  subset(data, POPULATION %in% c("Illinois",  "Missouri", "Georgia", "Mississippi", "Louisiana", "Pavia", "Bangkok", "Lockhart River", "Cairns", "Townsville", "Rockhampton", "Brisbane", "Sydney"))
)

# Set a fixed height and width for all plots
plot_height <- 6  # Adjust as needed
plot_width <- 8   # Adjust as needed

# Define a fixed legend width
legend_width <- 2  # Adjust as needed

# Define fixed x-axis limits
x_limits <- c(-0.4, 0.2)  # Adjust as needed


# Iterate through subsets and create/save plots
for (i in seq_along(subset_list)) {
  plot <- create_plot(subset_list[[i]])
  
  # Adjust the legend width in the plot
  plot <- plot + theme(legend.key.width = unit(legend_width, "cm"))
  
  # Set fixed x-axis limits
  plot <- plot + scale_x_continuous(limits = x_limits)
  
  ggsave(paste0("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/extra_data/filter1/wolbachia/pca/plot_repel_", i, ".png"), plot, height = plot_height, width = plot_width)
```