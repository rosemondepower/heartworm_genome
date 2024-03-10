## SNPs in Mitochondria and ChrX

## Select chrX

```bash
#!/bin/bash

# PBS directives 
#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N snps_chrX
#PBS -l select=1:ncpus=1:mem=4GB
#PBS -l walltime=00:10:00
#PBS -m e
#PBS -q defaultQ
#PBS -o snps_chrX.txt

# qsub ../snps_chrX.pbs

# load gatk
module load vcftools/0.1.14
module load bcftools/1.11

cd /scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping/filter

vcftools --vcf FINAL_SETS/nuclear_samples3x_missing0.9.recode.vcf \
--chr dirofilaria_immitis_chrX \
--recode --out FINAL_SETS/nuclear_samples3x_missing0.9.chrX
#After filtering, kept 81 out of 81 Individuals
#Outputting VCF file...
#After filtering, kept 64158 out of a possible 240259 Sites

vcftools --vcf FINAL_SETS/nuclear_samples3x_missing0.9.chrX.recode.vcf --remove-indels
#After filtering, kept 81 out of 81 Individuals
#After filtering, kept 64158 out of a possible 64158 Sites

vcftools --vcf FINAL_SETS/nuclear_samples3x_missing0.9.chrX.recode.vcf --keep-only-indels
#After filtering, kept 81 out of 81 Individuals
#After filtering, kept 0 out of a possible 64158 Sites

#Rename samples to more meaningful names:

cd FINAL_SETS

bcftools reheader -s rename.txt nuclear_samples3x_missing0.9.chrX.recode.vcf -o nuclear_samples3x_missing0.9.chrX.recode.RENAMED.vcf
```

Check that it worked:

```bash
# List all sample names in the original VCF
bcftools query -l /scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping/filter/FINAL_SETS/nuclear_samples3x_missing0.9.chrX.recode.RENAMED.vcf
```

Yep all the sample names are changed.



## Rename mito

```bash
#!/bin/bash

# PBS directives 
#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N vcf_rename_mt
#PBS -l select=1:ncpus=2:mem=20GB
#PBS -l walltime=00:20:00
#PBS -m abe
#PBS -q defaultQ
#PBS -o vcf_rename_mt.txt
#PBS -M rosemonde.power@sydney.edu.au

# qsub ../vcf_rename_mt.pbs

# load modules
module load bcftools/1.11

cd /scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping/filter/FINAL_SETS

bcftools reheader -s rename.txt mito_samples3x_missing1.recode.vcf -o mito_samples3x_missing1.recode.RENAMED.vcf
```

Check that it worked:

```bash
# List all sample names in the original VCF
bcftools query -l /scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping/filter/FINAL_SETS/mito_samples3x_missing1.recode.RENAMED.vcf
```

Yep all the sample names are changed.

## PCA

```R
######################################################################

# PCA for filt1 - Mitochondrial VCF & ChrX VCF

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
setwd("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/batch3/filter1/mt_chrX/input")

#Preparing the data for plotting
# Set colours for different cities

red_palette <- brewer.pal(n = 9, name = "YlOrRd")
blue_palette <- brewer.pal(n = 9, name = "Blues")
green_palette <- brewer.pal(n = 9, name = "Greens")
pink_palette <- brewer.pal(n=9, name = "PuRd")


scale_colour_pop <- function(...){
  ggplot2:::manual_scale(
    'colour', 
    values = setNames(
      c(red_palette[9],
        red_palette[8],
        red_palette[7],
        red_palette[6],
        red_palette[5],
        green_palette[9],
        green_palette[8],
        green_palette[6],
        green_palette[4],
        "purple4",
        "purple3",
        "mediumpurple1",
        "mediumpurple",
        pink_palette[5],
        blue_palette[9],
        blue_palette[8],
        blue_palette[7],
        blue_palette[6],
        blue_palette[5],
        blue_palette[4]),
      c('Illinois', 'Missouri', 'Georgia', 'Mississippi', 'Louisiana',
        'San Lorenzo', 'Puerto Armuelles', 'Boca Chica', 'San Jose',
        'Pavia', 'Bucharest', 'Giurgiu', 'Comana',
        'Bangkok',
        'Lockhart River', 'Cairns', 'Townsville', 'Rockhampton', 'Brisbane', 'Sydney')), 
    ...
  )
}



######################################################################

# Mito

######################################################################




#PCA on nuclear variants using genotypes
snpgdsClose(genofile) # you need this line to close any previous file open, otherwise it won't work if you want to re-run
vcf.in <- "C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/batch3/filter1/mt_chrX/input/mito_samples3x_missing1.recode.RENAMED.vcf"
gds<-snpgdsVCF2GDS(vcf.in, "mtDNA.gds", method="biallelic.only")
genofile <- snpgdsOpen(gds)


pca <-snpgdsPCA(genofile, num.thread=2, autosome.only = F)
samples <- as.data.frame(pca$sample.id)
colnames(samples) <- "name"

# Metadata file that describes where the samples come from
metadata_file <- "C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/batch3/location.csv"
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
                   HOST = metadata$host,
                   stringsAsFactors = FALSE)

#Generating levels for population and super-population variables
data$POPULATION <- factor(data$POPULATION, 
                          levels = c('Illinois', 'Missouri', 'Georgia', 'Mississippi', 'Louisiana',
                                     'San Lorenzo', 'Puerto Armuelles', 'Boca Chica', 'San Jose',
                                     'Pavia', 'Bucharest', 'Giurgiu', 'Comana',
                                     'Bangkok',
                                     'Lockhart River', 'Cairns', 'Townsville', 'Rockhampton', 'Brisbane', 'Sydney'))

data$HOST <- factor(data$HOST, 
                    levels = c('Dog', 'Fox', 'Leopard', 'Cat', 'Wildcat', 'Golden jackal'))


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
         y = paste0("PC2 variance: ",round(pca$varprop[2]*100,digits=2),"%"),
         title = "Mitochondrial") +
    theme(legend.key.width = unit(2, "cm"))
  return(p)
}


# List of data subsets (you can gradually add more samples)
subset_list <- list(
  subset(data, POPULATION %in% c("Illinois",  "Missouri", "Georgia", "Mississippi", "Louisiana")),
  subset(data, POPULATION %in% c("Illinois",  "Missouri", "Georgia", "Mississippi", "Louisiana", "San Lorenzo", "Puerto Armuelles", "Boca Chica", "San Jose")),
  subset(data, POPULATION %in% c("Illinois",  "Missouri", "Georgia", "Mississippi", "Louisiana", "San Lorenzo", "Puerto Armuelles", "Boca Chica", "San Jose","Pavia", "Bucharest", "Giurgiu", "Comana")),
  subset(data, POPULATION %in% c("Illinois",  "Missouri", "Georgia", "Mississippi", "Louisiana", "San Lorenzo", "Puerto Armuelles", "Boca Chica", "San Jose", "Pavia", "Bucharest", "Giurgiu", "Comana", "Bangkok")),
  subset(data, POPULATION %in% c("Illinois",  "Missouri", "Georgia", "Mississippi", "Louisiana", "San Lorenzo", "Puerto Armuelles", "Boca Chica", "San Jose", "Pavia", "Bucharest", "Giurgiu", "Comana", "Bangkok", "Lockhart River", "Cairns", "Townsville", "Rockhampton", "Brisbane", "Sydney"))
)



# Set a fixed height and width for all plots
plot_height <- 6  # Adjust as needed
plot_width <- 8   # Adjust as needed

# Define a fixed legend width
legend_width <- 2  # Adjust as needed

# Define fixed x-axis limits
x_limits <- c(-0.2, 0.6)  # Adjust as needed


# Iterate through subsets and create/save plots
for (i in seq_along(subset_list)) {
  plot <- create_plot(subset_list[[i]])
  
  # Adjust the legend width in the plot
  plot <- plot + theme(legend.key.width = unit(legend_width, "cm"))
  
  # Set fixed x-axis limits
  plot <- plot + scale_x_continuous(limits = x_limits)
  
  ggsave(paste0("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/batch3/filter1/mt_chrX/plot_mt_", i, ".png"), plot, height = plot_height, width = plot_width)
}





######################################################################

# Chr X

######################################################################



#PCA on nuclear variants using genotypes
snpgdsClose(genofile) # you need this line to close any previous file open, otherwise it won't work if you want to re-run
vcf.in <- "C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/batch3/filter1/mt_chrX/input/nuclear_samples3x_missing0.9.chrX.recode.RENAMED.vcf"
gds<-snpgdsVCF2GDS(vcf.in, "nucDNA.gds", method="biallelic.only")
genofile <- snpgdsOpen(gds)


pca <-snpgdsPCA(genofile, num.thread=2, autosome.only = F)
samples <- as.data.frame(pca$sample.id)
colnames(samples) <- "name"

# Metadata file that describes where the samples come from
metadata_file <- "C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/batch3/location.csv"
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
                   HOST = metadata$host,
                   stringsAsFactors = FALSE)

#Generating levels for population and super-population variables
data$POPULATION <- factor(data$POPULATION, 
                          levels = c('Illinois', 'Missouri', 'Georgia', 'Mississippi', 'Louisiana',
                                     'San Lorenzo', 'Puerto Armuelles', 'Boca Chica', 'San Jose',
                                     'Pavia', 'Bucharest', 'Giurgiu', 'Comana',
                                     'Bangkok',
                                     'Lockhart River', 'Cairns', 'Townsville', 'Rockhampton', 'Brisbane', 'Sydney'))

data$HOST <- factor(data$HOST, 
                    levels = c('Dog', 'Fox', 'Leopard', 'Cat', 'Wildcat', 'Golden jackal'))


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
         y = paste0("PC2 variance: ",round(pca$varprop[2]*100,digits=2),"%"),
         title = "Chromosome X") +
    theme(legend.key.width = unit(2, "cm"))
  return(p)
}


# List of data subsets (you can gradually add more samples)
subset_list <- list(
  subset(data, POPULATION %in% c("Illinois",  "Missouri", "Georgia", "Mississippi", "Louisiana")),
  subset(data, POPULATION %in% c("Illinois",  "Missouri", "Georgia", "Mississippi", "Louisiana", "San Lorenzo", "Puerto Armuelles", "Boca Chica", "San Jose")),
  subset(data, POPULATION %in% c("Illinois",  "Missouri", "Georgia", "Mississippi", "Louisiana", "San Lorenzo", "Puerto Armuelles", "Boca Chica", "San Jose","Pavia", "Bucharest", "Giurgiu", "Comana")),
  subset(data, POPULATION %in% c("Illinois",  "Missouri", "Georgia", "Mississippi", "Louisiana", "San Lorenzo", "Puerto Armuelles", "Boca Chica", "San Jose", "Pavia", "Bucharest", "Giurgiu", "Comana", "Bangkok")),
  subset(data, POPULATION %in% c("Illinois",  "Missouri", "Georgia", "Mississippi", "Louisiana", "San Lorenzo", "Puerto Armuelles", "Boca Chica", "San Jose", "Pavia", "Bucharest", "Giurgiu", "Comana", "Bangkok", "Lockhart River", "Cairns", "Townsville", "Rockhampton", "Brisbane", "Sydney"))
)



# Set a fixed height and width for all plots
plot_height <- 6  # Adjust as needed
plot_width <- 8   # Adjust as needed

# Define a fixed legend width
legend_width <- 2  # Adjust as needed

# Define fixed x-axis limits
x_limits <- c(-0.15, 0.4)  # Adjust as needed


# Iterate through subsets and create/save plots
for (i in seq_along(subset_list)) {
  plot <- create_plot(subset_list[[i]])
  
  # Adjust the legend width in the plot
  plot <- plot + theme(legend.key.width = unit(legend_width, "cm"))
  
  # Set fixed x-axis limits
  plot <- plot + scale_x_continuous(limits = x_limits)
  
  ggsave(paste0("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/batch3/filter1/mt_chrX/plot_chrX_", i, ".png"), plot, height = plot_height, width = plot_width)
}
```


# Mito tree

Make a mitochondrial tree using my data and data from doi: 10.1371/journal.pntd.0005028.

```bash
module load samtools/1.9
module load bcftools/1.11
module load tabix/0.2.6

DIR=/scratch/RDS-FSC-Heartworm_MLR-RW/batch3/fastq/Mito/my_data
REF_DIR=/project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_ALL/batch1/analysis/mapping
VCF_DIR=/scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping/filter/FINAL_SETS
cd $DIR

# make mitochondrial reference sequence
samtools faidx $REF_DIR/dimmitis_WSI_2.2.fa "dirofilaria_immitis_chrMtDNA" > $REF_DIR/dimmitis_WSI_2.2_chrMtDNA.fa
# index it
samtools faidx $REF_DIR/dimmitis_WSI_2.2_chrMtDNA.fa

# zip mito vcf file
bgzip -c $VCF_DIR/mito_samples3x_missing1.recode.RENAMED.vcf > $VCF_DIR/mito_samples3x_missing1.recode.RENAMED.vcf.gz
# index it
tabix -p vcf $VCF_DIR/mito_samples3x_missing1.recode.RENAMED.vcf.gz

# try 1 sample
bcftools consensus -f $REF_DIR/dimmitis_WSI_2.2_chrMtDNA.fa -s "JS6349" $VCF_DIR/mito_samples3x_missing1.recode.RENAMED.vcf.gz > JS6349_mito.fasta
# this worked, now loop through each sample
```

```bash
#!/bin/bash

# PBS directives 
#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N mito_seq
#PBS -l select=1:ncpus=2:mem=10GB
#PBS -l walltime=00:20:00
#PBS -m e
#PBS -q defaultQ
#PBS -o mito_seq.txt
#PBS -M rosemonde.power@sydney.edu.au

# Submit job
## qsub ../mito_seq.sh

module load bcftools/1.11

DIR=/scratch/RDS-FSC-Heartworm_MLR-RW/batch3/fastq/Mito/my_data
MITO_REF=/project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_ALL/batch1/analysis/mapping/dimmitis_WSI_2.2_chrMtDNA.fa
MITO_VCF=/scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping/filter/FINAL_SETS/mito_samples3x_missing1.recode.RENAMED.vcf.gz

cd ${DIR}

bcftools query -l ${MITO_VCF} | while read -r SAMPLE; do
echo ">${SAMPLE}" > ${DIR}/${SAMPLE}.fasta 
bcftools consensus -f ${MITO_REF} -s ${SAMPLE} ${MITO_VCF} | tail -n +2 >> ${DIR}/${SAMPLE}.fasta 
echo "FASTA sequence generated for ${SAMPLE}"
done

# write sample name to first line
# append sequence starting from 2nd line
```
