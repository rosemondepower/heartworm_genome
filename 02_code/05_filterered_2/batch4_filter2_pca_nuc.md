#  Principal component analysis: Filter 2 - More stringent

Make some PCA plots in R using the nuclear VCF file which was filtered with the 'filter2' set of parameters.

Summary of PCAs:

- Autosomes (chr1-4)
- Chr1 only
- Chr2 only
- Chr3 only
- Chr4 only
- ChrX


```R
######################################################################

# PCA for BATCH4 FILTER2

######################################################################
 # R version 4.4.1 (2024-06-14 ucrt) -- "Race for Your Life"

library(tidyverse) # v.2.0.0
library(gdsfmt) # v.1.40.0
library(SNPRelate) # v.1.38.0
library(ggsci) # v.3.2.0
library(ggpubr) # v.0.6.0
library(reshape2) # v.1.4.4
library(viridis) # v.0.6.5
library(vcfR) # v.1.15.0
library(factoextra) # v.1.0.7
library(ggrepel) # v.0.9.5
library(ggtree) # v.3.12.0
library(poppr) # v.2.9.6
library(adegenet) # v.2.1.10
library(ape) # v.5.8
library(RColorBrewer) # v.1.1-3
library(ggimage) # v.0.3.3


# Set working directory
setwd("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/batch4/FILTER2/NO_OUTGROUPS/pca_nuc/input")

#Preparing the data for plotting
# Set colours for different cities

red_palette1 <- brewer.pal(n = 9, name = "YlOrRd")
red_palette2 <- brewer.pal(n=9, name = "YlOrBr")
blue_palette <- brewer.pal(n = 9, name = "Blues")
green_palette1 <- brewer.pal(n = 9, name = "BuGn")
green_palette2 <- brewer.pal(n=9, name = "YlGn")


scale_colour_pop <- function(...){
  ggplot2:::manual_scale(
    'colour', 
    values = setNames(
      c(red_palette1[8],
        red_palette2[7],
        red_palette1[7],
        red_palette1[6],
        red_palette2[5],
        red_palette1[5],
        "purple4",
        "darkviolet",
        "blueviolet",
        "mediumpurple1",
        green_palette1[9], 
        green_palette2[8],
        green_palette1[7],
        green_palette2[7],
        green_palette2[5],
        green_palette1[5],
        "violetred1",
        "hotpink",
        blue_palette[9],
        blue_palette[8],
        blue_palette[7],
        blue_palette[6],
        blue_palette[5],
        blue_palette[4]),
      c("USA: Illinois", "USA: Missouri", "USA: Texas", "USA: Louisiana", "USA: Georgia", "USA: Florida",
        "CRI: San Jose", "PAN: Boca Chica", "PAN: Puerto Armuelles", "PAN: San Lorenzo",
        "ITA: Northeast", "ITA: Pavia", "ROU: Bucharest", "ROU: Comana", "ROU: Giurgiu", "GRC: Thessaloniki/Xanthi",
        "THA: Bangkok", "MYS: Selangor",
        "AUS: Lockhart River", "AUS: Cairns", "AUS: Townsville", "AUS: Rockhampton", "AUS: Brisbane", "AUS: Sydney")), 
    ...
  )
}








#PCA on nuclear variants using genotypes
snpgdsClose(genofile) # you need this line to close any previous file open, otherwise it won't work if you want to re-run
vcf.in <- "C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/batch4/FILTER2/NO_OUTGROUPS/pca_nuc/input/nuclear_samples3x_missing1.chr1to4.recode.vcf"
gds<-snpgdsVCF2GDS(vcf.in, "nucDNA.gds", method="biallelic.only")
genofile <- snpgdsOpen(gds)


pca <-snpgdsPCA(genofile, num.thread=2, autosome.only = F)
samples <- as.data.frame(pca$sample.id)
colnames(samples) <- "name"

# Metadata file that describes where the samples come from
metadata_file <- "C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/batch4/FILTER2/NO_OUTGROUPS/pca_nuc/location_nuclear.csv"
metadata <- read.csv(metadata_file, header = TRUE)

data <- data.frame(sample.id = pca$sample.id,
                   EV1 = pca$eigenvect[,1],  
                   EV2 = pca$eigenvect[,2],
                   EV3 = pca$eigenvect[,3],
                   EV4 = pca$eigenvect[,4],
                   EV5 = pca$eigenvect[,5],
                   EV6 = pca$eigenvect[,6],
                   POPULATION = metadata$location,
                   SAMPLEID = metadata$sample_name,
                   HOST = metadata$host,
                   stringsAsFactors = FALSE)

#Generating levels for population and super-population variables
data$POPULATION <- factor(data$POPULATION, 
                          levels = c("USA: Illinois", "USA: Missouri", "USA: Texas", "USA: Louisiana", "USA: Georgia", "USA: Florida",
                                     "CRI: San Jose", "PAN: Boca Chica", "PAN: Puerto Armuelles", "PAN: San Lorenzo",
                                     "ITA: Northeast", "ITA: Pavia", "ROU: Bucharest", "ROU: Comana", "ROU: Giurgiu", "GRC: Thessaloniki/Xanthi",
                                     "THA: Bangkok", "MYS: Selangor",
                                     "AUS: Lockhart River", "AUS: Cairns", "AUS: Townsville", "AUS: Rockhampton", "AUS: Brisbane", "AUS: Sydney"))

data$HOST <- factor(data$HOST, 
                          levels = c('Cat', 'Dog', 'Ferret', 'Fox', 'Golden jackal', 'Leopard', 'Unknown', 'Wildcat'))
# Note: 'Canine' samples were renamed as 'dog', unlabelled samples were assumed to be from dogs

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
    geom_point(data = subset(data_subset), aes(EV1, EV2, col = POPULATION), alpha = 0.8, size = 2) +
    xlim(x_limits) + ylim(y_limits) +  # Set fixed limits for x and y axes
    scale_colour_pop() +
    labs(x = paste0("PC1 variance: ",round(pca$varprop[1]*100,digits=2),"%"),
         y = paste0("PC2 variance: ",round(pca$varprop[2]*100,digits=2),"%"),
         title = "Nuclear",
         subtitle = "SNPS: 149,467") +
    theme(legend.key.width = unit(1, "cm"),
          legend.text = element_text(size = 6),
          legend.key.size = unit(0.4, "cm"),
          legend.position = "right") +
    guides (colour = guide_legend(ncol = 1))
  return(p)
}


# List of data subsets (you can gradually add more samples)
subset_list <- list(
  subset(data, POPULATION %in% c("USA: Illinois", "USA: Missouri", "USA: Texas", "USA: Louisiana", "USA: Georgia", "USA: Florida")),
  subset(data, POPULATION %in% c("USA: Illinois", "USA: Missouri", "USA: Texas", "USA: Louisiana", "USA: Georgia", "USA: Florida",
                                 "CRI: San Jose", "PAN: Boca Chica", "PAN: Puerto Armuelles", "PAN: San Lorenzo")),
  subset(data, POPULATION %in% c("USA: Illinois", "USA: Missouri", "USA: Texas", "USA: Mississippi", "USA: Louisiana", "USA: Georgia", "USA: Florida",
                                 "CRI: San Jose", "PAN: Boca Chica", "PAN: Puerto Armuelles", "PAN: San Lorenzo",
                                 "ITA: Northeast", "ITA: Pavia", "ROU: Bucharest", "ROU: Comana", "ROU: Giurgiu", "GRC: Thessaloniki/Xanthi")),
  subset(data, POPULATION %in% c("USA: Illinois", "USA: Missouri", "USA: Texas", "USA: Louisiana", "USA: Georgia", "USA: Florida",
                                 "CRI: San Jose", "PAN: Boca Chica", "PAN: Puerto Armuelles", "PAN: San Lorenzo",
                                 "ITA: Northeast", "ITA: Pavia", "ROU: Bucharest", "ROU: Comana", "ROU: Giurgiu", "GRC: Thessaloniki/Xanthi",
                                 "THA: Bangkok", "MYS: Selangor")),
  subset(data, POPULATION %in% c("USA: Illinois", "USA: Missouri", "USA: Texas", "USA: Louisiana", "USA: Georgia", "USA: Florida",
                                 "CRI: San Jose", "PAN: Boca Chica", "PAN: Puerto Armuelles", "PAN: San Lorenzo",
                                 "ITA: Northeast", "ITA: Pavia", "ROU: Bucharest", "ROU: Comana", "ROU: Giurgiu", "GRC: Thessaloniki/Xanthi",
                                 "THA: Bangkok", "MYS: Selangor",
                                 "AUS: Lockhart River", "AUS: Cairns", "AUS: Townsville", "AUS: Rockhampton", "AUS: Brisbane", "AUS: Sydney"))
)



# Set a fixed height and width for all plots
plot_height <- 6  # Adjust as needed
plot_width <- 8   # Adjust as needed

# Define a fixed legend width
legend_width <- 2  # Adjust as needed

# Define fixed x-axis limits
x_limits <- c(-0.15, 0.18)  # Adjust as needed
y_limits <- c(-0.15, 0.2)


# Iterate through subsets and create/save plots
for (i in seq_along(subset_list)) {
  plot <- create_plot(subset_list[[i]])
  
  # Adjust the legend width in the plot
  plot <- plot + theme(legend.key.width = unit(legend_width, "cm"))
  
  # Set fixed x-axis limits
  plot <- plot + scale_x_continuous(limits = x_limits)
  
  ggsave(paste0("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/batch4/FILTER2/NO_OUTGROUPS/pca_nuc/plot_nuc_", i, ".png"), plot, height = plot_height, width = plot_width)
}



# With ellipses
PC1_PC2_plot_ellipse <- ggplot(data, aes(x = EV1, y = EV2, color = POPULATION)) +
  theme_bw() +
  geom_point(alpha = 0.8, size = 2) +
  stat_ellipse(level = 0.95, geom = "path") +
  labs(x = paste0("PC1 variance: ", round(pca$varprop[1] * 100, digits = 2), "%"),
       y = paste0("PC2 variance: ", round(pca$varprop[2] * 100, digits = 2), "%"),
       title = "Nuclear",
       subtitle = "SNPS: 149,467") +
  scale_colour_pop() +
  theme(legend.key.width = unit(1, "cm"),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.4, "cm"),
        legend.position = "right") +
  guides(color = guide_legend(ncol=1))
PC1_PC2_plot_ellipse

ggsave(paste0("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/batch4/FILTER2/NO_OUTGROUPS/pca_nuc/plot_nuc_ellipse", ".png"), PC1_PC2_plot_ellipse, height = plot_height, width = plot_width)


  
#######################################################################################
# Label replicate samples to see if they show the same result and it's not an artefact
#######################################################################################

replicates <- c("AUS_BNE_AD_003", "AUS_BNE_AD_003_R", "AUS_SYD_AD_005", "AUS_SYD_AD_005_R", "AUS_SYD_AD_008", "AUS_SYD_AD_008_R", "PAN_PUE_AD_004", "PAN_PUE_AD_004_R")

rep_plot <- base_plot +
    geom_point(data = data, aes(EV1, EV2, col = POPULATION), alpha = 0.8, size = 2) +
    xlim(x_limits) + ylim(y_limits) +  # Set fixed limits for x and y axes
    scale_colour_pop() +
    labs(x = paste0("PC1 variance: ",round(pca$varprop[1]*100,digits=2),"%"),
         y = paste0("PC2 variance: ",round(pca$varprop[2]*100,digits=2),"%"),
         title = "Nuclear",
         subtitle = "SNPS: 149,467") +
    theme(legend.key.width = unit(1, "cm"),
          legend.text = element_text(size = 6),
          legend.key.size = unit(0.4, "cm"),
          legend.position = "right") +
    guides (colour = guide_legend(ncol = 1)) +
    geom_text_repel(data = subset(data, SAMPLEID %in% replicates), aes (EV1, EV2, label = SAMPLEID), size = 2)
rep_plot
  
ggsave(paste0("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/batch4/FILTER2/NO_OUTGROUPS/pca_nuc/plot_nuc_rep", ".png"), rep_plot, height = plot_height, width = plot_width)


#####################################################################################################################################
# Make a PCA using my other code (follows Javier's approach) to make sure I'm getting the same thing. Also get labels on the graph.
#####################################################################################################################################

# Make PCA plots

# PC1 vs PC2
PC1_PC2_plot <- ggplot(data, aes(x = EV1, y = EV2, color = POPULATION, label = SAMPLEID)) +
  theme_bw() +
  geom_point(alpha = 0.8, size = 2) +
  geom_text_repel(aes(label = SAMPLEID), box.padding = 0.05, point.padding = 0.01, segment.color = 'grey50', size = 3, hjust = 0, vjust = 0, max.overlaps = Inf, show.legend=FALSE) +
  labs(x = paste0("PC1 variance: ", round(pca$varprop[1] * 100, digits = 2), "%"),
       y = paste0("PC2 variance: ", round(pca$varprop[2] * 100, digits = 2), "%"),
       title = "Nuclear",
       subtitle = "SNPS: 149,467") +
  scale_colour_pop() +
  theme(legend.key.width = unit(1, "cm"),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.4, "cm"),
        legend.position = "right") +
  guides (color = guide_legend(ncol=1))

PC1_PC2_plot
ggsave(paste0("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/batch4/FILTER2/NO_OUTGROUPS/pca_nuc/PC1_PC2_plot_sampleid.png"), PC1_PC2_plot, height = plot_height, width = plot_width)

# PC3 vs PC4
PC3_PC4_plot <- ggplot(data, aes(x = EV3, y = EV4, color = POPULATION, label = SAMPLEID)) +
  theme_bw() +
  geom_point(alpha = 0.8, size = 2) +
 # geom_text_repel(aes(label = SAMPLEID), box.padding = 0.05, point.padding = 0.01, segment.color = 'grey50', size = 3, hjust = 0, vjust = 0, max.overlaps = Inf, show.legend=FALSE) +
  labs(x = paste0("PC3 variance: ", round(pca$varprop[3] * 100, digits = 2), "%"),
       y = paste0("PC4 variance: ", round(pca$varprop[4] * 100, digits = 2), "%"),
       title = "Nuclear",
       subtitle = "SNPS: 149,467") +
  scale_colour_pop() +
  theme(legend.key.width = unit(1, "cm"),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.4, "cm"),
        legend.position = "right") +
  guides (color = guide_legend(ncol=1))

PC3_PC4_plot
ggsave(paste0("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/batch4/FILTER2/NO_OUTGROUPS/pca_nuc/PC3_PC4_plot.png"), PC3_PC4_plot, height = plot_height, width = plot_width)

######################################################################

# PCs on linear plane

######################################################################

# PC1
constant_y <- 0
pc1_line_data <- data.frame(PC1 = data$EV1, Y = constant_y)
PC1_line <- ggplot(pc1_line_data, aes(x = PC1, y = Y)) +
  geom_point(data = data, aes(x = EV1, y = rep(constant_y, nrow(data)), color = POPULATION), alpha = 0.8, size = 2) +
  theme_bw() +
  labs(x = "PC1", y = NULL) +
  scale_colour_pop() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "none") +
  theme(panel.grid = element_blank())

# Print the plot
print(PC1_line)
# Save plot
ggsave("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/batch4/FILTER2/NO_OUTGROUPS/pca_nuc/PC1_line.png", height=2, width=8)

# PC2
constant_y <- 0
pc2_line_data <- data.frame(PC2 = data$EV2, Y = constant_y)
PC2_line <- ggplot(pc2_line_data, aes(x = PC2, y = Y)) +
  geom_point(data = data, aes(x = EV2, y = rep(constant_y, nrow(data)), color = POPULATION), alpha = 0.8, size = 2) +
  theme_bw() +
  labs(x = "PC2", y = NULL) +
  scale_colour_pop() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "none") +
  theme(panel.grid = element_blank())

# Print the plot
print(PC2_line)
# Save plot
ggsave("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/batch4/FILTER2/NO_OUTGROUPS/pca_nuc/PC2_line.png", height=2, width=8)



######################################################################

# Per chromosome

######################################################################



# Chromosome 1

snpgdsClose(genofile_chr1)
vcf.in_chr1 <- "C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/batch4/FILTER2/NO_OUTGROUPS/pca_nuc/input/nuclear_samples3x_missing1.chr1.recode.vcf"
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
                        POPULATION = metadata$location,
                        SAMPLEID = metadata$sample_name,
                        stringsAsFactors = FALSE)

#Generating levels for population and super-population variables
data_chr1$POPULATION <- factor(data_chr1$POPULATION, 
                               levels = c("USA: Illinois", "USA: Missouri", "USA: Texas", "USA: Louisiana", "USA: Georgia", "USA: Florida",
                                          "CRI: San Jose", "PAN: Boca Chica", "PAN: Puerto Armuelles", "PAN: San Lorenzo",
                                          "ITA: Northeast", "ITA: Pavia", "ROU: Bucharest", "ROU: Comana", "ROU: Giurgiu", "GRC: Thessaloniki/Xanthi",
                                          "THA: Bangkok", "MYS: Selangor",
                                          "AUS: Lockhart River", "AUS: Cairns", "AUS: Townsville", "AUS: Rockhampton", "AUS: Brisbane", "AUS: Sydney"))



# PC1 vs PC2
PC1_PC2_plot_chr1 <- ggplot(data_chr1, aes(x = EV1_chr1, y = EV2_chr1, color = POPULATION, label = SAMPLEID)) +
  theme_bw() +
  geom_point(alpha = 0.8, size = 2) +
  labs(x = paste0("PC1 variance: ", round(pca_chr1$varprop[1] * 100, digits = 2), "%"),
       y = paste0("PC2 variance: ", round(pca_chr1$varprop[2] * 100, digits = 2), "%"),
       title = "Nuclear Chr1",
       subtitle = "SNPS: 38,108") +
  scale_colour_pop() +
  theme(legend.key.width = unit(1, "cm"),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.4, "cm"),
        legend.position = "right") +
  guides(color = guide_legend(ncol=1))

PC1_PC2_plot_chr1

ggsave("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/batch4/FILTER2/NO_OUTGROUPS/pca_nuc/PC1_PC2_plot_chr1.png", PC1_PC2_plot_chr1, height = plot_height, width = plot_width)



# Chromosome 2

snpgdsClose(genofile_chr2)
vcf.in_chr2 <- "C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/batch4/FILTER2/NO_OUTGROUPS/pca_nuc/input/nuclear_samples3x_missing1.chr2.recode.vcf"
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
                        POPULATION = metadata$location,
                        SAMPLEID = metadata$sample_name,
                        stringsAsFactors = FALSE)

#Generating levels for population and super-population variables
data_chr2$POPULATION <- factor(data_chr2$POPULATION, 
                               levels = c("USA: Illinois", "USA: Missouri", "USA: Texas", "USA: Louisiana", "USA: Georgia", "USA: Florida",
                                          "CRI: San Jose", "PAN: Boca Chica", "PAN: Puerto Armuelles", "PAN: San Lorenzo",
                                          "ITA: Northeast", "ITA: Pavia", "ROU: Bucharest", "ROU: Comana", "ROU: Giurgiu", "GRC: Thessaloniki/Xanthi",
                                          "THA: Bangkok", "MYS: Selangor",
                                          "AUS: Lockhart River", "AUS: Cairns", "AUS: Townsville", "AUS: Rockhampton", "AUS: Brisbane", "AUS: Sydney"))



# PC1 vs PC2
PC1_PC2_plot_chr2 <- ggplot(data_chr2, aes(x = EV1_chr2, y = EV2_chr2, color = POPULATION, label = SAMPLEID)) +
  theme_bw() +
  geom_point(alpha = 0.8, size = 2) +
  labs(x = paste0("PC1 variance: ", round(pca_chr2$varprop[1] * 100, digits = 2), "%"),
       y = paste0("PC2 variance: ", round(pca_chr2$varprop[2] * 100, digits = 2), "%"),
       title = "Nuclear Chr2",
       subtitle = "SNPS: 33,649") +
  scale_colour_pop() +
  theme(legend.key.width = unit(1, "cm"),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.4, "cm"),
        legend.position = "right") +
  guides(color = guide_legend(ncol=1))

PC1_PC2_plot_chr2

ggsave("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/batch4/FILTER2/NO_OUTGROUPS/pca_nuc/PC1_PC2_plot_chr2.png", PC1_PC2_plot_chr2, height = plot_height, width = plot_width)




# Chromosome 3

snpgdsClose(genofile_chr3)
vcf.in_chr3 <- "C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/batch4/FILTER2/NO_OUTGROUPS/pca_nuc/input/nuclear_samples3x_missing1.chr3.recode.vcf"
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
                        POPULATION = metadata$location,
                        SAMPLEID = metadata$sample_name,
                        stringsAsFactors = FALSE)

#Generating levels for population and super-population variables
data_chr3$POPULATION <- factor(data_chr3$POPULATION, 
                               levels = c("USA: Illinois", "USA: Missouri", "USA: Texas", "USA: Louisiana", "USA: Georgia", "USA: Florida",
                                          "CRI: San Jose", "PAN: Boca Chica", "PAN: Puerto Armuelles", "PAN: San Lorenzo",
                                          "ITA: Northeast", "ITA: Pavia", "ROU: Bucharest", "ROU: Comana", "ROU: Giurgiu", "GRC: Thessaloniki/Xanthi",
                                          "THA: Bangkok", "MYS: Selangor",
                                          "AUS: Lockhart River", "AUS: Cairns", "AUS: Townsville", "AUS: Rockhampton", "AUS: Brisbane", "AUS: Sydney"))



# PC1 vs PC2
PC1_PC2_plot_chr3 <- ggplot(data_chr3, aes(x = EV1_chr3, y = EV2_chr3, color = POPULATION, label = SAMPLEID)) +
  theme_bw() +
  geom_point(alpha = 0.8, size = 2) +
  labs(x = paste0("PC1 variance: ", round(pca_chr3$varprop[1] * 100, digits = 2), "%"),
       y = paste0("PC2 variance: ", round(pca_chr3$varprop[2] * 100, digits = 2), "%"),
       title = "Nuclear Chr3",
       subtitle = "SNPS: 38,388") +
  scale_colour_pop() +
  theme(legend.key.width = unit(1, "cm"),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.4, "cm"),
        legend.position = "right") +
  guides(color = guide_legend(ncol=1))

PC1_PC2_plot_chr3

ggsave("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/batch4/FILTER2/NO_OUTGROUPS/pca_nuc/PC1_PC2_plot_chr3.png", PC1_PC2_plot_chr3, height = plot_height, width = plot_width)


# Chromosome 4

snpgdsClose(genofile_chr4)
vcf.in_chr4 <- "C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/batch4/FILTER2/NO_OUTGROUPS/pca_nuc/input/nuclear_samples3x_missing1.chr4.recode.vcf"
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
                        POPULATION = metadata$location,
                        SAMPLEID = metadata$sample_name,
                        stringsAsFactors = FALSE)

#Generating levels for population and super-population variables
data_chr4$POPULATION <- factor(data_chr4$POPULATION, 
                               levels = c("USA: Illinois", "USA: Missouri", "USA: Texas", "USA: Louisiana", "USA: Georgia", "USA: Florida",
                                          "CRI: San Jose", "PAN: Boca Chica", "PAN: Puerto Armuelles", "PAN: San Lorenzo",
                                          "ITA: Northeast", "ITA: Pavia", "ROU: Bucharest", "ROU: Comana", "ROU: Giurgiu", "GRC: Thessaloniki/Xanthi",
                                          "THA: Bangkok", "MYS: Selangor",
                                          "AUS: Lockhart River", "AUS: Cairns", "AUS: Townsville", "AUS: Rockhampton", "AUS: Brisbane", "AUS: Sydney"))



# PC1 vs PC2
PC1_PC2_plot_chr4 <- ggplot(data_chr4, aes(x = EV1_chr4, y = EV2_chr4, color = POPULATION, label = SAMPLEID)) +
  theme_bw() +
  geom_point(alpha = 0.8, size = 2) +
  labs(x = paste0("PC1 variance: ", round(pca_chr4$varprop[1] * 100, digits = 2), "%"),
       y = paste0("PC2 variance: ", round(pca_chr4$varprop[2] * 100, digits = 2), "%"),
       title = "Nuclear Chr4",
       subtitle = "SNPS: 39,322") +
  scale_colour_pop() +
  theme(legend.key.width = unit(1, "cm"),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.4, "cm"),
        legend.position = "right") +
  guides(color = guide_legend(ncol=1))

PC1_PC2_plot_chr4

ggsave("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/batch4/FILTER2/NO_OUTGROUPS/pca_nuc/PC1_PC2_plot_chr4.png", PC1_PC2_plot_chr4, height = plot_height, width = plot_width)



# Chromosome X

snpgdsClose(genofile_chrX)
vcf.in_chrX <- "C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/batch4/FILTER2/NO_OUTGROUPS/pca_nuc/input/nuclear_samples3x_missing1.chrX.recode.vcf"
gds_chrX<-snpgdsVCF2GDS(vcf.in_chrX, "nucDNA_chrX.gds", method="biallelic.only")
genofile_chrX <- snpgdsOpen(gds_chrX)


pca_chrX <-snpgdsPCA(genofile_chrX, num.thread=2, autosome.only = F)
samples_chrX <- as.data.frame(pca_chrX$sample.id)
colnames(samples_chrX) <- "name"

data_chrX <- data.frame(sample.id = pca_chrX$sample.id,
                        EV1_chrX = pca_chrX$eigenvect[,1],  
                        EV2_chrX = pca_chrX$eigenvect[,2],
                        EV3_chrX = pca_chrX$eigenvect[,3],
                        EV4_chrX = pca_chrX$eigenvect[,4],
                        EV5_chrX = pca_chrX$eigenvect[,5],
                        EV6_chrX = pca_chrX$eigenvect[,6],
                        POPULATION = metadata$location,
                        SAMPLEID = metadata$sample_name,
                        stringsAsFactors = FALSE)

#Generating levels for population and super-population variables
data_chrX$POPULATION <- factor(data_chrX$POPULATION, 
                               levels = c("USA: Illinois", "USA: Missouri", "USA: Texas", "USA: Louisiana", "USA: Georgia", "USA: Florida",
                                          "CRI: San Jose", "PAN: Boca Chica", "PAN: Puerto Armuelles", "PAN: San Lorenzo",
                                          "ITA: Northeast", "ITA: Pavia", "ROU: Bucharest", "ROU: Comana", "ROU: Giurgiu", "GRC: Thessaloniki/Xanthi",
                                          "THA: Bangkok", "MYS: Selangor",
                                          "AUS: Lockhart River", "AUS: Cairns", "AUS: Townsville", "AUS: Rockhampton", "AUS: Brisbane", "AUS: Sydney"))



# PC1 vs PC2
PC1_PC2_plot_chrX <- ggplot(data_chrX, aes(x = EV1_chrX, y = EV2_chrX, color = POPULATION, label = SAMPLEID)) +
  theme_bw() +
  geom_point(alpha = 0.8, size = 2) +
  labs(x = paste0("PC1 variance: ", round(pca_chrX$varprop[1] * 100, digits = 2), "%"),
       y = paste0("PC2 variance: ", round(pca_chrX$varprop[2] * 100, digits = 2), "%"),
       title = "Nuclear chrX",
       subtitle = "SNPS: 44,829") +
  scale_colour_pop() +
  theme(legend.position = "right") +  # Adjust legend position as needed
  guides(color = guide_legend(ncol=1))

PC1_PC2_plot_chrX

ggsave("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/batch4/FILTER2/NO_OUTGROUPS/pca_nuc/PC1_PC2_plot_chrX.png", PC1_PC2_plot_chrX, height = 8, width = 10)


######################################################################
# PCA based on the host
######################################################################




PC1_PC2_plot_host <- ggplot(data, aes(x = EV1, y = EV2, color = POPULATION, label = HOST)) +
  theme_bw() +
  geom_point(alpha = 0.8, size = 2) +
  labs(x = paste0("PC1 variance: ", round(pca$varprop[1] * 100, digits = 2), "%"),
       y = paste0("PC2 variance: ", round(pca$varprop[2] * 100, digits = 2), "%"),
       title = "Nuclear",
       subtitle = "SNPS: 149,467") +
  geom_text_repel(aes(label = HOST), box.padding = 0.05, point.padding = 0.01, segment.color = 'grey50', size = 3, hjust = 0, vjust = 0, max.overlaps = Inf, show.legend=FALSE) +
  scale_colour_pop() +
  theme(legend.key.width = unit(1, "cm"),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.4, "cm"),
        legend.position = "right") +
  guides (color = guide_legend(ncol=1))

PC1_PC2_plot_host

ggsave("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/batch4/FILTER2/NO_OUTGROUPS/pca_nuc/PC1_PC2_plot_host.png", PC1_PC2_plot_host, height = plot_height, width = plot_width)



#######################################################################################
# Make PCA focusing on Australian samples only
#######################################################################################

scale_colour_AUS <- function(...){
  ggplot2:::manual_scale(
    'colour', 
    values = setNames(
      c(blue_palette[9],
        blue_palette[8],
        blue_palette[7],
        blue_palette[6],
        blue_palette[5],
        blue_palette[4]),
      c("Lockhart River", "Cairns", "Townsville", "Rockhampton", "Brisbane", "Sydney")), 
    ...
  )
}


snpgdsClose(genofile_AUS)
vcf.in_AUS <- "C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/batch4/FILTER2/NO_OUTGROUPS/pca_nuc/input/nuclear_samples3x_missing1.chr1to4.AUS.recode.vcf"
gds_AUS<-snpgdsVCF2GDS(vcf.in_AUS, "nucDNA_AUS.gds", method="biallelic.only")
genofile_AUS <- snpgdsOpen(gds_AUS)


pca_AUS <-snpgdsPCA(genofile_AUS, num.thread=2, autosome.only = F)
samples_AUS <- as.data.frame(pca_AUS$sample.id)
colnames(samples_AUS) <- "name"

# Aus only metadata file
metadata_file_AUS <- "C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/batch4/FILTER2/NO_OUTGROUPS/pca_nuc/location_nuclear_AUS.csv"
metadata_AUS <- read.csv(metadata_file_AUS, header = TRUE)

data_AUS <- data.frame(sample.id = pca_AUS$sample.id,
                        EV1_AUS = pca_AUS$eigenvect[,1],  
                        EV2_AUS = pca_AUS$eigenvect[,2],
                        EV3_AUS = pca_AUS$eigenvect[,3],
                        EV4_AUS = pca_AUS$eigenvect[,4],
                        EV5_AUS = pca_AUS$eigenvect[,5],
                        EV6_AUS = pca_AUS$eigenvect[,6],
                        POPULATION = metadata_AUS$location,
                        SAMPLEID = metadata_AUS$sample_name,
                        stringsAsFactors = FALSE)

#Generating levels for population and super-population variables
data_AUS$POPULATION <- factor(data_AUS$POPULATION, 
                               levels = c("Lockhart River", "Cairns", "Townsville", "Rockhampton", "Brisbane", "Sydney"))



# PC1 vs PC2
x_limits <- c(-0.3, 0.3)  # Adjust as needed
y_limits <- c(-0.5, 0.3)

PC1_PC2_plot_AUS <- ggplot(data_AUS, aes(x = EV1_AUS, y = EV2_AUS, color = POPULATION)) +
  theme_bw() +
  geom_point(alpha = 0.8, size = 3) +
  xlim(x_limits) + ylim(y_limits) +  # Set fixed limits for x and y axes
  labs(x = paste0("PC1 variance: ", round(pca_AUS$varprop[1] * 100, digits = 2), "%"),
       y = paste0("PC2 variance: ", round(pca_AUS$varprop[2] * 100, digits = 2), "%"),
       title = "Nuclear AUS",
       subtitle = "SNPS: 88,316") +
  scale_colour_AUS() +
  theme(legend.key.width = unit(1, "cm"),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.4, "cm"),
        legend.position = "right") +
  guides(color = guide_legend(ncol=1))
  
PC1_PC2_plot_AUS

ggsave("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/batch4/FILTER2/NO_OUTGROUPS/pca_nuc/PC1_PC2_plot_AUS.png", PC1_PC2_plot_AUS, height = plot_height, width = plot_width)
```