#  PCA based on mitochondrial SNPs

## Principal component analysis: Filter 1 - Moderate stringency

```R
######################################################################

# PCA for filt1 - Mitochondrial VCF

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
setwd("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/R2_Sep25/pca_mito/input")

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
        "ITA: Northeast", "ITA: Pavia", "ROU: Bucharest", "ROU: Comana", "ROU: Giurgiu", "GRC: Thess/Xanthi",
        "THA: Bangkok", "MYS: Selangor",
        "AUS: Lockhart River", "AUS: Cairns", "AUS: Townsville", "AUS: Rockhampton", "AUS: Brisbane", "AUS: Sydney")), 
    ...
  )
}


######################################################################

# Mito

######################################################################




#PCA on nuclear variants using genotypes
snpgdsClose(genofile) # you need this line to close any previous file open, otherwise it won't work if you want to re-run
vcf.in <- "C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/R2_Sep25/pca_mito/input/dirofilaria_global.cohort.2025-06-18.mitoSNPs..keep_samples.dimmitis-only.recode.vcf"
gds<-snpgdsVCF2GDS(vcf.in, "mtDNA.gds", method="biallelic.only")
genofile <- snpgdsOpen(gds)


pca <-snpgdsPCA(genofile, num.thread=2, autosome.only = F)
samples <- as.data.frame(pca$sample.id)
colnames(samples) <- "name"

# Metadata file that describes where the samples come from
metadata_file <- "C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/R2_Sep25/pca_mito/location_mito.csv"
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
                                     "ITA: Northeast", "ITA: Pavia", "ROU: Bucharest", "ROU: Comana", "ROU: Giurgiu", "GRC: Thess/Xanthi",
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
         title = "Mitochondrial",
         subtitle = "SNPS: 57") +
    theme(legend.key.width = unit(0, "cm"),
          legend.text = element_text(size = 10),
          legend.key.size = unit(0.4, "cm"),
          legend.position = "right",
          axis.title = element_text(size = 16),
          axis.text = element_text(size = 14), 
          plot.title = element_text(size = 18), 
          plot.subtitle = element_text(size = 16)) +
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
                                 "ITA: Northeast", "ITA: Pavia", "ROU: Bucharest", "ROU: Comana", "ROU: Giurgiu", "GRC: Thess/Xanthi")),
  subset(data, POPULATION %in% c("USA: Illinois", "USA: Missouri", "USA: Texas", "USA: Louisiana", "USA: Georgia", "USA: Florida",
                                 "CRI: San Jose", "PAN: Boca Chica", "PAN: Puerto Armuelles", "PAN: San Lorenzo",
                                 "ITA: Northeast", "ITA: Pavia", "ROU: Bucharest", "ROU: Comana", "ROU: Giurgiu", "GRC: Thess/Xanthi",
                                 "THA: Bangkok", "MYS: Selangor")),
  subset(data, POPULATION %in% c("USA: Illinois", "USA: Missouri", "USA: Texas", "USA: Louisiana", "USA: Georgia", "USA: Florida",
                                 "CRI: San Jose", "PAN: Boca Chica", "PAN: Puerto Armuelles", "PAN: San Lorenzo",
                                 "ITA: Northeast", "ITA: Pavia", "ROU: Bucharest", "ROU: Comana", "ROU: Giurgiu", "GRC: Thess/Xanthi",
                                 "THA: Bangkok", "MYS: Selangor",
                                 "AUS: Lockhart River", "AUS: Cairns", "AUS: Townsville", "AUS: Rockhampton", "AUS: Brisbane", "AUS: Sydney"))
)



# Set a fixed height and width for all plots
plot_height <- 6  # Adjust as needed
plot_width <- 7.5   # Adjust as needed

# Define a fixed legend width
legend_width <- 0  # Adjust as needed

# Define fixed x-axis limits
x_limits <- c(-0.1, 0.4)  # Adjust as needed
y_limits <- c(-0.15, 0.23)


# Iterate through subsets and create/save plots
for (i in seq_along(subset_list)) {
  plot <- create_plot(subset_list[[i]])
  
  # Adjust the legend width in the plot
  plot <- plot + theme(legend.key.width = unit(legend_width, "cm"))
  
  # Set fixed x-axis limits
  plot <- plot + scale_x_continuous(limits = x_limits)
  
  ggsave(paste0("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/R2_Sep25/pca_mito/plot_mt_", i, ".png"), plot, height = plot_height, width = plot_width, dpi = 300)
  ggsave(paste0("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/R2_Sep25/pca_mito/plot_mt_", i, ".tif"), plot, height = plot_height, width = plot_width, dpi = 300)
  ggsave(paste0("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/R2_Sep25/pca_mito/plot_mt_", i, ".pdf"), plot, height = plot_height, width = plot_width, dpi = 300)
}


# Mito with replicates labelled

replicates <- c("AUS_BNE_AD_003", "AUS_BNE_AD_003_R", "AUS_SYD_AD_005", "AUS_SYD_AD_005_R", "AUS_SYD_AD_008", "AUS_SYD_AD_008_R", "MYS_SEL_AD_001", "MYS_SEL_AD_001_R", "PAN_PUE_AD_004", "PAN_PUE_AD_004_R")

rep_plot <- base_plot +
  geom_point(data = data, aes(EV1, EV2, col = POPULATION), alpha = 0.8, size = 2) +
  xlim(x_limits) + ylim(y_limits) +  # Set fixed limits for x and y axes
  scale_colour_pop() +
  labs(x = paste0("PC1 variance: ",round(pca$varprop[1]*100,digits=2),"%"),
       y = paste0("PC2 variance: ",round(pca$varprop[2]*100,digits=2),"%"),
       title = "Mitochondrial",
       subtitle = "SNPS: 57") +
  theme(legend.key.width = unit(0, "cm"),
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.4, "cm"),
        legend.position = "right",
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14), 
        plot.title = element_text(size = 18), 
        plot.subtitle = element_text(size = 16)) +
  guides (colour = guide_legend(ncol = 1)) +
geom_text_repel(data = subset(data, SAMPLEID %in% replicates), aes (EV1, EV2, label = SAMPLEID), size = 3.5, max.overlaps = 20)
rep_plot

ggsave(paste0("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/R2_Sep25/pca_mito/plot_mito_rep", ".tif"), rep_plot, height = plot_height, width = plot_width, dpi = 300)

ggsave(paste0("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/R2_Sep25/pca_mito/plot_mito_rep", ".png"), rep_plot, height = plot_height, width = plot_width, dpi = 300)

ggsave(paste0("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/R2_Sep25/pca_mito/plot_mito_rep", ".pdf"), rep_plot, height = plot_height, width = plot_width, dpi = 300)
```

