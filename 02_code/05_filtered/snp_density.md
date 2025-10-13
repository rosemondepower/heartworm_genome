# Dirofilaria immitis WGS Lab Book - SNP density

### Rose Power USYD 2023

## All samples

```R
library(vcfR)
library(GenomicRanges)
library(VariantAnnotation)


######################################################################################

# Set working directory
setwd("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/batch4/FILTER1/NO_OUTGROUPS/snp_density")


# How many SNPs per 100kb window for each chromosome?

# Set options to display numbers without scientific notation
options(scipen = 10)


# CHROM 1

## Read in vcf file
vcf_file_chr1  <-  "../pca_nuc/input/nuclear_samples3x_missing0.9.chr1.recode.vcf"
vcf_chr1 <- read.vcfR(vcf_file_chr1, verbose = FALSE)
# Extract SNP positions
snp_positions_chr1 <- getPOS(vcf_chr1)
# Define the window size
window_size <- 100000  # 100 kb
# Calculate the SNP density per 100 kb window
snp_density_chr1 <- tabulate(cut(snp_positions_chr1, breaks = seq(0, max(snp_positions_chr1), by = window_size)))
## seq(0, max(snp_positions), by = window_size): This creates a sequence of breakpoints for defining the windows. We start from 0 and increment by the window_size value. This sequence represents the start positions of each window.
## cut(snp_positions, breaks = seq(0, max(snp_positions), by = window_size)): The cut() function is used to bin the snp_positions into the defined windows. It assigns each position to the corresponding window based on its value. This creates a factor object with labels indicating the window to which each SNP position belongs.
## tabulate(): This function counts the occurrences of each factor level. In this case, it counts the number of SNP positions within each window.
## snp_density <-: Assigns the resulting counts to the variable snp_density, representing the SNP density per 100 kb window.


# Create the x-axis values for plotting
x_values_chr1 <- seq(0, max(snp_positions_chr1), by = window_size)[1:length(snp_density_chr1)] + window_size / 2
## seq(0, max(snp_positions), by = window_size): This creates a sequence of values representing the start positions of each window. The sequence starts from 0 and increments by the window_size value.
## [1:length(snp_density)]: This subset operation selects a subset of the sequence generated in the previous step. It ensures that the length of the subset matches the length of the snp_density vector.
## + window_size / 2: This adds half of the window_size value to each element in the subsetted sequence. This ensures that the x-values used for plotting represent the midpoint of each window rather than the start position.
## x_values <-: Assigns the resulting values to the variable x_values.
## By using the subset operation [1:length(snp_density)] on the sequence of window start positions, we ensure that the length of the subset matches the length of the snp_density vector. This ensures that the x-values align with the corresponding SNP density values for proper plotting.

# Define the MLR SNP positions
MLR_chr1 <- c(4935110,
              5521915,
              7142460,
              14630736) # NODE_12925, NODE_55751_B, NODE_42291, NODE_13063

# Plot the SNP density
tiff(file="snp_density_chr1.tif", width=800, height=800)

plot(x_values_chr1, snp_density_chr1, type = "l", main = "Chromosome 1", xlab = "Genomic Position", ylab = "SNP Density", ylim = c(0, 1000), col = "darkorchid", lwd = 2,
     cex.lab = 1.5,  # Increase label size
     cex.main = 2,   # Increase main title size
     cex.axis = 1.2, # Increase axis text size
     cex.sub = 1)

# Add dotted lines from SNP positions of interest
for (i in seq_along(MLR_chr1)) {
  lines(c(MLR_chr1[i], MLR_chr1[i]), c(0, 1000), col = "grey", lty = 2, lwd = 2)
}

# Add points for the resistance SNPs
points(MLR_chr1, rep(1000, length(MLR_chr1)), pch = 16, col = "blue", cex = 2)

dev.off()







#Comments from Jan: we can see that there are some peaks at the beginning and end. want to investigate these peaks further. Want to find out where they are (what region in the genome etc) and what's around it, and blast that region to see what genes come up. Do they have some importance? It's normal for SNPs to be in a U shape (more SNPs near the telomeres). Can put this plot into ggplot to make it look fancier.


# CHROM 2

## Read in vcf file
vcf_file_chr2  <-  "../pca_nuc/input/nuclear_samples3x_missing0.9.chr2.recode.vcf"
vcf_chr2 <- read.vcfR(vcf_file_chr2, verbose = FALSE)
# Extract SNP positions
snp_positions_chr2 <- getPOS(vcf_chr2)
# Calculate the SNP density per 100 kb window
snp_density_chr2 <- tabulate(cut(snp_positions_chr2, breaks = seq(0, max(snp_positions_chr2), by = window_size)))
# Create the x-axis values for plotting
x_values_chr2 <- seq(0, max(snp_positions_chr2), by = window_size)[1:length(snp_density_chr2)] + window_size / 2

# Define the MLR SNP positions
MLR_chr2 <- c(4869624,
              4978210,
              5195103,
              8719077,
              15137317) # NODE_5667, NODE_29168, NODE_17333, NODE_39492, NODE_58864


# Plot the SNP density
tiff(file="snp_density_chr2.tif", width=800, height=800)
plot(x_values_chr2, snp_density_chr2, type = "l", main = "Chromosome 2", xlab = "Genomic Position", ylab = "SNP Density", ylim = c(0, 1000), col = "deepskyblue3", lwd = 2,
     cex.lab = 1.5,  # Increase label size
     cex.main = 2,   # Increase main title size
     cex.axis = 1.2, # Increase axis text size
     cex.sub = 1)

# Add dotted lines from SNP positions of interest
for (i in seq_along(MLR_chr2)) {
  lines(c(MLR_chr2[i], MLR_chr2[i]), c(0, 1000), col = "grey", lty = 2)
}

# Add points for the resistance SNPs
points(MLR_chr2, rep(1000, length(MLR_chr2)), pch = 16, col = "blue", cex = 2)
dev.off()

# CHROM 3

## Read in vcf file
vcf_file_chr3  <-  "../pca_nuc/input/nuclear_samples3x_missing0.9.chr3.recode.vcf"
vcf_chr3 <- read.vcfR(vcf_file_chr3, verbose = FALSE)
# Extract SNP positions
snp_positions_chr3 <- getPOS(vcf_chr3)
# Calculate the SNP density per 100 kb window
snp_density_chr3 <- tabulate(cut(snp_positions_chr3, breaks = seq(0, max(snp_positions_chr3), by = window_size)))
# Create the x-axis values for plotting
x_values_chr3 <- seq(0, max(snp_positions_chr3), by = window_size)[1:length(snp_density_chr3)] + window_size / 2

# Define the MLR SNP positions
MLR_chr3 <- c(8166538,
              8167140,
              8273382,
              8326593,
              8471903,
              9197509,
              13617703,
              13669081,
              13966406,
              14388742) # NODE_48992_B, NODE_48992_A, NODE_30575, NODE_15709_A, NODE_21554, NODE_27461, NODE_45689, NODE_20587, NODE_29455, NODE_58162_B


# Plot the SNP density
tiff(file="snp_density_chr3.tif", width=800, height=800)
plot(x_values_chr3, snp_density_chr3, type = "l", main = "Chromosome 3", xlab = "Genomic Position", ylab = "SNP Density", ylim = c(0, 1000), col = "springgreen4", lwd = 2,
     cex.lab = 1.5,  # Increase label size
     cex.main = 2,   # Increase main title size
     cex.axis = 1.2, # Increase axis text size
     cex.sub = 1)

# Add dotted lines from SNP positions of interest
for (i in seq_along(MLR_chr3)) {
  lines(c(MLR_chr3[i], MLR_chr3[i]), c(0, 1000), col = "grey", lty = 2)
}

# Add points for the resistance SNPs
points(MLR_chr3, rep(1000, length(MLR_chr3)), pch = 16, col = "blue", cex = 2)
dev.off()

# CHROM 4

## Read in vcf file
vcf_file_chr4  <-  "../pca_nuc/input/nuclear_samples3x_missing0.9.chr4.recode.vcf"
vcf_chr4 <- read.vcfR(vcf_file_chr4, verbose = FALSE)
# Extract SNP positions
snp_positions_chr4 <- getPOS(vcf_chr4)
# Calculate the SNP density per 100 kb window
snp_density_chr4 <- tabulate(cut(snp_positions_chr4, breaks = seq(0, max(snp_positions_chr4), by = window_size)))
# Create the x-axis values for plotting
x_values_chr4 <- seq(0, max(snp_positions_chr4), by = window_size)[1:length(snp_density_chr4)] + window_size / 2

# Define the MLR SNP positions
MLR_chr4 <- c(1580345,
              1776502,
              1785791,
              4760324,
              9971192,
              11479915,
              11498385) # NODE_5266, NODE_46063, NODE_42003, NODE_4553, NODE_51661, NODE_48750_B, NODE_48750_C


# Plot the SNP density
tiff(file="snp_density_chr4.tif", width=800, height=800)
plot(x_values_chr4, snp_density_chr4, type = "l", main = "Chromosome 4", xlab = "Genomic Position", ylab = "SNP Density", ylim = c(0, 1000), col = "darkgoldenrod2", lwd = 2,
     cex.lab = 1.5,  # Increase label size
     cex.main = 2,   # Increase main title size
     cex.axis = 1.2, # Increase axis text size
     cex.sub = 1)

# Add dotted lines from SNP positions of interest
for (i in seq_along(MLR_chr4)) {
  lines(c(MLR_chr4[i], MLR_chr4[i]), c(0, 1000), col = "grey", lty = 2)
}

# Add points for the resistance SNPs
points(MLR_chr4, rep(1000, length(MLR_chr4)), pch = 16, col = "blue", cex = 2)
dev.off()



# CHROM X

## Read in vcf file
vcf_file_chrX  <-  "../pca_nuc/input/nuclear_samples3x_missing0.9.chrX.recode.vcf"
vcf_chrX <- read.vcfR(vcf_file_chrX, verbose = FALSE)
# Extract SNP positions
snp_positions_chrX <- getPOS(vcf_chrX)
# Calculate the SNP density per 100 kb window
snp_density_chrX <- tabulate(cut(snp_positions_chrX, breaks = seq(0, max(snp_positions_chrX), by = window_size)))
# Create the x-axis values for plotting
x_values_chrX <- seq(0, max(snp_positions_chrX), by = window_size)[1:length(snp_density_chrX)] + window_size / 2

# Define the MLR SNP positions
MLR_chrX <- c(4420284,
              4529339,
              4562363,
              4842006,
              4849468,
              5014917,
              5199737,
              5524823,
              5599492,
              15729163,
              15846722,
              16386485,
              21903706,
              26552841) # NODE_1514, NODE_38622_A, NODE_38622_D, NODE_617, NODE_5365, NODE_22259, NODE_12716, NODE_10349, NODE_9858, NODE_7986, NODE_29128, NODE_35336, NODE_47722_A, NODE_26225



# Plot the SNP density
tiff(file="snp_density_chrX.tif", width=800, height=800)
plot(x_values_chrX, snp_density_chrX, type = "l", main = "Chromosome X", xlab = "Genomic Position", ylab = "SNP Density", ylim = c(0, 1000), col = "salmon", lwd = 2,
     cex.lab = 1.5,  # Increase label size
     cex.main = 2,   # Increase main title size
     cex.axis = 1.2, # Increase axis text size
     cex.sub = 1)

# Add dotted lines from SNP positions of interest
for (i in seq_along(MLR_chrX)) {
  lines(c(MLR_chrX[i], MLR_chrX[i]), c(0, 1000), col = "grey", lty = 2)
}

# Add points for the resistance SNPs
points(MLR_chrX, rep(1000, length(MLR_chrX)), pch = 16, col = "blue", cex = 2)
dev.off()



# MITO

## Read in vcf file
vcf_file_mito  <-  "../pca_mito/input/mito_samples3x_missing0.9.recode.vcf"
vcf_mito <- read.vcfR(vcf_file_mito, verbose = FALSE)
# Extract SNP positions
snp_positions_mito <- getPOS(vcf_mito)
# Calculate the SNP density per 0.5 kb window
mito_window_size <- 100
snp_density_mito <- tabulate(cut(snp_positions_mito, breaks = seq(0, max(snp_positions_mito), by = mito_window_size)))
# Create the x-axis values for plotting
x_values_mito <- seq(0, max(snp_positions_mito), by = mito_window_size)[1:length(snp_density_mito)] + mito_window_size / 2

# Plot the SNP density
tiff(file="snp_density_mito.tif", width=800, height=800)
plot(x_values_mito, snp_density_mito, type = "l", main = "Mitochondria", xlab = "Genomic Position", ylab = "SNP Density", ylim = c(0, 10), col = "black", lwd = 2,
     cex.lab = 1.5,  # Increase label size
     cex.main = 2,   # Increase main title size
     cex.axis = 1.2, # Increase axis text size
     cex.sub = 1)

dev.off()

##########################################################################################


# Combine all 5 plots

#Layout graphs

tiff(file="snp_density.tif", width=800, height=1200)

par(mfrow = c(3, 2))
par(mar = c(5, 5, 5, 5) + 0.1)

plot(x_values_chr1, snp_density_chr1, type = "l", main = "Chromosome 1", xlab = "Genomic Position", ylab = "SNP Density", ylim = c(0, 1000), col = "darkorchid", lwd = 2, cex.lab = 2, cex.main = 3, cex.axis = 1.5)
for (i in seq_along(MLR_chr1)) {
  lines(c(MLR_chr1[i], MLR_chr1[i]), c(0, 1000), col = "grey", lty = 2)
}
points(MLR_chr1, rep(1000, length(MLR_chr1)), pch = 16, col = "black", cex = 2)

plot(x_values_chr2, snp_density_chr2, type = "l", main = "Chromosome 2", xlab = "Genomic Position", ylab = "SNP Density", ylim = c(0, 1000), col = "deepskyblue3", lwd = 2, cex.lab = 2, cex.main = 3, cex.axis = 1.5)
for (i in seq_along(MLR_chr2)) {
  lines(c(MLR_chr2[i], MLR_chr2[i]), c(0, 1000), col = "grey", lty = 2)
}
points(MLR_chr2, rep(1000, length(MLR_chr2)), pch = 16, col = "black", cex = 2)

plot(x_values_chr3, snp_density_chr3, type = "l", main = "Chromosome 3", xlab = "Genomic Position", ylab = "SNP Density", ylim = c(0, 1000), col = "springgreen4", lwd = 2, cex.lab = 2, cex.main = 3, cex.axis = 1.5)
for (i in seq_along(MLR_chr3)) {
  lines(c(MLR_chr3[i], MLR_chr3[i]), c(0, 1000), col = "grey", lty = 2)
}
points(MLR_chr3, rep(1000, length(MLR_chr3)), pch = 16, col = "black", cex = 2)

plot(x_values_chr4, snp_density_chr4, type = "l", main = "Chromosome 4", xlab = "Genomic Position", ylab = "SNP Density", ylim = c(0, 1000), col = "darkgoldenrod2", lwd = 2, cex.lab = 2, cex.main = 3, cex.axis = 1.5)
for (i in seq_along(MLR_chr4)) {
  lines(c(MLR_chr4[i], MLR_chr4[i]), c(0, 1000), col = "grey", lty = 2)
}
points(MLR_chr4, rep(1000, length(MLR_chr4)), pch = 16, col = "black", cex = 2)

plot(x_values_chrX, snp_density_chrX, type = "l", main = "Chromosome X", xlab = "Genomic Position", ylab = "SNP Density", ylim = c(0, 1000), col = "salmon", lwd = 2, cex.lab = 2, cex.main = 3, cex.axis = 1.5)
for (i in seq_along(MLR_chrX)) {
  lines(c(MLR_chrX[i], MLR_chrX[i]), c(0, 1000), col = "grey", lty = 2)
}
points(MLR_chrX, rep(1000, length(MLR_chrX)), pch = 16, col = "black", cex = 2)

plot(x_values_mito, snp_density_mito, type = "l", main = "Mitochondria", xlab = "Genomic Position", ylab = "SNP Density", ylim = c(0, 10), col = "black", lwd = 2, cex.lab = 2, cex.main = 3, cex.axis = 1.5)

dev.off()

```