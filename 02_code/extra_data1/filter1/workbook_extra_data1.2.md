# Dirofilaria immitis WGS Lab Book - Extra Data version 1.2

This is the same as version 1, but here I am consistently using the combined reference genome in the snps_qc steps.


## SNPs QC

### Querying SNP and INDEL QC profiles to determine thresholds for filters

Adopted from Javier's paper.

```bash
#!/bin/bash

# PBS directives 
#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N snps_qc
#PBS -l select=1:ncpus=1:mem=10GB
#PBS -l walltime=02:00:00
#PBS -m abe
#PBS -q defaultQ
#PBS -o snps_qc.txt

# qsub ../snps_qc.pbs


# load gatk
module load gatk/4.1.4.1

WORKING_DIR=/scratch/RDS-FSC-Heartworm_MLR-RW/extra_data1/analysis/mapping/filter/filter1.2

# set reference, vcf, and mitochondrial and Wb contig
REFERENCE=/project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_ALL/old_data/data/analysis/mapping/reference_di_wol_dog.fa
VCF=/scratch/RDS-FSC-Heartworm_MLR-RW/extra_data1/analysis/mapping/filter/filter1.2/Dirofilaria_immitis_Sep2023.vcf.gz
MIT_CONTIG=dirofilaria_immitis_chrMtDNA
WB_CONTIG=dirofilaria_immitis_chrWb


# Make .dict file for reference sequence
## gatk CreateSequenceDictionary -R ${REFERENCE} -O /project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_ALL/data/analysis/mapping/dimmitis_WSI_2.2.dict # already did this previously

cd ${WORKING_DIR}

# select nuclear SNPs
gatk SelectVariants \
--reference ${REFERENCE} \
--variant ${VCF} \
--select-type-to-include SNP \
--exclude-intervals ${MIT_CONTIG} \
--exclude-intervals ${WB_CONTIG} \
--output ${VCF%.vcf.gz}.nuclearSNPs.vcf

# select nuclear INDELs
gatk SelectVariants \
--reference ${REFERENCE} \
--variant ${VCF} \
--select-type-to-include INDEL \
--exclude-intervals ${MIT_CONTIG} \
--exclude-intervals ${WB_CONTIG} \
--output ${VCF%.vcf.gz}.nuclearINDELs.vcf

# select mitochondrial SNPs
gatk SelectVariants \
--reference ${REFERENCE} \
--variant ${VCF} \
--select-type-to-include SNP \
--intervals ${MIT_CONTIG} \
--output ${VCF%.vcf.gz}.mitoSNPs.vcf

# select mitochondrial INDELs
gatk SelectVariants \
--reference ${REFERENCE} \
--variant ${VCF} \
--select-type-to-include INDEL \
--intervals ${MIT_CONTIG} \
--output ${VCF%.vcf.gz}.mitoINDELs.vcf

# select WB SNPs
gatk SelectVariants \
--reference ${REFERENCE} \
--variant ${VCF} \
--select-type-to-include SNP \
--intervals ${WB_CONTIG} \
--output ${VCF%.vcf.gz}.WbSNPs.vcf

# select WB INDELs
gatk SelectVariants \
--reference ${REFERENCE} \
--variant ${VCF} \
--select-type-to-include INDEL \
--intervals ${WB_CONTIG} \
--output ${VCF%.vcf.gz}.WbINDELs.vcf

# make a table of nuclear SNP data
gatk VariantsToTable \
--reference ${REFERENCE} \
--variant Dirofilaria_immitis_Sep2023.nuclearSNPs.vcf \
--fields CHROM --fields POS --fields QUAL --fields QD --fields DP --fields MQ --fields MQRankSum --fields FS --fields ReadPosRankSum --fields SOR \
--output GVCFall_nuclearSNPs.table

# make a table of nuclear INDEL data
gatk VariantsToTable \
--reference ${REFERENCE} \
--variant Dirofilaria_immitis_Sep2023.nuclearINDELs.vcf \
--fields CHROM --fields POS --fields QUAL --fields QD --fields DP --fields MQ --fields MQRankSum --fields FS --fields ReadPosRankSum --fields SOR \
--output GVCFall_nuclearINDELs.table

# make a table of mito SNP data
gatk VariantsToTable \
--reference ${REFERENCE} \
--variant Dirofilaria_immitis_Sep2023.mitoSNPs.vcf \
--fields CHROM --fields POS --fields QUAL --fields QD --fields DP --fields MQ --fields MQRankSum --fields FS --fields ReadPosRankSum --fields SOR \
--output GVCFall_mitoSNPs.table

# make a table of mito INDEL data
gatk VariantsToTable \
--reference ${REFERENCE} \
--variant Dirofilaria_immitis_Sep2023.mitoINDELs.vcf \
--fields CHROM --fields POS --fields QUAL --fields QD --fields DP --fields MQ --fields MQRankSum --fields FS --fields ReadPosRankSum --fields SOR \
--output GVCFall_mitoINDELs.table

# make a table of Wb SNP data
gatk VariantsToTable \
--reference ${REFERENCE} \
--variant Dirofilaria_immitis_Sep2023.WbSNPs.vcf \
--fields CHROM --fields POS --fields QUAL --fields QD --fields DP --fields MQ --fields MQRankSum --fields FS --fields ReadPosRankSum --fields SOR \
--output GVCFall_WbSNPs.table

# make a table of Wb INDEL data
gatk VariantsToTable \
--reference ${REFERENCE} \
--variant Dirofilaria_immitis_Sep2023.WbINDELs.vcf \
--fields CHROM --fields POS --fields QUAL --fields QD --fields DP --fields MQ --fields MQRankSum --fields FS --fields ReadPosRankSum --fields SOR \
--output GVCFall_WbINDELs.table
```


### Make some density plots of the data and get quantiles in R

```R
library(patchwork)
require(data.table)
library(tidyverse)
library(gridExtra)

setwd("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/extra_data/filter1.2/snps_qc")

# Import data
## Nuclear SNPs & indels
VCF_nuclear_snps <- fread('GVCFall_nuclearSNPs.table', header = TRUE, fill=TRUE, na.strings=c("","NA"), sep = "\t")
VCF_nuclear_snps <- sample_frac(VCF_nuclear_snps, 0.2) # select fraction of rows
VCF_nuclear_indels <- fread('GVCFall_nuclearINDELs.table', header = TRUE, fill=TRUE, na.strings=c("","NA"), sep = "\t")
VCF_nuclear_indels <- sample_frac(VCF_nuclear_indels, 0.2)
dim(VCF_nuclear_snps) # returns the dimensions of the data frame
dim(VCF_nuclear_indels)
VCF_nuclear <- rbind(VCF_nuclear_snps, VCF_nuclear_indels) # joins multiple rows to form a single batch
VCF_nuclear$Variant <- factor(c(rep("SNPs", dim(VCF_nuclear_snps)[1]), rep("Indels", dim(VCF_nuclear_indels)[1]))) # make a new column saying whether it's a SNP or indel

## Mitochondrial SNPs & indels
VCF_mito_snps <- fread('GVCFall_mitoSNPs.table', header = TRUE, fill=TRUE, na.strings=c("","NA"), sep = "\t")
VCF_mito_indels <- fread('GVCFall_mitoINDELs.table', header = TRUE, fill=TRUE, na.strings=c("","NA"), sep = "\t")
dim(VCF_mito_snps)
dim(VCF_mito_indels)
VCF_mito <- rbind(VCF_mito_snps, VCF_mito_indels)
VCF_mito$Variant <- factor(c(rep("SNPs", dim(VCF_mito_snps)[1]), rep("Indels", dim(VCF_mito_indels)[1])))

## Wolbachia SNPs & indels
VCF_wb_snps <- fread('GVCFall_WbSNPs.table', header = TRUE, fill=TRUE, na.strings=c("","NA"), sep = "\t")
VCF_wb_indels <- fread('GVCFall_WbINDELs.table', header = TRUE, fill=TRUE, na.strings=c("","NA"), sep = "\t")
dim(VCF_wb_snps)
dim(VCF_wb_indels)
VCF_wb <- rbind(VCF_wb_snps, VCF_wb_indels)
VCF_wb$Variant <- factor(c(rep("SNPs", dim(VCF_wb_snps)[1]), rep("Indels", dim(VCF_wb_indels)[1])))

# Set colours for the plots
snps <- '#A9E2E4'
indels <- '#F4CCCA'

# make function which makes plots
fun_variant_summaries <- function(data, title){
  # gatk hardfilter: SNP & INDEL QUAL < 0
  QUAL_quant <- quantile(data$QUAL, c(.01,.99), na.rm=T)
  
  QUAL <-
    ggplot(data, aes(x=log10(QUAL), fill=Variant)) +
    geom_density(alpha=.3) +
    geom_vline(xintercept=0, size=0.7, col="red") +
    geom_vline(xintercept=c(log10(QUAL_quant[2]), log10(QUAL_quant[3])), size=0.7, col="blue") +
    #xlim(0,10000) +
    theme_bw() +
    labs(title=paste0(title,": QUAL"))
  
  
  # DP doesnt have a hardfilter
  DP_quant <- quantile(data$DP, c(.01,.99), na.rm=T)
  
  DP <-
    ggplot(data, aes(x=log10(DP), fill=Variant)) +
    geom_density(alpha=0.3) +
    geom_vline(xintercept=log10(DP_quant), col="blue") +
    theme_bw() +
    labs(title=paste0(title,": DP"))
  
  # gatk hardfilter: SNP & INDEL QD < 2
  QD_quant <- quantile(data$QD, c(.01,.99), na.rm=T)
  
  QD <-
    ggplot(data, aes(x=QD, fill=Variant)) +
    geom_density(alpha=.3) +
    geom_vline(xintercept=2, size=0.7, col="red") +
    geom_vline(xintercept=QD_quant, size=0.7, col="blue") +
    theme_bw() +
    labs(title=paste0(title,": QD"))
  
  # gatk hardfilter: SNP FS > 60, INDEL FS > 200
  FS_quant <- quantile(data$FS, c(.01,.99), na.rm=T)
  
  FS <-
    ggplot(data, aes(x=log10(FS), fill=Variant)) +
    geom_density(alpha=.3) +
    geom_vline(xintercept=c(log10(60), log10(200)), size=0.7, col="red") +
    geom_vline(xintercept=log10(FS_quant), size=0.7, col="blue") +
    #xlim(0,250) +
    theme_bw() +
    labs(title=paste0(title,": FS"))
  
  # gatk hardfilter: SNP & INDEL MQ < 30
  MQ_quant <- quantile(data$MQ, c(.01,.99), na.rm=T)
  
  MQ <-
    ggplot(data, aes(x=MQ, fill=Variant)) + geom_density(alpha=.3) +
    geom_vline(xintercept=40, size=0.7, col="red") +
    geom_vline(xintercept=MQ_quant, size=0.7, col="blue") +
    theme_bw() +
    labs(title=paste0(title,": MQ"))
  
  # gatk hardfilter: SNP MQRankSum < -20
  MQRankSum_quant <- quantile(data$MQRankSum, c(.01,.99), na.rm=T)
  
  MQRankSum <-
    ggplot(data, aes(x=log10(MQRankSum), fill=Variant)) + geom_density(alpha=.3) +
    geom_vline(xintercept=log10(-20), size=0.7, col="red") +
    geom_vline(xintercept=log10(MQRankSum_quant), size=0.7, col="blue") +
    theme_bw() +
    labs(title=paste0(title,": MQRankSum"))
  
  
  # gatk hardfilter: SNP SOR < 4 , INDEL SOR > 10
  SOR_quant <- quantile(data$SOR, c(.01, .99), na.rm=T)
  
  SOR <-
    ggplot(data, aes(x=SOR, fill=Variant)) +
    geom_density(alpha=.3) +
    geom_vline(xintercept=c(4, 10), size=1, colour = c(snps,indels)) +
    geom_vline(xintercept=SOR_quant, size=0.7, col="blue") +
    theme_bw() +
    labs(title=paste0(title,": SOR"))
  
  # gatk hardfilter: SNP ReadPosRankSum <-10 , INDEL ReadPosRankSum < -20
  ReadPosRankSum_quant <- quantile(data$ReadPosRankSum, c(.01,.99), na.rm=T)
  
  ReadPosRankSum <-
    ggplot(data, aes(x=ReadPosRankSum, fill=Variant)) +
    geom_density(alpha=.3) +
    geom_vline(xintercept=c(-10,10,-20,20), size=1, colour = c(snps,snps,indels,indels)) +
    xlim(-10, 10) +
    geom_vline(xintercept=ReadPosRankSum_quant, size=0.7, col="blue") +
    theme_bw() +
    labs(title=paste0(title,": ReadPosRankSum"))
  
  
  plot <- QUAL + DP + QD + FS + MQ + MQRankSum + SOR + ReadPosRankSum + plot_layout(ncol=2)
  
  print(plot)
  
  ggsave(paste0("plot_",title,"_variant_summaries.png"), height=20, width=15, type="cairo")
  
  
  # generate a table of quantiles for each variant feature
  QUAL_quant <- data %>% group_by(Variant) %>% summarise(quants = list(quantile(QUAL, probs = c(0.01,0.05,0.95,0.99),na.rm=T))) %>% unnest_wider(quants)
  QUAL_quant$name <- "QUAL"
  DP_quant <- data %>% group_by(Variant) %>% summarise(quants = list(quantile(DP, probs = c(0.01,0.05,0.95,0.99),na.rm=T))) %>% unnest_wider(quants)
  DP_quant$name <- "DP"
  QD_quant <- data %>% group_by(Variant) %>% summarise(quants = list(quantile(QD, probs = c(0.01,0.05,0.95,0.99),na.rm=T))) %>% unnest_wider(quants)
  QD_quant$name <- "QD"
  FS_quant <- data %>% group_by(Variant) %>% summarise(quants = list(quantile(FS, probs = c(0.01,0.05,0.95,0.99),na.rm=T))) %>% unnest_wider(quants)
  FS_quant$name <- "FS"
  MQ_quant <- data %>% group_by(Variant) %>% summarise(quants = list(quantile(MQ, probs = c(0.01,0.05,0.95,0.99),na.rm=T))) %>% unnest_wider(quants)
  MQ_quant$name <- "MQ"
  MQRankSum_quant <- data %>% group_by(Variant) %>% summarise(quants = list(quantile(MQRankSum, probs = c(0.01,0.05,0.95,0.99),na.rm=T))) %>% unnest_wider(quants)
  MQRankSum_quant$name <- "MQRankSum"
  SOR_quant <- data %>% group_by(Variant) %>% summarise(quants = list(quantile(SOR, probs = c(0.01,0.05,0.95,0.99),na.rm=T))) %>% unnest_wider(quants)
  SOR_quant$name <- "SOR"
  ReadPosRankSum_quant <- data %>% group_by(Variant) %>% summarise(quants = list(quantile(ReadPosRankSum, probs = c(0.01,0.05,0.95,0.99),na.rm=T))) %>% unnest_wider(quants)
  ReadPosRankSum_quant$name <- "ReadPosRankSum"
  
  quantiles <- bind_rows(QUAL_quant,DP_quant, QD_quant, FS_quant, MQ_quant, MQRankSum_quant, SOR_quant, ReadPosRankSum_quant)
  quantiles$name <- c("QUAL_Indels","QUAL_SNPs","DP_indels","DP_SNPs", "QD_indels","QD_SNPs", "FS_indels","FS_SNPs", "MQ_indels","MQ_SNPs", "MQRankSum_indels","MQRankSum_SNPs", "SOR_indels","SOR_SNPs","ReadPosRankSum_indels","ReadPosRankSum_SNPs")
  
  png(paste0("table_",title,"_variant_quantiles.png"), width=1000,height=500,bg = "white")
  print(quantiles)
  grid.table(quantiles)
  dev.off()
  
}

# run nuclear variants
fun_variant_summaries(VCF_nuclear,"nuclear")
# run mitochondrial variants
fun_variant_summaries(VCF_mito,"mitochondrial")
# run wolbachia variants
fun_variant_summaries(VCF_wb,"wolbachia")






# Plot for mtDNA
##The mitochondrial genome has fewer SNPs so we might want to plot the individual points instead of density. This will give us a better idea of how to filter this data.

  # gatk hardfilter: SNP & INDEL QUAL < 0
  QUAL <-
    ggplot(VCF_mito) +
    geom_point( aes(x=log10(QUAL), y=Variant, colour=Variant)) +
    #xlim(0,10000) +
    theme_bw() +
    labs(title=paste0("mitochondrial",": QUAL"))
  QUAL
  
  # DP doesnt have a hardfilter
  DP <-
    ggplot(VCF_mito) +
    geom_point( aes(x=log10(DP), y=Variant, colour=Variant)) +
    theme_bw() +
    labs(title=paste0("mitochondrial",": DP"))
  DP
  
  # gatk hardfilter: SNP & INDEL QD < 2
  QD <-
    ggplot(VCF_mito) +
    geom_point( aes(x=QD, y=Variant, colour=Variant)) +
    theme_bw() +
    labs(title=paste0("mitochondrial",": QD"))
  QD
  # Removed 1 rows containing missing values (`geom_point()`). 
  
  # gatk hardfilter: SNP FS > 60, INDEL FS > 200
  FS <-
    ggplot(VCF_mito) +
    geom_point( aes(x=log10(FS), y=Variant, colour=Variant)) +
    #xlim(0,250) +
    theme_bw() +
    labs(title=paste0("mitochondrial",": FS"))
  FS
  
  # gatk hardfilter: SNP & INDEL MQ < 30
  MQ <-
    ggplot(VCF_mito) +
    geom_point( aes(x=MQ, y=Variant, colour=Variant)) +
    theme_bw() +
    labs(title=paste0("mitochondrial",": MQ"))
  MQ
  
  # gatk hardfilter: SNP MQRankSum < -20
  MQRankSum <-
    ggplot(VCF_mito) +
    geom_point( aes(x=log10(MQRankSum), y=Variant, colour=Variant)) +
    theme_bw() +
    labs(title=paste0("mitochondrial",": MQRankSum"))
  MQRankSum
  # Removed 86 rows containing missing values (`geom_point()`). 
  
  
  # gatk hardfilter: SNP SOR < 4 , INDEL SOR > 10
  SOR <-
    ggplot(VCF_mito) +
    geom_point( aes(x=SOR, y=Variant, colour=Variant)) +
    theme_bw() +
    labs(title=paste0("mitochondrial",": SOR"))
  SOR
  
  # gatk hardfilter: SNP ReadPosRankSum <-10 , INDEL ReadPosRankSum < -20
  ReadPosRankSum <-
    ggplot(VCF_mito) +
    geom_point( aes(x=ReadPosRankSum, y=Variant, colour=Variant)) +
    theme_bw() +
    labs(title=paste0("mitochondrial",": ReadPosRankSum"))
  ReadPosRankSum
  # Removed 6 rows containing missing values (`geom_point()`). 
  
  
  plot <- QUAL + DP + QD + FS + MQ + MQRankSum + SOR + ReadPosRankSum + plot_layout(ncol=2)
  
  print(plot)
  
  ggsave(paste0("plot_mitochondrial_indiv_variant_summaries.png"), height=20, width=15, type="cairo")
  
  
  
  
  # Repeat the same thing as above but use geom_jitter instead of geom_point.
  # gatk hardfilter: SNP & INDEL QUAL < 0
  QUAL <-
    ggplot(VCF_mito) +
    geom_jitter( aes(x=log10(QUAL), y=Variant, colour=Variant)) +
    #xlim(0,10000) +
    theme_bw() +
    labs(title=paste0("mitochondrial",": QUAL"))
  QUAL
  
  # DP doesnt have a hardfilter
  DP <-
    ggplot(VCF_mito) +
    geom_jitter( aes(x=log10(DP), y=Variant, colour=Variant)) +
    theme_bw() +
    labs(title=paste0("mitochondrial",": DP"))
  DP
  
  # gatk hardfilter: SNP & INDEL QD < 2
  QD <-
    ggplot(VCF_mito) +
    geom_jitter( aes(x=QD, y=Variant, colour=Variant)) +
    theme_bw() +
    labs(title=paste0("mitochondrial",": QD"))
  QD
  # Removed 1 rows containing missing values (`geom_point()`). 
  
  # gatk hardfilter: SNP FS > 60, INDEL FS > 200
  FS <-
    ggplot(VCF_mito) +
    geom_jitter( aes(x=log10(FS), y=Variant, colour=Variant)) +
    #xlim(0,250) +
    theme_bw() +
    labs(title=paste0("mitochondrial",": FS"))
  FS
  
  # gatk hardfilter: SNP & INDEL MQ < 30
  MQ <-
    ggplot(VCF_mito) +
    geom_jitter( aes(x=MQ, y=Variant, colour=Variant)) +
    theme_bw() +
    labs(title=paste0("mitochondrial",": MQ"))
  MQ
  
  # gatk hardfilter: SNP MQRankSum < -20
  MQRankSum <-
    ggplot(VCF_mito) +
    geom_jitter( aes(x=log10(MQRankSum), y=Variant, colour=Variant)) +
    theme_bw() +
    labs(title=paste0("mitochondrial",": MQRankSum"))
  MQRankSum
  # Removed 86 rows containing missing values (`geom_point()`). 
  
  
  # gatk hardfilter: SNP SOR < 4 , INDEL SOR > 10
  SOR <-
    ggplot(VCF_mito) +
    geom_jitter( aes(x=SOR, y=Variant, colour=Variant)) +
    theme_bw() +
    labs(title=paste0("mitochondrial",": SOR"))
  SOR
  
  # gatk hardfilter: SNP ReadPosRankSum <-10 , INDEL ReadPosRankSum < -20
  ReadPosRankSum <-
    ggplot(VCF_mito) +
    geom_jitter( aes(x=ReadPosRankSum, y=Variant, colour=Variant)) +
    theme_bw() +
    labs(title=paste0("mitochondrial",": ReadPosRankSum"))
  ReadPosRankSum
  # Removed 6 rows containing missing values (`geom_point()`). 
  
  
  plot <- QUAL + DP + QD + FS + MQ + MQRankSum + SOR + ReadPosRankSum + plot_layout(ncol=2)
  
  print(plot)
  
  ggsave(paste0("plot_mitochondrial_indiv_variant_summaries_jitter.png"), height=20, width=15, type="cairo")
```


### Attending to the quantiles, thresholds for specific parameters are established

```bash
#!/bin/bash

# PBS directives 
#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N snps_qc2
#PBS -l select=1:ncpus=1:mem=10GB
#PBS -l walltime=00:10:00
#PBS -m e
#PBS -q defaultQ
#PBS -o snps_qc2.txt

# qsub ../snps_qc2.pbs


# load gatk
module load gatk/4.1.4.1

WORKING_DIR=/scratch/RDS-FSC-Heartworm_MLR-RW/extra_data1/analysis/mapping/filter/filter1.2

# set reference, vcf
REFERENCE=/project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_ALL/old_data/data/analysis/mapping/reference_di_wol_dog.fa

cd ${WORKING_DIR}

#Nuclear
gatk VariantFiltration \
--reference ${REFERENCE} \
--variant Dirofilaria_immitis_Sep2023.nuclearSNPs.vcf \
--filter-expression 'QUAL < 32 || DP < 232 || DP > 12119 || MQ < 29.80 || SOR > 8.439 || QD < 0.350 || FS > 110.401 || MQRankSum < -7.960 || ReadPosRankSum < -10.320 || ReadPosRankSum > 4.150' \
--filter-name "SNP_filtered" \
--output Dirofilaria_immitis_Sep2023.nuclearSNPs.filtered.vcf

gatk VariantFiltration \
--reference ${REFERENCE} \
--variant Dirofilaria_immitis_Sep2023.nuclearINDELs.vcf \
--filter-expression 'QUAL < 39 || DP < 448 || DP > 9240 || MQ < 36.54 || SOR > 8.155 || QD < 0.450 || FS > 102.212 || MQRankSum < -4.343 || ReadPosRankSum < -11.090 || ReadPosRankSum > 3.500' \
--filter-name "INDEL_filtered" \
--output Dirofilaria_immitis_Sep2023.nuclearINDELs.filtered.vcf

#Mitochondrial
gatk VariantFiltration \
--reference ${REFERENCE} \
--variant Dirofilaria_immitis_Sep2023.mitoSNPs.vcf \
--filter-expression ' QUAL < 73 || DP < 74770 || DP > 344706 || MQ < 44.90 || SOR > 12.837 || QD < 0.309 || FS > 118.826 || MQRankSum < -11.162 || ReadPosRankSum < -29.939 || ReadPosRankSum > 2.988 ' \
--filter-name "SNP_filtered" \
--output Dirofilaria_immitis_Sep2023.mitoSNPs.filtered.vcf

gatk VariantFiltration \
--reference ${REFERENCE} \
--variant Dirofilaria_immitis_Sep2023.mitoINDELs.vcf \
--filter-expression 'QUAL < 41 || DP < 97127 || DP > 302774 || MQ < 45.12 || SOR > 6.131 || QD < 0.024 || FS > 38.265 || MQRankSum < -9.223 || ReadPosRankSum < -19.980 || ReadPosRankSum > 7.935' \
--filter-name "INDEL_filtered" \
--output Dirofilaria_immitis_Sep2023.mitoINDELs.filtered.vcf

#Wolbachia
gatk VariantFiltration \
--reference ${REFERENCE} \
--variant Dirofilaria_immitis_Sep2023.WbSNPs.vcf \
--filter-expression ' QUAL < 30 || DP < 46205 || DP > 75065 || MQ < 40.00 || SOR > 8.439 || QD < 0.310 || FS > 140.016 || MQRankSum < -7.359 || ReadPosRankSum < -11.790 || ReadPosRankSum > 5.226' \
--filter-name "SNP_filtered" \
--output Dirofilaria_immitis_Sep2023.WbSNPs.filtered.vcf

gatk VariantFiltration \
--reference ${REFERENCE} \
--variant Dirofilaria_immitis_Sep2023.WbINDELs.vcf \
--filter-expression 'QUAL < 36 || DP < 53935 || DP > 76905 || MQ < 42.44 || SOR > 8.254 || QD < 0.370 || FS > 137.307 || MQRankSum < -4.466 || ReadPosRankSum < -11.908 || ReadPosRankSum > 3.043' \
--filter-name "INDEL_filtered" \
--output Dirofilaria_immitis_Sep2023.WbINDELs.filtered.vcf

# once done, count the filtered sites - funny use of "|" allows direct markdown table format
echo -e "| Filtered_VCF | Variants_PASS | Variants_FILTERED |\n| -- | -- | -- | " > filter.stats

for i in *filtered.vcf; do
     name=${i}; pass=$( grep -E 'PASS' ${i} | wc -l ); filter=$( grep -E 'filter' ${i} | wc -l );
     echo -e "| ${name} | ${pass} | ${filter} |" >> filter.stats
done
```


Filtering looks ok. Didn't lose too many SNPs.

Usually we would merge the SNP and INDEL files together to make a joined VCF. However, I want to disregard the indels moving forward and only focus on the SNPs. Some downstream tools don't like having indels in there.



### Filter genotypes based on x3 depth per genotype

```bash
#!/bin/bash

# PBS directives 
#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N snps_qc3
#PBS -l select=1:ncpus=1:mem=10GB
#PBS -l walltime=00:30:00
#PBS -m e
#PBS -q defaultQ
#PBS -o snps_qc3.txt

# qsub ../snps_qc3.pbs

# load gatk
module load gatk/4.1.4.1

WORKING_DIR=/scratch/RDS-FSC-Heartworm_MLR-RW/extra_data1/analysis/mapping/filter/filter1.2

# set reference, vcf
REFERENCE=/project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_ALL/old_data/data/analysis/mapping/reference_di_wol_dog.fa
VCF=/scratch/RDS-FSC-Heartworm_MLR-RW/extra_data1/analysis/mapping/filter/filter1.2/Dirofilaria_immitis_Sep2023.vcf.gz

cd ${WORKING_DIR}

#Nuclear 
gatk VariantFiltration \
--reference ${REFERENCE} \
--variant ${VCF%.vcf.gz}.nuclearSNPs.filtered.vcf \
--genotype-filter-expression ' DP < 3 '  \
--genotype-filter-name "DP_lt3" \
--output ${VCF%.vcf.gz}.nuclearSNPs.DPfiltered.vcf

gatk SelectVariants \
--reference ${REFERENCE} \
--variant ${VCF%.vcf.gz}.nuclearSNPs.DPfiltered.vcf \
--set-filtered-gt-to-nocall \
--output ${VCF%.vcf.gz}.nuclearSNPs.DPfilterNoCall.vcf

#Mito
gatk VariantFiltration \
--reference ${REFERENCE} \
--variant ${VCF%.vcf.gz}.mitoSNPs.filtered.vcf \
--genotype-filter-expression ' DP < 3 '  \
--genotype-filter-name "DP_lt3" \
--output ${VCF%.vcf.gz}.mitoSNPs.DPfiltered.vcf

gatk SelectVariants \
--reference ${REFERENCE} \
--variant ${VCF%.vcf.gz}.mitoSNPs.DPfiltered.vcf \
--set-filtered-gt-to-nocall \
--output ${VCF%.vcf.gz}.mitoSNPs.DPfilterNoCall.vcf

#wolbachia
gatk VariantFiltration \
--reference ${REFERENCE} \
--variant ${VCF%.vcf.gz}.WbSNPs.filtered.vcf \
--genotype-filter-expression ' DP < 3 '  \
--genotype-filter-name "DP_lt3" \
--output ${VCF%.vcf.gz}.WbSNPs.DPfiltered.vcf

gatk SelectVariants \
--reference ${REFERENCE} \
--variant ${VCF%.vcf.gz}.WbSNPs.DPfiltered.vcf \
--set-filtered-gt-to-nocall \
--output ${VCF%.vcf.gz}.WbSNPs.DPfilterNoCall.vcf
```

### Now we apply a set of standard filters for population genomics

```bash
#!/bin/bash

# PBS directives 
#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N snps_qc4
#PBS -l select=1:ncpus=1:mem=4GB
#PBS -l walltime=00:05:00
#PBS -m e
#PBS -q defaultQ
#PBS -o snps_qc4.txt

# qsub ../snps_qc4.pbs

# load gatk
module load gatk/4.1.4.1
module load vcftools/0.1.14

WORKING_DIR=/scratch/RDS-FSC-Heartworm_MLR-RW/extra_data1/analysis/mapping/filter/filter1.2

VCF=/scratch/RDS-FSC-Heartworm_MLR-RW/extra_data1/analysis/mapping/filter/filter1.2/Dirofilaria_immitis_Sep2023.vcf.gz


cd ${WORKING_DIR}

#Nuclear variants
vcftools \
--vcf ${VCF%.vcf.gz}.nuclearSNPs.DPfilterNoCall.vcf \
--remove-filtered-geno-all \
--remove-filtered-all \
--min-alleles 2 \
--max-alleles 2 \
--hwe 1e-06 \
--maf 0.02 \
--recode \
--recode-INFO-all \
--out ${VCF%.vcf.gz}.nuclear_SNPs.final
# After filtering, kept 61 out of 61 Individuals
#Outputting VCF file...
#After filtering, kept 250306 out of a possible 506734 Sites

# Changed --maf from 0.05 to 0.02 to follow Javier's code

#--- nuclear SNPs
vcftools --vcf Dirofilaria_immitis_Sep2023.nuclear_SNPs.final.recode.vcf --remove-indels
#After filtering, kept 61 out of 61 Individuals
#After filtering, kept 250306 out of a possible 250306 Sites

#--- nuclear  INDELs
vcftools --vcf Dirofilaria_immitis_Sep2023.nuclear_SNPs.final.recode.vcf --keep-only-indels
#After filtering, kept 61 out of 61 Individuals
#After filtering, kept 0 out of a possible 250306 Sites



#Mitochondrial variants
vcftools \
--vcf ${VCF%.vcf.gz}.mitoSNPs.DPfilterNoCall.vcf \
--remove-filtered-geno-all \
--remove-filtered-all \
--min-alleles 2 \
--max-alleles 2 \
--maf 0.02 \
--recode \
--recode-INFO-all \
--out ${VCF%.vcf.gz}.mito_SNPs.final
#After filtering, kept 61 out of 61 Individuals
#Outputting VCF file...
#After filtering, kept 32 out of a possible 99 Sites


#--- mito SNPs
vcftools --vcf Dirofilaria_immitis_Sep2023.mito_SNPs.final.recode.vcf --remove-indels
#After filtering, kept 61 out of 61 Individuals
#After filtering, kept 32 out of a possible 32 Sites

#--- mito INDELs
vcftools --vcf Dirofilaria_immitis_Sep2023.mito_SNPs.final.recode.vcf --keep-only-indels
#After filtering, kept 61 out of 61 Individuals
#After filtering, kept 0 out of a possible 32 Sites



#wolbachia variants
vcftools \
--vcf ${VCF%.vcf.gz}.WbSNPs.DPfilterNoCall.vcf \
--remove-filtered-geno-all \
--remove-filtered-all \
--min-alleles 2 \
--max-alleles 2 \
--maf 0.02 \
--recode \
--recode-INFO-all \
--out ${VCF%.vcf.gz}.Wb_SNPs.final
#After filtering, kept 61 out of 61 Individuals
#Outputting VCF file...
#After filtering, kept 480 out of a possible 3493 Sites


#--- Wb SNPs
vcftools --vcf Dirofilaria_immitis_Sep2023.Wb_SNPs.final.recode.vcf --remove-indels
# After filtering, kept 61 out of 61 Individuals
#After filtering, kept 480 out of a possible 480 Sites

#--- Wb INDELs
vcftools --vcf Dirofilaria_immitis_Sep2023.Wb_SNPs.final.recode.vcf --keep-only-indels
# After filtering, kept 61 out of 61 Individuals
#After filtering, kept 0 out of a possible 480 Sites
```

### Now, we are filtering by missingness

```bash
#!/bin/bash

# PBS directives 
#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N snps_qc5
#PBS -l select=1:ncpus=1:mem=10GB
#PBS -l walltime=00:10:00
#PBS -m e
#PBS -q defaultQ
#PBS -o snps_qc5.txt

# qsub ../snps_qc5.pbs

# load modules
module load gatk/4.1.4.1
module load vcftools/0.1.14

WORKING_DIR=/scratch/RDS-FSC-Heartworm_MLR-RW/extra_data1/analysis/mapping/filter/filter1.2

VCF=/scratch/RDS-FSC-Heartworm_MLR-RW/extra_data1/analysis/mapping/filter/filter1.2/Dirofilaria_immitis_Sep2023.vcf.gz

cd ${WORKING_DIR}

#determine missingness per individual
vcftools --vcf ${VCF%.vcf.gz}.nuclear_SNPs.final.recode.vcf --out nuclear --missing-indv
vcftools --vcf ${VCF%.vcf.gz}.mito_SNPs.final.recode.vcf --out mito --missing-indv
vcftools --vcf ${VCF%.vcf.gz}.Wb_SNPs.final.recode.vcf --out wb --missing-indv
```

### Check the missingess in R



### Let's check different thresholds for each dataset

```bash
#!/bin/bash

# PBS directives 
#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N snps_qc6
#PBS -l select=1:ncpus=1:mem=10GB
#PBS -l walltime=00:20:00
#PBS -m e
#PBS -q defaultQ
#PBS -o snps_qc6.txt

# qsub ../snps_qc6.pbs

# load gatk
module load gatk/4.1.4.1
module load vcftools/0.1.14

WORKING_DIR=/scratch/RDS-FSC-Heartworm_MLR-RW/extra_data1/analysis/mapping/filter/filter1.2

VCF=/scratch/RDS-FSC-Heartworm_MLR-RW/extra_data1/analysis/mapping/filter/filter1.2/Dirofilaria_immitis_Sep2023.vcf.gz

cd ${WORKING_DIR}

# For nuclear variants
for i in 0.7 0.8 0.9 1; do
     vcftools --vcf ${VCF%.vcf.gz}.nuclear_SNPs.final.recode.vcf --keep nuclear_samplelist.keep --max-missing ${i} ;
done

# max-missing = 0.7
# After filtering, kept 61 out of 61 Individuals
#After filtering, kept 245071 out of a possible 250306 Sites

# max-missing = 0.8
# After filtering, kept 61 out of 61 Individuals
#After filtering, kept 240560 out of a possible 250306 Sites

# max-missing = 0.9
# After filtering, kept 61 out of 61 Individuals
#After filtering, kept 205621 out of a possible 250306 Sites

# max-missing = 1
# After filtering, kept 61 out of 61 Individuals
#After filtering, kept 7029 out of a possible 250306 Sites


# For mito variants
for i in 0.7 0.8 0.9 1; do
     vcftools --vcf ${VCF%.vcf.gz}.mito_SNPs.final.recode.vcf --keep mito_samplelist.keep --max-missing ${i} ;
done

# max-missing = 0.7
# After filtering, kept 61 out of 61 Individuals
#After filtering, kept 32 out of a possible 32 Sites

# max-missing = 0.8
# After filtering, kept 61 out of 61 Individuals
#After filtering, kept 32 out of a possible 32 Sites

# max-missing = 0.9
# After filtering, kept 61 out of 61 Individuals
#After filtering, kept 32 out of a possible 32 Sites

# max-missing = 1
# After filtering, kept 61 out of 61 Individuals
#After filtering, kept 30 out of a possible 32 Sites


# For Wb variants
for i in 0.7 0.8 0.9 1; do
     vcftools --vcf ${VCF%.vcf.gz}.Wb_SNPs.final.recode.vcf --keep wb_samplelist.keep --max-missing ${i} ;
done

# max-missing = 0.7
# After filtering, kept 61 out of 61 Individuals
#After filtering, kept 480 out of a possible 480 Sites

# max-missing = 0.8
# After filtering, kept 61 out of 61 Individuals
#After filtering, kept 480 out of a possible 480 Sites

# max-missing = 0.9
# After filtering, kept 61 out of 61 Individuals
#After filtering, kept 473 out of a possible 480 Sites

# max-missing = 1
# After filtering, kept 61 out of 61 Individuals
#After filtering, kept 64 out of a possible 480 Sites
```

Selecting a max missingness of 0.9 for nuclear, 1 for mito and 0.9 for Wb is sensible.

```bash
#!/bin/bash

# PBS directives 
#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N snps_qc7
#PBS -l select=1:ncpus=1:mem=10GB
#PBS -l walltime=00:10:00
#PBS -m e
#PBS -q defaultQ
#PBS -o snps_qc7.txt

# qsub ../snps_qc7.pbs

# load gatk
module load gatk/4.1.4.1
module load vcftools/0.1.14

WORKING_DIR=/scratch/RDS-FSC-Heartworm_MLR-RW/extra_data1/analysis/mapping/filter/filter1.2

VCF=/scratch/RDS-FSC-Heartworm_MLR-RW/extra_data1/analysis/mapping/filter/filter1.2/Dirofilaria_immitis_Sep2023.vcf.gz

cd ${WORKING_DIR}

# For nuclear
vcftools --vcf ${VCF%.vcf.gz}.nuclear_SNPs.final.recode.vcf \
     --keep nuclear_samplelist.keep \
     --max-missing 0.9 \
     --recode --recode-INFO-all \
     --out FINAL_SETS/nuclear_samples3x_missing0.9

# For mito
vcftools --vcf ${VCF%.vcf.gz}.mito_SNPs.final.recode.vcf \
     --keep mito_samplelist.keep \
     --max-missing 1 \
     --recode --recode-INFO-all \
     --out FINAL_SETS/mito_samples3x_missing1

# For wb
vcftools --vcf ${VCF%.vcf.gz}.Wb_SNPs.final.recode.vcf \
     --keep wb_samplelist.keep \
     --max-missing 0.9 \
     --recode --recode-INFO-all \
     --out FINAL_SETS/wb_samples3x_missing0.9
```

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

cd /scratch/RDS-FSC-Heartworm_MLR-RW/extra_data1/analysis/mapping/filter/filter1.2

vcftools --vcf FINAL_SETS/nuclear_samples3x_missing0.9.recode.vcf \
--chr dirofilaria_immitis_chr1 \
--chr dirofilaria_immitis_chr2 \
--chr dirofilaria_immitis_chr3 \
--chr dirofilaria_immitis_chr4 \
--recode --out FINAL_SETS/nuclear_samples3x_missing0.9.chr1to4
# After filtering, kept 61 out of 61 Individuals
#Outputting VCF file...
#After filtering, kept 153596 out of a possible 205621 Sites


vcftools --vcf FINAL_SETS/nuclear_samples3x_missing0.9.chr1to4.recode.vcf --remove-indels
# After filtering, kept 61 out of 61 Individuals
#After filtering, kept 153596 out of a possible 153596 Sites


vcftools --vcf FINAL_SETS/nuclear_samples3x_missing0.9.chr1to4.recode.vcf --keep-only-indels
# After filtering, kept 61 out of 61 Individuals
#After filtering, kept 0 out of a possible 153596 Sites- this makes sense because I removed the indels earlier and only focused on the SNPs.
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
#PBS -M rosemonde.power@sydney.edu.au

# qsub ../vcf_rename.pbs

# load modules
module load bcftools/1.11

cd /scratch/RDS-FSC-Heartworm_MLR-RW/extra_data1/analysis/mapping/filter/filter1.2/FINAL_SETS

bcftools reheader -s rename.txt nuclear_samples3x_missing0.9.chr1to4.recode.vcf -o nuclear_samples3x_missing0.9.chr1to4.recode.RENAMED.vcf

bcftools reheader -s rename.txt nuclear_samples3x_missing0.9.recode.vcf -o nuclear_samples3x_missing0.9.recode.RENAMED.vcf
```


rename.txt:

ERR034940	ERR034940	
ERR034941	ERR034941	
ERR034942	ERR034942	
ERR034943	ERR034943	
JS6277	Dog2.1	
JS6278	Dog2.2	
JS6279	Dog2.3	
JS6280	Dog2.4	
JS6281	JS6281	
JS6342	JS6342	
JS6343	JS6343	
JS6344	JS6344	
JS6345	JS6345	
JS6346	Fox1	
JS6347	Fox2	
JS6349	JS6349	
JS6350	JS6350	
JS6351	JS6351	
JS6352	JS6352	
JS6353	JS6353	
JS6354	JS6354	
JS6355	JS6355	
JS6356	JS6356	
JS6357	JS6357	
JS6358	JS6358	
JS6359	JS6359	
JS6360	JS6360	
JS6368	JS6368	
JS6369	JS6369	
JS6370	JS6370	
JS6597_DKDN230025568-1A_HHMLHDSX7_L1	JS6597	
JS6598_DKDN230025569-1A_HHMLHDSX7_L1	JS6598	
JS6599_DKDN230025570-1A_HHMLHDSX7_L1	JS6599	
JS6600_DKDN230025571-1A_HHMLHDSX7_L1	JS6600	
JS6601_DKDN230025572-1A_HHMLHDSX7_L1	Dog3.1	
JS6602_DKDN230025573-1A_HHMLHDSX7_L1	Dog3.2	
JS6603_DKDN230025574-1A_HHMLHDSX7_L1	Dog3.3	
JS6604_DKDN230025575-1A_HHMLHDSX7_L1	JS6604	
JS6605_DKDN230025576-1A_HHMLHDSX7_L1	JS6605	
JS6606_DKDN230025577-1A_HHMLHDSX7_L1	JS6606	
JS6607_DKDN230025578-1A_HHMLHDSX7_L4	JS6344-REP	
JS6608_DKDN230025579-1A_HHMLHDSX7_L4	Dog2.3-REP	
JS6609_DKDN230025580-1A_HHMLHDSX7_L4	Dog1.5-REP	
JS6610_DKDN230025581-1A_HHMLHDSX7_L1	JS6610	
JS6611_DKDN230025582-1A_HHMLHDSX7_L1	JS6611	
JS6612_DKDN230025583-1A_HHMLHDSX7_L1	JS6612	
JS6613_DKDN230025584-1A_HHMLHDSX7_L1	JS6613	
JS6614_DKDN230025585-1A_HHMLHDSX7_L1	JS6614	
JS6615_DKDN230025586-1A_HHMLHDSX7_L1	JS6615	
JS6616_DKDN230025587-1A_HHMLHDSX7_L1	JS6616	
JS6617_DKDN230025588-1A_HHMLHDSX7_L1	JS6617	
SRR10533236	SRR10533236	
SRR10533237	SRR10533237	
SRR10533238	SRR10533238	
SRR10533239	SRR10533239	
SRR10533240	SRR10533240	
SRR13154013	Dog1.1	
SRR13154014	Dog1.2	
SRR13154015	Dog1.3	
SRR13154016	Dog1.4	
SRR13154017	Dog1.5	



Check that it worked:

```bash
# List all sample names in the original VCF
bcftools query -l nuclear_samples3x_missing0.9.chr1to4.recode.RENAMED.vcf
```

Yep all the sample names are changed.


## PCA

```R
######################################################################

# PCA for filt1

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
setwd("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/extra_data/filter1/vcf/input")

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
vcf.in <- "C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/extra_data/filter1/vcf/input/nuclear_samples3x_missing0.9.chr1to4.recode.RENAMED.vcf"
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
jpeg("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/extra_data/filter1/vcf/pca/barplot_eigenvalues.jpg", width = 800, height = 600)
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
# 16.23
PC2.variance
# 10.17
PC3.variance
# 9.173
PC4.variance
# 7.336


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
ggtitle("Filter1: PC1 vs PC2")
  
PC1_PC2_plot

ggsave("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/extra_data/filter1/vcf/pca/PC1_PC2_plot.png", PC1_PC2_plot, height = 6, width = 8)

# PC1 vs PC3
PC1_PC3_plot <- ggplot(data, aes(x = EV1, y = EV3, color = POPULATION, label = SAMPLEID)) +
  theme_bw() +
  geom_point(alpha = 0.8, size = 2) +
  geom_text_repel(aes(label = SAMPLEID), box.padding = 0.05, point.padding = 0.01, segment.color = 'grey50', size = 3, hjust = 0, vjust = 0, max.overlaps = Inf, show.legend=FALSE) +
  labs(x = paste0("PC1 variance: ", round(pca$varprop[1] * 100, digits = 2), "%"),
       y = paste0("PC3 variance: ", round(pca$varprop[3] * 100, digits = 2), "%")) +
  scale_colour_pop() +
  theme(legend.position = "right") +  # Adjust legend position as needed
  ggtitle("Filter1: PC1 vs PC3")

PC1_PC3_plot

ggsave("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/extra_data/filter1/vcf/pca/PC1_PC3_plot.png", PC1_PC3_plot, height = 6, width = 8)

# PC1 vs PC4
PC1_PC4_plot <- ggplot(data, aes(x = EV1, y = EV4, color = POPULATION, label = SAMPLEID)) +
  theme_bw() +
  geom_point(alpha = 0.8, size = 2) +
  geom_text_repel(aes(label = SAMPLEID), box.padding = 0.05, point.padding = 0.01, segment.color = 'grey50', size = 3, hjust = 0, vjust = 0, max.overlaps = Inf, show.legend=FALSE) +
  labs(x = paste0("PC1 variance: ", round(pca$varprop[1] * 100, digits = 2), "%"),
       y = paste0("PC4 variance: ", round(pca$varprop[4] * 100, digits = 2), "%")) +
  scale_colour_pop() +
  theme(legend.position = "right") +  # Adjust legend position as needed
  ggtitle("Filter1: PC1 vs PC4")

PC1_PC4_plot

ggsave("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/extra_data/filter1/vcf/pca/PC1_PC4_plot.png", PC1_PC4_plot, height = 6, width = 8)

# PC1 vs PC5
PC1_PC5_plot <- ggplot(data, aes(x = EV1, y = EV5, color = POPULATION, label = SAMPLEID)) +
  theme_bw() +
  geom_point(alpha = 0.8, size = 2) +
  geom_text_repel(aes(label = SAMPLEID), box.padding = 0.05, point.padding = 0.01, segment.color = 'grey50', size = 3, hjust = 0, vjust = 0, max.overlaps = Inf, show.legend=FALSE) +
  labs(x = paste0("PC1 variance: ", round(pca$varprop[1] * 100, digits = 2), "%"),
       y = paste0("PC5 variance: ", round(pca$varprop[5] * 100, digits = 2), "%")) +
  scale_colour_pop() +
  theme(legend.position = "right") +  # Adjust legend position as needed
  ggtitle("Filter1: PC1 vs PC5")

PC1_PC5_plot

ggsave("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/extra_data/filter1/vcf/pca/PC1_PC5_plot.png", PC1_PC5_plot, height = 6, width = 8)

# PC1 vs PC6
PC1_PC6_plot <- ggplot(data, aes(x = EV1, y = EV6, color = POPULATION, label = SAMPLEID)) +
  theme_bw() +
  geom_point(alpha = 0.8, size = 2) +
  geom_text_repel(aes(label = SAMPLEID), box.padding = 0.05, point.padding = 0.01, segment.color = 'grey50', size = 3, hjust = 0, vjust = 0, max.overlaps = Inf, show.legend=FALSE) +
  labs(x = paste0("PC1 variance: ", round(pca$varprop[1] * 100, digits = 2), "%"),
       y = paste0("PC6 variance: ", round(pca$varprop[6] * 100, digits = 2), "%")) +
  scale_colour_pop() +
  theme(legend.position = "right") +  # Adjust legend position as needed
  ggtitle("Filter1: PC1 vs PC6")

PC1_PC6_plot

ggsave("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/extra_data/filter1/vcf/pca/PC1_PC6_plot.png", PC1_PC6_plot, height = 6, width = 8)

# PC2 vs PC3
PC2_PC3_plot <- ggplot(data, aes(x = EV2, y = EV3, color = POPULATION, label = SAMPLEID)) +
  theme_bw() +
  geom_point(alpha = 0.8, size = 2) +
  geom_text_repel(aes(label = SAMPLEID), box.padding = 0.05, point.padding = 0.01, segment.color = 'grey50', size = 3, hjust = 0, vjust = 0, max.overlaps = Inf, show.legend=FALSE) +
  labs(x = paste0("PC2 variance: ", round(pca$varprop[2] * 100, digits = 2), "%"),
       y = paste0("PC3 variance: ", round(pca$varprop[3] * 100, digits = 2), "%")) +
  scale_colour_pop() +
  theme(legend.position = "right") +  # Adjust legend position as needed
  ggtitle("Filter1: PC2 vs PC3")

PC2_PC3_plot

ggsave("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/extra_data/filter1/vcf/pca/PC2_PC3_plot.png", PC2_PC3_plot, height = 6, width = 8)

# PC2 vs PC4
PC2_PC4_plot <- ggplot(data, aes(x = EV2, y = EV4, color = POPULATION, label = SAMPLEID)) +
  theme_bw() +
  geom_point(alpha = 0.8, size = 2) +
  geom_text_repel(aes(label = SAMPLEID), box.padding = 0.05, point.padding = 0.01, segment.color = 'grey50', size = 3, hjust = 0, vjust = 0, max.overlaps = Inf, show.legend=FALSE) +
  labs(x = paste0("PC2 variance: ", round(pca$varprop[2] * 100, digits = 2), "%"),
       y = paste0("PC4 variance: ", round(pca$varprop[4] * 100, digits = 2), "%")) +
  scale_colour_pop() +
  theme(legend.position = "right") +  # Adjust legend position as needed
  ggtitle("Filter1: PC2 vs PC4")

PC2_PC4_plot

ggsave("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/extra_data/filter1/vcf/pca/PC2_PC4_plot.png", PC2_PC4_plot, height = 6, width = 8)

# PC3 vs PC4
PC3_PC4_plot <- ggplot(data, aes(x = EV3, y = EV4, color = POPULATION, label = SAMPLEID)) +
  theme_bw() +
  geom_point(alpha = 0.8, size = 2) +
  geom_text_repel(aes(label = SAMPLEID), box.padding = 0.05, point.padding = 0.01, segment.color = 'grey50', size = 3, hjust = 0, vjust = 0, max.overlaps = Inf, show.legend=FALSE) +
  labs(x = paste0("PC3 variance: ", round(pca$varprop[3] * 100, digits = 2), "%"),
       y = paste0("PC4 variance: ", round(pca$varprop[4] * 100, digits = 2), "%")) +
  scale_colour_pop() +
  theme(legend.position = "right") +  # Adjust legend position as needed
  ggtitle("Filter1: PC3 vs PC4")

PC3_PC4_plot

ggsave("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/extra_data/filter1/vcf/pca/PC3_PC4_plot.png", PC3_PC4_plot, height = 6, width = 8)


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
ggsave("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/extra_data/filter1/vcf/pca/PC1_line.png", height=2, width=8)

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
ggsave("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/extra_data/filter1/vcf/pca/PC2_line.png", height=2, width=8)

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
ggsave("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/extra_data/filter1/vcf/pca/PC3_line.png", height=2, width=8)

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
ggsave("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/extra_data/filter1/vcf/pca/PC4_line.png", height=2, width=8)

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
ggsave("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/extra_data/filter1/vcf/pca/PC5_line.png", height=2, width=8)

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
ggsave("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/extra_data/filter1/vcf/pca/PC6_line.png", height=2, width=8)
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

cd /scratch/RDS-FSC-Heartworm_MLR-RW/mapping/extra_data/analysis/mapping/filter/filter1/FINAL_SETS

# chr1
vcftools --vcf nuclear_samples3x_missing0.9.chr1to4.recode.RENAMED.vcf \
--chr dirofilaria_immitis_chr1 \
--recode --out nuclear_samples3x_missing0.9.chr1
# After filtering, kept 61 out of 61 Individuals
#Outputting VCF file...
#After filtering, kept 40623 out of a possible 153555 Sites

# chr2
vcftools --vcf nuclear_samples3x_missing0.9.chr1to4.recode.RENAMED.vcf \
--chr dirofilaria_immitis_chr2 \
--recode --out nuclear_samples3x_missing0.9.chr2
# After filtering, kept 61 out of 61 Individuals
#Outputting VCF file...
#After filtering, kept 35706 out of a possible 153555 Sites

# chr3
vcftools --vcf nuclear_samples3x_missing0.9.chr1to4.recode.RENAMED.vcf \
--chr dirofilaria_immitis_chr3 \
--recode --out nuclear_samples3x_missing0.9.chr3
# After filtering, kept 61 out of 61 Individuals
#Outputting VCF file...
#After filtering, kept 39513 out of a possible 153555 Sites

# chr4
vcftools --vcf nuclear_samples3x_missing0.9.chr1to4.recode.RENAMED.vcf \
--chr dirofilaria_immitis_chr4 \
--recode --out nuclear_samples3x_missing0.9.chr4
# After filtering, kept 61 out of 61 Individuals
#Outputting VCF file...
#After filtering, kept 37713 out of a possible 153555 Sites
```

## PCA per chromosome

```R
######################################################################

# Per chromosome

######################################################################



# Chromosome 1

snpgdsClose(genofile_chr1)
vcf.in_chr1 <- "C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/extra_data/filter1/vcf/input/nuclear_samples3x_missing0.9.chr1.recode.vcf"
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
  ggtitle("Filter1: PC1 vs PC2 - Chromosome 1")

PC1_PC2_plot_chr1

ggsave("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/extra_data/filter1/vcf/pca/PC1_PC2_plot_chr1.png", PC1_PC2_plot_chr1, height = 6, width = 8)



# Chromosome 2

snpgdsClose(genofile_chr2)
vcf.in_chr2 <- "C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/extra_data/filter1/vcf/input/nuclear_samples3x_missing0.9.chr2.recode.vcf"
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
  ggtitle("Filter1: PC1 vs PC2 - Chromosome 2")

PC1_PC2_plot_chr2

ggsave("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/extra_data/filter1/vcf/pca/PC1_PC2_plot_chr2.png", PC1_PC2_plot_chr2, height = 6, width = 8)




# Chromosome 3

snpgdsClose(genofile_chr3)
vcf.in_chr3 <- "C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/extra_data/filter1/vcf/input/nuclear_samples3x_missing0.9.chr3.recode.vcf"
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
  ggtitle("Filter1: PC1 vs PC2 - Chromosome 3")

PC1_PC2_plot_chr3

ggsave("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/extra_data/filter1/vcf/pca/PC1_PC2_plot_chr3.png", PC1_PC2_plot_chr3, height = 6, width = 8)


# Chromosome 4

snpgdsClose(genofile_chr4)
vcf.in_chr4 <- "C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/extra_data/filter1/vcf/input/nuclear_samples3x_missing0.9.chr4.recode.vcf"
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
  ggtitle("Filter1: PC1 vs PC2 - Chromosome 4")

PC1_PC2_plot_chr4

ggsave("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/extra_data/filter1/vcf/pca/PC1_PC2_plot_chr4.png", PC1_PC2_plot_chr4, height = 6, width = 8)

```

# Concordance between technical replicates

I had technical replicates for 3 samples. Look at the concordance between these. If there is a small % that are discordant, where are they? I can use this to help QC some of my SNPs by potentially excluding certain parts of the genome which may be more error-prone.

```bash
#!/bin/bash

# PBS directives 
#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N replicates
#PBS -l select=1:ncpus=1:mem=2GB
#PBS -l walltime=00:10:00
#PBS -m e
#PBS -q defaultQ
#PBS -o replicates.txt

# Submit job
## qsub ../replicates.pbs

# load modules
module load bcftools/1.11
module load tabix/0.2.6

cd /scratch/RDS-FSC-Heartworm_MLR-RW/extra_data1/analysis/mapping/filter/filter1/REPLICATES

###################################################################################
# Extract variants for each technical replicate
###################################################################################

# Extract variants for JS6344
bcftools view -s JS6344 ../FINAL_SETS/nuclear_samples3x_missing0.9.chr1to4.recode.RENAMED.vcf -o JS6344.vcf
bgzip JS6344.vcf
tabix -p vcf JS6344.vcf.gz

# Extract variants for JS6344-REP
bcftools view -s JS6344-REP ../FINAL_SETS/nuclear_samples3x_missing0.9.chr1to4.recode.RENAMED.vcf -o JS6344-REP.vcf
bgzip JS6344-REP.vcf
tabix -p vcf JS6344-REP.vcf.gz

# Extract variants for Dog2.3
bcftools view -s Dog2.3 ../FINAL_SETS/nuclear_samples3x_missing0.9.chr1to4.recode.RENAMED.vcf -o Dog2.3.vcf
bgzip Dog2.3.vcf
tabix -p vcf Dog2.3.vcf.gz

# Extract variants for Dog2.3-REP
bcftools view -s Dog2.3-REP ../FINAL_SETS/nuclear_samples3x_missing0.9.chr1to4.recode.RENAMED.vcf -o Dog2.3-REP.vcf
bgzip Dog2.3-REP.vcf
tabix -p vcf Dog2.3-REP.vcf.gz

# Extract variants for Dog1.5
bcftools view -s Dog1.5 ../FINAL_SETS/nuclear_samples3x_missing0.9.chr1to4.recode.RENAMED.vcf -o Dog1.5.vcf
bgzip Dog1.5.vcf
tabix -p vcf Dog1.5.vcf.gz

# Extract variants for Dog1.5-REP
bcftools view -s Dog1.5-REP ../FINAL_SETS/nuclear_samples3x_missing0.9.chr1to4.recode.RENAMED.vcf -o Dog1.5-REP.vcf
bgzip Dog1.5-REP.vcf
tabix -p vcf Dog1.5-REP.vcf.gz


###################################################################################
# Extract the genotypes for each sample
###################################################################################

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\n' -s JS6344 JS6344.vcf.gz > JS6344_genotypes.txt
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\n' -s JS6344-REP JS6344-REP.vcf.gz > JS6344-REP_genotypes.txt

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\n' -s Dog2.3 Dog2.3.vcf.gz > Dog2.3_genotypes.txt
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\n' -s Dog2.3-REP Dog2.3-REP.vcf.gz > Dog2.3-REP_genotypes.txt

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\n' -s Dog1.5 Dog1.5.vcf.gz > Dog1.5_genotypes.txt
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\n' -s Dog1.5-REP Dog1.5-REP.vcf.gz > Dog1.5-REP_genotypes.txt
```


```R
# Replicates concordance

# Load modules
library(data.table)
library(ggplot2)

# Set working directory
setwd("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/extra_data/filter1/replicates")

########################################################
# Dog 1.5
########################################################

# Load genotypes from the two files
genotypes_Dog1.5 <- fread("Dog1.5_genotypes.txt", header = FALSE, col.names = c("CHROM", "POS", "REF", "ALT", "Dog1.5_Genotype"))
genotypes_Dog1.5_REP <- fread("Dog1.5-REP_genotypes.txt", header = FALSE, col.names = c("CHROM", "POS", "REF", "ALT", "Dog1.5_REP_Genotype"))

# Merge data.tables by CHROM and POS
merged_genotypes_Dog1.5 <- merge(genotypes_Dog1.5, genotypes_Dog1.5_REP, by = c("CHROM", "POS", "REF", "ALT"))

# Create a new column for concordance status
merged_genotypes_Dog1.5$Concordance <- ifelse(merged_genotypes_Dog1.5$Dog1.5_Genotype == merged_genotypes_Dog1.5$Dog1.5_REP_Genotype, "Concordant", "Discordant")

# Calculate discordance
discordant_positions_Dog1.5 <- subset(merged_genotypes_Dog1.5, Dog1.5_Genotype != Dog1.5_REP_Genotype)
discordant_percentage_Dog1.5 <- nrow(discordant_positions_Dog1.5) / nrow(merged_genotypes_Dog1.5) * 100

# Format the discordance percentage to two decimal places
discordant_percentage_Dog1.5_2dp <- sprintf("%.2f", discordant_percentage_Dog1.5)
discordant_percentage_Dog1.5_2dp

# Output discordant positions and percentage
write.table(discordant_positions_Dog1.5, "discordant_genotypes_Dog1.5.txt", sep = "\t", quote = FALSE, row.names = FALSE)
cat("Discordance Percentage:", discordant_percentage_Dog1.5, "%\n")

# Define custom labels for chromosomes
chr_labels <- c(
  "dirofilaria_immitis_chr1" = "Chr1",
  "dirofilaria_immitis_chr2" = "Chr2",
  "dirofilaria_immitis_chr3" = "Chr3",
  "dirofilaria_immitis_chr4" = "Chr4")

# Create a bar plot
ggplot(data = merged_genotypes_Dog1.5, aes(x = CHROM, fill = Concordance)) +
  geom_bar() +
  labs(title = "Discordance by Chromosome: Dog 1.5 vs Dog 1.5 REP", subtitle = paste("Discordance Percentage:", discordant_percentage_Dog1.5_2dp, "%"), y = "Count") +
  scale_fill_manual(values = c("Concordant" = "lightblue1", "Discordant" = "tomato2")) +
  theme_minimal() +
  scale_x_discrete(labels = chr_labels) +
  theme(plot.title = element_text(face = "bold"))

ggsave("discordance_plot_Dog1.5.png", plot = last_plot(), width = 8, height = 6, dpi = 600)



########################################################
# Dog 2.3
########################################################

# Load genotypes from the two files
genotypes_Dog2.3 <- fread("Dog2.3_genotypes.txt", header = FALSE, col.names = c("CHROM", "POS", "REF", "ALT", "Dog2.3_Genotype"))
genotypes_Dog2.3_REP <- fread("Dog2.3-REP_genotypes.txt", header = FALSE, col.names = c("CHROM", "POS", "REF", "ALT", "Dog2.3_REP_Genotype"))

# Merge data.tables by CHROM and POS
merged_genotypes_Dog2.3 <- merge(genotypes_Dog2.3, genotypes_Dog2.3_REP, by = c("CHROM", "POS", "REF", "ALT"))

# Create a new column for concordance status
merged_genotypes_Dog2.3$Concordance <- ifelse(merged_genotypes_Dog2.3$Dog2.3_Genotype == merged_genotypes_Dog2.3$Dog2.3_REP_Genotype, "Concordant", "Discordant")

# Calculate discordance
discordant_positions_Dog2.3 <- subset(merged_genotypes_Dog2.3, Dog2.3_Genotype != Dog2.3_REP_Genotype)
discordant_percentage_Dog2.3 <- nrow(discordant_positions_Dog2.3) / nrow(merged_genotypes_Dog2.3) * 100

# Format the discordance percentage to two decimal places
discordant_percentage_Dog2.3_2dp <- sprintf("%.2f", discordant_percentage_Dog2.3)
discordant_percentage_Dog2.3_2dp

# Output discordant positions and percentage
write.table(discordant_positions_Dog2.3, "discordant_genotypes_Dog2.3.txt", sep = "\t", quote = FALSE, row.names = FALSE)
cat("Discordance Percentage:", discordant_percentage_Dog2.3, "%\n")

# Define custom labels for chromosomes
chr_labels <- c(
  "dirofilaria_immitis_chr1" = "Chr1",
  "dirofilaria_immitis_chr2" = "Chr2",
  "dirofilaria_immitis_chr3" = "Chr3",
  "dirofilaria_immitis_chr4" = "Chr4")

# Create a bar plot
ggplot(data = merged_genotypes_Dog2.3, aes(x = CHROM, fill = Concordance)) +
  geom_bar() +
  labs(title = "Discordance by Chromosome: Dog 2.3 vs Dog 2.3 REP", subtitle = paste("Discordance Percentage:", discordant_percentage_Dog2.3_2dp, "%"), y = "Count") +
  scale_fill_manual(values = c("Concordant" = "darkseagreen3", "Discordant" = "tomato2")) +
  theme_minimal() +
  scale_x_discrete(labels = chr_labels) +
  theme(plot.title = element_text(face = "bold"))

ggsave("discordance_plot_Dog2.3.png", plot = last_plot(), width = 8, height = 6, dpi = 600)


########################################################
# JS6344
########################################################

# Load genotypes from the two files
genotypes_JS6344 <- fread("JS6344_genotypes.txt", header = FALSE, col.names = c("CHROM", "POS", "REF", "ALT", "JS6344_Genotype"))
genotypes_JS6344_REP <- fread("JS6344-REP_genotypes.txt", header = FALSE, col.names = c("CHROM", "POS", "REF", "ALT", "JS6344_REP_Genotype"))

# Merge data.tables by CHROM and POS
merged_genotypes_JS6344 <- merge(genotypes_JS6344, genotypes_JS6344_REP, by = c("CHROM", "POS", "REF", "ALT"))

# Create a new column for concordance status
merged_genotypes_JS6344$Concordance <- ifelse(merged_genotypes_JS6344$JS6344_Genotype == merged_genotypes_JS6344$JS6344_REP_Genotype, "Concordant", "Discordant")

# Calculate discordance
discordant_positions_JS6344 <- subset(merged_genotypes_JS6344, JS6344_Genotype != JS6344_REP_Genotype)
discordant_percentage_JS6344 <- nrow(discordant_positions_JS6344) / nrow(merged_genotypes_JS6344) * 100

# Format the discordance percentage to two decimal places
discordant_percentage_JS6344_2dp <- sprintf("%.2f", discordant_percentage_JS6344)
discordant_percentage_JS6344_2dp

# Output discordant positions and percentage
write.table(discordant_positions_JS6344, "discordant_genotypes_JS6344.txt", sep = "\t", quote = FALSE, row.names = FALSE)
cat("Discordance Percentage:", discordant_percentage_JS6344, "%\n")

# Define custom labels for chromosomes
chr_labels <- c(
  "dirofilaria_immitis_chr1" = "Chr1",
  "dirofilaria_immitis_chr2" = "Chr2",
  "dirofilaria_immitis_chr3" = "Chr3",
  "dirofilaria_immitis_chr4" = "Chr4")

# Create a bar plot
ggplot(data = merged_genotypes_JS6344, aes(x = CHROM, fill = Concordance)) +
  geom_bar() +
  labs(title = "Discordance by Chromosome: JS6344 vs JS6344 REP", subtitle = paste("Discordance Percentage:", discordant_percentage_JS6344_2dp, "%"), y = "Count") +
  scale_fill_manual(values = c("Concordant" = "lightpink", "Discordant" = "tomato2")) +
  theme_minimal() +
  scale_x_discrete(labels = chr_labels) +
  theme(plot.title = element_text(face = "bold"))

ggsave("discordance_plot_JS6344.png", plot = last_plot(), width = 8, height = 6, dpi = 600)





########################################################
# Are discordant variants shared between replicates?
########################################################

library(VennDiagram)

# Read the discordant genotypes files for all three pairs
Dog1.5_data <- read.table("discordant_genotypes_Dog1.5.txt", header = TRUE, sep = "\t")
Dog2.3_data <- read.table("discordant_genotypes_Dog2.3.txt", header = TRUE, sep = "\t")
JS6344_data <- read.table("discordant_genotypes_JS6344.txt", header = TRUE, sep = "\t")

# Extract the "CHROM" and "POS" columns for each pair
Dog1.5_variants <- with(Dog1.5_data, paste(CHROM, POS, sep = ":"))
Dog2.3_variants <- with(Dog2.3_data, paste(CHROM, POS, sep = ":"))
JS6344_variants <- with(JS6344_data, paste(CHROM, POS, sep = ":"))

# Create a Venn diagram with three circles
venn.plot <- venn.diagram(
  x = list(Dog1.5_variants, Dog2.3_variants, JS6344_variants),
  category.names = c("Dog1.5", "Dog2.3", "JS6344"),
  filename = NULL,
  output = TRUE,
  euler.d = TRUE, # For a more circular layout
  fontfamily = "sans", # Change font family to Arial (or your preferred font)
  category.cex = 2, # Increase category label size
  output = TRUE,
  lwd = 1, # Increase line width
  fill = c("lightblue1", "darkseagreen3", "lightpink"), # Customize fill colors
  alpha = 0.5, # Set transparency
  cat.pos = 0, # Position of category labels (0 = centered)
  margin = 0.05 # Adjust the margin
)

# Display the Venn diagram
grid.draw(venn.plot)

ggsave("discordance_venn_diagram.png", plot = venn.plot, width = 6, height = 6, dpi = 600)

```

Now I may want to exclude these discordant variant sites from my cohort VCF

```R

```
```bash
bcftools view -T ^discordant_variants.txt -O v -o filtered_cohort.vcf input_cohort.vcf
```





## Generate an ALL SITES variant set for running pixy properly

```bash
#!/bin/bash

# PBS directives 
#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N gatk_allsites.pbs
#PBS -l select=1:ncpus=1:mem=80GB
#PBS -l walltime=24:00:00
#PBS -m abe
#PBS -q defaultQ
#PBS -o gatk_allsites.txt
#PBS -M rosemonde.power@sydney.edu.au

# qsub ../gatk_allsites.pbs

# ALL SITES (including non-variant sites)

WORKING_DIR=/scratch/RDS-FSC-Heartworm_MLR-RW/mapping/extra_data/analysis/mapping/vcf
cd "${WORKING_DIR}"

cohort=Dirofilaria_immitis_Sep2023
ref=/project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_ALL/data/analysis/mapping/reference_di_wol_dog.fa

gvcf=/scratch/RDS-FSC-Heartworm_MLR-RW/mapping/extra_data/analysis/mapping/vcf/${cohort}.g.vcf.gz
vcf=/scratch/RDS-FSC-Heartworm_MLR-RW/mapping/extra_data/analysis/mapping/vcf/${cohort}_ALLSITES.vcf.gz

# Load modules
module load gatk/4.2.1.0

# Genotype cohort vcf
gatk --java-options "-Xmx28g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
        GenotypeGVCFs \
        -R ${ref} \
        -all-sites \
        -V ${gvcf} \
        -O ${vcf}
```

This is taking a while. Try running this on Sanger system to see which one finishes first.

```bash
# Transfer cohort.g.vcf.gz and associated tbi file into RDS
dt-script -P RDS-FSC-Heartworm_MLR-RW \
-m 50GB \
--from /scratch/RDS-FSC-Heartworm_MLR-RW/extra_data1/analysis/mapping/vcf/Dirofilaria_immitis_Sep2023.g.vcf.gz \
--to /rds/PRJ-Heartworm_MLR/scratch/RDS-FSC-Heartworm_MLR-RW/mapping/extra_data/analysis/mapping/vcf/Dirofilaria_immitis_Sep2023.g.vcf.gz
# done

dt-script -P RDS-FSC-Heartworm_MLR-RW \
-m 10GB \
--from /scratch/RDS-FSC-Heartworm_MLR-RW/extra_data1/analysis/mapping/vcf/Dirofilaria_immitis_Sep2023.g.vcf.gz.tbi \
--to /rds/PRJ-Heartworm_MLR/scratch/RDS-FSC-Heartworm_MLR-RW/mapping/extra_data/analysis/mapping/vcf/Dirofilaria_immitis_Sep2023.g.vcf.gz.tbi
# done

# Transfer data from RDS to Sanger
# reference
pscp -P 2227 Z:/PRJ-Heartworm_MLR/project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_ALL/data/analysis/mapping/reference_di_wol_dog.fa.gz rp24@localhost:/lustre/scratch125/pam/teams/team333/rp24/Diro
# done

# g.vcf.gz.tbi
pscp -P 2227 Z:/PRJ-Heartworm_MLR/scratch/RDS-FSC-Heartworm_MLR-RW/mapping/extra_data/analysis/mapping/vcf/Dirofilaria_immitis_Sep2023.g.vcf.gz.tbi rp24@localhost:/lustre/scratch125/pam/teams/team333/rp24/Diro


# ALL SITES (including non-variant sites)

WORKING_DIR=/project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_ALL/data/analysis
cd "${WORKING_DIR}"

cohort=Dirofilaria_immitis_June2023
config=/project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_ALL/data/analysis/info.txt
ref=/project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_ALL/data/analysis/mapping/reference_di_wol_dog.fa

gvcf=/scratch/RDS-FSC-Heartworm_MLR-RW/mapping/vcf/${cohort}_ALLSITES.g.vcf.gz
vcf=/scratch/RDS-FSC-Heartworm_MLR-RW/mapping/vcf/${cohort}_ALLSITES.vcf.gz

# Load modules
#module load gatk/4.2.1.0
module load gatk/4.1.4.1 #they only had this version available

# Genotype cohort vcf
bsub.py 10 gatk_allsites "gatk GenotypeGVCFs \
        -R reference_di_wol_dog.fa \
        -all-sites \
        -V Dirofilaria_immitis_Sep2023.g.vcf.gz \
        -O Dirofilaria_immitis_Sep2023_ALLSITES.vcf.gz"
```

This ran out of walltime and there was a fatal error. Try splitting it into jobs like Javier did.




## Split my cohort GVCF file into individual GVCFs for each chromosome (Sanger)

#bsub.py 10 SPLIT_GVCF "../SPLIT_GVCFS.sh"

```bash
#Creating environmental variables and loading modules
# working dir
WORKING_DIR=/lustre/scratch125/pam/teams/team333/rp24/Diro/data
# load gatk
module load gatk/4.1.4.1
#vcftools
module load vcftools/0.1.16-c4
# also need htslib for tabix
module load common-apps/htslib/1.9.229


REFERENCE=${WORKING_DIR}/REFERENCE/reference_di_wol_dog.fa
REF_DIR=${WORKING_DIR}/REFERENCE

# made a sequences list to allow splitting jobs per scaffold/contig
SEQUENCE_LIST=${WORKING_DIR}/VARIANTS/sequences.list
ulimit -c unlimited

COHORT_GVCF=${WORKING_DIR}/VARIANTS/Dirofilaria_immitis_Sep2023.g.vcf.gz
OUTPUT_DIR=${WORKING_DIR}/VARIANTS/


# Split my cohort GVCF file into individual GVCF files for each chromosome
while IFS= read -r CHR; do
    gatk SelectVariants \
        -R ${REFERENCE} \
        -V ${COHORT_GVCF} \
        -L ${CHR} \
        -O ${OUTPUT_DIR}/${CHR}.g.vcf.gz
done < "${SEQUENCE_LIST}"
```

That seemed to work.


## Run genotyping on each chromosome GVCF

bsub.py 1 genotype_jobs "../genotype_jobs.sh"

```bash
#Creating environmental variables and loading modules
# working dir
WORKING_DIR=/lustre/scratch125/pam/teams/team333/rp24/Diro/data
# load gatk
module load gatk/4.1.4.1
#vcftools
module load vcftools/0.1.16-c4
# also need htslib for tabix
module load common-apps/htslib/1.9.229

REFERENCE=${WORKING_DIR}/REFERENCE/reference_di_wol_dog.fa

# made a sequences list to allow splitting jobs per scaffold/contig
SEQUENCE_LIST=${WORKING_DIR}/VARIANTS/sequences.list
ulimit -c unlimited

OUTPUT_DIR=${WORKING_DIR}/VARIANTS/

# Split each chromosome up into separate jobs, and run genotyping on each individually.
n=1
while read SEQUENCE; do
     echo -e "gatk GenotypeGVCFs \
     -R ${REFERENCE} \
     -V ${OUTPUT_DIR}/${SEQUENCE}.g.vcf.gz \
     --intervals ${SEQUENCE} \
     --all-sites \
     --heterozygosity 0.015 \
     --indel-heterozygosity 0.01 \
     --annotation DepthPerAlleleBySample --annotation Coverage --annotation ExcessHet --annotation FisherStrand --annotation MappingQualityRankSumTest --annotation StrandOddsRatio --annotation RMSMappingQuality --annotation ReadPosRankSumTest --annotation DepthPerSampleHC --annotation QualByDepth \
     -O ${n}.${SEQUENCE}.cohort.vcf.gz" > run_hc_genotype.${SEQUENCE}.tmp.job_${n};
     let "n+=1";
done < ${SEQUENCE_LIST}
```

```bash
chmod a+x run_hc_genotype*

mkdir LOGFILES

# setup job conditions
JOBS=$( ls -1 run_hc_* | wc -l )
ID="U$(date +%s)"

# run
bsub -q long -R'span[hosts=1] select[mem>20000] rusage[mem=20000]' -n 4 -M20000 -J GATK_HC_GENOTYPE_${ID}_[1-$JOBS] -e LOGFILES/GATK_HC_GENOTYPE_${ID}_[1-$JOBS].e -o LOGFILES/GATK_HC_GENOTYPE_${ID}_[1-$JOBS].o "./run_hc_*\$LSB_JOBINDEX"
```

That seemed to work. But [6] (scaffold 1) is missing? I ran this twice and it was missing both times. Don't know what could've gone wrong there...
I moved all the all-site cohort vcf files into the data/VARIANTS/ALLSITES folder to keep things tidy. All the log files were moved to the scripts/logs folder.

## Bring the files together

bsub.py 4 genotype_concat "../genotype_concat.sh"

```bash
#vcftools
module load vcftools/0.1.16-c4
# also need htslib for tabix
module load common-apps/htslib/1.9.229

cd /lustre/scratch125/pam/teams/team333/rp24/Diro/data/VARIANTS/ALLSITES

# make list of vcfs
ls -1 *.cohort.vcf.gz | sort -n > vcf_files.list

# merge them
vcf-concat --files vcf_files.list > dirofilaria_immitis.cohort.allsites.vcf;
     bgzip dirofilaria_immitis.cohort.allsites.vcf;
     tabix -p vcf dirofilaria_immitis.cohort.allsites.vcf.gz
```

Do this later:
```bash
# clean up
rm run*
rm ^[0-9]* #remove files in the current directory whose filenames start with one or more digits (0-9)
rm *.g.vcf.gz*
```

## Now, let's filter nuclear variants and invariants

bsub.py 10 genotype_filter "../genotype_filter.sh"

```bash
#Creating environmental variables and loading modules
# working dir
WORKING_DIR=/lustre/scratch125/pam/teams/team333/rp24/Diro/data
cd ${WORKING_DIR}

# load gatk
module load gatk/4.1.4.1
#vcftools
module load vcftools/0.1.16-c4
# also need htslib for tabix
module load common-apps/htslib/1.9.229
# bsub.py
module load bsub.py

REFERENCE=${WORKING_DIR}/REFERENCE/reference_di_wol_dog.fa

VCF=${WORKING_DIR}/VARIANTS/ALLSITES/dirofilaria_immitis.cohort.allsites.vcf.gz
MIT_CONTIG=dirofilaria_immitis_chrMtDNA
WB_CONTIG=dirofilaria_immitis_chrWb

vcftools --gzvcf ${WORKING_DIR}/VARIANTS/ALLSITES/dirofilaria_immitis.cohort.allsites.vcf.gz --remove-indels
#After filtering, kept 61 out of 61 Individuals
#After filtering, kept 87974473 out of a possible 88674417 Sites
# Some weird warnings appeared - not sure if important?

#select nuclear invariants
gatk SelectVariants \
--reference ${REFERENCE} \
--variant ${VCF} \
--select-type-to-include NO_VARIATION \
--exclude-intervals ${MIT_CONTIG} \
--exclude-intervals ${WB_CONTIG} \
--output ${VCF%.vcf.gz}.nuclearINVARIANTs.vcf

vcftools --vcf ${WORKING_DIR}/VARIANTS/ALLSITES/dirofilaria_immitis.cohort.allsites.nuclearINVARIANTs.vcf --remove-indels
#After filtering, kept 61 out of 61 Individuals
#After filtering, kept 86562415 out of a possible 86601320 Sites
# Some weird warnings appeared - not sure if important?
```
Some guides recommend not to apply pop genomics filters for invariants


## In any way, lets merge the nuclear invariants with nuclear SNPs previously filtered

bsub.py 1 filter_nuclear_INVARITANTvcftools "../filter_nuclear_INVARITANTvcftools.sh"

```bash
#Creating environmental variables and loading modules
# working dir
WORKING_DIR=/lustre/scratch125/pam/teams/team333/rp24/Diro/data
cd ${WORKING_DIR}

# load gatk
module load gatk/4.1.4.1
#vcftools
module load vcftools/0.1.16-c4
# also need htslib for tabix
module load common-apps/htslib/1.9.229
# bsub.py
module load bsub.py

REFERENCE=${WORKING_DIR}/REFERENCE/reference_di_wol_dog.fa

VCF=${WORKING_DIR}/VARIANTS/ALLSITES/dirofilaria_immitis.cohort.allsites.vcf.gz

#Keep all the samples of the nuclear dataset and then merge with the SNPs-INDELs
vcftools --vcf ${VCF%.vcf.gz}.nuclearINVARIANTs.vcf \
--keep ${WORKING_DIR}/VARIANTS/nuclear_samplelist.keep \
--recode --out nuclearINVARIANTS_allsamples
#After filtering, kept 61 out of 61 Individuals
#After filtering, kept 86601320 out of a possible 86601320 Sites
# Weird warnings again?
```



Load default bsub.py
bsub.py 2 merge_nuclear_VARIANTsandINVARIANTs "../merge_nuclear_VARIANTsandINVARIANTs.sh"

```bash
#Creating environmental variables and loading modules
# working dir
WORKING_DIR=/lustre/scratch125/pam/teams/team333/rp24/Diro/data
cd ${WORKING_DIR}

# load gatk
module load gatk/4.1.4.1
#vcftools
module load vcftools/0.1.16-c4
# also need htslib for tabix
module load common-apps/htslib/1.9.229
# bsub.py
module load bsub.py

REFERENCE=${WORKING_DIR}/REFERENCE/reference_di_wol_dog.fa

VCF=${WORKING_DIR}/VARIANTS/ALLSITES/dirofilaria_immitis.cohort.allsites.vcf.gz


#And merge with the already filtered variants
gatk MergeVcfs \
--INPUT ${WORKING_DIR}/VARIANTS/FINAL_SETS/nuclear_samples3x_missing0.9.recode.vcf \
--INPUT ${WORKING_DIR}/VARIANTS/ALLSITES/nuclearINVARIANTS_allsamples.recode.vcf \
--OUTPUT ${WORKING_DIR}/VARIANTS/ALLSITES/nuclearVARIANTsandINVARIANTs_allsamples.recode.vcf
```
Successfully completed.



Now rename the output file, it still has the old long sample names.

bsub.py 1 rename_vcf "../rename_vcf.sh"

```bash
# load modules
module load bcftools/1.14--h88f3f91_0

cd /lustre/scratch125/pam/teams/team333/rp24/Diro/data/VARIANTS/ALLSITES

bcftools reheader -s rename.txt nuclearVARIANTsandINVARIANTs_allsamples.recode.vcf -o nuclearVARIANTsandINVARIANTs_allsamples.recode.RENAMED.vcf
```
Successfully completed.



bsub.py 1 merge_nuclear_VARIANTsandINVARIANTs_2 "../merge_nuclear_VARIANTsandINVARIANTs_2.sh"

```bash
#Creating environmental variables and loading modules
# working dir
WORKING_DIR=/lustre/scratch125/pam/teams/team333/rp24/Diro/data
cd ${WORKING_DIR}

# load gatk
module load gatk/4.1.4.1
#vcftools
module load vcftools/0.1.16-c4
# also need htslib for tabix
module load common-apps/htslib/1.9.229
# bsub.py
module load bsub.py

REFERENCE=${WORKING_DIR}/REFERENCE/reference_di_wol_dog.fa

#Also, we will select only the chrX to chr4, avoiding the scaffolds
vcftools \
--vcf ${WORKING_DIR}/VARIANTS/ALLSITES/nuclearVARIANTsandINVARIANTs_allsamples.recode.RENAMED.vcf \
--chr dirofilaria_immitis_chrX \
--chr dirofilaria_immitis_chr1 \
--chr dirofilaria_immitis_chr2 \
--chr dirofilaria_immitis_chr3 \
--chr dirofilaria_immitis_chr4 \
--remove-indels \
--recode --out ${WORKING_DIR}/VARIANTS/ALLSITES/nuclearSNPssandINVARIANTs.chrxto4
#After filtering, kept _ out of _ Individuals
#After filtering, kept _ out of a possible _ Sites

bgzip ${WORKING_DIR}/VARIANTS/ALLSITES/nuclearSNPssandINVARIANTs.chrxto4.recode.vcf;
tabix -p vcf ${WORKING_DIR}/VARIANTS/ALLSITES/nuclearSNPssandINVARIANTs.chrxto4.recode.vcf.gz
```
Successfully completed.
Then moved these output files into the FINAL_SETS folder. Now continue this work back in the workbook_extra_data1.md document.


## Run pixy using that file

```bash
#Install miniconda
cd /nfs/users/nfs_r/rp24
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
source /nfs/users/nfs_r/rp24/miniconda3/bin/activate

#Create and set the environment
conda create --name pixy
conda init --all
conda activate pixy

#Installing pixy
conda install -c conda-forge pixy
conda install -c bioconda htslib

#Few environemntal variables
WORKING_DIR=/lustre/scratch125/pam/teams/team333/rp24/Diro/data
VCF=${WORKING_DIR}/VARIANTS/FINAL_SETS/nuclearSNPssandINVARIANTs.chrxto4.recode.vcf.gz

cd ${WORKING_DIR}/PIXY
module load bsub.py

#Generate a population file
#country_pop.list:
ERR034940	USA
ERR034941	Italy
ERR034942	Italy
ERR034943	Italy
Dog2.1	Australia
Dog2.2	Australia
Dog2.3	Australia
Dog2.4	Australia
JS6281	Australia
JS6342	Australia
JS6343	Australia
JS6344	Australia
JS6345	Australia
Fox1	Australia
Fox2	Australia
JS6349	Australia
JS6350	Australia
JS6351	Australia
JS6352	Australia
JS6353	Australia
JS6354	Australia
JS6355	Australia
JS6356	Australia
JS6357	Australia
JS6358	Australia
JS6359	Australia
JS6360	Australia
JS6368	Australia
JS6369	Australia
JS6370	Australia
JS6597	USA
JS6598	USA
JS6599	USA
JS6600	USA
Dog3.1	Australia
Dog3.2	Australia
Dog3.3	Australia
JS6604	Thailand
JS6605	Thailand
JS6606	Australia
JS6344-REP	Australia
Dog2.3-REP	Australia
Dog1.5-REP	Australia
JS6610	Australia
JS6611	Australia
JS6612	Australia
JS6613	Australia
JS6614	Australia
JS6615	Australia
JS6616	Australia
JS6617	Australia
SRR10533236	USA
SRR10533237	USA
SRR10533238	USA
SRR10533239	USA
SRR10533240	USA
Dog1.1	Australia
Dog1.2	Australia
Dog1.3	Australia
Dog1.4	Australia
Dog1.5	Australia


# Submit a job to get the diversity per country:
bsub.py --queue long --threads 20 20 pixy_country \
"pixy --stats pi fst dxy \
--vcf ${VCF} \
--populations country_pop.list \
--window_size 100000 \
--n_cores 20 \
--output_prefix country"
```
Successfully completed.

Download the files and use it as input in R to estimate the values and generate plots.


## Nucleotide diversity (pi)

```R
library(tidyverse)
library(ggsci)
library(ggpubr)
library(patchwork)
library(ggridges)
library(dplyr)

setwd("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/extra_data/filter1.2/pixy")



##############################################################
# PI
##############################################################

# ChatGPT summary:

# Pi is used to measure genetic diversity within a population. It is the average pairwise nucleotide diversity (which is a measure of the extent of genetic variation at the nucleotide level).

# get nucleotide diversity (pi) data from pixy output
pi_data <- read.table("country_pi.txt", header=T)
pi_data$chromosome <- str_remove(pi_data$chromosome, 'dirofilaria_immitis_')

# subset
pi_data_Australia <- pi_data %>%
  filter(pop=="Australia")
pi_data_USA <- pi_data %>%
  filter(pop=="USA")
pi_data_Italy <- pi_data %>%
  filter(pop=="Italy")
pi_data_Thailand <- pi_data %>%
  filter(pop=="Thailand")



# filter pi data to remove small scaffolds not in the linkage groups, and to number the rows per group to help with plotting
pi_data <- pi_data %>%
  group_by(pop) %>%
  mutate(position = 1:n())

# get position for vertical lines used in plot to delineate the linkage groups
pi_data %>%
  group_by(chromosome) %>%
  summarise(max = max(position, na.rm = TRUE))

'
# A tibble: 5 × 2
  chromosome   max
  <chr>      <int>
1 chr1         440
2 chr2         592
3 chr3         744
4 chr4         885
5 chrX         283
'

#Let's add the chr type variable
pi_data <- pi_data %>%
  mutate(chr_type = ifelse(str_detect(chromosome, "X"), "sexchr", "autosome"))

# calculate the median Pi and checking the ratio of sex-to-autosome diversity. Should be about 0.75, as Trichuris is XX/XY
pi_data_sex_median <- pi_data %>%
  group_by(chr_type) %>%
  summarise(median = median(avg_pi, na.rm = TRUE))

# sex/autosome ratio is 0.0001615095 / 0.0003148628 = 0.512952 (far off 0.75 expected of diversity on sex chromosome relative to autosome)


pi_data_pop_sex_median <- pi_data %>%
  group_by(pop, chr_type) %>%
  summarise(median = median(avg_pi, na.rm = TRUE))

# plot 1 - genome wide plots per population
plot_1_pi <- ggplot(pi_data, aes(position*100000, avg_pi, col=chromosome)) +
  geom_point(size=1) +
  facet_grid(pop~.) +
  scale_color_tron() +
  theme_bw() +
  theme(legend.position="top", 
        legend.title = element_blank()) +
  labs(x="Genomic Position", y="Nucleotide Diversity (Pi)")
plot_1_pi

# plot 2 - density plots of pi per group
plot_2_pi <- ggplot(pi_data, aes(avg_pi, chr_type, fill=chr_type), guide="none") +
  geom_density_ridges(quantile_lines=TRUE, quantile_fun=function(x,...)median(x), size=0.5) +
  theme_bw() + theme(legend.position = "none", axis.text.y=element_blank()) +
  facet_grid(pop~.) +
  xlim(0, 0.005) +
  scale_x_continuous(breaks = c(0, 0.0025, 0.005)) +
  scale_fill_npg() +
  labs(x="Nucleotide Diversity (Pi)", y="Density")
plot_2_pi

# bring it together
plot_1_pi + plot_2_pi +  plot_layout(widths = c(5, 1))
ggsave("plots_genomewide_and_density_Pi.png", width=9, height=6)

#Now a boxplot of the pi value per population

scale_colour_pop <- function(...){
  ggplot2:::manual_scale(
    'colour', 
    values = setNames(
      c('dodgerblue4',
        'mediumpurple2',
        'springgreen4',
        'orangered3'), 
      c('Australia', 'Italy', 
        'Thailand',
        'USA')), 
    ...
  )
}

boxplot_pi <- ggplot(pi_data, aes(pop, avg_pi, col=pop)) +
  geom_jitter(size = 1, alpha = 0.5) +
  geom_boxplot(fill=NA, col="black") +
  labs(x = "Population" , y = "Nucleotide diversity (Pi)", colour = "Population") +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_colour_pop ()+
  ylim(0, 0.003)

boxplot_pi
# Removed 4 rows containing non-finite values (`stat_boxplot()`). 
# Removed 138 rows containing missing values (`geom_point()`). 

ggsave("plots_boxplot_pop_Pi.png", width=4, height=4)

# Removed 4 rows containing non-finite values (`stat_boxplot()`). 
# Removed 171 rows containing missing values (`geom_point()`). 




#let's check the statistical significance of the differences between pi values
#Let's explore the normality
pi_data %>% filter(chromosome != 'chrX') %>%
  ggplot(., aes(x = avg_pi, colour = pop)) +
  geom_histogram() +
  theme_bw() +
  facet_grid(. ~ pop) +
  scale_colour_pop ()
ggsave("plots_normality_Pi.png", width=8, height=4) # this generates a histogram plot of the avg_pi values

#Confirm the absence of normality by shapiro
pi_data_Australia <- pi_data_Australia %>% filter(chromosome != 'chrX') # filter the data frame to exclude rows with ChrX
Australia_shapiro <- shapiro.test(pi_data_Australia$avg_pi) # perform Shapiro-Wilk normality test on avg_pi. The shapiro-Wilk test is used to assess whether a dataset follows a normal distribution
print(Australia_shapiro)
#Shapiro-Wilk normality test
#data:  pi_data_Australia$avg_pi
#W = 0.84701, p-value < 2.2e-16'  ##p is <0.05 so data significantly deviates from a normal distribution



pi_data_Italy <- pi_data_Italy %>% filter(chromosome != 'chrX')
Italy_shapiro <- shapiro.test(pi_data_Italy$avg_pi)
print(Italy_shapiro)
#Shapiro-Wilk normality test
#data:  pi_data_Italy$avg_pi
#W = 0.59426, p-value < 2.2e-16

pi_data_Thailand <- pi_data_Thailand %>% filter(chromosome != 'chrX')
Thailand_shapiro <- shapiro.test(pi_data_Thailand$avg_pi)
print(Thailand_shapiro)
# Shapiro-Wilk normality test
#data:  pi_data_Thailand$avg_pi
#W = 0.79255, p-value < 2.2e-16

pi_data_USA <- pi_data_USA %>% filter(chromosome != 'chrX')
USA_shapiro <- shapiro.test(pi_data_USA$avg_pi)
print(USA_shapiro)
# Shapiro-Wilk normality test
#data:  pi_data_USA$avg_pi
#W = 0.80547, p-value < 2.2e-16



#Wilcoxson test for everyone
# The Wilcoxson test determines whether the distribution of paired data significantly differs from a hypothetical median/reference value. It's a non-parametric alternative to the paired t-test (which doesn't assume that the data follows a normal distribution)
wilcox.test(pi_data_Australia$avg_pi, pi_data_Italy$avg_pi)
# W = 331814, p-value < 2.2e-16
wilcox.test(pi_data_Australia$avg_pi, pi_data_Thailand$avg_pi)
# W = 265525, p-value < 2.2e-16
wilcox.test(pi_data_Australia$avg_pi, pi_data_USA$avg_pi)
# W = 189143, p-value = 0.1881
wilcox.test(pi_data_Italy$avg_pi, pi_data_Thailand$avg_pi)
# W = 96563, p-value < 2.2e-16
wilcox.test(pi_data_Italy$avg_pi, pi_data_USA$avg_pi)
# W = 36040, p-value < 2.2e-16
wilcox.test(pi_data_Thailand$avg_pi, pi_data_USA$avg_pi)
# W = 104414, p-value < 2.2e-16


#let's generate some dataframe for the estatistics
# subset
pi_data_Australia <- pi_data %>%
  filter(pop=="Australia") %>%
  filter(chromosome != 'chrX')

pi_data_Italy <- pi_data %>%
  filter(pop=="Italy") %>%
  filter(chromosome != 'chrX') 

pi_data_Thailand <- pi_data %>%
  filter(pop=="Thailand") %>%
  filter(chromosome != 'chrX')

pi_data_USA <- pi_data %>%
  filter(pop=="USA") %>%
  filter(chromosome != 'chrX')

```


## Dxy and Fst

```R
##############################################################
# Dxy and Fst
##############################################################

# ChatGPT summary:

# Dxy (Nei's genetic distance) is the average number of nucleotide differences per site between two populations. It provides a measure of genetic divergence between populations. Dxy ranges from 0 (no difference) to a max value that depends on the length of the sequences. Higher Dxy = greater genetic differentiation.

# Fst (fixation index) is a measure of population genetic structure that quantifies the proportion of genetic variation within populations compared to the total genetic variation in the entire set of populations. Fst values range from 0 to 1 (0 = no genetic differentiation, so all pops have identical genetic structure. 1 = complete genetic differentiation, so each population has unique genetic variants and there's no gene flow/shared variation).


# load data
dxy_data <- read.table("country_dxy.txt", header=T)
dxy_data_renamed <- dxy_data %>% mutate(pop1 = case_when(pop1 == 'Australia' ~ 'AUS', pop1 == 'Italy' ~ 'ITL', pop1 == 'Thailand' ~ 'THAI', TRUE ~ pop1))
dxy_data_renamed <- dxy_data_renamed %>% mutate(pop2 = case_when(pop2 == 'Australia' ~ 'AUS', pop2 == 'Italy' ~ 'ITL', pop2 == 'Thailand' ~ 'THAI', TRUE ~ pop2))

fst_data <- read.table("country_fst.txt", header=T)
fst_data_renamed <- fst_data %>% mutate(pop1 = case_when(pop1 == 'Australia' ~ 'AUS', pop1 == 'Italy' ~ 'ITL', pop1 == 'Thailand' ~ 'THAI', TRUE ~ pop1))
fst_data_renamed <- fst_data_renamed %>% mutate(pop2 = case_when(pop2 == 'Australia' ~ 'AUS', pop2 == 'Italy' ~ 'ITL', pop2 == 'Thailand' ~ 'THAI', TRUE ~ pop2))

# add some columns
dxy_data_renamed$data_type <- "Dxy"
dxy_data_renamed <- mutate(dxy_data_renamed,
                   comparison = paste(pop1, pop2, sep = '_v_'))
fst_data_renamed$data_type <- "Fst"
fst_data_renamed <- mutate(fst_data_renamed,
                   comparison = paste(pop1, pop2, sep = '_v_'))
# subset and merge dataframes
dxy_data_sub <- dxy_data_renamed %>% select(comparison,data_type,chromosome,window_pos_1,window_pos_2,avg_dxy)
colnames(dxy_data_sub) <- c("comparison","data_type","chromosome","window_pos_1","window_pos_2","value")

fst_data_sub <- fst_data_renamed %>% select(comparison,data_type,chromosome,window_pos_1,window_pos_2,avg_wc_fst)
colnames(fst_data_sub) <- c("comparison","data_type","chromosome","window_pos_1","window_pos_2","value")

# cheeky fix to get matching rows from both datasets
tmp_dxy <- semi_join(dxy_data_sub, fst_data_sub, by=c("comparison", "chromosome", "window_pos_1", "window_pos_2"))
tmp_fst <- semi_join(fst_data_sub, dxy_data_sub,  by=c("comparison", "chromosome", "window_pos_1", "window_pos_2"))

# join the datasets together to create a single dataframe
data <- full_join(tmp_dxy, tmp_fst)
data$chromosome <- str_remove(data$chromosome, 'dirofilaria_immitis_')

# add numbering to help with plotting
data <- data %>%
  group_by(comparison, data_type) %>%
  mutate(position = 1:n())

# add sex chromosome information
data <- data %>%
  mutate(chr_type = ifelse(str_detect(chromosome, "X"), "sexchr", "autosome"))


# summarise median data for Fst and Dxy
data %>%
  group_by(comparison,data_type) %>%
  summarise(median = median(value, na.rm = TRUE))
'
# A tibble: 12 × 3
# Groups:   comparison [6]
   comparison data_type   median
   <chr>      <chr>        <dbl>
 1 AUS_v_THAI Dxy       0.000518
 2 AUS_v_THAI Fst       0.187   
 3 ITL_v_AUS  Dxy       0.000528
 4 ITL_v_AUS  Fst       0.338   
 5 ITL_v_THAI Dxy       0.000588
 6 ITL_v_THAI Fst       0.728   
 7 USA_v_AUS  Dxy       0.000583
 8 USA_v_AUS  Fst       0.206   
 9 USA_v_ITL  Dxy       0.000521
10 USA_v_ITL  Fst       0.333   
11 USA_v_THAI Dxy       0.000523
12 USA_v_THAI Fst       0.220  
'


#Plotting Dxy

# plot 1 - genome wide plots per comparison
plot_1_dxy <- data %>%
  filter(data_type =='Dxy')%>%
  ggplot(., aes(position*100000, value, col=chromosome)) +
  geom_point(size=1) +
  facet_grid(comparison~.) +
  scale_color_tron() +
  theme_bw() +
  theme(legend.position="top", 
        legend.title = element_blank()) +
  labs(x="Genomic Position", y="Dxy")
plot_1_dxy

# plot 2 - density plots of dxy per group
plot_2_dxy <- data %>%
  filter(data_type=='Dxy')%>%
  ggplot(aes(value, chr_type, fill=chr_type), guide="none") +
  geom_density_ridges(quantile_lines=TRUE, quantile_fun=function(x,...)median(x), size=0.5) +
  theme_bw() + theme(legend.position = "none", axis.text.y=element_blank()) +
  facet_grid(comparison~.) +
  scale_fill_npg() +
  labs(x="Dxy", y="Density")
plot_2_dxy

# bring it together
dxy_plot <- plot_1_dxy + plot_2_dxy +  plot_layout(widths = c(5, 1))
dxy_plot
ggsave("plots_genomewide_and_density_dxy.png", width=9, height=6)

#some additional plots
boxplot_dxy <- data %>%
  filter(data_type =='Dxy') %>% 
  ggplot(., aes(comparison, value, col=comparison)) +
  geom_jitter(size = 1, alpha = 0.5) +
  geom_boxplot(fill=NA, col="black") +
  labs(x = "Population" , y = "Dxy", colour = "Population") +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
boxplot_dxy

density_dxy <- data %>%
  filter(data_type =='Dxy') %>% 
  ggplot(., aes(x=value, group = comparison, fill=comparison)) +
  geom_density(adjust=1.5, alpha=.4) +
  labs(x = "Dxy" , y = "density", colour = "Comparison") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
density_dxy

# AUS_v_THAI & ITL_v_AUS
scatter_dxy_a <-  data %>%
  filter(data_type == 'Dxy') %>%
  as.data.frame() %>%
  select('position', 'comparison', 'value', 'chromosome')%>%
  as_tibble() %>%
  pivot_wider(names_from = comparison, values_from = value) %>%
  ggplot(., aes(x=AUS_v_THAI, y = ITL_v_AUS, col=chromosome)) +
  geom_point(alpha=.4) +
  labs(x = "AUS_v_THAI" , y = "ITL_v_AUS", colour = "Comparison") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
scatter_dxy_a

# AUS_v_THAI & USA_v_AUS
scatter_dxy_b <-  data %>%
  filter(data_type == 'Dxy') %>%
  as.data.frame() %>%
  select('position', 'comparison', 'value', 'chromosome')%>%
  as_tibble() %>%
  pivot_wider(names_from = comparison, values_from = value) %>%
  ggplot(., aes(x=AUS_v_THAI, y = USA_v_AUS, col=chromosome)) +
  geom_point(alpha=.4) +
  labs(x = "AUS_v_THAI" , y = "USA_v_AUS", colour = "Comparison") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
scatter_dxy_b

# ITL_v_AUS & USA_v_AUS
scatter_dxy_c <-  data %>%
  filter(data_type == 'Dxy') %>%
  as.data.frame() %>%
  select('position', 'comparison', 'value', 'chromosome')%>%
  as_tibble() %>%
  pivot_wider(names_from = comparison, values_from = value) %>%
  ggplot(., aes(x=ITL_v_AUS, y = USA_v_AUS, col=chromosome)) +
  geom_point(alpha=.4) +
  labs(x = "ITL_v_AUS" , y = "USA_v_AUS", colour = "Comparison") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
scatter_dxy_c

genome_pos_dxy_a <-  data %>%
  filter(data_type == 'Dxy') %>%
  as.data.frame() %>%
  select('position', 'comparison', 'value', 'chromosome')%>%
  as_tibble() %>%
  pivot_wider(names_from = comparison, values_from = value) %>%
  mutate("xx" = ITL_v_AUS - USA_v_AUS,
         "yy" = ITL_v_AUS - USA_v_ITL) %>%
  ggplot(., aes(position*100000, xx, col=chromosome)) +
  geom_point(size=1) +
  scale_color_tron() +
  theme_bw() +
  theme(legend.position="top", 
        legend.title = element_blank()) +
  labs(x="Position", y="AUS vs ITL - AUS vs USA")
genome_pos_dxy_a

genome_pos_dxy_b <-  data %>%
  filter(data_type == 'Dxy') %>%
  as.data.frame() %>%
  select('position', 'comparison', 'value', 'chromosome')%>%
  as_tibble() %>%
  pivot_wider(names_from = comparison, values_from = value) %>%
  mutate("xx" = ITL_v_AUS - USA_v_AUS,
         "yy" = ITL_v_AUS - USA_v_ITL) %>%
  ggplot(., aes(position*100000, yy, col=chromosome)) +
  geom_point(size=1) +
  scale_color_tron() +
  theme_bw() +
  theme(legend.position="top", 
        legend.title = element_blank()) +
  labs(x="Position", y="AUS vs ITL - ITL vs USA")
genome_pos_dxy_b

genome_pos_dxy_aa <-  data %>%
  filter(data_type == 'Dxy') %>%
  as.data.frame() %>%
  select('position', 'comparison', 'value', 'chromosome')%>%
  as_tibble() %>%
  pivot_wider(names_from = comparison, values_from = value) %>%
  mutate("xx" = ITL_v_AUS / USA_v_AUS,
         "yy" = ITL_v_AUS / USA_v_ITL) %>%
  ggplot(., aes(position*100000, xx, col=chromosome)) +
  geom_point(size=1) +
  scale_color_tron() +
  theme_bw() +
  theme(legend.position="top", 
        legend.title = element_blank()) +
  labs(x="Position", y="AUS vs ITL / AUS vs USA")
genome_pos_dxy_aa

genome_pos_dxy_bb <-  data %>%
  filter(data_type == 'Dxy') %>%
  as.data.frame() %>%
  select('position', 'comparison', 'value', 'chromosome')%>%
  as_tibble() %>%
  pivot_wider(names_from = comparison, values_from = value) %>%
  mutate("xx" = ITL_v_AUS / USA_v_AUS,
         "yy" = ITL_v_AUS / USA_v_ITL) %>%
  ggplot(., aes(position*100000, yy, col=chromosome)) +
  geom_point(size=1) +
  scale_color_tron() +
  theme_bw() +
  theme(legend.position="top", 
        legend.title = element_blank()) +
  labs(x="Position", y="AUS vs ITL / ITL vs USA")
genome_pos_dxy_bb

genome_pos_dxy_cc <-  data %>%
  filter(data_type == 'Dxy') %>%
  as.data.frame() %>%
  select('position', 'comparison', 'value', 'chromosome')%>%
  as_tibble() %>%
  pivot_wider(names_from = comparison, values_from = value) %>%
  mutate("xx" = ITL_v_AUS / USA_v_AUS,
         "yy" = ITL_v_AUS / USA_v_ITL,
         "zz" = USA_v_AUS / USA_v_ITL) %>%
  ggplot(., aes(position*100000, zz, col=chromosome)) +
  geom_point(size=1) +
  scale_color_tron() +
  theme_bw() +
  theme(legend.position="top", 
        legend.title = element_blank()) +
  labs(x="Position", y="AUS_v_USA / ITL vs USA")
genome_pos_dxy_cc

ggarrange(scatter_dxy_a, scatter_dxy_b, common.legend = T)
ggsave("scatter_dxy_ab.png", width=9, height=6)

ggarrange(genome_pos_dxy_a, genome_pos_dxy_b, common.legend = T, ncol = 1)
ggsave("genome_pos_dxy_ab.png", width=9, height=6)

ggarrange(genome_pos_dxy_aa, genome_pos_dxy_bb, genome_pos_dxy_cc,
          common.legend = T, ncol = 1)
ggsave("genome_pos_dxy_aabbcc.png", width=9, height=6)

#Trying to plot lines
#getting chr positions
chr <- c('X', '1', '2', '3', '4')

for (i in chr) {
  x <- data %>%
    filter(chromosome == paste0('chr', i))
  print(quantile(x$position))
}

data$comparison <- str_replace_all(data$comparison, '_', ' ')
data$comparison <- str_replace_all(data$comparison, 'v', 'vs')

dxy_lineplot <- data %>%
  filter(data_type=='Dxy')%>%
  ggplot(., aes(position*100000, value, col=comparison)) +
  geom_point(size=0.2, alpha = 0.1) +
  geom_vline(xintercept=c(281*100000,
                          438*100000,
                          590*100000,
                          740*100000),
             size=1, linetype="dashed", col='grey41')+
  annotate(geom="text", x=140*100000, y=0.004, label="ChrX")+
  annotate(geom="text", x=359.5*100000, y=0.004, label="Chr1")+
  annotate(geom="text", x=545*100000, y=0.004, label="Chr2")+
  annotate(geom="text", x=665*100000, y=0.004, label="Chr3")+
  annotate(geom="text", x=810.5*100000, y=0.004, label="Chr4")+
  geom_line() +
  scale_color_brewer(type = 'div', palette = 'Accent') +
  theme_bw() +
  theme(legend.position="top", 
        legend.title = element_blank()) +
  labs(x="Genomic Position", y="Dxy")

dxy_lineplot
ggsave("dxy_lineplot.png", width=9, height=6)

# Plotting Fst

# plot 1 - genome wide plots per population
plot_1_fst <- data %>%
  filter(data_type=='Fst')%>%
  ggplot(., aes(position*100000, value, col=chromosome)) +
  geom_point(size=1) +
  facet_grid(comparison~.) +
  scale_color_tron() +
  theme_bw() +
  theme(legend.position="top", 
        legend.title = element_blank()) +
  labs(x="Genomic Position", y="Fst")
plot_1_fst


# plot 2 - density plots of Fst per group
plot_2_fst <- data %>%
  filter(data_type=='Fst')%>%
  ggplot(aes(value, chr_type, fill=chr_type), guide="none") +
  geom_density_ridges(quantile_lines=TRUE, quantile_fun=function(x,...)median(x), size=0.5) +
  theme_bw() + theme(legend.position = "none", axis.text.y=element_blank()) +
  facet_grid(comparison~.) +
  scale_fill_npg() +
  labs(x="Fst", y="Density")
plot_2_fst

# bring it together
plot_1_fst + plot_2_fst +  plot_layout(widths = c(5, 1))
ggsave("plots_genomewide_and_density_fst.png", width=9, height=6)
```


## Explore the Australian populations

bsub.py 1 select_aus "../select_aus.sh"

```bash
# load gatk
module load gatk/4.1.4.1
#vcftools
module load vcftools/0.1.16-c4
# also need htslib for tabix
module load common-apps/htslib/1.9.229
# bsub.py
module load bsub.py

#Let's select only the Australian samples
WORKING_DIR=/lustre/scratch125/pam/teams/team333/rp24/Diro/data
cd ${WORKING_DIR}/VARIANTS/FINAL_SETS/

vcftools --gzvcf nuclearSNPssandINVARIANTs.chrxto4.recode.vcf.gz \
--keep ../AUS_samplelist.keep \
--recode --out AUS_nuclearSNPssandINVARIANTs.chrxto4

bgzip AUS_nuclearSNPssandINVARIANTs.chrxto4.recode.vcf;
tabix -p vcf AUS_nuclearSNPssandINVARIANTs.chrxto4.recode.vcf.gz
```
Successfully completed

```bash
#Now, let's run pixy
VCF=${WORKING_DIR}/VARIANTS/FINAL_SETS/AUS_nuclearSNPssandINVARIANTs.chrxto4.recode.vcf.gz

# AUS_pop.list
Dog2.1	Sydney
Dog2.2	Sydney
Dog2.3	Sydney
Dog2.4	Sydney
JS6281	Lockhart River Cooktown
JS6342	Brisbane
JS6343	Brisbane
JS6344	Brisbane
JS6345	Brisbane
Fox1	Brisbane
Fox2	Brisbane
JS6349	Cairns
JS6350	Cairns
JS6351	Townsville
JS6352	Townsville
JS6353	Townsville
JS6354	Townsville
JS6355	Townsville
JS6356	Townsville
JS6357	Townsville
JS6358	Townsville
JS6359	Townsville
JS6360	Townsville
JS6368	Townsville
JS6369	Townsville
JS6370	Townsville
Dog3.1	Sydney
Dog3.2	Sydney
Dog3.3	Sydney
JS6606	Rockhampton
JS6344-REP	Brisbane
Dog2.3-REP	Sydney
Dog1.5-REP	Sydney
JS6610	Brisbane
JS6611	Brisbane
JS6612	Townsville
JS6613	Townsville
JS6614	Townsville
JS6615	Townsville
JS6616	Townsville
JS6617	Townsville
Dog1.1	Sydney
Dog1.2	Sydney
Dog1.3	Sydney
Dog1.4	Sydney
Dog1.5	Sydney

conda activate pixy
module load bsub.py

cd /lustre/scratch125/pam/teams/team333/rp24/Diro/data/PIXY

#submitting jobs
bsub.py --queue long --threads 20 20 aus_pop \
"pixy --stats pi fst dxy \
--vcf ${VCF} \
--populations AUS_pop.list \
--window_size 100000 \
--n_cores 20 \
--output_prefix AUS_pop"
```
Successfully completed


## Nucleotide diversity (pi)

```R
library(tidyverse)
library(ggsci)
library(ggpubr)
library(patchwork)
library(ggridges)
library(dplyr)

setwd("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/extra_data/filter1.2/pixy_aus")

##############################################################
# PI
##############################################################

# Per population
# get nucleotide diversity (pi) data from pixy output
pi_data <- read.table("aus_pop_pi.txt", header=T)
pi_data$chromosome <- str_remove(pi_data$chromosome, 'dirofilaria_immitis_')

# filter pi data to remove small scaffolds not in the linkage groups, and to number the rows per group to help with plotting
pi_data <- pi_data %>%
  group_by(pop) %>%
  mutate(position = 1:n())

#LEt's add the chr type variable
pi_data <- pi_data %>%
  mutate(chr_type = ifelse(str_detect(chromosome, "X"), "sexchr", "autosome"))

# calculate the median Pi and checking the ratio of sex-to-autosome diversity. Should be about 0.75, as Trichuris is XX/XY
pi_data_sex_median <- pi_data %>%
  group_by(chr_type) %>%
  summarise(median = median(avg_pi, na.rm = TRUE))
'

'
# calculation

pi_data_pop_sex_median <- pi_data %>%
  group_by(pop, chr_type) %>%
  summarise(median = median(avg_pi, na.rm = TRUE))
'

'

#Now a boxplot of the pi value per population
pi_data$pop <- factor(pi_data$pop, 
                           levels = c('Lockhart River Cooktown', 'Cairns', 'Townsville', 'Rockhampton', 
                                      'Brisbane', 'Sydney'))

scale_colour_javier_PCA <- function(...){
  ggplot2:::manual_scale(
    'colour', 
    values = setNames(
      c('royalblue1', 'blue3',
        'turquoise3',
        'green',
        'lightgoldenrod1', 'darkgoldenrod1', 'orange', 'orange2',
        'orange3', 'tomato1', 'red2', 'darkred'), 
      c('QUE', 'NSW', 
        'PAV',
        'SIC',
        'MCH', 'ILL', 'TEN', 'ARK', 
        'GEO', 'MIP', 'TEX', 'LOU')), 
    ...
  )
}

boxplot_pi_usa_pop <- ggplot(pi_data, aes(pop, avg_pi, col=pop)) +
  geom_jitter(size = 1, alpha = 0.5) +
  geom_boxplot(fill=NA, col="grey15", linewidth = 1) +
  labs(x = "Population" , y = "Nucleotide diversity (Pi)") +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_color_npg() +
  ylim(0, 0.003) + scale_colour_javier_PCA()

#Per sample type
# get nucleotide diversity (pi) data from pixy output
pi_data <- read.table("data/usa_sampletype_pi.txt", header=T)
pi_data$chromosome <- str_remove(pi_data$chromosome, 'dirofilaria_immitis_')

# filter pi data to remove small scaffolds not in the linkage groups, and to number the rows per group to help with plotting
pi_data <- pi_data %>%
  group_by(pop) %>%
  mutate(position = 1:n())

#LEt's add the chr type variable
pi_data <- pi_data %>%
  mutate(chr_type = ifelse(str_detect(chromosome, "X"), "sexchr", "autosome"))
# calculate the median Pi and checking the ratio of sex-to-autosome diversity. Should be about 0.75, as Trichuris is XX/XY
pi_data_sex_median <- pi_data %>%
  group_by(chr_type) %>%
  summarise(median = median(avg_pi, na.rm = TRUE))
'

'
# 0.000378 / 0.000704 = 0.5369318 (far off 0.75 expected of diversity on sex chromosome relative to autosome)

pi_data_pop_sex_median <- pi_data %>%
  group_by(pop, chr_type) %>%
  summarise(median = median(avg_pi, na.rm = TRUE))
'

'
#Now a boxplot of the pi value per population
boxplot_pi_usa_sampletype <- ggplot(pi_data, aes(pop, avg_pi, col=pop, shape = pop)) +
  geom_jitter(size = 2, alpha = 0.5) +
  geom_boxplot(fill=NA, col="grey15", linewidth = 1) +
  labs(x = "Population" , y = "Nucleotide diversity (Pi)") +
  theme_bw() +
  scale_color_tron()+
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylim(0, 0.003)

#and arranging
ggarrange(boxplot_pi_usa_pop, boxplot_pi_usa_sampletype, labels = c('a', 'b'),
          ncol = 2, widths =  c(1.5, 1))
ggsave("FigS8_plots_diversity_USA.png", width=9, height=5.5)
```






