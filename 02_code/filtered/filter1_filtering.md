# Filter 1 - Moderate stringency

Filter the VCF file with moderate stringency.

## SNPs QC

```bash
# Load modules
module load vcftools/0.1.16-c4
module load bsub.py/0.42.1
module load common-apps/htslib/1.9.229

# Copy VCF into my own folder
cp /lustre/scratch125/pam/teams/team333/sd21/dirofilaria_immitis/POPGEN/NEWDATA_2024/VARIANTS/gatk_hc_DIMMITIS_POPGEN/GATK_HC_MERGED/DIMMITIS_POPGEN.cohort.2024-06-27.* /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/04_VARIANTS/FILTER1

VCF=/lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/04_VARIANTS/FILTER1/DIMMITIS_POPGEN.cohort.2024-06-27.vcf.gz

cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/04_VARIANTS/FILTER1/NO_OUTGROUPS

# Remove outgroup samples so I have a VCF of D. immitis samples only
bsub.py 2 vcf_no_outgroups "vcftools --gzvcf ${VCF} --keep dimmitis_samplelist.keep --recode --recode-INFO-all --out DIMMITIS_POPGEN.cohort.2024-06-27.NO_OUTGROUPS"
bsub.py --done "vcf_no_outgroups" 1 vcf_no_outgroups_bgzip "bgzip DIMMITIS_POPGEN.cohort.2024-06-27.NO_OUTGROUPS.recode.vcf"
bsub.py --done "vcf_no_outgroups_bgzip" 1 vcf_no_outgroups_index "tabix -p vcf DIMMITIS_POPGEN.cohort.2024-06-27.NO_OUTGROUPS.recode.vcf.gz"
# Update VCF environmental variable
VCF=/lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/04_VARIANTS/FILTER1/NO_OUTGROUPS/DIMMITIS_POPGEN.cohort.2024-06-27.NO_OUTGROUPS.recode.vcf.gz
```

### Querying SNP and INDEL QC profiles to determine thresholds for filters

Adopted from Javier's paper.

```bash
bsub.py --done "vcf_no_outgroups_index" 10 run_snps_qc "run_snps_qc.sh"
```

```bash
# Load modules
module load gatk/4.1.4.1

WORKING_DIR=/lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/04_VARIANTS/FILTER1/NO_OUTGROUPS

# set reference, vcf, and mitochondrial and Wb contig
REFERENCE=/lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/01_REF/dimmitis_WSI_2.2.fa
VCF=/lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/04_VARIANTS/FILTER1/NO_OUTGROUPS/DIMMITIS_POPGEN.cohort.2024-06-27.NO_OUTGROUPS.recode.vcf.gz
MIT_CONTIG=dirofilaria_immitis_chrMtDNA
WB_CONTIG=dirofilaria_immitis_chrWb
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
--variant DIMMITIS_POPGEN.cohort.2024-06-27.NO_OUTGROUPS.recode.nuclearSNPs.vcf \
--fields CHROM --fields POS --fields QUAL --fields QD --fields DP --fields MQ --fields MQRankSum --fields FS --fields ReadPosRankSum --fields SOR \
--output GVCFall_nuclearSNPs.table

# make a table of nuclear INDEL data
gatk VariantsToTable \
--reference ${REFERENCE} \
--variant DIMMITIS_POPGEN.cohort.2024-06-27.NO_OUTGROUPS.recode.nuclearINDELs.vcf \
--fields CHROM --fields POS --fields QUAL --fields QD --fields DP --fields MQ --fields MQRankSum --fields FS --fields ReadPosRankSum --fields SOR \
--output GVCFall_nuclearINDELs.table

# make a table of mito SNP data
gatk VariantsToTable \
--reference ${REFERENCE} \
--variant DIMMITIS_POPGEN.cohort.2024-06-27.NO_OUTGROUPS.recode.mitoSNPs.vcf \
--fields CHROM --fields POS --fields QUAL --fields QD --fields DP --fields MQ --fields MQRankSum --fields FS --fields ReadPosRankSum --fields SOR \
--output GVCFall_mitoSNPs.table

# make a table of mito INDEL data
gatk VariantsToTable \
--reference ${REFERENCE} \
--variant DIMMITIS_POPGEN.cohort.2024-06-27.NO_OUTGROUPS.recode.mitoINDELs.vcf \
--fields CHROM --fields POS --fields QUAL --fields QD --fields DP --fields MQ --fields MQRankSum --fields FS --fields ReadPosRankSum --fields SOR \
--output GVCFall_mitoINDELs.table

# make a table of Wb SNP data
gatk VariantsToTable \
--reference ${REFERENCE} \
--variant DIMMITIS_POPGEN.cohort.2024-06-27.NO_OUTGROUPS.recode.WbSNPs.vcf \
--fields CHROM --fields POS --fields QUAL --fields QD --fields DP --fields MQ --fields MQRankSum --fields FS --fields ReadPosRankSum --fields SOR \
--output GVCFall_WbSNPs.table

# make a table of Wb INDEL data
gatk VariantsToTable \
--reference ${REFERENCE} \
--variant DIMMITIS_POPGEN.cohort.2024-06-27.NO_OUTGROUPS.recode.WbINDELs.vcf \
--fields CHROM --fields POS --fields QUAL --fields QD --fields DP --fields MQ --fields MQRankSum --fields FS --fields ReadPosRankSum --fields SOR \
--output GVCFall_WbINDELs.table
```


### Make some density plots of the data and get quantiles in R

```R
# Make some density plots of the data and get quantiles in R

## RUN 06/07/2024 R version 4.4.1 (2024-06-14 ucrt) -- "Race for Your Life"

# load libraries
library(patchwork) # v.1.2.0
require(data.table) # v.1.15.4
library(tidyverse) # v.2.0.0
library(gridExtra) # v.2.3

sessionInfo()

setwd("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/batch4/FILTER1/snps_qc")

set.seed(123)

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
  # Removed 2 rows containing missing values or values outside the scale range (`geom_point()`). 
  
  # gatk hardfilter: SNP MQRankSum < -20
  MQRankSum <-
    ggplot(VCF_mito) +
    geom_point( aes(x=log10(MQRankSum), y=Variant, colour=Variant)) +
    theme_bw() +
    labs(title=paste0("mitochondrial",": MQRankSum"))
  MQRankSum
  # Removed 1621 rows containing missing values (`geom_point()`).  
  
  
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
  # Removed 633 rows containing missing values (`geom_point()`). 
  
  
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
  # Removed 2 rows containing missing values or values outside the scale range (`geom_point()`). 
  
  # gatk hardfilter: SNP MQRankSum < -20
  MQRankSum <-
    ggplot(VCF_mito) +
    geom_jitter( aes(x=log10(MQRankSum), y=Variant, colour=Variant)) +
    theme_bw() +
    labs(title=paste0("mitochondrial",": MQRankSum"))
  MQRankSum
  # Removed 1621 rows containing missing values (`geom_point()`). 
  
  
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
  # Removed 633 rows containing missing values (`geom_point()`). 
  
  
  plot <- QUAL + DP + QD + FS + MQ + MQRankSum + SOR + ReadPosRankSum + plot_layout(ncol=2)
  
  print(plot)
  
  ggsave(paste0("plot_mitochondrial_indiv_variant_summaries_jitter.png"), height=20, width=15, type="cairo")
```


### Attending to the quantiles, thresholds for specific parameters are established

```bash
# load modules
module load bsub.py/0.42.1
module load gatk/4.1.4.1

WORKING_DIR=/lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/04_VARIANTS/FILTER1/NO_OUTGROUPS
cd ${WORKING_DIR}

# set reference
REFERENCE=/lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/01_REF/dimmitis_WSI_2.2.fa

#Nuclear
bsub.py 1 filter_nuclearSNPs "gatk VariantFiltration \
--reference ${REFERENCE} \
--variant DIMMITIS_POPGEN.cohort.2024-06-27.NO_OUTGROUPS.recode.nuclearSNPs.vcf \
--filter-expression 'QUAL < 40 || DP < 1837 || DP > 11972 || MQ < 21.52 || SOR > 7.112 || QD < 0.780 || FS > 8.275 || MQRankSum < -5.870 || ReadPosRankSum < -2.792 || ReadPosRankSum > 2.240' \
--filter-name "SNP_filtered" \
--output DIMMITIS_POPGEN.cohort.2024-06-27.NO_OUTGROUPS.nuclearSNPs.filtered.vcf"

bsub.py 1 filter_nuclearINDELS "gatk VariantFiltration \
--reference ${REFERENCE} \
--variant DIMMITIS_POPGEN.cohort.2024-06-27.NO_OUTGROUPS.recode.nuclearINDELs.vcf \
--filter-expression 'QUAL < 46 || DP < 787 || DP > 12600 || MQ < 22.83 || SOR > 6.736 || QD < 0.990 || FS > 7.324 || MQRankSum < -4.797 || ReadPosRankSum < -4.178 || ReadPosRankSum > 1.960' \
--filter-name "INDEL_filtered" \
--output DIMMITIS_POPGEN.cohort.2024-06-27.NO_OUTGROUPS.nuclearINDELs.filtered.vcf"

#Mitochondrial
bsub.py 1 filter_mitoSNPs "gatk VariantFiltration \
--reference ${REFERENCE} \
--variant DIMMITIS_POPGEN.cohort.2024-06-27.NO_OUTGROUPS.recode.mitoSNPs.vcf \
--filter-expression 'QUAL < 77 || DP < 199695 || DP > 659990 || MQ < 20.92 || SOR > 12.402 || QD < 3.320 || FS > 41.869 || MQRankSum < -5.622 || ReadPosRankSum < -3.003 || ReadPosRankSum > 5.452' \
--filter-name "SNP_filtered" \
--output DIMMITIS_POPGEN.cohort.2024-06-27.NO_OUTGROUPS.mitoSNPs.filtered.vcf"

bsub.py 1 filter_mitoINDELS "gatk VariantFiltration \
--reference ${REFERENCE} \
--variant DIMMITIS_POPGEN.cohort.2024-06-27.NO_OUTGROUPS.recode.mitoINDELs.vcf \
--filter-expression 'QUAL < 61 || DP < 200280 || DP > 655589 || MQ < 21.36 || SOR > 14.646 || QD < 0.288 || FS > 82.622 || MQRankSum < -7.928 || ReadPosRankSum < -7.416 || ReadPosRankSum > 6.170' \
--filter-name "INDEL_filtered" \
--output DIMMITIS_POPGEN.cohort.2024-06-27.NO_OUTGROUPS.mitoINDELs.filtered.vcf"

#Wolbachia
bsub.py 1 filter_WbSNPs "gatk VariantFiltration \
--reference ${REFERENCE} \
--variant DIMMITIS_POPGEN.cohort.2024-06-27.NO_OUTGROUPS.recode.WbSNPs.vcf \
--filter-expression 'QUAL < 51 || DP < 93526 || DP > 140940 || MQ < 22.00 || SOR > 11.661 || QD < 1.790 || FS > 46.426 || MQRankSum < -15.506 || ReadPosRankSum < -6.165 || ReadPosRankSum > 13.108' \
--filter-name "SNP_filtered" \
--output DIMMITIS_POPGEN.cohort.2024-06-27.NO_OUTGROUPS.WbSNPs.filtered.vcf"

bsub.py 1 filter_WbINDELS "gatk VariantFiltration \
--reference ${REFERENCE} \
--variant DIMMITIS_POPGEN.cohort.2024-06-27.NO_OUTGROUPS.recode.WbINDELs.vcf \
--filter-expression 'QUAL < 42 || DP < 93456 || DP > 145612 || MQ < 22.017 || SOR > 11.268 || QD < 1.210 || FS > 59.220 || MQRankSum < -6.716 || ReadPosRankSum < -6.413 || ReadPosRankSum > 6.014' \
--filter-name "INDEL_filtered" \
--output DIMMITIS_POPGEN.cohort.2024-06-27.NO_OUTGROUPS.WbINDELs.filtered.vcf"

# once done, count the filtered sites - funny use of "|" allows direct markdown table format
echo -e "| Filtered_VCF | Variants_PASS | Variants_FILTERED |\n| -- | -- | -- | " > filter.stats

for i in *filtered.vcf; do
     name=${i}; pass=$( grep -E 'PASS' ${i} | wc -l ); filter=$( grep -E 'filter' ${i} | wc -l );
     echo -e "| ${name} | ${pass} | ${filter} |" >> filter.stats
done

mv *.o *.e LOGS
```
This is the summary of the filtered variants ('filter.stats'):

| Filtered_VCF | Variants_PASS | Variants_FILTERED |
| -- | -- | -- | 
| DIMMITIS_POPGEN.cohort.2024-06-27.NO_OUTGROUPS.mitoINDELs.filtered.vcf | 240 | 30 |
| DIMMITIS_POPGEN.cohort.2024-06-27.NO_OUTGROUPS.mitoSNPs.filtered.vcf | 1511 | 141 |
| DIMMITIS_POPGEN.cohort.2024-06-27.NO_OUTGROUPS.nuclearINDELs.filtered.vcf | 748192 | 85599 |
| DIMMITIS_POPGEN.cohort.2024-06-27.NO_OUTGROUPS.nuclearSNPs.filtered.vcf | 1775848 | 189001 |
| DIMMITIS_POPGEN.cohort.2024-06-27.NO_OUTGROUPS.WbINDELs.filtered.vcf | 8897 | 731 |
| DIMMITIS_POPGEN.cohort.2024-06-27.NO_OUTGROUPS.WbSNPs.filtered.vcf | 69952 | 5110 |


Filtering looks ok.

Usually we would merge the SNP and INDEL files together to make a joined VCF. However, I want to disregard the indels moving forward and only focus on the SNPs. Some downstream tools don't like having indels in there.


### Filter genotypes based on x3 depth per genotype

```bash
#Nuclear 
bsub.py 1 filter_nuclear_GT \
"gatk VariantFiltration \
--reference ${REFERENCE} \
--variant ${VCF%.recode.vcf.gz}.nuclearSNPs.filtered.vcf \
--genotype-filter-expression ' DP < 3 '  \
--genotype-filter-name "DP_lt3" \
--output ${VCF%.vcf.gz}.nuclearSNPs.DPfiltered.vcf"

bsub.py --done "filter_nuclear_GT" 1 filter_nuclear_GT2 \
"gatk SelectVariants \
--reference ${REFERENCE} \
--variant ${VCF%.vcf.gz}.nuclearSNPs.DPfiltered.vcf \
--set-filtered-gt-to-nocall \
--output ${VCF%.vcf.gz}.nuclearSNPs.DPfilterNoCall.vcf"

#Mito
bsub.py 1 filter_mito_GT \
"gatk VariantFiltration \
--reference ${REFERENCE} \
--variant ${VCF%.recode.vcf.gz}.mitoSNPs.filtered.vcf \
--genotype-filter-expression ' DP < 3 '  \
--genotype-filter-name "DP_lt3" \
--output ${VCF%.vcf.gz}.mitoSNPs.DPfiltered.vcf"

bsub.py --done "filter_mito_GT" 1 filter_mito_GT2 \
"gatk SelectVariants \
--reference ${REFERENCE} \
--variant ${VCF%.vcf.gz}.mitoSNPs.DPfiltered.vcf \
--set-filtered-gt-to-nocall \
--output ${VCF%.vcf.gz}.mitoSNPs.DPfilterNoCall.vcf"

#wolbachia
bsub.py 1 filter_Wb_GT \
"gatk VariantFiltration \
--reference ${REFERENCE} \
--variant ${VCF%.recode.vcf.gz}.WbSNPs.filtered.vcf \
--genotype-filter-expression ' DP < 3 '  \
--genotype-filter-name "DP_lt3" \
--output ${VCF%.vcf.gz}.WbSNPs.DPfiltered.vcf"

bsub.py --done "filter_Wb_GT" 1 filter_Wb_GT2 \
"gatk SelectVariants \
--reference ${REFERENCE} \
--variant ${VCF%.vcf.gz}.WbSNPs.DPfiltered.vcf \
--set-filtered-gt-to-nocall \
--output ${VCF%.vcf.gz}.WbSNPs.DPfilterNoCall.vcf"
```


### Now we apply a set of standard filters for population genomics

bsub.py 10 run_standard_filt "run_standard_filt.sh"

```bash
# load modules
module load vcftools/0.1.16-c4

cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/04_VARIANTS/FILTER1/NO_OUTGROUPS
VCF=/lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/04_VARIANTS/FILTER1/NO_OUTGROUPS/DIMMITIS_POPGEN.cohort.2024-06-27.NO_OUTGROUPS.recode.vcf.gz

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
#After filtering, kept 138 out of 138 Individuals
#Outputting VCF file...
#After filtering, kept 267840 out of a possible 1964841 Sites

#--- nuclear SNPs
vcftools --vcf DIMMITIS_POPGEN.cohort.2024-06-27.NO_OUTGROUPS.recode.nuclear_SNPs.final.recode.vcf --remove-indels
#After filtering, kept 138 out of 138 Individuals
#After filtering, kept 267840 out of a possible 267840 Sites

#--- nuclear  INDELs
vcftools --vcf DIMMITIS_POPGEN.cohort.2024-06-27.NO_OUTGROUPS.recode.nuclear_SNPs.final.recode.vcf --keep-only-indels
#After filtering, kept 138 out of 138 Individuals
#After filtering, kept 0 out of a possible 267840 Sites


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
#After filtering, kept 138 out of 138 Individuals
#Outputting VCF file...
#After filtering, kept 32 out of a possible 1644 Sites

#--- mito SNPs
vcftools --vcf DIMMITIS_POPGEN.cohort.2024-06-27.NO_OUTGROUPS.recode.mito_SNPs.final.recode.vcf --remove-indels
#After filtering, kept 138 out of 138 Individuals
#After filtering, kept 32 out of a possible 32 Sites

#--- mito INDELs
vcftools --vcf DIMMITIS_POPGEN.cohort.2024-06-27.NO_OUTGROUPS.recode.mito_SNPs.final.recode.vcf --keep-only-indels
#After filtering, kept 138 out of 138 Individuals
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
#After filtering, kept 138 out of 138 Individuals
#Outputting VCF file...
#After filtering, kept 280 out of a possible 75054 Sites

#--- Wb SNPs
vcftools --vcf DIMMITIS_POPGEN.cohort.2024-06-27.NO_OUTGROUPS.recode.Wb_SNPs.final.recode.vcf --remove-indels
#After filtering, kept 138 out of 138 Individuals
#After filtering, kept 280 out of a possible 280 Sites

#--- Wb INDELs
vcftools --vcf DIMMITIS_POPGEN.cohort.2024-06-27.NO_OUTGROUPS.recode.Wb_SNPs.final.recode.vcf --keep-only-indels
#After filtering, kept 138 out of 138 Individuals
#After filtering, kept 0 out of a possible 280 Sites
```

### Now, we are filtering by missingness

```bash
#determine missingness per individual
vcftools --vcf ${VCF%.vcf.gz}.nuclear_SNPs.final.recode.vcf --out nuclear --missing-indv
vcftools --vcf ${VCF%.vcf.gz}.mito_SNPs.final.recode.vcf --out mito --missing-indv
vcftools --vcf ${VCF%.vcf.gz}.Wb_SNPs.final.recode.vcf --out wb --missing-indv
```

### Check the missingess in R

```R
 # Check missingness
  data_nuclear <- read.delim("nuclear.imiss", header=T)
  data_mito <- read.delim("mito.imiss", header=T)
  data_wb <- read.delim("wb.imiss", header=T)
  
  #creating the function - per sample
  fun_plot_missingness <- function(data,title) {
    
    plot <- ggplot(data, aes(INDV, 1-F_MISS)) +
      geom_boxplot(color = 'brown') +
      geom_point(size = 1, color = 'brown4') +
      theme_bw() +
      labs(x="Sample ID", y="Proportion of total variants present (1-missingness)")+
      ggtitle(title) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    print(plot)
    ggsave(paste0("plot_missingness_figure_",title,".png"))
  }
  
  # plotting for each dataset
  fun_plot_missingness(data_nuclear, "nuclear_variants")

  fun_plot_missingness(data_mito,"mitochondrial_variants")

  fun_plot_missingness(data_wb, "wb_variants")
```

```bash
# For nuclear (n=128) - nuclear_samplelist.keep
## remove samples with missingness less than 0.5 (AUS_BNE_AD_005, AUS_SYD_AD_016, AUS_TVS_AD_009, MYS_SEL_AD_001_R, ROU_GIU_AD_001_R, USA_GEO_MF_001, USA_GEO_MF_002, USA_ILL_MF_001, USA_LOU_MF_001, USA_MIP_MF_001)
AUS_BNE_AD_001
AUS_BNE_AD_002
AUS_BNE_AD_003
AUS_BNE_AD_003_R
AUS_BNE_AD_004
AUS_BNE_AD_006
AUS_BNE_AD_008
AUS_BNE_AD_009
AUS_CNS_AD_001
AUS_CNS_AD_002
AUS_LHR_AD_001
AUS_ROK_AD_001
AUS_SYD_AD_001
AUS_SYD_AD_002
AUS_SYD_AD_003
AUS_SYD_AD_004
AUS_SYD_AD_005
AUS_SYD_AD_005_R
AUS_SYD_AD_006
AUS_SYD_AD_007
AUS_SYD_AD_008
AUS_SYD_AD_008_R
AUS_SYD_AD_009
AUS_SYD_AD_010
AUS_SYD_AD_011
AUS_SYD_AD_012
AUS_SYD_AD_013
AUS_SYD_AD_014
AUS_SYD_AD_015
AUS_SYD_AD_017
AUS_TVS_AD_001
AUS_TVS_AD_002
AUS_TVS_AD_003
AUS_TVS_AD_004
AUS_TVS_AD_005
AUS_TVS_AD_006
AUS_TVS_AD_007
AUS_TVS_AD_008
AUS_TVS_AD_010
AUS_TVS_AD_011
AUS_TVS_AD_012
AUS_TVS_AD_013
AUS_TVS_AD_014
AUS_TVS_AD_015
AUS_TVS_AD_016
AUS_TVS_AD_017
AUS_TVS_AD_018
AUS_TVS_AD_019
CRI_SJO_AD_001
GRC_XAN_AD_001
GRC_XAN_AD_002
GRC_XAN_AD_003
GRC_XAN_AD_004
GRC_XAN_AD_005
GRC_XAN_AD_006
GRC_XAN_AD_007
GRC_XAN_AD_008
GRC_XAN_AD_009
GRC_XAN_AD_010
GRC_XAN_AD_011
GRC_XAN_AD_012
ITA_NEA_AD_001
ITA_NEA_AD_002
ITA_NEA_AD_003
ITA_PAV_AD_001
MYS_SEL_AD_001
PAN_BOC_AD_001
PAN_BOC_AD_002
PAN_BOC_AD_003
PAN_BOC_AD_004
PAN_BOC_AD_005
PAN_PUE_AD_001
PAN_PUE_AD_002
PAN_PUE_AD_003
PAN_PUE_AD_004
PAN_PUE_AD_004_R
PAN_PUE_AD_005
PAN_PUE_AD_006
PAN_SLO_AD_001
PAN_SLO_AD_002
PAN_SLO_AD_003
ROU_BUC_AD_001
ROU_COM_AD_001
ROU_GIU_AD_001
THA_BKK_AD_001
THA_BKK_AD_002
THA_BKK_AD_003
THA_BKK_AD_004
THA_BKK_AD_005
THA_BKK_AD_006
USA_FLO_AD_001
USA_FLO_AD_002
USA_FLO_AD_003
USA_FLO_AD_004
USA_FLO_AD_005
USA_FLO_AD_006
USA_FLO_AD_007
USA_FLO_AD_008
USA_FLO_AD_009
USA_FLO_AD_010
USA_FLO_AD_011
USA_FLO_AD_012
USA_FLO_AD_013
USA_FLO_AD_014
USA_FLO_AD_015
USA_FLO_AD_016
USA_GEO_AD_001
USA_ILL_AD_001
USA_ILL_AD_002
USA_LOU_AD_001
USA_LOU_AD_002
USA_MIS_AD_001
USA_MIS_AD_002
USA_TEX_AD_001
USA_TEX_AD_002
USA_TEX_AD_003
USA_TEX_AD_004
USA_TEX_AD_005
USA_TEX_AD_006
USA_TEX_AD_007
USA_TEX_AD_008
USA_TEX_AD_009
USA_TEX_AD_010
USA_TEX_AD_011
USA_TEX_AD_012
USA_TEX_AD_013
USA_TEX_AD_014
USA_TEX_AD_015


# For mithochondiral (n=132) - mito_samplelist.keep
# removed: (ROU_GIU_AD_001_R, USA_GEO_MF_001, USA_GEO_MF_002, USA_ILL_MF_001, USA_LOU_MF_001, USA_MIP_MF_001)
AUS_BNE_AD_001
AUS_BNE_AD_002
AUS_BNE_AD_003
AUS_BNE_AD_003_R
AUS_BNE_AD_004
AUS_BNE_AD_005
AUS_BNE_AD_006
AUS_BNE_AD_008
AUS_BNE_AD_009
AUS_CNS_AD_001
AUS_CNS_AD_002
AUS_LHR_AD_001
AUS_ROK_AD_001
AUS_SYD_AD_001
AUS_SYD_AD_002
AUS_SYD_AD_003
AUS_SYD_AD_004
AUS_SYD_AD_005
AUS_SYD_AD_005_R
AUS_SYD_AD_006
AUS_SYD_AD_007
AUS_SYD_AD_008
AUS_SYD_AD_008_R
AUS_SYD_AD_009
AUS_SYD_AD_010
AUS_SYD_AD_011
AUS_SYD_AD_012
AUS_SYD_AD_013
AUS_SYD_AD_014
AUS_SYD_AD_015
AUS_SYD_AD_016
AUS_SYD_AD_017
AUS_TVS_AD_001
AUS_TVS_AD_002
AUS_TVS_AD_003
AUS_TVS_AD_004
AUS_TVS_AD_005
AUS_TVS_AD_006
AUS_TVS_AD_007
AUS_TVS_AD_008
AUS_TVS_AD_009
AUS_TVS_AD_010
AUS_TVS_AD_011
AUS_TVS_AD_012
AUS_TVS_AD_013
AUS_TVS_AD_014
AUS_TVS_AD_015
AUS_TVS_AD_016
AUS_TVS_AD_017
AUS_TVS_AD_018
AUS_TVS_AD_019
CRI_SJO_AD_001
GRC_XAN_AD_001
GRC_XAN_AD_002
GRC_XAN_AD_003
GRC_XAN_AD_004
GRC_XAN_AD_005
GRC_XAN_AD_006
GRC_XAN_AD_007
GRC_XAN_AD_008
GRC_XAN_AD_009
GRC_XAN_AD_010
GRC_XAN_AD_011
GRC_XAN_AD_012
ITA_NEA_AD_001
ITA_NEA_AD_002
ITA_NEA_AD_003
ITA_PAV_AD_001
MYS_SEL_AD_001
MYS_SEL_AD_001_R
PAN_BOC_AD_001
PAN_BOC_AD_002
PAN_BOC_AD_003
PAN_BOC_AD_004
PAN_BOC_AD_005
PAN_PUE_AD_001
PAN_PUE_AD_002
PAN_PUE_AD_003
PAN_PUE_AD_004
PAN_PUE_AD_004_R
PAN_PUE_AD_005
PAN_PUE_AD_006
PAN_SLO_AD_001
PAN_SLO_AD_002
PAN_SLO_AD_003
ROU_BUC_AD_001
ROU_COM_AD_001
ROU_GIU_AD_001
THA_BKK_AD_001
THA_BKK_AD_002
THA_BKK_AD_003
THA_BKK_AD_004
THA_BKK_AD_005
THA_BKK_AD_006
USA_FLO_AD_001
USA_FLO_AD_002
USA_FLO_AD_003
USA_FLO_AD_004
USA_FLO_AD_005
USA_FLO_AD_006
USA_FLO_AD_007
USA_FLO_AD_008
USA_FLO_AD_009
USA_FLO_AD_010
USA_FLO_AD_011
USA_FLO_AD_012
USA_FLO_AD_013
USA_FLO_AD_014
USA_FLO_AD_015
USA_FLO_AD_016
USA_GEO_AD_001
USA_ILL_AD_001
USA_ILL_AD_002
USA_LOU_AD_001
USA_LOU_AD_002
USA_MIS_AD_001
USA_MIS_AD_002
USA_TEX_AD_001
USA_TEX_AD_002
USA_TEX_AD_003
USA_TEX_AD_004
USA_TEX_AD_005
USA_TEX_AD_006
USA_TEX_AD_007
USA_TEX_AD_008
USA_TEX_AD_009
USA_TEX_AD_010
USA_TEX_AD_011
USA_TEX_AD_012
USA_TEX_AD_013
USA_TEX_AD_014
USA_TEX_AD_015


# For wb (n=125) - wb_samplelist.keep
# remove: (AUS_SYD_AD_002, GRC_XAN_AD_011, MYS_SEL_AD_001, MYS_SEL_AD_001_R, PAN_BOC_AD_001, PAN_BOC_AD_003, ROU_COM_AD_001, ROU_GIU_AD_001_R, USA_GEO_MF_001, USA_GEO_MF_002, USA_ILL_MF_001, USA_LOU_MF_001, USA_MIP_MF_001)
AUS_BNE_AD_001
AUS_BNE_AD_002
AUS_BNE_AD_003
AUS_BNE_AD_003_R
AUS_BNE_AD_004
AUS_BNE_AD_005
AUS_BNE_AD_006
AUS_BNE_AD_008
AUS_BNE_AD_009
AUS_CNS_AD_001
AUS_CNS_AD_002
AUS_LHR_AD_001
AUS_ROK_AD_001
AUS_SYD_AD_001
AUS_SYD_AD_003
AUS_SYD_AD_004
AUS_SYD_AD_005
AUS_SYD_AD_005_R
AUS_SYD_AD_006
AUS_SYD_AD_007
AUS_SYD_AD_008
AUS_SYD_AD_008_R
AUS_SYD_AD_009
AUS_SYD_AD_010
AUS_SYD_AD_011
AUS_SYD_AD_012
AUS_SYD_AD_013
AUS_SYD_AD_014
AUS_SYD_AD_015
AUS_SYD_AD_016
AUS_SYD_AD_017
AUS_TVS_AD_001
AUS_TVS_AD_002
AUS_TVS_AD_003
AUS_TVS_AD_004
AUS_TVS_AD_005
AUS_TVS_AD_006
AUS_TVS_AD_007
AUS_TVS_AD_008
AUS_TVS_AD_009
AUS_TVS_AD_010
AUS_TVS_AD_011
AUS_TVS_AD_012
AUS_TVS_AD_013
AUS_TVS_AD_014
AUS_TVS_AD_015
AUS_TVS_AD_016
AUS_TVS_AD_017
AUS_TVS_AD_018
AUS_TVS_AD_019
CRI_SJO_AD_001
GRC_XAN_AD_001
GRC_XAN_AD_002
GRC_XAN_AD_003
GRC_XAN_AD_004
GRC_XAN_AD_005
GRC_XAN_AD_006
GRC_XAN_AD_007
GRC_XAN_AD_008
GRC_XAN_AD_009
GRC_XAN_AD_010
GRC_XAN_AD_012
ITA_NEA_AD_001
ITA_NEA_AD_002
ITA_NEA_AD_003
ITA_PAV_AD_001
PAN_BOC_AD_002
PAN_BOC_AD_004
PAN_BOC_AD_005
PAN_PUE_AD_001
PAN_PUE_AD_002
PAN_PUE_AD_003
PAN_PUE_AD_004
PAN_PUE_AD_004_R
PAN_PUE_AD_005
PAN_PUE_AD_006
PAN_SLO_AD_001
PAN_SLO_AD_002
PAN_SLO_AD_003
ROU_BUC_AD_001
ROU_GIU_AD_001
THA_BKK_AD_001
THA_BKK_AD_002
THA_BKK_AD_003
THA_BKK_AD_004
THA_BKK_AD_005
THA_BKK_AD_006
USA_FLO_AD_001
USA_FLO_AD_002
USA_FLO_AD_003
USA_FLO_AD_004
USA_FLO_AD_005
USA_FLO_AD_006
USA_FLO_AD_007
USA_FLO_AD_008
USA_FLO_AD_009
USA_FLO_AD_010
USA_FLO_AD_011
USA_FLO_AD_012
USA_FLO_AD_013
USA_FLO_AD_014
USA_FLO_AD_015
USA_FLO_AD_016
USA_GEO_AD_001
USA_ILL_AD_001
USA_ILL_AD_002
USA_LOU_AD_001
USA_LOU_AD_002
USA_MIS_AD_001
USA_MIS_AD_002
USA_TEX_AD_001
USA_TEX_AD_002
USA_TEX_AD_003
USA_TEX_AD_004
USA_TEX_AD_005
USA_TEX_AD_006
USA_TEX_AD_007
USA_TEX_AD_008
USA_TEX_AD_009
USA_TEX_AD_010
USA_TEX_AD_011
USA_TEX_AD_012
USA_TEX_AD_013
USA_TEX_AD_014
USA_TEX_AD_015
```


### Let's check different thresholds for each dataset

bsub.py 2 run_check_thresh "run_check_thresh.sh"

```bash
# load modules
module load vcftools/0.1.16-c4

cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/04_VARIANTS/FILTER1/NO_OUTGROUPS
VCF=/lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/04_VARIANTS/FILTER1/NO_OUTGROUPS/DIMMITIS_POPGEN.cohort.2024-06-27.NO_OUTGROUPS.recode.vcf.gz

# For nuclear variants
for i in 0.7 0.8 0.9 1; do
     vcftools --vcf ${VCF%.vcf.gz}.nuclear_SNPs.final.recode.vcf --keep nuclear_samplelist.keep --max-missing ${i} ;
done

# max-missing = 0.7
#After filtering, kept 128 out of 138 Individuals
#After filtering, kept 267246 out of a possible 267840 Sites

# max-missing = 0.8
#After filtering, kept 128 out of 138 Individuals
#After filtering, kept 266451 out of a possible 267840 Sites

# max-missing = 0.9
#After filtering, kept 128 out of 138 Individuals
#After filtering, kept 264558 out of a possible 267840 Sites

# max-missing = 1
#After filtering, kept 128 out of 138 Individuals
#After filtering, kept 194553 out of a possible 267840 Sites


# For mito variants
for i in 0.7 0.8 0.9 1; do
     vcftools --vcf ${VCF%.vcf.gz}.mito_SNPs.final.recode.vcf --keep mito_samplelist.keep --max-missing ${i} ;
done

# max-missing = 0.7
#After filtering, kept 132 out of 138 Individuals
#After filtering, kept 32 out of a possible 32 Sites

# max-missing = 0.8
#After filtering, kept 132 out of 138 Individuals
#After filtering, kept 32 out of a possible 32 Sites

# max-missing = 0.9
#After filtering, kept 132 out of 138 Individuals
#After filtering, kept 32 out of a possible 32 Sites

# max-missing = 1
#After filtering, kept 132 out of 138 Individuals
#After filtering, kept 30 out of a possible 32 Sites


# For Wb variants
for i in 0.7 0.8 0.9 1; do
     vcftools --vcf ${VCF%.vcf.gz}.Wb_SNPs.final.recode.vcf --keep wb_samplelist.keep --max-missing ${i} ;
done

# max-missing = 0.7
#After filtering, kept 125 out of 138 Individuals
#After filtering, kept 280 out of a possible 280 Sites

# max-missing = 0.8
#After filtering, kept 125 out of 138 Individuals
#After filtering, kept 280 out of a possible 280 Sites

# max-missing = 0.9
#After filtering, kept 125 out of 138 Individuals
#After filtering, kept 280 out of a possible 280 Sites

# max-missing = 1
#After filtering, kept 125 out of 138 Individuals
#After filtering, kept 58 out of a possible 280 Sites
```

Selecting a max missingness of 0.9 for nuclear, 0.9 for mito and 0.9 for Wb is sensible.

bsub.py 1 run_thresh "run_thresh.sh"

```bash
# load modules
module load vcftools/0.1.16-c4

cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/04_VARIANTS/FILTER1/NO_OUTGROUPS
VCF=/lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/04_VARIANTS/FILTER1/NO_OUTGROUPS/DIMMITIS_POPGEN.cohort.2024-06-27.NO_OUTGROUPS.recode.vcf.gz
mkdir FINAL_SETS

# For nuclear
vcftools --vcf ${VCF%.vcf.gz}.nuclear_SNPs.final.recode.vcf \
     --keep nuclear_samplelist.keep \
     --max-missing 0.9 \
     --recode --recode-INFO-all \
     --out FINAL_SETS/nuclear_samples3x_missing0.9

# For mito
vcftools --vcf ${VCF%.vcf.gz}.mito_SNPs.final.recode.vcf \
     --keep mito_samplelist.keep \
     --max-missing 0.9 \
     --recode --recode-INFO-all \
     --out FINAL_SETS/mito_samples3x_missing0.9

# For wb
vcftools --vcf ${VCF%.vcf.gz}.Wb_SNPs.final.recode.vcf \
     --keep wb_samplelist.keep \
     --max-missing 0.9 \
     --recode --recode-INFO-all \
     --out FINAL_SETS/wb_samples3x_missing0.9
```


### Select only the variants in the chr 1 to chr4, avoiding the chrX and the scaffolds. Select chr1-4 and chrX separately.

bsub.py 4 run_select_chr "run_select_chr.sh"

```bash
# load modules
module load vcftools/0.1.16-c4

cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/04_VARIANTS/FILTER1/NO_OUTGROUPS/FINAL_SETS

# chr1-4
vcftools --vcf nuclear_samples3x_missing0.9.recode.vcf \
--chr dirofilaria_immitis_chr1 \
--chr dirofilaria_immitis_chr2 \
--chr dirofilaria_immitis_chr3 \
--chr dirofilaria_immitis_chr4 \
--recode --out nuclear_samples3x_missing0.9.chr1to4
#After filtering, kept 128 out of 128 Individuals
#After filtering, kept 196550 out of a possible 264558 Sites


vcftools --vcf nuclear_samples3x_missing0.9.chr1to4.recode.vcf --remove-indels
#After filtering, kept 128 out of 128 Individuals
#After filtering, kept 196550 out of a possible 196550 Sites

vcftools --vcf nuclear_samples3x_missing0.9.chr1to4.recode.vcf --keep-only-indels
#After filtering, kept 128 out of 128 Individuals
#After filtering, kept 0 out of a possible 196550 Sites - this makes sense because I removed the indels earlier and only focused on the SNPs.


# chr1
vcftools --vcf nuclear_samples3x_missing0.9.recode.vcf \
--chr dirofilaria_immitis_chr1 \
--recode --out nuclear_samples3x_missing0.9.chr1
#After filtering, kept 128 out of 128 Individuals
#After filtering, kept 50075 out of a possible 264558 Sites

# chr2
vcftools --vcf nuclear_samples3x_missing0.9.recode.vcf \
--chr dirofilaria_immitis_chr2 \
--recode --out nuclear_samples3x_missing0.9.chr2
#After filtering, kept 128 out of 128 Individuals
#After filtering, kept 44937 out of a possible 264558 Sites

# chr3
vcftools --vcf nuclear_samples3x_missing0.9.recode.vcf \
--chr dirofilaria_immitis_chr3 \
--recode --out nuclear_samples3x_missing0.9.chr3
#After filtering, kept 128 out of 128 Individuals
#After filtering, kept 50855 out of a possible 264558 Sites

# chr4
vcftools --vcf nuclear_samples3x_missing0.9.recode.vcf \
--chr dirofilaria_immitis_chr4 \
--recode --out nuclear_samples3x_missing0.9.chr4
#After filtering, kept 128 out of 128 Individuals
#After filtering, kept 50683 out of a possible 264558 Sites

# chrX
vcftools --vcf nuclear_samples3x_missing0.9.recode.vcf \
--chr dirofilaria_immitis_chrX \
--recode --out nuclear_samples3x_missing0.9.chrX
#After filtering, kept 128 out of 128 Individuals
#After filtering, kept 67561 out of a possible 264558 Sites
```

```bash
# Australian samples only
grep '^AUS' ../nuclear_samplelist.keep > ../AUS_samplelist.keep

bsub.py 4 select_aus "vcftools --vcf nuclear_samples3x_missing0.9.chr1to4.recode.vcf \
--keep ../AUS_samplelist.keep \
--recode --out nuclear_samples3x_missing0.9.chr1to4.AUS"
```

**PCA**

![](images/plot_nuc_5.png)

![](images/PC3_PC4_plot.png)

![](images/plot_mt_5.png)

![](images/plot_wb_5.png)



