# Treemix

## Treemix

Prep:

```bash
# Install
conda create --name py-popgen
conda activate py-popgen
conda install -c bioconda treemix # version 1.12
conda activate py-popgen

# load modules
module load bcftools/1.14--h88f3f91_0
module load vcftools/0.1.16-c4
module load common-apps/htslib/1.9.229
module load bsub.py
module load plink/1.90b6.18--h516909a_0

cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/05_ANALYSIS/TREEMIX

# Make .bed file for autosome only
awk '$1 ~ /^dirofilaria_immitis_chr[1-4]$/ {print}' /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/01_REF/dimmitis_WSI_2.2.bed > dimmitis_WSI_2.2_autosome.bed

export LD_LIBRARY_PATH=/nfs/users/nfs_r/rp24/miniconda3/envs/py-popgen/lib

ln -s /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/04_VARIANTS/FILTER1/OUTGROUPS/DIMMITIS_POPGEN.cohort.2024-06-27.nuclear_SNPs.final.recode.vcf

vcftools \
--vcf DIMMITIS_POPGEN.cohort.2024-06-27.nuclear_SNPs.final.recode.vcf \
--chr dirofilaria_immitis_chr1 \
--chr dirofilaria_immitis_chr2 \
--chr dirofilaria_immitis_chr3 \
--chr dirofilaria_immitis_chr4 \
--max-missing 1 \
--keep treemix.keep \
--recode --stdout | gzip > treemix.vcf.gz
# After filtering, kept 129 out of 143 Individuals
# After filtering, kept 19009 out of a possible 268135 Sites



./ldPruning.sh treemix.vcf.gz
# working...
#finished, new file treemix.LDpruned filtered for LD in 10 kb windows, shifting by 10 kb with LD threshold 0.1

# before pruning
vcftools --gzvcf treemix.vcf.gz
# After filtering, kept 129 out of 129 Individuals
# After filtering, kept 19009 out of a possible 19009 Sites

# after pruning
vcftools --gzvcf treemix.LDpruned.vcf.gz
# After filtering, kept 129 out of 129 Individuals
# After filtering, kept 10242 out of a possible 10242 Sites

bcftools query -l treemix.vcf.gz | awk '{split($1,pop,"."); print $1"\t"$1"\t"pop[2]}' > data.clust
# Edited my own data.clust file. Exported the above output into Excel, added in my own metadata columns manually, then transferred back into the farm.

# this "data.clust" needs to be edited to put the population in, eg:

"
AUS_BNE_AD_001	AUS_BNE_AD_001	AUS_BNE
AUS_BNE_AD_002	AUS_BNE_AD_002	AUS_BNE
AUS_BNE_AD_003	AUS_BNE_AD_003	AUS_BNE
AUS_BNE_AD_003_R	AUS_BNE_AD_003_R	AUS_BNE
AUS_BNE_AD_004	AUS_BNE_AD_004	AUS_BNE
AUS_BNE_AD_006	AUS_BNE_AD_006	AUS_BNE
AUS_BNE_AD_008	AUS_BNE_AD_008	AUS_BNE
AUS_BNE_AD_009	AUS_BNE_AD_009	AUS_BNE
AUS_CNS_AD_001	AUS_CNS_AD_001	AUS_CNS
AUS_CNS_AD_002	AUS_CNS_AD_002	AUS_CNS
AUS_LHR_AD_001	AUS_LHR_AD_001	AUS_LHR

"


cat data.clust | cut -f3 | sort | uniq > populations.list

./vcf2treemix.sh treemix.LDpruned.vcf.gz data.clust

# run treemix, across a range of migration edges, and across a range of seeds
#--- the range of seeds is used by optM below to estimate the optimal number of edges.
bsub.py 10 run_treemix "./treemix.sh"
```

```bash
cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/05_ANALYSIS/TREEMIX

for i in {0..5}; do
     for j in {0..10}; do
          treemix -i treemix.LDpruned.treemix.frq.gz -seed $j -m $i -o treemix.m_$i.s_$j -bootstrap -k 500  > treemix_${i}_log &
     done;
done

# I have D. ursi as an outgroup so it can help root the tree
# -m is migration events
```

Successfully completed.


## Plotting treemix data

```R
# Plotting treemix data

# load libraries
library(RColorBrewer)
library(R.utils)
library(ggplot2)

setwd("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/batch4/FILTER1/OUTGROUPS/treemix")

source("plotting_funcs.R")
prefix="treemix"

# plot trees across range of migration edges
par(mfrow=c(1,2))

# generate pdf file containing multiple lots of Treemix trees
pdf("treemix_edges_2-6.pdf", width = 10, height = 10)
par(mar = c(5, 4, 4, 8))
for(edge in 0:5){
  plot_tree(cex=0.8,paste0(prefix,".m_",edge,".s_4"))
  title(paste(edge,"edges"))
}
dev.off()


# generate series of pdf files - each containing a treemix plot and a residual plot for a diff number of edges
for(edge in 0:5){
  pdf(paste0("plot_treemix_tree_m-",edge,".pdf"))
  plot_tree(cex=0.8,paste0(prefix,".m_",edge,".s_3"))
  plot_resid(stem=paste0(prefix,".m_",edge,".s_3"),pop_order="populations.list")
  title(paste(edge,"edges"))
  dev.off()
}


# plot residuals across range of migration edges - positive values suggest admixture
for(edge in 0:5){
  pdf(paste0("plot_treemix_residuals_m-",edge,".pdf"))
  plot_resid(stem=paste0(prefix,".m_",edge,".s_4"),pop_order="populations.list")
  title(paste(edge,"edges"))
        dev.off()
}
for(edge in 0:5){
  tiff(paste0("plot_treemix_residuals_m-",edge,".tiff"))
  plot_resid(stem=paste0(prefix,".m_",edge,".s_4"),pop_order="populations.list")
  title(paste(edge,"edges"))
  dev.off()
}

#####################################################################################################

# Estimating the optimal number of migration edges

# run in the folder with the treemix output files
#install.packages("OptM")

library(OptM)

optM(folder="./")

#The maximum value for delta m was 3.4477 at m = 1 edges.

# remake the plots, using 1 migration edge
prefix="treemix"

par(mfrow=c(1,1)) # single plot per page

plot_tree(cex=0.8,paste0(prefix,".m_0.s_4"))
title(paste(0,"edges"))
ggsave("plot_treemix_m-0.tif", dpi = 300)

plot_tree(cex=0.8,paste0(prefix,".m_1.s_4"))
title(paste(1,"edges"))

# save a few other edges
plot_tree(cex=0.8,paste0(prefix,".m_2.s_4"))
title(paste(2,"edges"))

plot_tree(cex=0.8,paste0(prefix,".m_3.s_4"))
title(paste(3,"edges"))

# plot_resid(stem=paste0(prefix,".",2),pop_order="populations.list")
# title(paste(1,"edges")) ## didn't work?


prefix="treemix"
tiff("plot_treemix_tree_m-1.tiff", width = 1200, height = 1200, res = 250)
par(mfrow=c(1,1)) # single plot per page
plot_tree(cex=0.7,paste0(prefix,".m_1.s_4"))
title(paste(1,"edges"))
dev.off()
```


# havent done the below

## Calculate the variance explained by the data

```bash

Rscript treemixVarianceExplained.R treemix.m_1.s_2

#Standard error for all entries in the covariance matrix estimated from the data 0.00192256580473373
# Variance of relatedness between populations explained by the model 0.879974704806534
```