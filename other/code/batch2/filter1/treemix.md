# Treemix

## Treemix

Prep:

```bash
# Install
conda create --name py-popgen
conda activate py-popgen
conda install -c bioconda treemix

# Make .bed file for autosome only
awk '$1 ~ /^dirofilaria_immitis_chr[1-4]$/ {print}' /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/01_REF/dimmitis_WSI_2.2.bed > dimmitis_WSI_2.2_autosome.bed
```



bsub.py 1 treemix "../treemix.sh"

```bash
conda activate py-popgen

# load modules
module load bcftools/1.14--h88f3f91_0
#vcftools
module load vcftools/0.1.16-c4
# also need htslib for tabix
module load common-apps/htslib/1.9.229
#bsub.py
module load bsub.py
# plink
module load plink/1.90b6.18--h516909a_0

export LD_LIBRARY_PATH=/nfs/users/nfs_r/rp24/miniconda3/envs/py-popgen/lib

cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/02_VARIANTS/TREEMIX

vcftools \
--gzvcf ../FINAL_SETS/nuclear_samples3x_missing0.9.chr1to4.recode.RENAMED.vcf.gz \
--chr dirofilaria_immitis_chr1 \
--chr dirofilaria_immitis_chr2 \
--chr dirofilaria_immitis_chr3 \
--chr dirofilaria_immitis_chr4 \
--max-missing 1 \
--keep nuclear_samples3x_RENAMED.list \
--recode --stdout | gzip > treemix.vcf.gz
# using the bed file kept cutting out the chr1 variant for some reason?
# After filtering, kept 60 out of 61 Individuals
# After filtering, kept 5412 out of a possible 153555 Sites


./ldPruning.sh treemix.vcf.gz
# took ldPruning.sh from Steve's folder
# finished, new file treemix.LDpruned filtered for LD in 10 kb windows, shifting by 10 kb with LD threshold 0.1

# before pruning
vcftools --gzvcf treemix.vcf.gz
# After filtering, kept 60 out of 60 Individuals
# After filtering, kept 5412 out of a possible 5412 Sites

# after pruning
vcftools --gzvcf treemix.LDpruned.vcf.gz
# After filtering, kept 60 out of 60 Individuals
# After filtering, kept 2262 out of a possible 2262 Sites

bcftools query -l treemix.vcf.gz > data.clust
# Edited my own data.clust file. Exported the above output into Excel, added in my own metadata columns manually, then transferred back into the farm.

# this "data.clust" needs to be edited to put the population in, eg.

"
AN_DNK_COG_EN_0012	AN_DNK_COG_EN_0012	ANCIENT
AN_NLD_KAM_EN_0034	AN_NLD_KAM_EN_0034	ANCIENT
MN_CHN_GUA_HS_001	MN_CHN_GUA_HS_001	CHN
MN_CHN_GUA_HS_002	MN_CHN_GUA_HS_002	CHN
MN_CHN_GUA_HS_003	MN_CHN_GUA_HS_003	CHN
"


cat data.clust | cut -f3 | sort | uniq > populations.list

./vcf2treemix.sh treemix.LDpruned.vcf.gz data.clust
# took vcf2treemix.sh from Steve's folder

# run treemix, across a range of migration edges, and across a range of seeds
#--- the range of seeds is used by optM below to estimate the optimal number of edges.

for i in {0..5}; do
     for j in {0..10}; do
          treemix -i treemix.LDpruned.treemix.frq.gz -seed $j -m $i -o treemix.m_$i.s_$j -bootstrap -k 500  > treemix_${i}_log &
     done;
done

# Steve used colobus monkey as root - this helps determine the direction of the branches. I don't really have a root... run without for now (but think about including later - maybe D. repens or Onchocerca).
# -m is migration events
```


## Plotting treemix data

```R
# Plotting treemix data

# load libraries
library(RColorBrewer)
library(R.utils)

setwd("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/extra_data/filter1/treemix")

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
```


## Estimating the optimal number of migration edges

```R
# run in the folder with the treemix output files
#install.packages("OptM")

library(OptM)

optM(folder="./")

#The maximum value for delta m was 0.9236 at m = 3 edges.

# remake the plots, using 3 migration edges
prefix="treemix"

par(mfrow=c(1,1))

plot_tree(cex=0.8,paste0(prefix,".m_3.s_4"))
title(paste(3,"edges"))

plot_resid(stem=paste0(prefix,".",3),pop_order="populations.list")
title(paste(3,"edges"))
# this last part didn't work...?

```

## Calculate the variance explained by the data

```bash

Rscript treemixVarianceExplained.R treemix.m_1.s_2

#Standard error for all entries in the covariance matrix estimated from the data 0.00192256580473373
# Variance of relatedness between populations explained by the model 0.879974704806534
```