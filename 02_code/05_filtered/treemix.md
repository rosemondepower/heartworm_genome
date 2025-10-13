# Treemix

## Treemix

Prep:

```bash
# Install
module load ISG/conda
conda create --name py-popgen
conda activate py-popgen
conda install -c bioconda treemix # NOW version 1.13
conda activate py-popgen

# load modules
module load bcftools/1.14--h88f3f91_0
module load vcftools/0.1.16-c4
module load htslib-1.19/perl-5.38.0
module load bsub.py
module load plink/1.90b6.18--h516909a_0

cd /lustre/scratch127/pam/teams/team333/rp24/dirofilaria_immitis/R2/treemix
REF=/lustre/scratch127/pam/teams/team333/sd21/dirofilaria_immitis/project_global-diversity/01_references/dimmitis_WSI_2.2.fa

# make bed file for reference
# Get the scaffolds/positions.
head ${REF}.fai
# Column 1 is the chromosome/scaffold, column 2 is how long it is, then there's some other info.

# Get chromosome, then start and end positions
awk '{print $1, "1", $2}' OFS="\t" ${REF}.fai | head

# Save this info as a bed file
awk '{print $1, "1", $2}' OFS="\t" ${REF}.fai > dimmitis_WSI_2.2.bed
# Now we have a nice bed file that has info telling us where things are

# Make .bed file for autosome only
awk '$1 ~ /^dirofilaria_immitis_chr[1-4]$/ {print}' dimmitis_WSI_2.2.bed > dimmitis_WSI_2.2_autosome.bed

export LD_LIBRARY_PATH=/software/conda/users/rp24/py-popgen/lib

ln -s /lustre/scratch127/pam/teams/team333/sd21/dirofilaria_immitis/project_global-diversity/04_analysis/04_pca/dirofilaria_global.cohort.2025-06-18.nuclearSNPs..final.autosomal.recode.vcf .

# get DI samples, no reps - D. ursi as outgroup
vcftools \
--vcf dirofilaria_global.cohort.2025-06-18.nuclearSNPs..final.autosomal.recode.vcf \
--chr dirofilaria_immitis_chr1 \
--chr dirofilaria_immitis_chr2 \
--chr dirofilaria_immitis_chr3 \
--chr dirofilaria_immitis_chr4 \
--max-missing 1 \
--keep keep_samples.dimmitis-only.nuclear_snps.no-reps_ursi.list \
--recode --stdout | gzip > treemix.vcf.gz
#After filtering, kept 125 out of 138 Individuals
#After filtering, kept 20429 out of a possible 218393 Sites

chmod a+x ldPruning.sh
./ldPruning.sh treemix.vcf.gz
#working...
#finished, new file treemix.LDpruned filtered for LD in 10 kb windows, shifting by 10 kb with LD threshold 0.1

# before pruning
vcftools --gzvcf treemix.vcf.gz
#After filtering, kept 125 out of 125 Individuals
#After filtering, kept 20429 out of a possible 20429 Sites

# after pruning
vcftools --gzvcf treemix.LDpruned.vcf.gz
#After filtering, kept 125 out of 125 Individuals
#After filtering, kept 10764 out of a possible 10764 Sites

bcftools query -l treemix.vcf.gz | awk '{split($1,pop,"."); print $1"\t"$1"\t"pop[2]}' > data.clust
# Edited my own data.clust file. Exported the above output into Excel, added in my own metadata columns manually, then transferred back into the farm.

# this "data.clust" needs to be edited to put the population in, eg:

"
AUS_BNE_AD_001	AUS_BNE_AD_001	AUS_BNE
AUS_BNE_AD_002	AUS_BNE_AD_002	AUS_BNE
AUS_BNE_AD_003	AUS_BNE_AD_003	AUS_BNE
AUS_BNE_AD_004	AUS_BNE_AD_004	AUS_BNE
AUS_BNE_AD_006	AUS_BNE_AD_006	AUS_BNE
AUS_BNE_AD_008	AUS_BNE_AD_008	AUS_BNE
AUS_BNE_AD_009	AUS_BNE_AD_009	AUS_BNE
AUS_CNS_AD_001	AUS_CNS_AD_001	AUS_CNS
AUS_CNS_AD_002	AUS_CNS_AD_002	AUS_CNS
AUS_LHR_AD_001	AUS_LHR_AD_001	AUS_LHR

"


cat data.clust | cut -f3 | sort | uniq > populations.list

# the below needed python2 to run, so make a new environment to install and run this
conda deactivate
conda create --name python2
conda install python=2.7
which python # added this to the header of plink2treemix.py

./vcf2treemix.sh treemix.LDpruned.vcf.gz data.clust

# run treemix, across a range of migration edges, and across a range of seeds
#--- the range of seeds is used by optM below to estimate the optimal number of edges.
bsub.py 10 run_treemix "./treemix.sh"
```

```bash
cd /lustre/scratch127/pam/teams/team333/rp24/dirofilaria_immitis/R2/treemix

for i in {0..5}; do
     for j in {0..10}; do
          treemix -i treemix.LDpruned.treemix.frq.gz -seed $j -m $i -o treemix.m_$i.s_$j -bootstrap -k 500  > treemix_${i}_log &
     done;
done

# I have D. ursi as an outgroup so it can help root the tree
# -m is migration events
```

Successfully completed.


```bash
for m in {0..5}; do
    echo "Best seed for m=$m:"
    grep "Exiting ln(likelihood) with $m migration events" treemix.m_${m}.s_*.llik | \
        awk -F':' '{print $1, $NF}' | sort -k2 -n -r | head -n 1
done

'
Best seed for m=0:
treemix.m_0.s_10.llik  -177.452
Best seed for m=1:
treemix.m_1.s_6.llik  287.717
Best seed for m=2:
treemix.m_2.s_6.llik  366.105
Best seed for m=3:
treemix.m_3.s_6.llik  399.677
Best seed for m=4:
treemix.m_4.s_6.llik  462.024
Best seed for m=5:
treemix.m_5.s_6.llik  504.541
'
# Best seed in general seems to be 6
```


## Plotting treemix data

```R
# Plotting treemix data

# load libraries
library(RColorBrewer)
library(R.utils)
library(ggplot2)

setwd("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/R2_Sep25/treemix")

source("plotting_funcs.R")
prefix="treemix"

# plot trees across range of migration edges
par(mfrow=c(1,2))

# generate pdf file containing multiple lots of Treemix trees
pdf("treemix_edges_2-6.pdf", width = 10, height = 10)
par(mar = c(5, 4, 4, 8))
for(edge in 0:5){
  plot_tree(cex=0.8,paste0(prefix,".m_",edge,".s_6"))
  title(paste(edge,"edges"))
}
dev.off()


# generate series of pdf files - each containing a treemix plot and a residual plot for a diff number of edges
for(edge in 0:5){
  pdf(paste0("plot_treemix_tree_m-",edge,".pdf"))
  plot_tree(cex=0.7,paste0(prefix,".m_",edge,".s_6"))
  plot_resid(stem=paste0(prefix,".m_",edge,".s_6"),pop_order="populations.list")
  title(paste(edge,"edges"))
  dev.off()
}


# plot residuals across range of migration edges - positive values suggest admixture
for(edge in 0:5){
  pdf(paste0("plot_treemix_residuals_m-",edge,".pdf"))
  plot_resid(stem=paste0(prefix,".m_",edge,".s_6"),pop_order="populations.list")
  title(paste(edge,"edges"))
        dev.off()
}
for(edge in 0:5){
  tiff(paste0("plot_treemix_residuals_m-",edge,".tiff"))
  plot_resid(stem=paste0(prefix,".m_",edge,".s_6"),pop_order="populations.list")
  title(paste(edge,"edges"))
  dev.off()
}

#####################################################################################################

# Estimating the optimal number of migration edges

# run in the folder with the treemix output files
#install.packages("OptM")

library(OptM)

optM(folder="./")

#The maximum value for delta m was 4.1012 at m = 1 edges.

# remake the plots, using 1 migration edge

prefix="treemix"
tiff("plot_treemix_tree_m-1.tiff", width = 1200, height = 1200, res = 250)
par(mfrow=c(1,1)) # single plot per page
plot_tree(cex=0.7,paste0(prefix,".m_1.s_6"))
title(paste(1,"edges"))
dev.off()


for(edge in 0:5){
  pdf(paste0("plot_treemix_tree_m-",edge,".pdf"))
  par(mar = c(5, 4, 4, 8))
  plot_tree(cex=1,paste0(prefix,".m_",edge,".s_6"))
  title(paste(edge,"edges"))
  dev.off()
}
```