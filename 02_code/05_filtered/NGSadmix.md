# NGSadmix

## Run NGSADmix

Use filter1 vcf file, autosomes only, no outgroups

```bash
# load modules
module load bsub.py/0.42.1
module load vcftools/0.1.16-c4
module load common-apps/htslib/1.9.229

cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/05_ANALYSIS/NGSADMIX

bsub.py 10 run_ngsadmix "run_ngsadmix.sh"
```

```bash
cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/05_ANALYSIS/NGSADMIX

mkdir CHROMOSOMES_PL

ln -s /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/04_VARIANTS/FILTER1/NO_OUTGROUPS/FINAL_SETS/nuclear_samples3x_missing0.9.chr1to4.recode.vcf

cat /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/01_REF/dimmitis_WSI_2.2.fa.fai| cut -f1 | grep -v "chrX\|scaffold\|Wb\|MtDNA" | while read -r CHR; do
     vcftools --vcf nuclear_samples3x_missing0.9.chr1to4.recode.vcf --max-missing 1 --out CHROMOSOMES_PL/${CHR} --BEAGLE-PL --chr ${CHR};
done

# merge the data from individual chromosomes into a single dataset
cd CHROMOSOMES_PL

cat $(ls -1 *PL | head -n1 ) | head -n1 > merged.PL

for i in *BEAGLE.PL; do
     cat ${i} | grep -v "marker" >> merged.PL;
done

# chromosomes=$(cat ../../01_REF/trichuris_trichiura.fa.fai | cut -f1 | grep -v "MITO" | while read -r CHROMOSOME; do printf "$CHROMOSOME,"; done | sed 's/,$//g')
# vcftools --gzvcf ../../04_VARIANTS/GATK_HC_MERGED/nuclear_samples3x_missing0.8_animalPhonly.recode.vcf.gz  --out CHROMOSOMES_PL/all_chromosomes --BEAGLE-PL --chr ${chromosomes}

head -n1 dirofilaria_immitis_chr1.BEAGLE.PL > chromosome.PL

cat dirofilaria_immitis_chr1.BEAGLE.PL dirofilaria_immitis_chr2.BEAGLE.PL dirofilaria_immitis_chr3.BEAGLE.PL dirofilaria_immitis_chr4.BEAGLE.PL | grep -v "marker" | sort -t ":" -k1,1 -k2,2n >> chromosome.PL
```
Successfully completed.

```bash
# run admixture for multiple values of K
for j in 1 2 3 4 5; do
     for i in 2 3 4 5 6 7 8 9 10; do
          bsub.py --queue long --threads 10 3 run_ngsadmix_multiK "/nfs/users/nfs_r/rp24/software/NGSadmix -likes CHROMOSOMES_PL/chromosome.PL -K ${i} -P 10 -seed ${j} -minMaf 0.02 -misTol 0.9 -o k_${i}_s_${j}_out" ;
     done;
done

# -likes is the beagle format filename with genotype likelihoods
# -K is the number of clusters
# -P is the number of threads
# -seed is the seed for initial guess in EM algorithm
# -minMaf is the minimum minor allele frequency
# -misTol is the tolerance for considering a site as missing. Default = 0.05. To include high quality genotypes only, increase this value (e.g. to 0.9 like Steve has done).
```
All 45 jobs successfully completed.




## Admixture plots

```R
# load libraries
library(ggsci)
library(patchwork)
library(tidyverse)
library(reshape2)
library(dplyr)

setwd("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/batch4/FILTER1/NO_OUTGROUPS/ngsadmix")

# make a function

plot_admixture <- function(data,title) {
  
  # get metadata
  samples <- read.delim("samples.list", header=F)
  names(samples) <- c("sample_ID") #assign column name
  metadata <- read.delim("metadata.txt", header=F)
  names(metadata) <- c("sampleID", "country", "country_city", "city") #assign column name

  desired_city_order <- c("SYD", "BNE", "ROK", "TSV", "CNS", "LHR", "BKK", "SEL",  "GEO", "ILL", "LOU", "MIS", "TEX", "FLO", "NEA", "PAV", "BUC", "GIU", "COM", "XAN", "SJO", "PUE", "SLO",  "BOC") 
  metadata$city <- factor(metadata$city, levels = desired_city_order)
  
  
  
  
  
  # read data
  data <- read.delim(data, sep=" ", header=F)
  #names(data) <- paste("ancestral", 1:ncol(data), sep="")
  
  # bring metadata and data together
  data <- cbind(samples, metadata, data)
  data <- melt(data, id.vars=c("sample_ID", "sampleID", "country","country_city", "city"))
  
  # make plot
  ggplot(data,aes(sample_ID,value,fill=variable)) +
    geom_col(color = "gray", size = 0.1)+
    facet_grid(~ city, switch = "x", scales = "free", space = "free") +
    theme_minimal() + labs(title = title, y = "Ancestry", x = NULL) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_discrete(expand = expansion(add = 0.7)) +
    theme(panel.spacing.x = unit(0.1, "lines"),
          axis.text.x = element_blank(),
          # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          strip.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          panel.grid = element_blank()) +
    scale_fill_npg(guide = "none")
}

s = 3
# run function for each value of K
k_2_plot <- plot_admixture(paste0("input/k_2_s_",s,"_out.qopt"), "K = 2")
k_3_plot <- plot_admixture(paste0("input/k_3_s_",s,"_out.qopt"), "K = 3")
k_4_plot <- plot_admixture(paste0("input/k_4_s_",s,"_out.qopt"), "K = 4")
k_5_plot <- plot_admixture(paste0("input/k_5_s_",s,"_out.qopt"), "K = 5")
k_6_plot <- plot_admixture(paste0("input/k_6_s_",s,"_out.qopt"), "K = 6")
k_7_plot <- plot_admixture(paste0("input/k_7_s_",s,"_out.qopt"), "K = 7")
k_8_plot <- plot_admixture(paste0("input/k_8_s_",s,"_out.qopt"), "K = 8")
k_9_plot <- plot_admixture(paste0("input/k_9_s_",s,"_out.qopt"), "K = 9")
k_10_plot <- plot_admixture(paste0("input/k_10_s_",s,"_out.qopt"), "K = 10")

# bring the plots together
k_2_plot + k_3_plot + k_4_plot + k_5_plot + k_6_plot + k_7_plot + k_8_plot + k_9_plot + k_10_plot + plot_layout(ncol=1)

# save it
ggsave("admixture_plots_k2-10.tif", height=15, width=10)
ggsave("admixture_plots_k2-10.pdf", height=15, width=10)
```

- need to determine the optimal K, at least from what the data suggests.
- usually there is a cross validation approach for tools like STRUCTURE and ADMIXTURE, but there doesnt seem to be one for NGSadmix


## Clumpak to determine optimal K

Website: https://clumpak.tau.ac.il/bestK.html

In command line
```bash
(for log in `ls k*.log`; do
     grep -Po 'like=\K[^ ]+' $log;
done) > logfile
```

Back to R
```R
logs <- as.data.frame(read.table("logfile"))

logs$K <- c(rep("10", 5), rep("2", 5), rep("3", 5), rep("4", 5), rep("5", 5), rep("6", 5), rep("7", 5), rep("8", 5), rep("9", 5))

write.table(logs[, c(2, 1)], "logfile_formatted", row.names = F,
    col.names = F, quote = F)
```

Used Clumpak to determine the optimal K:
Optimal K by Evanno is: 3

But k=6 appears to illustrate the continental separation described in the PCA.

```R
# plot k = 3 (optimal k)
k_3_plot
ggsave("admixture_plots_k3.tif", height=1.5, width=10)

# plot k = 6
k_6_plot
ggsave("admixture_plots_k6.tif", height=1.5, width=10)
```