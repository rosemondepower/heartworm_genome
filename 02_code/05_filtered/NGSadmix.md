# NGSadmix

## Run NGSADmix

Use filter1 vcf file, autosomes only, no outgroups

```bash
# load modules
module load PaM/environment
module load bsub.py/0.42.1
module load vcftools/0.1.16-c4
module load htslib-1.19/perl-5.38.0

cd /lustre/scratch127/pam/teams/team333/rp24/dirofilaria_immitis/R2/ngsadmix

bsub.py 10 run_ngsadmix "run_ngsadmix.sh"
```

```bash
cd /lustre/scratch127/pam/teams/team333/rp24/dirofilaria_immitis/R2/ngsadmix

mkdir CHROMOSOMES_PL

ln -s ../vcf/dirofilaria_global.cohort.2025-06-18.nuclearSNPs..final.autosomal.recode.vcf.keep_samples.dimmitis-only.recode.vcf

cat /lustre/scratch127/pam/teams/team333/sd21/dirofilaria_immitis/project_global-diversity/01_references/dimmitis_WSI_2.2.fa.fai | cut -f1 | grep -v "chrX\|scaffold\|Wb\|MtDNA" | while read -r CHR; do
     vcftools --vcf dirofilaria_global.cohort.2025-06-18.nuclearSNPs..final.autosomal.recode.vcf.keep_samples.dimmitis-only.recode.vcf --max-missing 1 --out CHROMOSOMES_PL/${CHR} --BEAGLE-PL --chr ${CHR};
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

# check if it ran ok
grep -i "Successfully completed" run_ngsadmix_multiK.o | wc -l
grep -i "error" run_ngsadmix_multiK.o
grep -i "error" run_ngsadmix_multiK.e
grep -i "exit" run_ngsadmix_multiK.o
grep -i "exit" run_ngsadmix_multiK.e

```
All 45 jobs successfully completed.

## Best convergence

```bash
# We want to assess convergence and find the run with the best log likelihood
# In the log file, the log likelihood of the estimates is called "best like"

for i in 2 3 4 5 6 7 8 9 10; do
     echo "K=${i}"
     for j in 1 2 3 4 5; do
          cat k_${i}_s_${j}_out.log | grep "best like" | awk -F"[ =]" '{print $3}' >> allK$i.likes
     done
     cat -n allK$i.likes | sort -rhk2
done
# -r reverses it so highest values first, -h makes human readable sorting, -k2 is sort by 2nd column
'
K=2
     2  -11517320.455859
     3  -11517320.522817
     4  -11639201.678445
     5  -11639201.680981
     1  -11639201.699417
K=3
     3  -10248651.254286
     4  -10248651.360505
     1  -10248651.363292
     2  -10248651.400245
     5  -10248651.443510
K=4
     1  -9641370.673188
     3  -9641370.712490
     5  -9722592.219498
     4  -9722592.296260
     2  -9769603.771627
K=5
     2  -9115978.132956
     5  -9115980.104525
     4  -9162652.447205
     1  -9243240.316138
     3  -9243240.340050
K=6
     2  -8636945.748991
     5  -8636945.799871
     3  -8636945.821184
     4  -8895687.690766
     1  -8942518.666818
K=7
     4  -8295019.266199
     3  -8420648.687927
     5  -8420648.752086
     1  -8420649.403551
     2  -8486599.219966
K=8
     1  -8069630.731791
     4  -8069630.820926
     3  -8078768.742214
     5  -8078769.721609
     2  -8078770.339767
K=9
     5  -7811514.573990
     2  -7918575.249979
     4  -7924909.211166
     3  -7962364.344344
     1  -8063219.490502
K=10
     3  -7596359.300553
     5  -7596359.364086
     1  -7596359.423381
     2  -7649532.717928
     4  -7775280.242831
'
# the best overall seed for all the Ks seems to be 2.
```




## Admixture plots

```R
# load libraries
library(ggsci)
library(patchwork)
library(tidyverse)
library(reshape2)
library(dplyr)
library(RColorBrewer) # v.1.1-3

setwd("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/R2_Sep25/ngsadmix")

red_palette1 <- brewer.pal(n = 9, name = "YlOrRd")
red_palette2 <- brewer.pal(n=9, name = "YlOrBr")
blue_palette <- brewer.pal(n = 9, name = "Blues")
green_palette1 <- brewer.pal(n = 9, name = "BuGn")
green_palette2 <- brewer.pal(n=9, name = "YlGn")

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
  data <- read.table(data, header=FALSE)
  #names(data) <- paste("ancestral", 1:ncol(data), sep="")
  
  # bring metadata and data together
  data <- cbind(samples, metadata, data)
  data <- melt(data, id.vars=c("sample_ID", "sampleID", "country","country_city", "city"))
  
  # make plot
  ggplot(data,aes(sample_ID,value,fill=variable)) +
    geom_col(color = "gray", size = 0.1)+
    facet_grid(~ city, switch = "x", scales = "free_x", space = "free_x") +
    theme_minimal() + labs(title = title, y = "Ancestry", x = NULL) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_discrete(expand = expansion(add = 0)) +
    theme(plot.title = element_text(),
          axis.title.x = element_text(),
          panel.spacing.x = unit(0.1, "lines"),
          axis.text.x = element_blank(),
          axis.text.y = element_text(),
          # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          strip.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          panel.grid = element_blank()) +
    scale_fill_manual(values = c("V1"=green_palette1[6],
                                 "V2"="tomato3",
                                 "V3"="darkmagenta",
                                 "V4"=blue_palette[6],
                                 "V5"="palevioletred1",
                                 "V6"="mediumpurple3",
                                 "V7"=blue_palette[4],
                                 "V8"="darkseagreen",
                                 "V9"="khaki1",
                                 "V10"= "lightsalmon"),
                      guide = "none")
}


s = 2
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
ggsave("admixture_plots_k2-10.tif", dpi = 600, height=15, width=10)
ggsave("admixture_plots_k2-10.pdf", dpi = 600, height=15, width=10)
ggsave("admixture_plots_k2-10.eps", dpi = 600, height=15, width=10)
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
Optimal K by Evanno is: 3 # not updated

But k=7 appears to illustrate the continental separation described in the PCA.

```R
# plot k = 3 (optimal k)
k_3_plot
ggsave("admixture_plots_k3.tif", dpi = 600, height=1.5, width=10)
ggsave("admixture_plots_k3.pdf", dpi = 600, height=1.5, width=10)
ggsave("admixture_plots_k3.eps", dpi = 600, height=1.5, width=10)

# plot k = 6
k_6_plot
ggsave("admixture_plots_k6.tif", dpi = 600, height=1.5, width=10)
ggsave("admixture_plots_k6.pdf", dpi = 600, height=1.5, width=10)
ggsave("admixture_plots_k6.eps", dpi = 600, height=1.5, width=10)

# plot k = 7
k_7_plot
ggsave("admixture_plots_k7.tif", dpi = 600, height=1.5, width=10)
ggsave("admixture_plots_k7.pdf", dpi = 600, height=1.5, width=10)
ggsave("admixture_plots_k7.eps", dpi = 600, height=1.5, width=10)

###############################################################################################


# k = 7 with sample labels


plot_admixture <- function(data,title) {
  
  # get metadata
  samples <- read.delim("samples.list", header=F)
  names(samples) <- c("sample_ID") #assign column name
  metadata <- read.delim("metadata.txt", header=F)
  names(metadata) <- c("sampleID", "country", "country_city", "city") #assign column name
  
  desired_city_order <- c("SYD", "BNE", "ROK", "TSV", "CNS", "LHR", "BKK", "SEL",  "GEO", "ILL", "LOU", "MIS", "TEX", "FLO", "NEA", "PAV", "BUC", "GIU", "COM", "XAN", "SJO", "PUE", "SLO",  "BOC") 
  metadata$city <- factor(metadata$city, levels = desired_city_order)
  
  
  
  
  
  # read data
  data <- read.table(data, header=FALSE)
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
          # axis.text.x = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 4),
          strip.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          panel.grid = element_blank()) +
    scale_fill_manual(values = c("V1"=green_palette1[6],
                                 "V2"="tomato3",
                                 "V3"="darkmagenta",
                                 "V4"=blue_palette[6],
                                 "V5"="palevioletred1",
                                 "V6"="mediumpurple3",
                                 "V7"=blue_palette[4],
                                 "V8"="darkseagreen",
                                 "V9"="khaki1",
                                 "V10"= "lightsalmon"),
                      guide = "none")
}

s = 2
# run function for each value of K
k_7_plot_ID <- plot_admixture(paste0("input/k_7_s_",s,"_out.qopt"), "K = 7")

k_7_plot_ID
ggsave("admixture_plots_k7_ID.tif", dpi = 600, height=3, width=10)
ggsave("admixture_plots_k7_ID.pdf", dpi = 600, height=3, width=10)
```