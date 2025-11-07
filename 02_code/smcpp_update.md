# Dirofilaria immitis - global diversity

## SMC++

### input data
```bash
cd /nfs/users/nfs_s/sd21/lustre_link/dirofilaria_immitis/project_global-diversity/04_analysis/05_smc++ 

ln -s ../../01_references/dimmitis_WSI_2.2.fai
ln -s ../03_variants-filtered/dirofilaria_global.cohort.2025-06-18.nuclearSNPs..final.autosomal.recode.vcf
ln -s /lustre/scratch127/pam/users/sd21/dirofilaria_immitis/project_global-diversity/03_metadata/keep_samples.dimmitis-only.nuclear_snps.no-reps.list
```



### generate VCFs for each population
```bash
# 
cat ../03_variants-filtered/dirofilaria_global.cohort.2025-06-18.nuclearSNPs..final.autosomal.recode.vcf | bgzip > dirofilaria_global.cohort.2025-06-18.nuclearSNPs..final.autosomal.recode.vcf.gz
tabix dirofilaria_global.cohort.2025-06-18.nuclearSNPs..final.autosomal.recode.vcf.gz

# Get sample names for each population
ASIA=$(awk -F'_' '$1 == "MYS" || $1 == "THA" {print}' keep_samples.dimmitis-only.nuclear_snps.no-reps.list | paste -sd ',')
AUS=$(awk -F'_' '$1 == "AUS" {print}' keep_samples.dimmitis-only.nuclear_snps.no-reps.list | paste -sd ',')
CENAM=$(awk -F'_' '$1 == "PAN" || $1 == "CRI" {print}' keep_samples.dimmitis-only.nuclear_snps.no-reps.list | paste -sd ',')
EUR=$(awk -F'_' '$1 == "GRC" || $1 == "ITA" || $1 == "ROU" {print}' keep_samples.dimmitis-only.nuclear_snps.no-reps.list | paste -sd ',')
USA=$(awk -F'_' '$1 == "USA" {print}' keep_samples.dimmitis-only.nuclear_snps.no-reps.list | paste -sd ',')
declare -A populations
populations=(
  [ASIA]=$ASIA
  [AUS]=$AUS
  [CENAM]=$CENAM
  [EUR]=$EUR
  [USA]=$USA
)

awk -F'_' '$1 == "MYS" || $1 == "THA" {print}' keep_samples.dimmitis-only.nuclear_snps.no-reps.list > ASIA.keep
awk -F'_' '$1 == "AUS" {print}' keep_samples.dimmitis-only.nuclear_snps.no-reps.list > AUS.keep
awk -F'_' '$1 == "PAN" || $1 == "CRI" {print}' keep_samples.dimmitis-only.nuclear_snps.no-reps.list > CENAM.keep
awk -F'_' '$1 == "GRC" || $1 == "ITA" || $1 == "ROU" {print}' keep_samples.dimmitis-only.nuclear_snps.no-reps.list  > EUR.keep
awk -F'_' '$1 == "USA" {print}' keep_samples.dimmitis-only.nuclear_snps.no-reps.list > USA.keep

# generate VCFs for each population, based on keep lists
for i in *keep; do 
    vcftools --gzvcf dirofilaria_global.cohort.2025-06-18.nuclearSNPs..final.autosomal.recode.vcf.gz \
    --max-missing 1 \
    --keep ${i} \
    --recode  \
    --stdout | bgzip > ${i%.keep}.vcf.gz;
    tabix ${i%.keep}.vcf.gz;
    done
```



### generate mask files for smc++

- Approach 1: mask everything apart from variable sites
```bash
# Generate mask file
VCF=dirofilaria_global.cohort.2025-06-18.nuclearSNPs..final.autosomal.recode.vcf
FAI=dimmitis_WSI_2.2.fai

# Step 1: Get VCF sites as BED
bcftools query -f '%CHROM\t%POS0\t%POS\n' cleaned.vcf.gz | sort -k1,1 -k2,2n > vcf_sites.bed

# Step 2: Create genome BED
awk -v OFS="\t" '{print $1, 0, $2}' $FAI > genome_full.bed

# Step 3: Subtract to get missing regions
bedtools subtract -a genome_full.bed -b vcf_sites.bed > missing_sites.bed
bgzip missing_sites.bed
tabix missing_sites.bed.gz
```

- Approach 2: mask regions in which reads do not map uniquely, based on Heng Li's SNPable scripts
```bash
cd /lustre/scratch127/pam/users/sd21/dirofilaria_immitis/project_global-diversity/04_analysis/05_smc++/mappability_mask

ln -s ../../../01_references/dimmitis_WSI_2.2.fa

/nfs/users/nfs_s/sd21/lustre_link/software/seqbility-20091110/splitfa dimmitis_WSI_2.2.fa 35 | split -l 20000000

bwa index dimmitis_WSI_2.2.fa

echo -e 'for i in x*; do bwa aln -R 1000000 -O 3 -E 3 dimmitis_WSI_2.2.fa ${i} > ${i}.sai; done' > run_bwa

chmod a+x run_bwa

bsub.py 10 bwaaln ./run_bwa

for i in *sai; do bwa samse -f ${i%.sai}.sam dimmitis_WSI_2.2.fa ${i} ${i%.sai}; done

# once mapping is completed, compress to save space
gzip *.sam

# make the raw mask
gzip -dc x??.sam.gz | /nfs/users/nfs_s/sd21/lustre_link/software/seqbility-20091110/gen_raw_mask.pl > rawMask_35.fa

# make the final mask
/nfs/users/nfs_s/sd21/lustre_link/software/seqbility-20091110/gen_mask -l 35 -r 0.5 rawMask_35.fa > mask_35_50.fa

# make bed files per chromosome of the category 3 positions, tool from here: https://github.com/stschiff/msmc-tools 
python makeMappabilityMask.py

# position categories
# c=3: the majortiy of overlapping 35-mers are mapped uniquely and without 1-mismatch (or 1-difference, depending on the BWA command line) hits.
# c=2: the majority of overlapping 35-mers are unique and c!=3.
# c=1: the majority of overlapping 35-mers are non-unique.
# c=0: all the 35-mers overlapping x cannot be mapped due to excessive ambiguous bases.

# count how many positions for each position in the genome
for i in 0 1 2 3; do
     echo -e "SNPtype: ${i}";
     cat mask_35_50.fa | grep -v ">" | grep -o "${i}" | wc -l;
done

#SNPtype: 0
39671
#SNPtype: 1
7125343
#SNPtype: 2
1476461
#SNPtype: 3
80834336

zcat *mask.bed.gz > dimmitis.mapability-mask.bed

samtools faidx dimmitis_WSI_2.2.fa
cat dimmitis_WSI_2.2.fa.fai | awk '{print $1,$2}' OFS="\t" > dimmitis_WSI_2.2.genome

bedtools complement -i dimmitis.mapability-mask.bed -g dimmitis_WSI_2.2.genome > dimmitis.mapability-mask.bad-regions.bed

cat dimmitis.mapability-mask.bad-regions.bed | bgzip > dimmitis.mapability-mask.bad-regions.bed.gz; tabix dimmitis.mapability-mask.bad-regions.bed.gz

```









```bash

get_samples() {
  local pattern="$1"
  awk -F'_' "$pattern {print}" keep_samples.dimmitis-only.nuclear_snps.no-reps.list | paste -sd ','
}

# Extract population-specific sample lists
ASIA=$(get_samples '$1 == "MYS" || $1 == "THA"')
AUS=$(get_samples '$1 == "AUS"')
CENAM=$(get_samples '$1 == "PAN" || $1 == "CRI"')
EUR=$(get_samples '$1 == "GRC" || $1 == "ITA" || $1 == "ROU"')
USA=$(get_samples '$1 == "USA"')

# Declare associative array
declare -A populations=(
  [ASIA]="$ASIA"
  [AUS]="$AUS"
  [CENAM]="$CENAM"
  [EUR]="$EUR"
  [USA]="$USA"
)

mkdir DATA_MASK
mkdir DATA_NO-MASK
mkdir DATA_MASK_MAP

# masked - all
for pop in "${!populations[@]}"; do
    samples="${populations[$pop]}"
    vcf_file="${pop}.vcf.gz"  # e.g. ASIA.vcf.gz

    # Step 1: Convert VCF to SMC
    for chr in {1..4}; do
        smc++ vcf2smc \
            "$vcf_file" \
            --mask missing_sites.bed.gz \
            "DATA_MASK/${pop}.chr${chr}.smc.gz" \
            "dirofilaria_immitis_chr${chr}" \
            "${pop}:${samples}"
    done

    # Step 2: Estimate demographic history
    smc++ estimate --timepoints 1 1000000 -o "RESULTS_MASK/${pop}/" 2.7e-9 \
        DATA_MASK/${pop}.chr*.smc.gz

    # Step 3: Plot results with different generation times
    for g in 1 2.5 4; do
        smc++ plot -g "$g" -c "SMCPP_${pop}_g${g}_mask.pdf" "RESULTS_MASK/${pop}/model.final.json"
    done
done



# no mask
for pop in "${!populations[@]}"; do
    samples="${populations[$pop]}"
    vcf_file="${pop}.vcf.gz"  # e.g. ASIA.vcf.gz

    # Step 1: Convert VCF to SMC
    for chr in {1..4}; do
        smc++ vcf2smc \
            "$vcf_file" \
            "DATA_NO-MASK/${pop}.chr${chr}.smc.gz" \
            "dirofilaria_immitis_chr${chr}" \
            "${pop}:${samples}"
    done

    # Step 2: Estimate demographic history
    smc++ estimate --timepoints 1 1000000 -o "RESULTS_NO-MASK/${pop}/" 2.7e-9 \
        DATA_NO-MASK/${pop}.chr*.smc.gz

    # Step 3: Plot results with different generation times
    for g in 1 2.5 4; do
        smc++ plot -g "$g" -c "SMCPP_${pop}_g${g}_no-mask.pdf" "RESULTS_NO-MASK/${pop}/model.final.json"
    done
done



# masked 
for pop in "${!populations[@]}"; do
    samples="${populations[$pop]}"
    vcf_file="${pop}.vcf.gz"  # e.g. ASIA.vcf.gz

    # Step 1: Convert VCF to SMC
    for chr in {1..4}; do
        smc++ vcf2smc \
            "$vcf_file" \
            --mask mappability_mask/dimmitis.mapability-mask.bad-regions.bed.gz \
            "DATA_MASK_MAP/${pop}.chr${chr}.smc.gz" \
            "dirofilaria_immitis_chr${chr}" \
            "${pop}:${samples}"
    done

    # Step 2: Estimate demographic history
    smc++ estimate --timepoints 1 1000000 -o "RESULTS_MASK_MAP/${pop}/" 2.7e-9 \
        DATA_MASK_MAP/${pop}.chr*.smc.gz

    # Step 3: Plot results with different generation times
    for g in 1 2.5 4; do
        smc++ plot -g "$g" -c "SMCPP_${pop}_g${g}_mask-map.pdf" "RESULTS_MASK_MAP/${pop}/model.final.json"
    done
done

```



### compare the no mask, all mask, and SNPable masks

```R
plot_ne_region_per_generation <- function(region, generation){

require(tidyverse)

data_mask <- read.table(paste0("SMCPP_",region,"_g",generation,"_mask.csv"), header=T, sep=",")
data_mask$group <- "full_mask"

data_map_mask <- read.table(paste0("SMCPP_",region,"_g",generation,"_mask-map.csv"), header=T, sep=",")
data_map_mask$group <- "mapping_mask"

data_no_mask <- read.table(paste0("SMCPP_",region,"_g",generation,"_no-mask.csv"), header=T, sep=",")
data_no_mask$group <- "no_mask"

data <- bind_rows(data_mask, data_map_mask, data_no_mask)

ggplot(data, aes(x, y, col=group)) +
  geom_rect(aes(xmin=5000, xmax=12000, ymin=1, ymax=1e7), fill="grey", col=NA) +
  geom_rect(aes(xmin=14000, xmax=40000, ymin=1, ymax=1e7), fill="grey", col=NA) +
  geom_rect(aes(xmin=116000, xmax=130000, ymin=1, ymax=1e7), fill="grey", col=NA) +
  geom_vline(xintercept=500) +
  geom_line(linewidth=1.5) + 
  labs(title = paste0("Region = ", region, ", Generation time (yrs) =", generation ), y="Ne", x="Years before present") + 
  scale_y_log10(limits = c(1,1e7)) + 
  scale_x_log10(limits = c(5,1e7)) +
  theme_bw()

}

plot_Ne_ASIA <- plot_ne_region_per_generation("ASIA", 1)
plot_Ne_AUS <- plot_ne_region_per_generation("AUS", 1)
plot_Ne_USA <- plot_ne_region_per_generation("USA", 1)
plot_Ne_CENAM <- plot_ne_region_per_generation("CENAM", 1)
plot_Ne_EUR <- plot_ne_region_per_generation("EUR", 1)

library(patchwork)

plot_Ne_ASIA + plot_Ne_AUS + plot_Ne_USA + plot_Ne_CENAM + plot_Ne_EUR + plot_layout(guides = "collect")

plot_Ne_ASIA <- plot_ne_region_per_generation("ASIA", 2.5)
plot_Ne_AUS <- plot_ne_region_per_generation("AUS", 2.5)
plot_Ne_USA <- plot_ne_region_per_generation("USA", 2.5)
plot_Ne_CENAM <- plot_ne_region_per_generation("CENAM", 2.5)
plot_Ne_EUR <- plot_ne_region_per_generation("EUR", 2.5)

library(patchwork)

plot_Ne_ASIA + plot_Ne_AUS + plot_Ne_USA + plot_Ne_CENAM + plot_Ne_EUR + plot_layout(guides = "collect")



plot_Ne_ASIA <- plot_ne_region_per_generation("ASIA", 4)
plot_Ne_AUS <- plot_ne_region_per_generation("AUS", 4)
plot_Ne_USA <- plot_ne_region_per_generation("USA", 4)
plot_Ne_CENAM <- plot_ne_region_per_generation("CENAM", 4)
plot_Ne_EUR <- plot_ne_region_per_generation("EUR", 4)

library(patchwork)

plot_Ne_ASIA + plot_Ne_AUS + plot_Ne_USA + plot_Ne_CENAM + plot_Ne_EUR + plot_layout(guides = "collect")


```

















### Run jackknifing
- want to generate some confidence intervals on the smc++ estimates
- will do it by rerunning smc++ multiple times - the first with all chromosomes, but then subsequently, dropping one chromosome at a time.

```bash
module load smcpp/1.15.4--py39hac1eaed_0
module load htslib-1.19/perl-5.38.0


bsub.py --queue long 10 smc_jackknife "bash ./run_smcpp_jackknife.sh
```
- where: run_smcpp_jackknife.sh

```bash


get_samples() {
  local pattern="$1"
  awk -F'_' "$pattern {print}" keep_samples.dimmitis-only.nuclear_snps.no-reps.list | paste -sd ','
}

# Extract population-specific sample lists
ASIA=$(get_samples '$1 == "MYS" || $1 == "THA"')
AUS=$(get_samples '$1 == "AUS"')
CENAM=$(get_samples '$1 == "PAN" || $1 == "CRI"')
EUR=$(get_samples '$1 == "GRC" || $1 == "ITA" || $1 == "ROU"')
USA=$(get_samples '$1 == "USA"')

# Declare associative array
declare -A populations=(
  [ASIA]="$ASIA"
  [AUS]="$AUS"
  [CENAM]="$CENAM"
  [EUR]="$EUR"
  [USA]="$USA"
)


# masked
for pop in "${!populations[@]}"; do
    samples="${populations[$pop]}"

    mkdir -p "DATA_${pop}"

    # Step 1: Convert VCF to SMC
    for chr in {1..4}; do
        chr_name="dirofilaria_immitis_chr${chr}"
        out_vcf="${pop}.chr${chr}.vcf.gz"

        # Subset VCF by chromosome
        vcftools --gzvcf "${pop}.vcf.gz" \
                 --chr "$chr_name" \
                 --recode --stdout | bgzip -c > "DATA_${pop}"/"$out_vcf"
        tabix -p vcf "DATA_${pop}"/"$out_vcf"

        # Convert to SMC++
        smc++ vcf2smc \
            "DATA_${pop}"/"$out_vcf" \
            --mask missing_sites.bed.gz \
            "DATA_${pop}/${pop}.chr${chr}.smc.gz" \
            "$chr_name" \
            "${pop}:${samples}"
        done
done

for pop in AUS CENAM EUR USA ASIA; do 

  mkdir -p "RESULT_${pop}"

  bsub.py 10 --threads 4 smc_estimate_${pop}_chr_all  "smc++ estimate -o RESULT_${pop}/jackknife_${pop}_all_chr --timepoints 1 1000000 --cores 4 2.7e-9 DATA_${pop}/${pop}.chr[1234].smc.gz"
  bsub.py 10 --threads 4  smc_estimate_${pop}_chr-minus-c1  "smc++ estimate -o RESULT_${pop}/jackknife_${pop}_minus_chr1 --timepoints 1 1000000 --cores 4 2.7e-9 DATA_${pop}/${pop}.chr[234].smc.gz"
  bsub.py 10 --threads 4  smc_estimate_${pop}_chr-minus-c2  "smc++ estimate -o RESULT_${pop}/jackknife_${pop}_minus_chr2 --timepoints 1 1000000 --cores 4 2.7e-9 DATA_${pop}/${pop}.chr[134].smc.gz"
  bsub.py 10 --threads 4  smc_estimate_${pop}_chr-minus-c3  "smc++ estimate -o RESULT_${pop}/jackknife_${pop}_minus_chr3 --timepoints 1 1000000 --cores 4 2.7e-9 DATA_${pop}/${pop}.chr[124].smc.gz"
  bsub.py 10 --threads 4  smc_estimate_${pop}_chr-minus-c4  "smc++ estimate -o RESULT_${pop}/jackknife_${pop}_minus_chr4 --timepoints 1 1000000 --cores 4 2.7e-9 DATA_${pop}/${pop}.chr[123].smc.gz";

done
```






# plot and generate CSVs
```bash

for pop in AUS CENAM EUR USA ASIA; do 
  for g in 1 2.5 4 ; do
    for chromosome in all_chr minus_chr1 minus_chr2 minus_chr3 minus_chr4; do

    smc++ plot -g ${g} -c RESULT_${pop}/SMCPP_${pop}_g${g}_${chromosome}.pdf RESULT_${pop}/jackknife_${pop}_${chromosome}/model.final.json;

    done
  done
done

```







### testing the variation in jackknife data
```R
plot_ne_region_per_generation <- function(region, generation){

require(tidyverse)

data_all_chr <- read.table(paste0("RESULT_",region,"/SMCPP_",region,"_g",generation,"_all_chr.csv"), header=T, sep=",")
data_all_chr$group <- "all_chr"

data_minus_chr1 <- read.table(paste0("RESULT_",region,"/SMCPP_",region,"_g",generation,"_minus_chr1.csv"), header=T, sep=",")
data_minus_chr1$group <- "minus_chr1"

data_minus_chr2 <- read.table(paste0("RESULT_",region,"/SMCPP_",region,"_g",generation,"_minus_chr2.csv"), header=T, sep=",")
data_minus_chr2$group <- "minus_chr2"

data_minus_chr3 <- read.table(paste0("RESULT_",region,"/SMCPP_",region,"_g",generation,"_minus_chr3.csv"), header=T, sep=",")
data_minus_chr3$group <- "minus_chr3"

data_minus_chr4 <- read.table(paste0("RESULT_",region,"/SMCPP_",region,"_g",generation,"_minus_chr4.csv"), header=T, sep=",")
data_minus_chr4$group <- "minus_chr4"


data <- bind_rows(data_all_chr, data_minus_chr1, data_minus_chr2, data_minus_chr3, data_minus_chr4)

ggplot(data, aes(x, y, col=group)) +
  geom_rect(aes(xmin=5000, xmax=12000, ymin=1, ymax=1e7), fill="grey", col=NA) +
  geom_rect(aes(xmin=14000, xmax=40000, ymin=1, ymax=1e7), fill="grey", col=NA) +
  geom_rect(aes(xmin=116000, xmax=130000, ymin=1, ymax=1e7), fill="grey", col=NA) +
  geom_vline(xintercept=500) +
  geom_line(linewidth=1.5) + 
  labs(title = paste0("Region = ", region, ", Generation time (yrs) =", generation ), y="Ne", x="Years before present") + 
  scale_y_log10(limits = c(1,1e7)) + 
  scale_x_log10(limits = c(5,1e7)) +
  theme_bw()

}

plot_AUS_jk <- plot_ne_region_per_generation("AUS", 2.5)

plot_CENAM_jk <- plot_ne_region_per_generation("CENAM", 2.5)



plot_USA_jk <-  plot_ne_region_per_generation("USA", 2.5)


plot_EUR_jk <- plot_ne_region_per_generation("EUR", 2.5)

plot_ASIA_jk <- plot_ne_region_per_generation("ASIA", 2.5)

plot_AUS_jk + plot_CENAM_jk + plot_USA_jk + plot_EUR_jk + plot_ASIA_jk + plot_layout(guides = "collect")




```



```R
library(tidyverse)
library(scales)
library(patchwork)


data_files <- list.files(path = Sys.glob("RESULT*"), pattern = "\\.csv$", full.names = TRUE)


# load data using file names, and make a formatted data frame
data <-
        purrr::map_df(data_files, function(x) {
	      
        fname <- basename(x)
        fname_no_ext <- str_remove(fname, "\\.csv$")
        parts <- str_split(fname_no_ext, "_", simplify = TRUE)
        

        read.delim(x, header = T, sep=",") %>% 
          mutate(filename = fname_no_ext) %>% 
          bind_cols(as_tibble(as.data.frame(parts)))
      	} )

colnames(data) <- c("region", "time", "Ne", "plot_type", "plot_num", "sample_name", "tool", "region.2", "generation", "jk", "chr")

all_data <- data %>% select(sample_name, region, generation, chr, time, Ne, jk)

# geom_smooth, all populations, different generations
# ggplot(all_data, aes(x = time, y = Ne, col = region, fill = region)) + geom_smooth(span = 0.1, alpha=0.1) + scale_x_log10() + scale_y_log10() + facet_grid(.~generation) + theme_bw()

#pop_colours <- c("ASIA" = "hotpink", "AUS" = "cornflowerblue", "CENAM" = "purple", "EUR" = "forestgreen", "USA" = "tomato2")

# ggplot(data=subset(all_data, jk=="all"), aes(x = time, y = Ne, col = region, fill = region)) + geom_line() + scale_x_log10() + scale_y_log10() + facet_grid(generation~.) + theme_bw() +
#   scale_colour_manual(values = pop_colours) +
#   scale_fill_manual(values = pop_colours)


data_g2.5 <- filter(all_data, generation == "g2.5")
data_g1 <- filter(all_data, generation == "g1")
data_g4 <- filter(all_data, generation == "g4")

pop_colours <- c("ASIA" = "hotpink", "AUS" = "cornflowerblue", "CENAM" = "purple", "EUR" = "forestgreen", "USA" = "tomato2")

# plot_allpop_g2.5_smooth <- 
#   ggplot(data_g2.5, aes(x = time, y = Ne, col = region, fill = region)) + 
#   geom_rect(aes(xmin=14000, xmax=40000, ymin=1e1, ymax=1e7), fill="grey", col=NA) +
#   scale_x_log10() + 
#   stat_smooth(span = 0.05, alpha=0.1, linetype="dashed") + 
#   scale_y_log10() + 
#   theme_bw() +
#   scale_colour_manual(values = pop_colours)+
#   scale_fill_manual(values = pop_colours)

# plot_allpop_g2.5_allchr <- 
#   ggplot(data_g2.5, aes(x = time, y = Ne, col = region, fill = region)) +
#   geom_rect(aes(xmin=14000, xmax=40000, ymin=1e1, ymax=1e7), fill="grey", col=NA) +
#   scale_x_log10() + 
#     geom_line(data=subset(data_g2.5,jk=="all"), aes(x = time, y = Ne, col = region), linewidth=1) +
#   scale_y_log10() + 
#   theme_bw() +
#   scale_colour_manual(values = pop_colours) +
#   scale_fill_manual(values = pop_colours)


# plot_allpop_g2.5_smooth + plot_allpop_g2.5_allchr





plot_g2.5 <-  ggplot(data_g2.5, aes(x = time, y = Ne, group = sample_name)) +  
  geom_rect(aes(xmin=14000, xmax=40000, ymin=1e1, ymax=1e7), fill="lightgrey", col=NA) +
  geom_rect(aes(xmin=116000, xmax=130000, ymin=1e1, ymax=1e7), fill="black", col=NA) +
  geom_rect(aes(xmin=58000, xmax=72000, ymin=1e1, ymax=1e7), fill="darkgrey", col=NA) +
  geom_line(aes(col=region), alpha=0.3, linewidth=0.5) + 
  geom_line(data=subset(data_g2.5,jk=="all"), aes(x = time, y = Ne, col = region), linewidth=1.5, alpha=0.75)+
  scale_y_log10(labels = comma) + 
  scale_x_log10(labels = comma) +
  theme_bw() + 
  scale_colour_manual(values = pop_colours)+
  labs(x="Years before present", y="Effective population size (Ne)")+
  theme(legend.position="bottom")

plot_g2.5

ggsave("figure_Ne_all-populations.pdf", height=100, width=100, units="mm")

plot_g1 <-  ggplot(data_g1, aes(x = time, y = Ne, group = sample_name)) +  
  geom_rect(aes(xmin=14000, xmax=40000, ymin=1e1, ymax=1e7), fill="lightgrey", col=NA) +
  geom_rect(aes(xmin=116000, xmax=130000, ymin=1e1, ymax=1e7), fill="black", col=NA) +
  geom_rect(aes(xmin=58000, xmax=72000, ymin=1e1, ymax=1e7), fill="darkgrey", col=NA) +
  geom_line(aes(col=region), alpha=0.3, linewidth=0.5) + 
  geom_line(data=subset(data_g1,jk=="all"), aes(x = time, y = Ne, col = region), linewidth=1.5, alpha=0.75)+
  scale_y_log10(labels = comma) + 
  scale_x_log10(labels = comma) +
  theme_bw() +
  scale_colour_manual(values = pop_colours)+
  labs(x="Years before present", y="Effective population size (Ne)")+
  theme(legend.position="bottom")

plot_g4 <-  ggplot(data_g4, aes(x = time, y = Ne, group = sample_name)) +  
  geom_rect(aes(xmin=14000, xmax=40000, ymin=1e1, ymax=1e7), fill="lightgrey", col=NA) +
  geom_rect(aes(xmin=116000, xmax=130000, ymin=1e1, ymax=1e7), fill="black", col=NA) +
  geom_rect(aes(xmin=58000, xmax=72000, ymin=1e1, ymax=1e7), fill="darkgrey", col=NA) +
  geom_line(aes(col=region), alpha=0.3, linewidth=0.5) + 
  geom_line(data=subset(data_g4,jk=="all"), aes(x = time, y = Ne, col = region), linewidth=1.5, alpha=0.75)+
  scale_y_log10(labels = comma) + 
  scale_x_log10(labels = comma) +
  theme_bw() + 
  scale_colour_manual(values = pop_colours)+
  labs(x="Years before present", y="Effective population size (Ne)")+
  theme(legend.position="bottom")

plot_g1 + plot_g2.5 + plot_g4 + plot_layout(ncol=1, guides = "collect") & theme(legend.position = 'bottom')










  # ggplot(data_g2.5, aes(x = time, y = Ne)) + 
  # geom_smooth(data=subset(data_g2.5,jk=="minus"), aes(x = time, y = Ne, col = region, fill = region), span = 0.05, linewidth=1, alpha=0.1)+
  # scale_y_log10(labels = comma) + 
  # scale_x_log10(labels = comma) +
  # theme_void() +
  # scale_colour_manual(values = pop_colours)+
  # scale_fill_manual(values = pop_colours) + 
  # theme(legend.position="none") +
  # labs(x="Years before present", y="Effective population size (Ne)")













data_CENAM_EUR  <- filter(data_g2.5, region == "CENAM" | region == "EUR")

# plot_CENAM_EUR_g2.5_smooth <- 
#   ggplot(data_CENAM_EUR, aes(x = time, y = Ne, col = region, fill = region)) + 
#   geom_rect(aes(xmin=5000, xmax=12000, ymin=1e1, ymax=1e7), fill="grey", col=NA) +
#   geom_rect(aes(xmin=14000, xmax=40000, ymin=1e1, ymax=1e7), fill="grey", col=NA) +
#   geom_rect(aes(xmin=116000, xmax=130000, ymin=1e1, ymax=1e7), fill="grey", col=NA) +
#   geom_vline(xintercept=500) +
#   scale_x_log10() + 
#   geom_smooth(span = 0.1, alpha=0.1, linetype="dashed") + 
#   scale_y_log10() + 
#   facet_grid(generation~.) + 
#   theme_bw() +
#   scale_colour_manual(values = pop_colours)+
#   scale_fill_manual(values = pop_colours)

# plot_CENAM_EUR_g2.5_allchr <- 
#   ggplot(data_CENAM_EUR, aes(x = time, y = Ne, col = region, fill = region)) +
#   geom_rect(aes(xmin=5000, xmax=12000, ymin=1e1, ymax=1e7), fill="grey", col=NA) +
#   geom_rect(aes(xmin=14000, xmax=40000, ymin=1e1, ymax=1e7), fill="grey", col=NA) +
#   geom_rect(aes(xmin=116000, xmax=130000, ymin=1e1, ymax=1e7), fill="grey", col=NA) +
#   geom_vline(xintercept=500) + 
#   scale_x_log10() + 
#     geom_line(data=subset(data_CENAM_EUR,jk=="all"), aes(x = time, y = Ne, col = region), linewidth=1) +
#   scale_y_log10() + 
#   facet_grid(generation~.) + 
#   theme_bw() +
#   scale_colour_manual(values = pop_colours) +
#   scale_fill_manual(values = pop_colours)

# library(patchwork)

# plot_CENAM_EUR_g2.5_smooth + plot_CENAM_EUR_g2.5_allchr + plot_layout(guides = "collect")

plot_CENAM_EUR <- 
  ggplot(data_CENAM_EUR, aes(x = time, y = Ne, group = sample_name)) +  
  geom_rect(aes(xmin=14000, xmax=40000, ymin=1e1, ymax=1e7), fill="lightgrey", col=NA) +
  geom_line(aes(col=region), alpha=0.4, linewidth=0.5) + 
  geom_line(data=subset(data_CENAM_EUR,jk=="all"), aes(x = time, y = Ne, col = region), linewidth=1.5, alpha=0.75)+
  geom_vline(xintercept=500, linetype="dashed") + 
  scale_y_log10(labels = comma) + 
  scale_x_log10(labels = comma) +
  theme_bw() +
  scale_colour_manual(values = pop_colours)+
  scale_fill_manual(values = pop_colours) +
  labs(x="Years before present", y="Effective population size (Ne)")+
  theme(legend.position="bottom")








data_AUS_ASIA  <- filter(data_g2.5, region == "AUS" | region == "ASIA")

# plot_AUS_ASIA_g2.5_smooth <- 
#   ggplot(data_AUS_ASIA, aes(x = time, y = Ne, col = region, fill = region)) + 
#   geom_rect(aes(xmin=5000, xmax=12000, ymin=1e1, ymax=1e7), fill="grey", col=NA) +
#   geom_rect(aes(xmin=14000, xmax=40000, ymin=1e1, ymax=1e7), fill="grey", col=NA) +
#   geom_rect(aes(xmin=116000, xmax=130000, ymin=1e1, ymax=1e7), fill="grey", col=NA) +
#   geom_vline(xintercept=500) +
#   scale_x_log10() + 
#   geom_smooth(span = 0.1, alpha=0.1, linetype="dashed") + 
#   scale_y_log10() + 
#   facet_grid(generation~.) + 
#   theme_bw() +
#   scale_colour_manual(values = pop_colours)+
#   scale_fill_manual(values = pop_colours)

# plot_AUS_ASIA_g2.5_allchr <- 
#   ggplot(data_AUS_ASIA, aes(x = time, y = Ne, col = region, fill = region)) +
#   geom_rect(aes(xmin=5000, xmax=12000, ymin=1e1, ymax=1e7), fill="grey", col=NA) +
#   geom_rect(aes(xmin=14000, xmax=40000, ymin=1e1, ymax=1e7), fill="grey", col=NA) +
#   geom_rect(aes(xmin=116000, xmax=130000, ymin=1e1, ymax=1e7), fill="grey", col=NA) +
#   geom_vline(xintercept=500) + 
#   scale_x_log10() + 
#     geom_line(data=subset(data_AUS_ASIA,jk=="all"), aes(x = time, y = Ne, col = region), linewidth=1) +
#   scale_y_log10() + 
#   facet_grid(generation~.) + 
#   theme_bw() +
#   scale_colour_manual(values = pop_colours) +
#   scale_fill_manual(values = pop_colours)

# library(patchwork)

# plot_AUS_ASIA_g2.5_smooth + plot_AUS_ASIA_g2.5_allchr + plot_layout(guides = "collect")


plot_AUS_ASIA <-
  ggplot(data_AUS_ASIA, aes(x = time, y = Ne, group = sample_name)) +  
  geom_rect(aes(xmin=14000, xmax=40000, ymin=1e1, ymax=1e7), fill="lightgrey", col=NA) +
  geom_rect(aes(xmin=5000, xmax=12000, ymin=1e1, ymax=1e7), fill="grey", col=NA) +
  geom_line(aes(col=region), alpha=0.4, linewidth=0.5) + 
  geom_line(data=subset(data_AUS_ASIA,jk=="all"), aes(x = time, y = Ne, col = region), linewidth=1.5, alpha=0.75)+
  scale_y_log10(labels = comma) + 
  scale_x_log10(labels = comma) +
  theme_bw() +
  scale_colour_manual(values = pop_colours)+
  scale_fill_manual(values = pop_colours) +
  labs(x="Years before present", y="Effective population size (Ne)")+
  theme(legend.position="bottom")







data_USA_ASIA  <- filter(data_g2.5, region == "USA" | region == "ASIA")

# plot_USA_ASIA_g2.5_smooth <- 
#   ggplot(data_USA_ASIA, aes(x = time, y = Ne, col = region, fill = region)) + 
#   geom_rect(aes(xmin=5000, xmax=12000, ymin=1e1, ymax=1e7), fill="grey", col=NA) +
#   geom_rect(aes(xmin=14000, xmax=40000, ymin=1e1, ymax=1e7), fill="grey", col=NA) +
#   geom_rect(aes(xmin=116000, xmax=130000, ymin=1e1, ymax=1e7), fill="grey", col=NA) +
#   geom_vline(xintercept=500) +
#   scale_x_log10() + 
#   geom_smooth(span = 0.1, alpha=0.1, linetype="dashed") + 
#   scale_y_log10() + 
#   facet_grid(generation~.) + 
#   theme_bw() +
#   scale_colour_manual(values = pop_colours)+
#   scale_fill_manual(values = pop_colours)

# plot_USA_ASIA_g2.5_allchr <- 
#   ggplot(data_USA_ASIA, aes(x = time, y = Ne, col = region, fill = region)) +
#   geom_rect(aes(xmin=5000, xmax=12000, ymin=1e1, ymax=1e7), fill="grey", col=NA) +
#   geom_rect(aes(xmin=14000, xmax=40000, ymin=1e1, ymax=1e7), fill="grey", col=NA) +
#   geom_rect(aes(xmin=116000, xmax=130000, ymin=1e1, ymax=1e7), fill="grey", col=NA) +
#   geom_vline(xintercept=500) + 
#   scale_x_log10() + 
#     geom_line(data=subset(data_USA_ASIA,jk=="all"), aes(x = time, y = Ne, col = region), linewidth=1) +
#   scale_y_log10() + 
#   facet_grid(generation~.) + 
#   theme_bw() +
#   scale_colour_manual(values = pop_colours) +
#   scale_fill_manual(values = pop_colours)

# library(patchwork)

# plot_USA_ASIA_g2.5_smooth + plot_USA_ASIA_g2.5_allchr + plot_layout(guides = "collect")


plot_USA_ASIA <-
  ggplot(data_USA_ASIA, aes(x = time, y = Ne, group = sample_name)) +  
  geom_rect(aes(xmin=14000, xmax=40000, ymin=1e1, ymax=1e7), fill="grey90", col=NA, alpha=0.1) +
  geom_rect(aes(xmin=11000, xmax=30000, ymin=1e1, ymax=1e7), fill="grey", col=NA, alpha=0.1) +
  geom_rect(aes(xmin=60000, xmax=70000, ymin=1e1, ymax=1e7), fill="grey", col=NA, alpha=0.1) +
  geom_line(aes(col=region), alpha=0.4, linewidth=0.5) + 
  geom_line(data=subset(data_USA_ASIA,jk=="all"), aes(x = time, y = Ne, col = region), linewidth=1.5, alpha=0.75)+
  scale_y_log10(labels = comma) + 
  scale_x_log10(labels = comma) +
  theme_bw() +
  scale_colour_manual(values = pop_colours)+
  scale_fill_manual(values = pop_colours) +
  labs(x="Years before present", y="Effective population size (Ne)") +
  theme(legend.position="bottom")




plot_CENAM_EUR + plot_AUS_ASIA + plot_USA_ASIA

ggsave("figure_Ne_pairwise_hypothesis-testing.pdf", height=100, width=300, units="mm")

```









```R
summary_data <- all_data %>%
  group_by(region, time, generation) %>%
  summarise(
    Ne_median = median(Ne),
    Ne_lower = quantile(Ne, 0.025),
    Ne_upper = quantile(Ne, 0.975),
    .groups = "drop"
  )

pop_colours <- c("ASIA" = "hotpink", "AUS" = "cornflowerblue", "CENAM" = "purple", "EUR" = "forestgreen", "USA" = "tomato2")

test <- filter(summary_data, generation == "g2.5")

ggplot(summary_data, aes(x = time, y = Ne_median, col = region, fill = region)) +
  geom_rect(aes(xmin=5000, xmax=12000, ymin=1e1, ymax=1e7), fill="grey", col=NA) +
  geom_rect(aes(xmin=14000, xmax=40000, ymin=1e1, ymax=1e7), fill="grey", col=NA) +
  geom_rect(aes(xmin=116000, xmax=130000, ymin=1e1, ymax=1e7), fill="grey", col=NA) +
  geom_vline(xintercept=500) +
  geom_ribbon(aes(ymin = Ne_lower, ymax = Ne_upper), alpha = 0.1) +
  #geom_line(linewidth=1.5) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    title = "SMC++ Bootstrap: Effective Population Size by Population",
    x = "Generations ago (log scale)",
    y = "Effective population size (log scale)"
  ) +
  theme_minimal() +
  theme(legend.position = "top") + 
  facet_grid(region ~ generation) +
  scale_colour_manual(values = pop_colours) +
  scale_fill_manual(values = pop_colours)


ggplot(test, aes(x = time, y = Ne_median, col = region, fill = region)) +
  geom_rect(aes(xmin=5000, xmax=12000, ymin=1e1, ymax=1e7), fill="grey", col=NA) +
  geom_rect(aes(xmin=14000, xmax=40000, ymin=1e1, ymax=1e7), fill="grey", col=NA) +
  geom_rect(aes(xmin=116000, xmax=130000, ymin=1e1, ymax=1e7), fill="grey", col=NA) +
  geom_vline(xintercept=500) +
  geom_ribbon(aes(ymin = Ne_lower, ymax = Ne_upper), alpha = 0.1) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    title = "SMC++ Bootstrap: Effective Population Size by Population",
    x = "Generations ago (log scale)",
    y = "Effective population size (log scale)"
  ) +
  theme_minimal() +
  theme(legend.position = "top") + 
  facet_grid(region ~ .) +
  scale_colour_manual(values = pop_colours)








all_chr <- filter(all_data, jk == "all")

pop_colours <- c("ASIA" = "hotpink", "AUS" = "cornflowerblue", "CENAM" = "purple", "EUR" = "forestgreen", "USA" = "tomato2")

ggplot(all_chr, aes(x = time, y = Ne, col = region)) +
  geom_line(linewidth=1.5) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    title = "SMC++ Bootstrap: Effective Population Size by Population",
    x = "Generations ago (log scale)",
    y = "Effective population size (log scale)"
  ) +
  theme_minimal() +
  theme(legend.position = "top") + 
  facet_grid(. ~ generation) +
  scale_colour_manual(values = pop_colours)



minus_chr <- filter(all_data, jk == "minus")

pop_colours <- c("ASIA" = "hotpink", "AUS" = "cornflowerblue", "CENAM" = "purple", "EUR" = "forestgreen", "USA" = "tomato2")

ggplot(minus_chr, aes(x = time, y = Ne, col = region)) +
  geom_line(linewidth=1.5) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    title = "SMC++ Bootstrap: Effective Population Size by Population",
    x = "Years before present (log scale)",
    y = "Effective population size (log scale)"
  ) +
  theme_minimal() +
  theme(legend.position = "top") + 
  facet_grid(chr ~ generation) +
  scale_colour_manual(values = pop_colours)

  ```



  ### split times

```bash
mkdir SPLIT
cd SPLIT

# Create necessary directories
mkdir DATA
mkdir g1 g2.5 g4
mkdir DATA/ASIA_AUS DATA/ASIA_CENAM DATA/ASIA_EUR DATA/ASIA_USA
mkdir DATA/AUS_CENAM DATA/AUS_EUR DATA/AUS_USA DATA/AUS_ASIA
mkdir DATA/CENAM_EUR DATA/CENAM_USA DATA/CENAM_ASIA DATA/CENAM_AUS
mkdir DATA/EUR_USA DATA/EUR_CENAM DATA/EUR_ASIA DATA/EUR_AUS
mkdir DATA/USA_EUR DATA/USA_CENAM DATA/USA_AUS DATA/USA_ASIA




# Get sample names for each population
ASIA=$(awk -F'_' '$1 == "MYS" || $1 == "THA" {print $0}' ../keep_samples.dimmitis-only.nuclear_snps.no-reps.list | paste -sd ',')
AUS=$(awk -F'_' '$1 == "AUS" {print $0}' ../keep_samples.dimmitis-only.nuclear_snps.no-reps.list | paste -sd ',')
CENAM=$(awk -F'_' '$1 == "PAN" || $1 == "CRI" {print $0}' ../keep_samples.dimmitis-only.nuclear_snps.no-reps.list | paste -sd ',')
EUR=$(awk -F'_' '$1 == "GRC" || $1 == "ITA" || $1 == "ROU" {print $0}' ../keep_samples.dimmitis-only.nuclear_snps.no-reps.list | paste -sd ',')
USA=$(awk -F'_' '$1 == "USA" {print $0}' ../keep_samples.dimmitis-only.nuclear_snps.no-reps.list | paste -sd ',')

declare -A populations
populations=(
  [ASIA]="$ASIA"
  [AUS]="$AUS"
  [CENAM]="$CENAM"
  [EUR]="$EUR"
  [USA]="$USA"
)

module load smcpp/1.15.4--py39hac1eaed_0
module load htslib-1.19/perl-5.38.0

# Function to create dataset and run smc++ commands
process_pair() {
  local pop1=$1
  local pop2=$2
  local out_prefix="./${pop1}_${pop2}"
  local dir_prefix="DATA/${pop1}_${pop2}"
  
  # Generate VCF for the pair
  vcftools --gzvcf ../dirofilaria_global.cohort.2025-06-18.nuclearSNPs..final.autosomal.recode.vcf.gz \
    $(echo ${populations[$pop1]} | tr ',' '\n' | awk '{print "--indv "$0}') \
    $(echo ${populations[$pop2]} | tr ',' '\n' | awk '{print "--indv "$0}') \
    --max-missing 1 --recode --out $out_prefix

  bgzip -f ${out_prefix}.recode.vcf
  tabix ${out_prefix}.recode.vcf.gz

  # create datasets containing the joint frequency spectrum for both pops
  for chr in {1..4}; do
    smc++ vcf2smc --mask ../missing_sites.bed.gz ${out_prefix}.recode.vcf.gz ${dir_prefix}/${pop1}_${pop2}.chr${chr}.smc.gz dirofilaria_immitis_chr${chr} ${pop1}:${populations[$pop1]} ${pop2}:${populations[$pop2]}
    smc++ vcf2smc --mask ../missing_sites.bed.gz ${out_prefix}.recode.vcf.gz ${dir_prefix}/${pop2}_${pop1}.chr${chr}.smc.gz dirofilaria_immitis_chr${chr} ${pop2}:${populations[$pop2]} ${pop1}:${populations[$pop1]}
  done

  # copy smc++ data for individual pops
  cp ../DATA_${pop1}/${pop1}*.smc.gz ${dir_prefix}
  cp ../DATA_${pop2}/${pop2}*.smc.gz ${dir_prefix}

  # run split
  smc++ split --timepoints 1 1000000 -o ${pop1}_${pop2} ../RESULT_${pop1}/jackknife_${pop1}_all_chr/model.final.json ../RESULT_${pop2}/jackknife_${pop2}_all_chr/model.final.json ${dir_prefix}/*.smc.gz
  smc++ plot -g 1 -c g1/SMCPP_${pop1}_${pop2}_t5m_g1.pdf ${pop1}_${pop2}/model.final.json
  smc++ plot -g 2.5 -c g2.5/SMCPP_${pop1}_${pop2}_t5m_g2.5.pdf ${pop1}_${pop2}/model.final.json
  smc++ plot -g 4 -c g4/SMCPP_${pop1}_${pop2}_t5m_g4.pdf ${pop1}_${pop2}/model.final.json
}

# Process each population pair
process_pair ASIA AUS
process_pair ASIA CENAM
process_pair ASIA EUR
process_pair ASIA USA

process_pair AUS CENAM
process_pair AUS EUR
process_pair AUS USA
process_pair AUS ASIA

process_pair CENAM EUR
process_pair CENAM USA
process_pair CENAM ASIA
process_pair CENAM AUS

process_pair EUR USA
process_pair EUR ASIA
process_pair EUR CENAM
process_pair EUR AUS

process_pair USA EUR
process_pair USA ASIA
process_pair USA CENAM
process_pair USA AUS

```




# get split times
```R

get_split_times <- function(generation_time){

require(tidyverse)

ASIA_AUS_data <- read.delim(paste0("g",generation_time,"/SMCPP_ASIA_AUS_t5m_g",generation_time,".csv"), header = TRUE, sep = ",")
ASIA_CENAM_data <- read.delim(paste0("g",generation_time,"/SMCPP_ASIA_CENAM_t5m_g",generation_time,".csv"), header = TRUE, sep = ",")
ASIA_EUR_data <- read.delim(paste0("g",generation_time,"/SMCPP_ASIA_EUR_t5m_g",generation_time,".csv"), header = TRUE, sep = ",")
ASIA_USA_data <- read.delim(paste0("g",generation_time,"/SMCPP_ASIA_USA_t5m_g",generation_time,".csv"), header = TRUE, sep = ",")

AUS_CENAM_data <- read.delim(paste0("g",generation_time,"/SMCPP_AUS_CENAM_t5m_g",generation_time,".csv"), header = TRUE, sep = ",")
AUS_EUR_data <- read.delim(paste0("g",generation_time,"/SMCPP_AUS_EUR_t5m_g",generation_time,".csv"), header = TRUE, sep = ",")
AUS_USA_data <- read.delim(paste0("g",generation_time,"/SMCPP_AUS_USA_t5m_g",generation_time,".csv"), header = TRUE, sep = ",")
AUS_ASIA_data <- read.delim(paste0("g",generation_time,"/SMCPP_AUS_ASIA_t5m_g",generation_time,".csv"), header = TRUE, sep = ",")

CENAM_EUR_data <- read.delim(paste0("g",generation_time,"/SMCPP_CENAM_EUR_t5m_g",generation_time,".csv"), header = TRUE, sep = ",")
CENAM_USA_data <- read.delim(paste0("g",generation_time,"/SMCPP_CENAM_USA_t5m_g",generation_time,".csv"), header = TRUE, sep = ",")
CENAM_ASIA_data <- read.delim(paste0("g",generation_time,"/SMCPP_CENAM_ASIA_t5m_g",generation_time,".csv"), header = TRUE, sep = ",")
CENAM_AUS_data <- read.delim(paste0("g",generation_time,"/SMCPP_CENAM_AUS_t5m_g",generation_time,".csv"), header = TRUE, sep = ",")

EUR_USA_data <- read.delim(paste0("g",generation_time,"/SMCPP_EUR_USA_t5m_g",generation_time,".csv"), header = TRUE, sep = ",")
EUR_ASIA_data <- read.delim(paste0("g",generation_time,"/SMCPP_EUR_ASIA_t5m_g",generation_time,".csv"), header = TRUE, sep = ",")
EUR_CENAM_data <- read.delim(paste0("g",generation_time,"/SMCPP_EUR_CENAM_t5m_g",generation_time,".csv"), header = TRUE, sep = ",")
EUR_AUS_data <- read.delim(paste0("g",generation_time,"/SMCPP_EUR_AUS_t5m_g",generation_time,".csv"), header = TRUE, sep = ",")

USA_AUS_data <- read.delim(paste0("g",generation_time,"/SMCPP_USA_AUS_t5m_g",generation_time,".csv"), header = TRUE, sep = ",")
USA_ASIA_data <- read.delim(paste0("g",generation_time,"/SMCPP_USA_ASIA_t5m_g",generation_time,".csv"), header = TRUE, sep = ",")
USA_EUR_data <- read.delim(paste0("g",generation_time,"/SMCPP_USA_EUR_t5m_g",generation_time,".csv"), header = TRUE, sep = ",")
USA_CENAM_data <- read.delim(paste0("g",generation_time,"/SMCPP_USA_CENAM_t5m_g",generation_time,".csv"), header = TRUE, sep = ",")

# Add pair identifiers and subset IDs
datasets <- list(
  ASIA_AUS = ASIA_AUS_data,
  ASIA_CENAM = ASIA_CENAM_data,
  ASIA_EUR = ASIA_EUR_data,
  ASIA_USA = ASIA_USA_data,
  AUS_CENAM = AUS_CENAM_data,
  AUS_EUR = AUS_EUR_data,
  AUS_USA = AUS_USA_data,
  AUS_ASIA = AUS_ASIA_data,
  CENAM_EUR = CENAM_EUR_data,
  CENAM_USA = CENAM_USA_data,
  CENAM_ASIA = CENAM_ASIA_data,
  CENAM_AUS = CENAM_AUS_data,
  EUR_USA = EUR_USA_data,
  EUR_ASIA = EUR_ASIA_data,
  EUR_CENAM = EUR_CENAM_data,
  EUR_AUS = EUR_AUS_data,
  USA_AUS = USA_AUS_data,
  USA_EUR = USA_EUR_data,
  USA_CENAM = USA_CENAM_data,
  USA_ASIA = USA_ASIA_data
)

for (name in names(datasets)) {
  datasets[[name]] <- datasets[[name]] %>%
    mutate(Pair = name, subset_id = name)
}

# Combine all datasets
all_data <- bind_rows(datasets)

get_second_pop_end_point <- function(data) {
  switch_point <- which(data$label != data$label[1])[1]
  second_pop_data <- data[switch_point:nrow(data), ]
  end_x <- max(second_pop_data$x, na.rm = TRUE)
  return(end_x)
}

# Apply the function to each dataset and store the results
end_points <- data.frame(
  subset_id = names(datasets),
  end_x = sapply(datasets, get_second_pop_end_point)
)

print(end_points)

write.csv(end_points, paste0("split_dates_g",generation_time,".csv"), row.names=FALSE)

}

get_split_times("1")
get_split_times("2.5")
get_split_times("4")

```



### testing different splines in smc++ estimate
```bash
cd /lustre/scratch127/pam/users/sd21/dirofilaria_immitis/project_global-diversity/04_analysis/05_smc++

# Get sample names for each population
ASIA=$(awk -F'_' '$1 == "MYS" || $1 == "THA" {print $0}' keep_samples.dimmitis-only.nuclear_snps.no-reps.list | paste -sd ',')
AUS=$(awk -F'_' '$1 == "AUS" {print $0}' keep_samples.dimmitis-only.nuclear_snps.no-reps.list | paste -sd ',')
CENAM=$(awk -F'_' '$1 == "PAN" || $1 == "CRI" {print $0}' keep_samples.dimmitis-only.nuclear_snps.no-reps.list | paste -sd ',')
EUR=$(awk -F'_' '$1 == "GRC" || $1 == "ITA" || $1 == "ROU" {print $0}' keep_samples.dimmitis-only.nuclear_snps.no-reps.list | paste -sd ',')
USA=$(awk -F'_' '$1 == "USA" {print $0}' keep_samples.dimmitis-only.nuclear_snps.no-reps.list | paste -sd ',')

declare -A populations
populations=(
  [ASIA]="$ASIA"
  [AUS]="$AUS"
  [CENAM]="$CENAM"
  [EUR]="$EUR"
  [USA]="$USA"
)

mkdir DATA_SPLINES
mkdir RESULTS_SPLINES

# 
for pop in "${!populations[@]}"; do
    samples="${populations[$pop]}"
    vcf_file="${pop}.vcf.gz"  # e.g. ASIA.vcf.gz

    # Step 1: Convert VCF to SMC
    for chr in {1..4}; do
        smc++ vcf2smc \
            $vcf_file \
            --mask missing_sites.bed.gz \
            DATA_SPLINES/${pop}.chr${chr}.smc.gz \
            "dirofilaria_immitis_chr${chr}" \
            "${pop}:${samples}"
    done

    # Step 2: Estimate demographic history
    smc++ estimate --spline piecewise --timepoints 1 500000 -o "RESULTS_SPLINES/${pop}_PIECEWISE/" 2.7e-9 \
        DATA_SPLINES/${pop}.chr*.smc.gz

    smc++ estimate --spline cubic --timepoints 1 500000 -o "RESULTS_SPLINES/${pop}_CUBIC/" 2.7e-9 \
        DATA_SPLINES/${pop}.chr*.smc.gz

     smc++ estimate --spline pchip --timepoints 1 500000 -o "RESULTS_SPLINES/${pop}_PCHIP/" 2.7e-9 \
        DATA_SPLINES/${pop}.chr*.smc.gz


    # Step 3: Plot results with different generation times
    for g in 1 2.5 4; do
        smc++ plot -g "$g" -c "SMCPP_${pop}_g${g}_piecewise.pdf" "RESULTS_SPLINES/${pop}_PIECEWISE/model.final.json"
        smc++ plot -g "$g" -c "SMCPP_${pop}_g${g}_cubic.pdf" "RESULTS_SPLINES/${pop}_CUBIC/model.final.json"
        smc++ plot -g "$g" -c "SMCPP_${pop}_g${g}_pchip.pdf" "RESULTS_SPLINES/${pop}_PCHIP/model.final.json"
    done
done
```
