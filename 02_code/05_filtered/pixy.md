# Dirofilaria immitis WGS Lab Book - Pixy

Ideas for running pixy:
- Populations = region (AUS, ASIA, USA, EUR, CENAM) (code below - adopted from Javi's paper)
- Populations = host (dog, cat, fox, ferret etc) (will adapt below code for this if everything looks ok)
- Populations = city (all samples worldwide) or perhaps just focus on specific countries, e.g. run pixy for Australia only and USA only.                                            

## Generate an ALL SITES variant set for running pixy

### Create cohort gvcf

```bash
module load bsub.py/0.42.1
module load gatk/4.1.4.1

cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/04_VARIANTS/ALLSITES
                   
# collect all sample g.vcfs (from all batches) into a list, to make input for CombineGVCFs
find /lustre/scratch125/pam/teams/team333/sd21/dirofilaria_immitis/POPGEN/NEWDATA_2024/VARIANTS/gatk_hc_DIMMITIS_POPGEN/GATK_HC_GVCFs -type f -name "*.g.vcf.gz" > /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/04_VARIANTS/ALLSITES/gvcf.list

args=$(xargs -I {} echo -n "-V {} " < gvcf.list)

# Remove trailing space
args=$(echo "$args" | sed 's/ *$//')

echo $args

SPLIT="split1.list split2.list split3.list split4.list split5.list split6.list"

# merge gvcfs (split per chromosome)
for file in $SPLIT; do
 base=$(basename $file .txt)
  bsub.py --queue basement 50 gatk_combinegvcfs_${base} \
  "gatk CombineGVCFs -R /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/01_REF/dimmitis_WSI_2.2.fa --intervals ${file} ${args} -O /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/04_VARIANTS/ALLSITES/DIMMITIS_POPGEN_ALLSITES_${base}.g.vcf.gz" ;
done
## successfully completed

# try running another job with threads to be faster
for file in $SPLIT; do
 base=$(basename $file .txt)
  bsub.py --queue basement --threads 8 20 gatk_combinegvcfs_${base}_NEW \
  "gatk CombineGVCFs -R /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/01_REF/dimmitis_WSI_2.2.fa --intervals ${file} ${args} -O /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/04_VARIANTS/ALLSITES_NEW/DIMMITIS_POPGEN_ALLSITES_${base}.g.vcf.gz" ;
done
## still took a while, disregard this
```

### Run genotyping on each chromosome separately

```bash
module load PaM/environment
module load bsub.py/0.42.1
module load gatk/4.1.4.1

for file in $SPLIT; do
 base=$(basename $file .txt)
  bsub.py --queue long 20 gatk_genotypevcfs_${base} \
  "gatk GenotypeGVCFs -R /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/01_REF/dimmitis_WSI_2.2.fa -V /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/04_VARIANTS/ALLSITES/DIMMITIS_POPGEN_ALLSITES_${base}.g.vcf.gz --intervals ${file} --all-sites -O DIMMITIS_POPGEN_ALLSITES_${base}.vcf.gz -G StandardAnnotation -G AS_StandardAnnotation" ;
done
# Successfully completed
```

### Bring vcf files together

```bash
module load htslib-1.19/perl-5.38.0 # this is a different version to what I've previously used - check what version Steve used for variant calling and how I can access that with the new farm layout
module load vcftools/0.1.16-c4

cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/04_VARIANTS/ALLSITES

# make list of vcfs
ls -1 *.list.vcf.gz | sort -n > vcf_files.list
# .list.vcf ensures that the g.vcfs aren't included

# merge vcfs
vcf-concat --files vcf_files.list > DIMMITIS_POPGEN_ALLSITES.vcf;
  bgzip DIMMITIS_POPGEN_ALLSITES.vcf;
  tabix -p vcf DIMMITIS_POPGEN_ALLSITES.vcf.gz
```

### Filter nuclear variants and invariants

```bash
VCF=/lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/04_VARIANTS/ALLSITES/DIMMITIS_POPGEN_ALLSITES.vcf.gz
MIT_CONTIG=dirofilaria_immitis_chrMtDNA
WB_CONTIG=dirofilaria_immitis_chrWb
REFERENCE=/lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/01_REF/dimmitis_WSI_2.2.fa

vcftools --gzvcf DIMMITIS_POPGEN_ALLSITES.vcf.gz --remove-indels
#After filtering, kept 143 out of 143 Individuals
#After filtering, kept 87117573 out of a possible 88130674 Sites

#select nuclear invariants
bsub.py 1 select_nuclearINVARIANTs "gatk SelectVariants \
--reference ${REFERENCE} \
--variant ${VCF} \
--select-type-to-include NO_VARIATION \
--exclude-intervals ${MIT_CONTIG} \
--exclude-intervals ${WB_CONTIG} \
--output ${VCF%.vcf.gz}.nuclearINVARIANTs.vcf"
# successfully completed

vcftools --vcf DIMMITIS_POPGEN_ALLSITES.nuclearINVARIANTs.vcf --remove-indels
#After filtering, kept 143 out of 143 Individuals
#After filtering, kept 84245264 out of a possible 84322763 Sites
```

### Merge these nuclear invariants with the nuclear SNPs I previously filtered

```bash
#Keep n=128 samples in the nuclear dataset (nuclear_samplelist.keep) and merge these with the SNPs
bsub.py 1 filter_nuclear_INVARIANTvcftools "vcftools --vcf ${VCF%.vcf.gz}.nuclearINVARIANTs.vcf \
--keep /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/04_VARIANTS/FILTER1/NO_OUTGROUPS/nuclear_samplelist.keep \
--recode --out nuclearINVARIANTS_128samples"
#After filtering, kept 128 out of 143 Individuals
#After filtering, kept 84322763 out of a possible 84322763 Sites

# Merge with the already filtered variants
bsub.py --done "filter_nuclear_INVARIANTvcftools" 1 merge_nuclear_VARIANTsandINVARIANTs "gatk MergeVcfs \
--INPUT /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/04_VARIANTS/FILTER1/NO_OUTGROUPS/FINAL_SETS/nuclear_samples3x_missing0.9.recode.vcf \
--INPUT nuclearINVARIANTS_128samples.recode.vcf \
--OUTPUT nuclearVARIANTsandINVARIANTs_128samples.recode.vcf"

# Select only the chrX to ch4 (avoid the scaffolds) and remove replicate samples
bsub.py 1 filter_CHR "vcftools \
--vcf nuclearVARIANTsandINVARIANTs_128samples.recode.vcf \
--chr dirofilaria_immitis_chrX \
--chr dirofilaria_immitis_chr1 \
--chr dirofilaria_immitis_chr2 \
--chr dirofilaria_immitis_chr3 \
--chr dirofilaria_immitis_chr4 \
--keep nuclear_samplelist_norep.keep \
--remove-indels \
--recode --out nuclearSNPssandINVARIANTs_124samples.chrXto4"
#After filtering, kept 124 out of 128 Individuals
#After filtering, kept 84210689 out of a possible 84587321 Sites


bgzip nuclearSNPssandINVARIANTs_124samples.chrXto4.recode.vcf;
tabix -p vcf nuclearSNPssandINVARIANTs_124samples.chrXto4.recode.vcf.gz
mv nuclearSNPssandINVARIANTs_124samples.chrXto4.recode.vcf.* ../FILTER1/NO_OUTGROUPS/FINAL_SETS/
```

### Run pixy on nuclear variants & invariants

```bash
# Create pixy environment
module load conda/24.7.1-2
conda create --name pixy
conda init --all
conda activate pixy

# Install pixy & htslib
conda install -c conda-forge pixy
conda install -c bioconda htslib

# Environmental variables
WORKING_DIR=/lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA
VCF=${WORKING_DIR}/03_ANALYSIS/04_VARIANTS/FILTER1/NO_OUTGROUPS/FINAL_SETS/nuclearSNPssandINVARIANTs_124samples.chrXto4.recode.vcf.gz

cd ${WORKING_DIR}/03_ANALYSIS/05_ANALYSIS/PIXY

# Population files per region and host
# region_pop.list
'
AUS_BNE_AD_001	AUS
AUS_BNE_AD_002	AUS
AUS_BNE_AD_003	AUS
AUS_BNE_AD_004	AUS
AUS_BNE_AD_006	AUS
AUS_BNE_AD_008	AUS
AUS_BNE_AD_009	AUS
AUS_CNS_AD_001	AUS
AUS_CNS_AD_002	AUS
AUS_LHR_AD_001	AUS
AUS_ROK_AD_001	AUS
AUS_SYD_AD_001	AUS
AUS_SYD_AD_002	AUS
AUS_SYD_AD_003	AUS
AUS_SYD_AD_004	AUS
AUS_SYD_AD_005	AUS
AUS_SYD_AD_006	AUS
AUS_SYD_AD_007	AUS
AUS_SYD_AD_008	AUS
AUS_SYD_AD_009	AUS
AUS_SYD_AD_010	AUS
AUS_SYD_AD_011	AUS
AUS_SYD_AD_012	AUS
AUS_SYD_AD_013	AUS
AUS_SYD_AD_014	AUS
AUS_SYD_AD_015	AUS
AUS_SYD_AD_017	AUS
AUS_TVS_AD_001	AUS
AUS_TVS_AD_002	AUS
AUS_TVS_AD_003	AUS
AUS_TVS_AD_004	AUS
AUS_TVS_AD_005	AUS
AUS_TVS_AD_006	AUS
AUS_TVS_AD_007	AUS
AUS_TVS_AD_008	AUS
AUS_TVS_AD_010	AUS
AUS_TVS_AD_011	AUS
AUS_TVS_AD_012	AUS
AUS_TVS_AD_013	AUS
AUS_TVS_AD_014	AUS
AUS_TVS_AD_015	AUS
AUS_TVS_AD_016	AUS
AUS_TVS_AD_017	AUS
AUS_TVS_AD_018	AUS
AUS_TVS_AD_019	AUS
CRI_SJO_AD_001	CENAM
GRC_XAN_AD_001	EUR
GRC_XAN_AD_002	EUR
GRC_XAN_AD_003	EUR
GRC_XAN_AD_004	EUR
GRC_XAN_AD_005	EUR
GRC_XAN_AD_006	EUR
GRC_XAN_AD_007	EUR
GRC_XAN_AD_008	EUR
GRC_XAN_AD_009	EUR
GRC_XAN_AD_010	EUR
GRC_XAN_AD_011	EUR
GRC_XAN_AD_012	EUR
ITA_NEA_AD_001	EUR
ITA_NEA_AD_002	EUR
ITA_NEA_AD_003	EUR
ITA_PAV_AD_001	EUR
MYS_SEL_AD_001	ASIA
PAN_BOC_AD_001	CENAM
PAN_BOC_AD_002	CENAM
PAN_BOC_AD_003	CENAM
PAN_BOC_AD_004	CENAM
PAN_BOC_AD_005	CENAM
PAN_PUE_AD_001	CENAM
PAN_PUE_AD_002	CENAM
PAN_PUE_AD_003	CENAM
PAN_PUE_AD_004	CENAM
PAN_PUE_AD_005	CENAM
PAN_PUE_AD_006	CENAM
PAN_SLO_AD_001	CENAM
PAN_SLO_AD_002	CENAM
PAN_SLO_AD_003	CENAM
ROU_BUC_AD_001	EUR
ROU_COM_AD_001	EUR
ROU_GIU_AD_001	EUR
THA_BKK_AD_001	ASIA
THA_BKK_AD_002	ASIA
THA_BKK_AD_003	ASIA
THA_BKK_AD_004	ASIA
THA_BKK_AD_005	ASIA
THA_BKK_AD_006	ASIA
USA_FLO_AD_001	USA
USA_FLO_AD_002	USA
USA_FLO_AD_003	USA
USA_FLO_AD_004	USA
USA_FLO_AD_005	USA
USA_FLO_AD_006	USA
USA_FLO_AD_007	USA
USA_FLO_AD_008	USA
USA_FLO_AD_009	USA
USA_FLO_AD_010	USA
USA_FLO_AD_011	USA
USA_FLO_AD_012	USA
USA_FLO_AD_013	USA
USA_FLO_AD_014	USA
USA_FLO_AD_015	USA
USA_FLO_AD_016	USA
USA_GEO_AD_001	USA
USA_ILL_AD_001	USA
USA_ILL_AD_002	USA
USA_LOU_AD_001	USA
USA_LOU_AD_002	USA
USA_MIS_AD_001	USA
USA_MIS_AD_002	USA
USA_TEX_AD_001	USA
USA_TEX_AD_002	USA
USA_TEX_AD_003	USA
USA_TEX_AD_004	USA
USA_TEX_AD_005	USA
USA_TEX_AD_006	USA
USA_TEX_AD_007	USA
USA_TEX_AD_008	USA
USA_TEX_AD_009	USA
USA_TEX_AD_010	USA
USA_TEX_AD_011	USA
USA_TEX_AD_012	USA
USA_TEX_AD_013	USA
USA_TEX_AD_014	USA
USA_TEX_AD_015	USA
'

# host_pop.list
'
AUS_BNE_AD_001	Dog
AUS_BNE_AD_002	Dog
AUS_BNE_AD_003	Dog
AUS_BNE_AD_004	Dog
AUS_BNE_AD_006	Fox
AUS_BNE_AD_008	Dog
AUS_BNE_AD_009	Dog
AUS_CNS_AD_001	Dog
AUS_CNS_AD_002	Dog
AUS_LHR_AD_001	Dog
AUS_ROK_AD_001	Dog
AUS_SYD_AD_001	Dog
AUS_SYD_AD_002	Dog
AUS_SYD_AD_003	Dog
AUS_SYD_AD_004	Dog
AUS_SYD_AD_005	Dog
AUS_SYD_AD_006	Dog
AUS_SYD_AD_007	Dog
AUS_SYD_AD_008	Dog
AUS_SYD_AD_009	Dog
AUS_SYD_AD_010	Dog
AUS_SYD_AD_011	Dog
AUS_SYD_AD_012	Dog
AUS_SYD_AD_013	Dog
AUS_SYD_AD_014	Dog
AUS_SYD_AD_015	Dog
AUS_SYD_AD_017	Dog
AUS_TVS_AD_001	Dog
AUS_TVS_AD_002	Dog
AUS_TVS_AD_003	Dog
AUS_TVS_AD_004	Dog
AUS_TVS_AD_005	Dog
AUS_TVS_AD_006	Dog
AUS_TVS_AD_007	Dog
AUS_TVS_AD_008	Dog
AUS_TVS_AD_010	Dog
AUS_TVS_AD_011	Dog
AUS_TVS_AD_012	Dog
AUS_TVS_AD_013	Dog
AUS_TVS_AD_014	Dog
AUS_TVS_AD_015	Dog
AUS_TVS_AD_016	Dog
AUS_TVS_AD_017	Dog
AUS_TVS_AD_018	Dog
AUS_TVS_AD_019	Dog
CRI_SJO_AD_001	Cat
GRC_XAN_AD_001	Dog
GRC_XAN_AD_002	Dog
GRC_XAN_AD_003	Dog
GRC_XAN_AD_004	Dog
GRC_XAN_AD_005	Dog
GRC_XAN_AD_006	Dog
GRC_XAN_AD_007	Dog
GRC_XAN_AD_008	Dog
GRC_XAN_AD_009	Dog
GRC_XAN_AD_010	Dog
GRC_XAN_AD_011	Dog
GRC_XAN_AD_012	Dog
ITA_NEA_AD_001	Dog
ITA_NEA_AD_002	Fox
ITA_NEA_AD_003	Fox
ITA_PAV_AD_001	Dog
MYS_SEL_AD_001	Dog
PAN_BOC_AD_001	Dog
PAN_BOC_AD_002	Dog
PAN_BOC_AD_003	Dog
PAN_BOC_AD_004	Dog
PAN_BOC_AD_005	Dog
PAN_PUE_AD_001	Dog
PAN_PUE_AD_002	Dog
PAN_PUE_AD_003	Dog
PAN_PUE_AD_004	Dog
PAN_PUE_AD_005	Dog
PAN_PUE_AD_006	Dog
PAN_SLO_AD_001	Dog
PAN_SLO_AD_002	Dog
PAN_SLO_AD_003	Dog
ROU_BUC_AD_001	Leopard
ROU_COM_AD_001	Wildcat
ROU_GIU_AD_001	Golden_jackal
THA_BKK_AD_001	Dog
THA_BKK_AD_002	Dog
THA_BKK_AD_003	Dog
THA_BKK_AD_004	Dog
THA_BKK_AD_005	Dog
THA_BKK_AD_006	Dog
USA_FLO_AD_001	Dog
USA_FLO_AD_002	Dog
USA_FLO_AD_003	Dog
USA_FLO_AD_004	Dog
USA_FLO_AD_005	Dog
USA_FLO_AD_006	Dog
USA_FLO_AD_007	Dog
USA_FLO_AD_008	Dog
USA_FLO_AD_009	Dog
USA_FLO_AD_010	Dog
USA_FLO_AD_011	Dog
USA_FLO_AD_012	Dog
USA_FLO_AD_013	Dog
USA_FLO_AD_014	Dog
USA_FLO_AD_015	Dog
USA_FLO_AD_016	Dog
USA_GEO_AD_001	Dog
USA_ILL_AD_001	Dog
USA_ILL_AD_002	Dog
USA_LOU_AD_001	Dog
USA_LOU_AD_002	Dog
USA_MIS_AD_001	Dog
USA_MIS_AD_002	Dog
USA_TEX_AD_001	Ferret
USA_TEX_AD_002	Dog
USA_TEX_AD_003	Dog
USA_TEX_AD_004	Dog
USA_TEX_AD_005	Dog
USA_TEX_AD_006	Dog
USA_TEX_AD_007	Dog
USA_TEX_AD_008	Dog
USA_TEX_AD_009	Dog
USA_TEX_AD_010	Dog
USA_TEX_AD_011	Cat
USA_TEX_AD_012	Cat
USA_TEX_AD_013	Dog
USA_TEX_AD_014	Dog
USA_TEX_AD_015	Cat
'


# Run pixy per region
bsub.py --queue long --threads 10 40 pixy_region \
"pixy --stats pi fst dxy \
--vcf ${VCF} \
--populations region_pop.list \
--window_size 100000 \
--n_cores 10 \
--output_prefix region"
# Successfully completed

# Run pixy per host
bsub.py --queue long --threads 10 40 pixy_host \
"pixy --stats pi fst dxy \
--vcf ${VCF} \
--populations host_pop.list \
--window_size 100000 \
--n_cores 10 \
--output_prefix host"
# Successfully completed
## Keep in mind that the sample sizes for non-dog hosts are quite low, but will still be interesting to explore
```

Download the file and use it as input in R to estimate the values of pi, fst and dxy and generate plots.

### Pixy by region

```R
# Batch 4 - Pixy for Pi, Fst and Dxy by REGION

library(tidyverse)
library(ggsci)
library(ggpubr)
library(patchwork)
library(ggridges)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(purrr)
library(stats)
library(graphics)
library(grDevices)
library(utils)
library(datasets)
library(methods)
library(base)
require(maps)
library(mapdata)
library(readxl)
library(ozmaps) 
library(grid)
library(gridExtra)
library(ggrepel)
library(ggnewscale)
library(reshape)
library(wesanderson)

setwd("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/batch4/FILTER1/NO_OUTGROUPS/pixy/region")

##### Nucleotide diversity (Pi) #####

# get nucleotide diversity (pi) data from pixy output
pi_data <- read.table("input/region_pi.txt", header=T)
pi_data$chromosome <- str_remove(pi_data$chromosome, 'dirofilaria_immitis_')

# subset
pi_data_AUS <- pi_data %>%
  filter(pop=="AUS")
pi_data_ASIA <- pi_data %>%
  filter(pop=="ASIA")
pi_data_USA <- pi_data %>%
  filter(pop=="USA")
pi_data_EUR <- pi_data %>%
  filter(pop=="EUR")
pi_data_CENAM <- pi_data %>%
  filter(pop=="CENAM")

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

# calculate the median Pi and check the ratio of sex-to-autosome diversity. Should be about 0.75, as D. immitis is XX/XY
pi_data_sex_median <- pi_data %>%
  group_by(chr_type) %>%
  summarise(median = median(avg_pi, na.rm = TRUE))
pi_data_sex_median

'
# A tibble: 2 × 2
  chr_type   median
  <chr>       <dbl>
1 autosome 0.000467
2 sexchr   0.000181
'
# 0.000181 / 0.000467 = 0.3875803 (way off 0.75)

pi_data_pop_sex_median <- pi_data %>%
  group_by(pop, chr_type) %>%
  summarise(median = median(avg_pi, na.rm = TRUE))
pi_data_pop_sex_median

'
# A tibble: 10 × 3
# Groups:   pop [5]
   pop   chr_type   median
   <chr> <chr>       <dbl>
 1 ASIA  autosome 0.000289
 2 ASIA  sexchr   0.000109
 3 AUS   autosome 0.000563
 4 AUS   sexchr   0.000192
 5 CENAM autosome 0.000372
 6 CENAM sexchr   0.000143
 7 EUR   autosome 0.000448
 8 EUR   sexchr   0.000200
 9 USA   autosome 0.000736
10 USA   sexchr   0.000319
'

# plot 1 - genome wide plots per population
plot_1_pi <- ggplot(pi_data, aes(position*100000, avg_pi, col=chromosome)) +
  geom_point(size=0.8) +
  facet_grid(pop~.) +
  scale_color_npg(guide = guide_legend(override.aes = list(size = 4))) +
  theme_bw() +
  theme(legend.position="top", 
        legend.title = element_blank(),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 18),
        legend.text = element_text(size = 16)) +
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
  theme(legend.position="top",
        legend.title = element_blank(),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 18),
        legend.text = element_text(size = 16)) +
  labs(x="Nucleotide Diversity (Pi)", y="Density")
plot_2_pi

# combine plots
plot_1_pi + plot_2_pi +  plot_layout(widths = c(5, 1))
ggsave("genomewide_and_density_Pi_region.tif", width=14, height=8, dpi = 300)
ggsave("genomewide_and_density_Pi_region.png", width=14, height=8, dpi = 300)

#Now a boxplot of the pi value per population

red_palette1 <- brewer.pal(n = 9, name = "YlOrRd")
red_palette2 <- brewer.pal(n=9, name = "YlOrBr")
blue_palette <- brewer.pal(n = 9, name = "Blues")
green_palette1 <- brewer.pal(n = 9, name = "BuGn")
green_palette2 <- brewer.pal(n=9, name = "YlGn")

scale_colour_region <- function(...){
  ggplot2:::manual_scale(
    'colour', 
    values = setNames(
      c(blue_palette[5],
        "violetred1",
        red_palette1[7], 
        "blueviolet",
        green_palette1[7]),
    c("AUS", "ASIA", "USA", "EUR", "CENAM")), 
  ...
)
}    

boxplot_pi <- ggplot(pi_data, aes(pop, avg_pi, col=pop)) +
  geom_jitter(size = 1, alpha = 0.5) +
  geom_boxplot(fill=NA, col="black") +
  labs(x = "Population" , y = "Nucleotide diversity (Pi)", colour = "Population") +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 12)
        ) +
  scale_colour_region ()+
  ylim(0, 0.003)

boxplot_pi

ggsave("boxplot_Pi_region.tif", width=6, height=6, dpi = 300)
ggsave("boxplot_Pi_region.png", width=6, height=6, dpi = 300)
# higher pi value = more diverse population
# lower pi value = less diverse population
## Looks like USA and AUS have higher diversity. Asia has the lowest diversity but this is probably because we have a small sample size all from the same country (Thailand).


#let's check the statistical significance of the differences between pi values
#Let's explore the normality
pi_data %>% filter(chromosome != 'chrX') %>%
  ggplot(., aes(x = avg_pi, colour = pop)) +
  geom_histogram() +
  theme_bw() +
  facet_grid(. ~ pop)

#Confirm the absence of normality by shapiro
## If p > 0.05 the data is normal, if p < 0.05 the data is NOT normal
pi_data_AUS <- pi_data_AUS %>% filter(chromosome != 'chrX')
AUS_shapiro <- shapiro.test(pi_data_AUS$avg_pi)
print(AUS_shapiro)
# W = 0.91924, p-value < 2.2e-16

pi_data_ASIA <- pi_data_ASIA %>% filter(chromosome != 'chrX')
ASIA_shapiro <- shapiro.test(pi_data_ASIA$avg_pi)
print(ASIA_shapiro)
# W = 0.87735, p-value < 2.2e-16

pi_data_USA <- pi_data_USA %>% filter(chromosome != 'chrX')
USA_shapiro <- shapiro.test(pi_data_USA$avg_pi)
print(USA_shapiro)
# W = 0.91893, p-value < 2.2e-16

pi_data_EUR <- pi_data_EUR %>% filter(chromosome != 'chrX')
EUR_shapiro <- shapiro.test(pi_data_EUR$avg_pi)
print(EUR_shapiro)
# W = 0.91477, p-value < 2.2e-16

pi_data_CENAM <- pi_data_CENAM %>% filter(chromosome != 'chrX')
CENAM_shapiro <- shapiro.test(pi_data_CENAM$avg_pi)
print(CENAM_shapiro)
# W = 0.89779, p-value < 2.2e-16

#Wilcoxson test for all pairs of pops
## This tests whether the mean values of 2 dependent groups differ significantly from each other
## If p < 0.05 they are sig different
wilcox.test(pi_data_AUS$avg_pi, pi_data_ASIA$avg_pi)
# W = 252617, p-value < 2.2e-16
wilcox.test(pi_data_AUS$avg_pi, pi_data_USA$avg_pi)
# W = 144348, p-value = 1.001e-09
wilcox.test(pi_data_AUS$avg_pi, pi_data_EUR$avg_pi)
# W = 212831, p-value = 1.581e-07
wilcox.test(pi_data_AUS$avg_pi, pi_data_CENAM$avg_pi)
# W = 230641, p-value = 2.499e-16
wilcox.test(pi_data_ASIA$avg_pi, pi_data_USA$avg_pi)
# W = 84717, p-value < 2.2e-16
wilcox.test(pi_data_ASIA$avg_pi, pi_data_EUR$avg_pi)
# W = 138529, p-value = 1.507e-12
wilcox.test(pi_data_ASIA$avg_pi, pi_data_CENAM$avg_pi)
# W = 154134, p-value = 7.226e-06
wilcox.test(pi_data_USA$avg_pi, pi_data_EUR$avg_pi)
# W = 244966, p-value < 2.2e-16
wilcox.test(pi_data_USA$avg_pi, pi_data_CENAM$avg_pi)
# W = 260336, p-value < 2.2e-16
wilcox.test(pi_data_EUR$avg_pi, pi_data_CENAM$avg_pi)
# W = 198948, p-value = 0.003265

## They're all sig different


#let's generate some dataframes for the estatistics
# subset
pi_data_AUS <- pi_data %>%
  filter(pop=="AUS") %>%
  filter(chromosome != 'chrX')

pi_data_ASIA <- pi_data %>%
  filter(pop=="ASIA") %>%
  filter(chromosome != 'chrX') 

pi_data_USA <- pi_data %>%
  filter(pop=="USA") %>%
  filter(chromosome != 'chrX')

pi_data_EUR <- pi_data %>%
  filter(pop=="EUR") %>%
  filter(chromosome != 'chrX')

pi_data_CENAM <- pi_data %>%
  filter(pop=="CENAM") %>%
  filter(chromosome != 'chrX')




#### Dxy and Fst ####

# load data
dxy_data <- read.table("input/region_dxy.txt", header=T)
fst_data <- read.table("input/region_fst.txt", header=T)

# add some columns
dxy_data$data_type <- "Dxy"
dxy_data <- mutate(dxy_data,
                   comparison = paste(pop1, pop2, sep = '_v_'))
fst_data$data_type <- "Fst"
fst_data <- mutate(fst_data,
                   comparison = paste(pop1, pop2, sep = '_v_'))

# subset and merge dataframes
dxy_data_sub <- dxy_data %>% select(comparison,data_type,chromosome,window_pos_1,window_pos_2,avg_dxy)
colnames(dxy_data_sub) <- c("comparison","data_type","chromosome","window_pos_1","window_pos_2","value")

fst_data_sub <- fst_data %>% select(comparison,data_type,chromosome,window_pos_1,window_pos_2,avg_wc_fst)
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
# A tibble: 20 × 3
# Groups:   comparison [10]
   comparison   data_type   median
   <chr>        <chr>        <dbl>
 1 ASIA_v_USA   Dxy       0.000599
 2 ASIA_v_USA   Fst       0.185   
 3 AUS_v_ASIA   Dxy       0.000512
 4 AUS_v_ASIA   Fst       0.221   
 5 AUS_v_CENAM  Dxy       0.000570
 6 AUS_v_CENAM  Fst       0.287   
 7 AUS_v_EUR    Dxy       0.000573
 8 AUS_v_EUR    Fst       0.263   
 9 AUS_v_USA    Dxy       0.000629
10 AUS_v_USA    Fst       0.163   
11 CENAM_v_ASIA Dxy       0.000607
12 CENAM_v_ASIA Fst       0.467   
13 CENAM_v_EUR  Dxy       0.000492
14 CENAM_v_EUR  Fst       0.242   
15 CENAM_v_USA  Dxy       0.000628
16 CENAM_v_USA  Fst       0.209   
17 EUR_v_ASIA   Dxy       0.000618
18 EUR_v_ASIA   Fst       0.399   
19 EUR_v_USA    Dxy       0.000632
20 EUR_v_USA    Fst       0.183   
'


#Plotting Dxy

# plot 1 - genome wide plots per comparison
plot_1_dxy <- data %>%
  filter(data_type =='Dxy')%>%
  ggplot(., aes(position*100000, value, col=chromosome)) +
  geom_point(size=1.5) +
  facet_grid(comparison~.) +
  scale_color_npg(guide = guide_legend(override.aes = list(size = 4))) +
  theme_bw() +
  theme(legend.position="top", 
        legend.title = element_blank(),
        axis.title = element_text(size = 32),
        axis.text = element_text(size = 16),
        strip.text = element_text(size = 14),
        legend.text = element_text(size = 20)) +
  labs(x="Genomic Position", y="Dxy")
plot_1_dxy

# plot 2 - density plots of dxy per group
plot_2_dxy <- data %>%
  filter(data_type=='Dxy')%>%
  ggplot(aes(value, chr_type, fill=chr_type), guide="none") +
  geom_density_ridges(quantile_lines=TRUE, quantile_fun=function(x,...)median(x), size=0.5) +
  theme_bw() + theme(legend.position = "top", 
                     legend.title = element_blank(),
                     axis.text.y=element_blank(),
                     axis.title = element_text(size = 32),
                     axis.text = element_text(size = 16),
                     strip.text = element_text(size = 14),
                     legend.text = element_text(size = 20)) +
  facet_grid(comparison~.) +
  scale_fill_npg() +
  labs(x="Dxy", y="Density")
plot_2_dxy

# combine plots
dxy_plot <- plot_1_dxy + plot_2_dxy +  plot_layout(widths = c(5, 1))
ggsave("genomewide_and_density_dxy_region.tif", width=16, height=18, dpi = 300)
ggsave("genomewide_and_density_dxy_region.png", width=16, height=18, dpi = 300)
# higher Dxy = more divergent

#some additional plots
boxplot_dxy <- data %>%
  filter(data_type =='Dxy') %>% 
  ggplot(., aes(comparison, value, col=comparison)) +
  geom_jitter(size = 1, alpha = 0.5) +
  geom_boxplot(fill=NA, col="black") +
  labs(x = "Population" , y = "Dxy", colour = "Population") +
  theme_bw() +
  scale_color_npg() +
  theme(legend.position = "none",
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 18),
        legend.text = element_text(size = 16)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
boxplot_dxy
ggsave("boxplot_dxy_region.tif", width=14, height = 6, dpi = 300)
ggsave("boxplot_dxy_region.png", width=14, height = 6, dpi = 300)

density_dxy <- data %>%
  filter(data_type =='Dxy') %>% 
  ggplot(., aes(x=value, group = comparison, fill=comparison)) +
  geom_density(adjust=1.5, alpha=.4) +
  labs(x = "Dxy" , y = "Density", colour = "Comparison") +
  theme_bw() +
  scale_fill_npg() +
  theme(legend.position = c(0.88,0.72),
        legend.title = element_blank(),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
density_dxy
ggsave("density_dxy_region.tif", width = 8, height = 6, dpi = 300)
ggsave("density_dxy_region.png", width = 8, height = 6, dpi = 300)

# AUS_v_ASIA & AUS_v_EUR
scatter_dxy_a <-  data %>%
  filter(data_type == 'Dxy') %>%
  as.data.frame() %>%
  select('position', 'comparison', 'value', 'chromosome')%>%
  as_tibble() %>%
  pivot_wider(names_from = comparison, values_from = value) %>%
  ggplot(., aes(x=AUS_v_ASIA, y = AUS_v_EUR, col=chromosome)) +
  geom_point(alpha=.4) +
  labs(x = "AUS_v_ASIA" , y = "AUS_v_EUR", colour = "Comparison") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
scatter_dxy_a
# can repeat for other pairs

genome_pos_dxy_a <-  data %>%
  filter(data_type == 'Dxy') %>%
  as.data.frame() %>%
  select('position', 'comparison', 'value', 'chromosome')%>%
  as_tibble() %>%
  pivot_wider(names_from = comparison, values_from = value) %>%
  mutate("xx" = AUS_v_ASIA - AUS_v_EUR,
         "yy" = AUS_v_ASIA - EUR_v_ASIA) %>%
  ggplot(., aes(position*100000, xx, col=chromosome)) +
  geom_point(size=1) +
  scale_color_tron() +
  theme_bw() +
  theme(legend.position="top", 
        legend.title = element_blank()) +
  labs(x="Position", y="AUS vs ASIA - AUS vs EUR")
genome_pos_dxy_a
# can repeat for other pairs

genome_pos_dxy_aa <-  data %>%
  filter(data_type == 'Dxy') %>%
  as.data.frame() %>%
  select('position', 'comparison', 'value', 'chromosome')%>%
  as_tibble() %>%
  pivot_wider(names_from = comparison, values_from = value) %>%
  mutate("xx" = AUS_v_ASIA / AUS_v_EUR,
         "yy" = AUS_v_ASIA / EUR_v_ASIA) %>%
  ggplot(., aes(position*100000, xx, col=chromosome)) +
  geom_point(size=1) +
  scale_color_tron() +
  theme_bw() +
  theme(legend.position="top", 
        legend.title = element_blank()) +
  labs(x="Position", y="AUS vs ASIA / AUS vs EUR")
genome_pos_dxy_aa
# can repeat for other pairs

# ggarrange above plots if using them


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
  annotate(geom="text", x=140*100000, y=0.004, label="ChrX", size = 6)+
  annotate(geom="text", x=359.5*100000, y=0.004, label="Chr1", size = 6)+
  annotate(geom="text", x=520*100000, y=0.004, label="Chr2", size = 6)+
  annotate(geom="text", x=665*100000, y=0.004, label="Chr3", size = 6)+
  annotate(geom="text", x=810.5*100000, y=0.004, label="Chr4", size = 6)+
  geom_line() +
  scale_color_npg() +
  theme_bw() +
  theme(legend.position="top", 
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 14)) +
  labs(x="Genomic Position", y="Dxy") +
  guides(col = guide_legend(override.aes = list(linewidth = 2)))
dxy_lineplot
ggsave("lineplot_dxy_region.tif", width=14, height=8, dpi = 300)
ggsave("lineplot_dxy_region.png", width=14, height=8, dpi = 300)





# Plotting Fst

# plot 1 - genome wide plots per population
plot_1_fst <- data %>%
  filter(data_type=='Fst')%>%
  ggplot(., aes(position*100000, value, col=chromosome)) +
  geom_point(size=1.5) +
  facet_grid(comparison~.) +
  scale_color_npg(guide = guide_legend(override.aes = list(size = 4))) +
  theme_bw() +
  theme(legend.position="top", 
        legend.title = element_blank(),
        axis.title = element_text(size = 32),
        axis.text = element_text(size = 16),
        strip.text = element_text(size = 14),
        legend.text = element_text(size = 20)) +
  labs(x="Genomic Position", y="Fst")
plot_1_fst


# plot 2 - density plots of Fst per group
plot_2_fst <- data %>%
  filter(data_type=='Fst')%>%
  ggplot(aes(value, chr_type, fill=chr_type), guide="none") +
  geom_density_ridges(quantile_lines=TRUE, quantile_fun=function(x,...)median(x), size=0.5) +
  theme_bw() + 
  theme(legend.position = "none", 
        axis.text.y=element_blank(),
        axis.title = element_text(size = 32),
        axis.text = element_text(size = 16),
        strip.text = element_text(size = 14),
        legend.text = element_text(size = 20)) +
  facet_grid(comparison~.) +
  scale_fill_npg() +
  labs(x="Fst", y="Density")
plot_2_fst

# combine plots
plot_1_fst + plot_2_fst +  plot_layout(widths = c(5, 1))
ggsave("genomewide_and_density_fst_region.tif", width=16, height=18, dpi = 300)
ggsave("genomewide_and_density_fst_region.png", width=16, height=18, dpi = 300)
# Fst closer to 1 = more genetically distinct
# Fst closer to 0 = more genetically similar


# Region map metadata. Just have one random point representing each region.
location <- read.csv("input/region_metadata.csv", header = TRUE)

## Make world map data
world_map <- map_data("world")

# Set colors for the points
red_palette1 <- brewer.pal(n = 9, name = "YlOrRd")
red_palette2 <- brewer.pal(n=9, name = "YlOrBr")
blue_palette <- brewer.pal(n = 9, name = "Blues")
green_palette1 <- brewer.pal(n = 9, name = "BuGn")
green_palette2 <- brewer.pal(n=9, name = "YlGn")

scale_colour_region <- 
  c('USA' = red_palette2[7],
    'CENAM' = "blueviolet",
    'EUR' = green_palette1[7],
    'ASIA' = "violetred1",
    'AUS' = blue_palette[5])
# chose one random location/point for each region

# Get median fst for each comparison
fst_median <- data %>%
  filter(data_type == "Fst") %>%
  group_by(comparison) %>%
  summarise(FST_MEDIAN = median(value, na.rm = TRUE))
# Save as csv file and then add required metadata such as population longitudes/latitudes
write.csv(fst_median, "fst_median.csv", row.names = TRUE)

# read new csv file with map metadata and median fst between groups
fst_map <- read.csv("fst_map.csv", header = TRUE)

# Fst by region on a map
plot_3_fst <- ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = "grey90") +
  geom_segment(data=fst_map, aes(x=POP1_LONG, y=POP1_LAT, xend=POP2_LONG, yend=POP2_LAT, color = FST_MEDIAN), size = 0.7) +
  scale_color_gradient(low = "gold", high = "red3", name = expression(F[ST])) +
  geom_point(data = location, aes(x = Longitude, y = Latitude, fill = Region), size = 6, shape = 21) +
  scale_fill_manual(values = scale_colour_region, limits = c("USA", "CENAM", "EUR", "ASIA", "AUS"), name = "Region") +
  geom_text_repel(data = location, aes(x = Longitude, y = Latitude, label = Region), size = 4, fontface = "bold", nudge_y = 7) +
  theme_void() +
  theme(
    legend.position = c(0.1,0.3),) +
  guides(fill = "none") +
  coord_fixed()
plot_3_fst

# Save plot
ggsave("fst_map_region.tif", bg = "white", height=5, width=10, dpi = 300)
ggsave("fst_map_region.png", bg = "white", height=5, width=10, dpi = 300)


# Fst by region on a heatmap
fst_heatmap <- read.csv("fst_heatmap.csv", header = TRUE)

plot_4_fst <- ggplot(fst_heatmap, aes(x = POP1, y=POP2, fill = FST_MEDIAN)) + 
  geom_tile(color = "white", lwd = 1.5, linetype = 1) +
  coord_fixed() +
  scale_fill_material("indigo") +
  guides(fill = guide_colourbar(title = expression(F[ST]), barwidth = 1, barheight = 5)) +
  theme_minimal() +
  labs(x = "Population", y = "Population") +
  theme(axis.text = element_text(size = 12),
    axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 10),
        panel.grid.major = element_blank())
plot_4_fst

ggsave("fst_heatmap_region.tif", bg = "white", height=5, width = 8, dpi = 300)
ggsave("fst_heatmap_region.png", bg = "white", height=5, width = 8, dpi = 300)
```

### Pixy by host

```R
# Batch 4 - Pixy for Pi, Fst and Dxy by HOST

library(tidyverse)
library(ggsci)
library(ggpubr)
library(patchwork)
library(ggridges)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(purrr)
library(stats)
library(graphics)
library(grDevices)
library(utils)
library(datasets)
library(methods)
library(base)
require(maps)
library(mapdata)
library(readxl)
library(ozmaps) 
library(grid)
library(gridExtra)
library(ggrepel)
library(ggnewscale)
library(reshape)
library(wesanderson)

setwd("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/batch4/FILTER1/NO_OUTGROUPS/pixy/host")

##### Nucleotide diversity (Pi) #####

# get nucleotide diversity (pi) data from pixy output
pi_data <- read.table("input/host_pi.txt", header=T)
pi_data$chromosome <- str_remove(pi_data$chromosome, 'dirofilaria_immitis_')

# subset
pi_data_Ferret <- pi_data %>%
  filter(pop=="Ferret")
pi_data_Golden_jackal <- pi_data %>%
  filter(pop=="Golden_jackal")
pi_data_Wildcat <- pi_data %>%
  filter(pop=="Wildcat")
pi_data_Leopard <- pi_data %>%
  filter(pop=="Leopard")
pi_data_Cat <- pi_data %>%
  filter(pop=="Cat")
pi_data_Fox <- pi_data %>%
  filter(pop=="Fox")
pi_data_Dog <- pi_data %>%
  filter(pop=="Dog")

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

# calculate the median Pi and check the ratio of sex-to-autosome diversity. Should be about 0.75, as D. immitis is XX/XY
pi_data_sex_median <- pi_data %>%
  group_by(chr_type) %>%
  summarise(median = median(avg_pi, na.rm = TRUE))
pi_data_sex_median

'
# A tibble: 2 × 2
  chr_type   median
  <chr>       <dbl>
1 autosome 0.000524
2 sexchr   0.000183
'
# 0.000183 / 0.000524 = 0.3492366 (way off 0.75)

pi_data_pop_sex_median <- pi_data %>%
  group_by(pop, chr_type) %>%
  summarise(median = median(avg_pi, na.rm = TRUE))
pi_data_pop_sex_median

'
# A tibble: 14 × 3
# Groups:   pop [7]
   pop           chr_type    median
   <chr>         <chr>        <dbl>
 1 Cat           autosome 0.000685 
 2 Cat           sexchr   0.000295 
 3 Dog           autosome 0.000701 
 4 Dog           sexchr   0.000274 
 5 Ferret        autosome 0.00121  
 6 Ferret        sexchr   0.000523 
 7 Fox           autosome 0.000614 
 8 Fox           sexchr   0.000212 
 9 Golden_jackal autosome 0.000177 
10 Golden_jackal sexchr   0.0000108
11 Leopard       autosome 0.000110 
12 Leopard       sexchr   0.0000320
13 Wildcat       autosome 0.000189 
14 Wildcat       sexchr   0.0000308
'

# plot 1 - genome wide plots per population
plot_1_pi <- ggplot(pi_data, aes(position*100000, avg_pi, col=chromosome)) +
  geom_point(size=0.8) +
  facet_grid(pop~.) +
  scale_color_npg(guide = guide_legend(override.aes = list(size = 4))) +
  theme_bw() +
  theme(legend.position="top", 
        legend.title = element_blank(),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 10),
        legend.text = element_text(size = 16)) +
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
  theme(legend.position="top",
        legend.title = element_blank(),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 10),
        legend.text = element_text(size = 16)) +
  labs(x="Nucleotide Diversity (Pi)", y="Density")
plot_2_pi

# combine plots
plot_1_pi + plot_2_pi +  plot_layout(widths = c(5, 1))
ggsave("genomewide_and_density_Pi_host.tif", width=14, height=8, dpi = 300)
ggsave("genomewide_and_density_Pi_host.png", width=14, height=8, dpi = 300)

#Now a boxplot of the pi value per population
boxplot_pi <- ggplot(pi_data, aes(pop, avg_pi, col=pop)) +
  geom_jitter(size = 1, alpha = 0.5) +
  geom_boxplot(fill=NA, col="black") +
  labs(x = "Population" , y = "Nucleotide diversity (Pi)", colour = "Population") +
  scale_color_npg() +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 12)
        )

boxplot_pi

ggsave("boxplot_Pi_host.tif", width=8, height=6, dpi = 300)
ggsave("boxplot_Pi_host.png", width=8, height=6, dpi = 300)
# higher pi value = more diverse population
# lower pi value = less diverse population
## Looks like USA and AUS have higher diversity. Asia has the lowest diversity but this is probably because we have a small sample size all from the same country (Thailand).


#let's check the statistical significance of the differences between pi values
#Let's explore the normality
pi_data %>% filter(chromosome != 'chrX') %>%
  ggplot(., aes(x = avg_pi, colour = pop)) +
  geom_histogram() +
  theme_bw() +
  facet_grid(. ~ pop)

#Confirm the absence of normality by shapiro
## If p > 0.05 the data is normal, if p < 0.05 the data is NOT normal
pi_data_Cat <- pi_data_Cat %>% filter(chromosome != 'chrX')
Cat_shapiro <- shapiro.test(pi_data_Cat$avg_pi)
print(Cat_shapiro)
# W = 0.92623, p-value < 2.2e-16

pi_data_Dog <- pi_data_Dog %>% filter(chromosome != 'chrX')
Dog_shapiro <- shapiro.test(pi_data_Dog$avg_pi)
print(Dog_shapiro)
# W = 0.92208, p-value < 2.2e-16

pi_data_Ferret <- pi_data_Ferret %>% filter(chromosome != 'chrX')
Ferret_shapiro <- shapiro.test(pi_data_Ferret$avg_pi)
print(Ferret_shapiro)
# W = 0.93423, p-value = 1.269e-15

pi_data_Fox <- pi_data_Fox %>% filter(chromosome != 'chrX')
Fox_shapiro <- shapiro.test(pi_data_Fox$avg_pi)
print(Fox_shapiro)
# W = 0.93245, p-value = 7.521e-16

pi_data_Golden_jackal <- pi_data_Golden_jackal %>% filter(chromosome != 'chrX')
Golden_jackal_shapiro <- shapiro.test(pi_data_Golden_jackal$avg_pi)
print(Golden_jackal_shapiro)
# W = 0.79614, p-value < 2.2e-16

pi_data_Leopard <- pi_data_Leopard %>% filter(chromosome != 'chrX')
Leopard_shapiro <- shapiro.test(pi_data_Leopard$avg_pi)
print(Leopard_shapiro)
# W = 0.74801, p-value < 2.2e-16

pi_data_Wildcat <- pi_data_Wildcat %>% filter(chromosome != 'chrX')
Wildcat_shapiro <- shapiro.test(pi_data_Wildcat$avg_pi)
print(Wildcat_shapiro)
# W = 0.77864, p-value < 2.2e-16

#Wilcoxson test for all pairs of pops
## This tests whether the mean values of 2 dependent groups differ significantly from each other
## If p < 0.05 they are sig different
wilcox.test(pi_data_Cat$avg_pi, pi_data_Dog$avg_pi)
# W = 173979, p-value = 0.2312
wilcox.test(pi_data_Cat$avg_pi, pi_data_Ferret$avg_pi)
# W = 115875, p-value < 2.2e-16
wilcox.test(pi_data_Cat$avg_pi, pi_data_Fox$avg_pi)
# W = 195866, p-value = 0.01507
wilcox.test(pi_data_Cat$avg_pi, pi_data_Golden_jackal$avg_pi)
# W = 270425, p-value < 2.2e-16
wilcox.test(pi_data_Cat$avg_pi, pi_data_Leopard$avg_pi)
# W = 283042, p-value < 2.2e-16
wilcox.test(pi_data_Cat$avg_pi, pi_data_Wildcat$avg_pi)
# W = 269034, p-value < 2.2e-16
wilcox.test(pi_data_Dog$avg_pi, pi_data_Ferret$avg_pi)
# W = 119780, p-value < 2.2e-16
wilcox.test(pi_data_Dog$avg_pi, pi_data_Fox$avg_pi)
# W = 202948, p-value = 0.0003125
wilcox.test(pi_data_Dog$avg_pi, pi_data_Golden_jackal$avg_pi)
# W = 274052, p-value < 2.2e-16
wilcox.test(pi_data_Dog$avg_pi, pi_data_Leopard$avg_pi)
# W = 287047, p-value < 2.2e-16
wilcox.test(pi_data_Dog$avg_pi, pi_data_Wildcat$avg_pi)
# W = 272505, p-value < 2.2e-16
wilcox.test(pi_data_Ferret$avg_pi, pi_data_Fox$avg_pi)
# W = 257601, p-value < 2.2e-16
wilcox.test(pi_data_Ferret$avg_pi, pi_data_Golden_jackal$avg_pi)
# W = 304815, p-value < 2.2e-16
wilcox.test(pi_data_Ferret$avg_pi, pi_data_Leopard$avg_pi)
# W = 313123, p-value < 2.2e-16
wilcox.test(pi_data_Ferret$avg_pi, pi_data_Wildcat$avg_pi)
# W = 303776, p-value < 2.2e-16
wilcox.test(pi_data_Fox$avg_pi, pi_data_Golden_jackal$avg_pi)
# W = 261257, p-value < 2.2e-16
wilcox.test(pi_data_Fox$avg_pi, pi_data_Leopard$avg_pi)
# W = 274564, p-value < 2.2e-16
wilcox.test(pi_data_Fox$avg_pi, pi_data_Wildcat$avg_pi)
# W = 260233, p-value < 2.2e-16
wilcox.test(pi_data_Golden_jackal$avg_pi, pi_data_Leopard$avg_pi)
# W = 192939, p-value = 0.05047
wilcox.test(pi_data_Golden_jackal$avg_pi, pi_data_Wildcat$avg_pi)
# W = 186263, p-value = 0.3988
wilcox.test(pi_data_Leopard$avg_pi, pi_data_Wildcat$avg_pi)
# W = 174470, p-value = 0.2607


#let's generate some dataframes for the estatistics
# subset
pi_data_Cat <- pi_data %>%
  filter(pop=="Cat") %>%
  filter(chromosome != 'chrX')

pi_data_Dog <- pi_data %>%
  filter(pop=="Dog") %>%
  filter(chromosome != 'chrX') 

pi_data_Ferret <- pi_data %>%
  filter(pop=="Ferret") %>%
  filter(chromosome != 'chrX')

pi_data_Fox <- pi_data %>%
  filter(pop=="Fox") %>%
  filter(chromosome != 'chrX')

pi_data_Golden_jackal <- pi_data %>%
  filter(pop=="Golden_jackal") %>%
  filter(chromosome != 'chrX')

pi_data_Leopard <- pi_data %>%
  filter(pop=="Leopard") %>%
  filter(chromosome != 'chrX')

pi_data_Wildcat <- pi_data %>%
  filter(pop=="Wildcat") %>%
  filter(chromosome != 'chrX')


#### Dxy and Fst ####

# load data
dxy_data <- read.table("input/host_dxy.txt", header=T)
fst_data <- read.table("input/host_fst.txt", header=T)

# add some columns
dxy_data$data_type <- "Dxy"
dxy_data <- mutate(dxy_data,
                   comparison = paste(pop1, pop2, sep = '_v_'))
fst_data$data_type <- "Fst"
fst_data <- mutate(fst_data,
                   comparison = paste(pop1, pop2, sep = '_v_'))

# subset and merge dataframes
dxy_data_sub <- dxy_data %>% select(comparison,data_type,chromosome,window_pos_1,window_pos_2,avg_dxy)
colnames(dxy_data_sub) <- c("comparison","data_type","chromosome","window_pos_1","window_pos_2","value")

fst_data_sub <- fst_data %>% select(comparison,data_type,chromosome,window_pos_1,window_pos_2,avg_wc_fst)
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
  summarise(median = median(value, na.rm = TRUE)) %>%
  print(n = 42)
'
# A tibble: 42 × 3
# Groups:   comparison [21]
   comparison              data_type    median
   <chr>                   <chr>         <dbl>
 1 Cat_v_Ferret            Dxy        0.000615
 2 Cat_v_Ferret            Fst       -0.105   
 3 Cat_v_Golden_jackal     Dxy        0.000601
 4 Cat_v_Golden_jackal     Fst        0.0318  
 5 Cat_v_Leopard           Dxy        0.000567
 6 Cat_v_Leopard           Fst        0.0470  
 7 Cat_v_Wildcat           Dxy        0.000577
 8 Cat_v_Wildcat           Fst        0.0429  
 9 Dog_v_Cat               Dxy        0.000584
10 Dog_v_Cat               Fst        0.0226  
11 Dog_v_Ferret            Dxy        0.000659
12 Dog_v_Ferret            Fst       -0.0435  
13 Dog_v_Fox               Dxy        0.000549
14 Dog_v_Fox               Fst       -0.00135 
15 Dog_v_Golden_jackal     Dxy        0.000549
16 Dog_v_Golden_jackal     Fst        0.0126  
17 Dog_v_Leopard           Dxy        0.000534
18 Dog_v_Leopard           Fst        0.0249  
19 Dog_v_Wildcat           Dxy        0.000540
20 Dog_v_Wildcat           Fst        0.0174  
21 Fox_v_Cat               Dxy        0.000573
22 Fox_v_Cat               Fst        0.0407  
23 Fox_v_Ferret            Dxy        0.000645
24 Fox_v_Ferret            Fst        0.0571  
25 Fox_v_Golden_jackal     Dxy        0.000493
26 Fox_v_Golden_jackal     Fst        0.105   
27 Fox_v_Leopard           Dxy        0.000504
28 Fox_v_Leopard           Fst        0.109   
29 Fox_v_Wildcat           Dxy        0.000494
30 Fox_v_Wildcat           Fst        0.0944  
31 Golden_jackal_v_Ferret  Dxy        0.000662
32 Golden_jackal_v_Ferret  Fst        0       
33 Leopard_v_Ferret        Dxy        0.000651
34 Leopard_v_Ferret        Fst        0       
35 Leopard_v_Golden_jackal Dxy        0.000257
36 Leopard_v_Golden_jackal Fst        0       
37 Leopard_v_Wildcat       Dxy        0.000241
38 Leopard_v_Wildcat       Fst        0       
39 Wildcat_v_Ferret        Dxy        0.000652
40 Wildcat_v_Ferret        Fst        0       
41 Wildcat_v_Golden_jackal Dxy        0.000261
42 Wildcat_v_Golden_jackal Fst 
'


#Plotting Dxy

# plot 1 - genome wide plots per comparison
plot_1_dxy <- data %>%
  filter(data_type =='Dxy')%>%
  ggplot(., aes(position*100000, value, col=chromosome)) +
  geom_point(size=1.5) +
  facet_grid(comparison~.) +
  scale_color_npg(guide = guide_legend(override.aes = list(size = 4))) +
  theme_bw() +
  theme(legend.position="top", 
        legend.title = element_blank(),
        axis.title = element_text(size = 32),
        axis.text = element_text(size = 16),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 20)) +
  labs(x="Genomic Position", y="Dxy")
plot_1_dxy

# plot 2 - density plots of dxy per group
plot_2_dxy <- data %>%
  filter(data_type=='Dxy')%>%
  ggplot(aes(value, chr_type, fill=chr_type), guide="none") +
  geom_density_ridges(quantile_lines=TRUE, quantile_fun=function(x,...)median(x), size=0.5) +
  theme_bw() + theme(legend.position = "top", 
                     legend.title = element_blank(),
                     axis.text.y=element_blank(),
                     axis.title = element_text(size = 32),
                     axis.text = element_text(size = 16),
                     strip.text = element_text(size = 8),
                     legend.text = element_text(size = 20)) +
  facet_grid(comparison~.) +
  scale_fill_npg() +
  labs(x="Dxy", y="Density")
plot_2_dxy

# combine plots
dxy_plot <- plot_1_dxy + plot_2_dxy +  plot_layout(widths = c(5, 1))
ggsave("genomewide_and_density_dxy_host.tif", width=20, height=30, dpi = 300)
ggsave("genomewide_and_density_dxy_host.png", width=20, height=30, dpi = 300)
# higher Dxy = more divergent

#some additional plots
boxplot_dxy <- data %>%
  filter(data_type =='Dxy') %>% 
  ggplot(., aes(comparison, value, col=comparison)) +
  geom_jitter(size = 1, alpha = 0.5) +
  geom_boxplot(fill=NA, col="black") +
  labs(x = "Host" , y = "Dxy", colour = "Population") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 18),
        legend.text = element_text(size = 16)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
boxplot_dxy
ggsave("boxplot_dxy_host.tif", width=20, height = 6, dpi = 300)
ggsave("boxplot_dxy_host.png", width=20, height = 6, dpi = 300)

density_dxy <- data %>%
  filter(data_type =='Dxy') %>% 
  ggplot(., aes(x=value, group = comparison, fill=comparison)) +
  geom_density(adjust=1.5, alpha=.4) +
  labs(x = "Dxy" , y = "Density", colour = "Comparison") +
  theme_bw() +
  theme(legend.position = c(0.73,0.77),
        legend.title = element_blank(),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
density_dxy
ggsave("density_dxy_host.tif", width = 8, height = 8, dpi = 300)
ggsave("density_dxy_host.png", width = 8, height = 8, dpi = 300)


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
  annotate(geom="text", x=140*100000, y=0.004, label="ChrX", size = 6)+
  annotate(geom="text", x=359.5*100000, y=0.004, label="Chr1", size = 6)+
  annotate(geom="text", x=520*100000, y=0.004, label="Chr2", size = 6)+
  annotate(geom="text", x=665*100000, y=0.004, label="Chr3", size = 6)+
  annotate(geom="text", x=810.5*100000, y=0.004, label="Chr4", size = 6)+
  geom_line() +
  theme_bw() +
  theme(legend.position="top", 
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 14)) +
  labs(x="Genomic Position", y="Dxy") +
  guides(col = guide_legend(override.aes = list(linewidth = 2)))
dxy_lineplot
ggsave("lineplot_dxy_host.tif", width=14, height=8, dpi = 300)
ggsave("lineplot_dxy_host.png", width=14, height=8, dpi = 300)


# Plotting Fst

# plot 1 - genome wide plots per population
plot_1_fst <- data %>%
  filter(data_type=='Fst')%>%
  ggplot(., aes(position*100000, value, col=chromosome)) +
  geom_point(size=1) +
  facet_grid(comparison~.) +
  scale_color_npg(guide = guide_legend(override.aes = list(size = 4))) +
  theme_bw() +
  theme(legend.position="top", 
        legend.title = element_blank(),
        axis.title = element_text(size = 32),
        axis.text = element_text(size = 16),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 20)) +
  labs(x="Genomic Position", y="Fst")
plot_1_fst
# Removed 468 rows containing missing values or values outside the scale range (`geom_point()`). 

# plot 2 - density plots of Fst per group
plot_2_fst <- data %>%
  filter(data_type=='Fst')%>%
  ggplot(aes(value, chr_type, fill=chr_type), guide="none") +
  geom_density_ridges(quantile_lines=TRUE, quantile_fun=function(x,...)median(x), size=0.5) +
  theme_bw() + 
  theme(legend.position = "none", 
        axis.text.y=element_blank(),
        axis.title = element_text(size = 32),
        axis.text = element_text(size = 16),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 20)) +
  facet_grid(comparison~.) +
  scale_fill_npg() +
  labs(x="Fst", y="Density")
plot_2_fst
# Removed 468 rows containing non-finite outside the scale range (`stat_density_ridges()`). 

# combine plots
plot_1_fst + plot_2_fst +  plot_layout(widths = c(5, 1))
## Fst values should only range between 0 and 1. Set anything <0 to 0 and re-do the plots.

# Data cleaning
data_fst_clean <- data %>% filter(data_type == 'Fst') %>% mutate(value = ifelse(value < 0, 0, value)) # this takes my original data, filters it to only include rows where the data_type is Fst. Then mutate function is used to create/modify columns. It modifies the value column, checks if the value is less than 0. If the value is <0, it replaces it with 0. If the value is >0, it keeps the original value.

# Re-do plots with cleaned data
# plot 1 - genome wide plots per population
plot_1_fst_clean <- data_fst_clean %>%
  filter(data_type=='Fst')%>%
  ggplot(., aes(position*100000, value, col=chromosome)) +
  geom_point(size=1) +
  facet_grid(comparison~.) +
  scale_color_tron() +
  theme_bw() +
  theme(legend.position="top", 
        legend.title = element_blank(),
        axis.title = element_text(size = 32),
        axis.text = element_text(size = 16),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 20)) +
  labs(x="Genomic Position", y="Fst")
plot_1_fst_clean


# plot 2 - density plots of Fst per group
plot_2_fst_clean <- data_fst_clean %>%
  filter(data_type=='Fst')%>%
  ggplot(aes(value, chr_type, fill=chr_type), guide="none") +
  geom_density_ridges(quantile_lines=TRUE, quantile_fun=function(x,...)median(x), size=0.5) +
  theme_bw() + theme(legend.position = "none", axis.text.y=element_blank()) +
  facet_grid(comparison~.) +
  scale_fill_npg() +
  labs(x="Fst", y="Density")
plot_2_fst_clean

# combine plots
plot_1_fst_clean + plot_2_fst_clean +  plot_layout(widths = c(5, 1))
ggsave("genomewide_and_density_fst_host.tif", width=20, height=30, dpi = 300)
ggsave("genomewide_and_density_fst_host.png", width=20, height=30, dpi = 300)
# Fst closer to 1 = more genetically distinct
# Fst closer to 0 = more genetically similar
# Doesn't look too great because some hosts only have 1 sample, and there's a bunch of NA data and negative values I had to clean.

# Fst by host on a heatmap
# Get median fst for each comparison
fst_median <- data_fst_clean %>%
  filter(data_type == "Fst") %>%
  group_by(comparison) %>%
  summarise(FST_MEDIAN = median(value, na.rm = TRUE))
# Save as csv file and then add required metadata such as population longitudes/latitudes
write.csv(fst_median, "fst_median.csv", row.names = TRUE)

fst_heatmap <- read.csv("fst_heatmap.csv", header = TRUE)

plot_4_fst <- ggplot(fst_heatmap, aes(x = POP1, y=POP2, fill = FST_MEDIAN)) + 
  geom_tile(color = "white", lwd = 1.5, linetype = 1) +
  coord_fixed() +
  scale_fill_material("indigo") +
  guides(fill = guide_colourbar(title = expression(F[ST]), barwidth = 1, barheight = 5)) +
  theme_minimal() +
  labs(x = "Host", y = "Host") +
  theme(axis.text = element_text(size = 12),
    axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 10),
        panel.grid.major = element_blank())
plot_4_fst

ggsave("fst_heatmap_host.tif", bg = "white", height=5, width = 8, dpi = 300)
ggsave("fst_heatmap_host.png", bg = "white", height=5, width = 8, dpi = 300)
```
When we break it down by host, we see some differences which could be inflating the region results we got. Re-do pixy by region for dog samples only.

### Pixy by region - DOG ONLY

```bash
conda activate pixy

module load PaM/environment
module load bsub.py/0.42.1
module load vcftools/0.1.16-c4
module load htslib-1.19/perl-5.38.0

WORKING_DIR=/lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA

cd ${WORKING_DIR}/03_ANALYSIS/04_VARIANTS/FILTER1/NO_OUTGROUPS/FINAL_SETS

# Select dog samples only in vcf
vcftools --gzvcf nuclearSNPssandINVARIANTs_124samples.chrXto4.recode.vcf.gz \
--keep ../../../../05_ANALYSIS/PIXY/DOG_ONLY_samplelist.keep \
--recode --out DOG_ONLY_nuclearSNPsandINVARIANTs.chrXto4
# After filtering, kept 113 out of 124 Individuals
# After filtering, kept 84210689 out of a possible 84210689 Sites

##up to here
bgzip DOG_ONLY_nuclearSNPsandINVARIANTs.chrXto4.recode.vcf;
tabix -p vcf DOG_ONLY_nuclearSNPsandINVARIANTs.chrXto4.recode.vcf.gz

# run pixy
VCF=${WORKING_DIR}/03_ANALYSIS/04_VARIANTS/FILTER1/NO_OUTGROUPS/FINAL_SETS/DOG_ONLY_nuclearSNPsandINVARIANTs.chrXto4.recode.vcf.gz

cd ../../../../05_ANALYSIS/PIXY

bsub.py --queue long --threads 10 20 pixy_region_DOG_ONLY \
"pixy --stats pi fst dxy \
--vcf ${VCF} \
--populations region_DOG_ONLY_pop.list \
--window_size 100000 \
--n_cores 10 \
--output_prefix region_DOG_ONLY"
# Successfully completed
```

```R
# Batch 4 - Pixy for Pi, Fst and Dxy by REGION for dog samples only

library(tidyverse)
library(ggsci)
library(ggpubr)
library(patchwork)
library(ggridges)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(purrr)
library(stats)
library(graphics)
library(grDevices)
library(utils)
library(datasets)
library(methods)
library(base)
require(maps)
library(mapdata)
library(readxl)
library(ozmaps) 
library(grid)
library(gridExtra)
library(ggrepel)
library(ggnewscale)
library(reshape)
library(wesanderson)

setwd("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/batch4/FILTER1/NO_OUTGROUPS/pixy/region_DOG_ONLY")

##### Nucleotide diversity (Pi) #####

# get nucleotide diversity (pi) data from pixy output
pi_data <- read.table("input/region_DOG_ONLY_pi.txt", header=T)
pi_data$chromosome <- str_remove(pi_data$chromosome, 'dirofilaria_immitis_')

# subset
pi_data_AUS <- pi_data %>%
  filter(pop=="AUS")
pi_data_ASIA <- pi_data %>%
  filter(pop=="ASIA")
pi_data_USA <- pi_data %>%
  filter(pop=="USA")
pi_data_EUR <- pi_data %>%
  filter(pop=="EUR")
pi_data_CENAM <- pi_data %>%
  filter(pop=="CENAM")

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

# calculate the median Pi and check the ratio of sex-to-autosome diversity. Should be about 0.75, as D. immitis is XX/XY
pi_data_sex_median <- pi_data %>%
  group_by(chr_type) %>%
  summarise(median = median(avg_pi, na.rm = TRUE))
pi_data_sex_median

'
# A tibble: 2 × 2
  chr_type   median
  <chr>       <dbl>
1 autosome 0.000459
2 sexchr   0.000177
'
# 0.000177 / 0.000459 = 0.3856209 (way off 0.75)

pi_data_pop_sex_median <- pi_data %>%
  group_by(pop, chr_type) %>%
  summarise(median = median(avg_pi, na.rm = TRUE))
pi_data_pop_sex_median

'
# A tibble: 10 × 3
# Groups:   pop [5]
   pop   chr_type   median
   <chr> <chr>       <dbl>
 1 ASIA  autosome 0.000289
 2 ASIA  sexchr   0.000109
 3 AUS   autosome 0.000561
 4 AUS   sexchr   0.000190
 5 CENAM autosome 0.000352
 6 CENAM sexchr   0.000139
 7 EUR   autosome 0.000429
 8 EUR   sexchr   0.000186
 9 USA   autosome 0.000739
10 USA   sexchr   0.000323
'

# plot 1 - genome wide plots per population
plot_1_pi <- ggplot(pi_data, aes(position*100000, avg_pi, col=chromosome)) +
  geom_point(size=0.8) +
  facet_grid(pop~.) +
  scale_color_npg(guide = guide_legend(override.aes = list(size = 4))) +
  theme_bw() +
  theme(legend.position="top", 
        legend.title = element_blank(),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 18),
        legend.text = element_text(size = 16)) +
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
  theme(legend.position="top",
        legend.title = element_blank(),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 18),
        legend.text = element_text(size = 16)) +
  labs(x="Nucleotide Diversity (Pi)", y="Density")
plot_2_pi

# combine plots
plot_1_pi + plot_2_pi +  plot_layout(widths = c(5, 1))
ggsave("genomewide_and_density_Pi_region_DOG_ONLY.tif", width=14, height=8, dpi = 300)
ggsave("genomewide_and_density_Pi_region_DOG_ONLY.png", width=14, height=8, dpi = 300)

#Now a boxplot of the pi value per population

red_palette1 <- brewer.pal(n = 9, name = "YlOrRd")
red_palette2 <- brewer.pal(n=9, name = "YlOrBr")
blue_palette <- brewer.pal(n = 9, name = "Blues")
green_palette1 <- brewer.pal(n = 9, name = "BuGn")
green_palette2 <- brewer.pal(n=9, name = "YlGn")

scale_colour_region <- function(...){
  ggplot2:::manual_scale(
    'colour', 
    values = setNames(
      c(blue_palette[5],
        "violetred1",
        red_palette1[7], 
        "blueviolet",
        green_palette1[7]),
    c("AUS", "ASIA", "USA", "EUR", "CENAM")), 
  ...
)
}    

boxplot_pi <- ggplot(pi_data, aes(pop, avg_pi, col=pop)) +
  geom_jitter(size = 1, alpha = 0.5) +
  geom_boxplot(fill=NA, col="black") +
  labs(x = "Population" , y = "Nucleotide diversity (Pi)", colour = "Population") +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 18)
        ) +
  scale_colour_region ()

boxplot_pi

ggsave("boxplot_Pi_region_DOG_ONLY.tif", width=6, height=6, dpi = 300)
ggsave("boxplot_Pi_region_DOG_ONLY.png", width=6, height=6, dpi = 300)
# higher pi value = more diverse population
# lower pi value = less diverse population
## Looks like USA and AUS have higher diversity.


#let's check the statistical significance of the differences between pi values
#Let's explore the normality
pi_data %>% filter(chromosome != 'chrX') %>%
  ggplot(., aes(x = avg_pi, colour = pop)) +
  geom_histogram() +
  theme_bw() +
  facet_grid(. ~ pop)

#Confirm the absence of normality by shapiro
## If p > 0.05 the data is normal, if p < 0.05 the data is NOT normal
pi_data_AUS <- pi_data_AUS %>% filter(chromosome != 'chrX')
AUS_shapiro <- shapiro.test(pi_data_AUS$avg_pi)
print(AUS_shapiro)
# W = 0.91953, p-value < 2.2e-16

pi_data_ASIA <- pi_data_ASIA %>% filter(chromosome != 'chrX')
ASIA_shapiro <- shapiro.test(pi_data_ASIA$avg_pi)
print(ASIA_shapiro)
# W = 0.87735, p-value < 2.2e-16

pi_data_USA <- pi_data_USA %>% filter(chromosome != 'chrX')
USA_shapiro <- shapiro.test(pi_data_USA$avg_pi)
print(USA_shapiro)
# W = 0.91881, p-value < 2.2e-16

pi_data_EUR <- pi_data_EUR %>% filter(chromosome != 'chrX')
EUR_shapiro <- shapiro.test(pi_data_EUR$avg_pi)
print(EUR_shapiro)
# W = 0.91837, p-value < 2.2e-16

pi_data_CENAM <- pi_data_CENAM %>% filter(chromosome != 'chrX')
CENAM_shapiro <- shapiro.test(pi_data_CENAM$avg_pi)
print(CENAM_shapiro)
# W = 0.89723, p-value < 2.2e-16

#Wilcoxson test for all pairs of pops
## This tests whether the mean values of 2 dependent groups differ significantly from each other
## If p < 0.05 they are sig different
wilcox.test(pi_data_AUS$avg_pi, pi_data_ASIA$avg_pi)
# W = 252451, p-value < 2.2e-16
wilcox.test(pi_data_AUS$avg_pi, pi_data_USA$avg_pi)
# W = 143196, p-value = 2.973e-10
wilcox.test(pi_data_AUS$avg_pi, pi_data_EUR$avg_pi)
# W = 218254, p-value = 8.156e-10
wilcox.test(pi_data_AUS$avg_pi, pi_data_CENAM$avg_pi)
# W = 234515, p-value < 2.2e-16
wilcox.test(pi_data_ASIA$avg_pi, pi_data_USA$avg_pi)
# W = 84065, p-value < 2.2e-16
wilcox.test(pi_data_ASIA$avg_pi, pi_data_EUR$avg_pi)
# W = 144658, p-value = 1.38e-09
wilcox.test(pi_data_ASIA$avg_pi, pi_data_CENAM$avg_pi)
# W = 158052, p-value = 0.0001243
wilcox.test(pi_data_USA$avg_pi, pi_data_EUR$avg_pi)
# W = 250455, p-value < 2.2e-16
wilcox.test(pi_data_USA$avg_pi, pi_data_CENAM$avg_pi)
# W = 264557, p-value < 2.2e-16
wilcox.test(pi_data_EUR$avg_pi, pi_data_CENAM$avg_pi)
# W = 197085, p-value = 0.008468

## They're all sig different


#let's generate some dataframes for the estatistics
# subset
pi_data_AUS <- pi_data %>%
  filter(pop=="AUS") %>%
  filter(chromosome != 'chrX')

pi_data_ASIA <- pi_data %>%
  filter(pop=="ASIA") %>%
  filter(chromosome != 'chrX') 

pi_data_USA <- pi_data %>%
  filter(pop=="USA") %>%
  filter(chromosome != 'chrX')

pi_data_EUR <- pi_data %>%
  filter(pop=="EUR") %>%
  filter(chromosome != 'chrX')

pi_data_CENAM <- pi_data %>%
  filter(pop=="CENAM") %>%
  filter(chromosome != 'chrX')




#### Dxy and Fst ####

# load data
dxy_data <- read.table("input/region_DOG_ONLY_dxy.txt", header=T)
fst_data <- read.table("input/region_DOG_ONLY_fst.txt", header=T)

# add some columns
dxy_data$data_type <- "Dxy"
dxy_data <- mutate(dxy_data,
                   comparison = paste(pop1, pop2, sep = '_v_'))
fst_data$data_type <- "Fst"
fst_data <- mutate(fst_data,
                   comparison = paste(pop1, pop2, sep = '_v_'))

# subset and merge dataframes
dxy_data_sub <- dxy_data %>% select(comparison,data_type,chromosome,window_pos_1,window_pos_2,avg_dxy)
colnames(dxy_data_sub) <- c("comparison","data_type","chromosome","window_pos_1","window_pos_2","value")

fst_data_sub <- fst_data %>% select(comparison,data_type,chromosome,window_pos_1,window_pos_2,avg_wc_fst)
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
# A tibble: 20 × 3
# Groups:   comparison [10]
   comparison   data_type   median
   <chr>        <chr>        <dbl>
 1 ASIA_v_CENAM Dxy       0.000613
 2 ASIA_v_CENAM Fst       0.473   
 3 ASIA_v_USA   Dxy       0.000600
 4 ASIA_v_USA   Fst       0.185   
 5 AUS_v_ASIA   Dxy       0.000513
 6 AUS_v_ASIA   Fst       0.223   
 7 AUS_v_CENAM  Dxy       0.000571
 8 AUS_v_CENAM  Fst       0.290   
 9 AUS_v_EUR    Dxy       0.000582
10 AUS_v_EUR    Fst       0.278   
11 AUS_v_USA    Dxy       0.000630
12 AUS_v_USA    Fst       0.165   
13 CENAM_v_USA  Dxy       0.000634
14 CENAM_v_USA  Fst       0.211   
15 EUR_v_ASIA   Dxy       0.000624
16 EUR_v_ASIA   Fst       0.433   
17 EUR_v_CENAM  Dxy       0.000495
18 EUR_v_CENAM  Fst       0.268   
19 EUR_v_USA    Dxy       0.000633
20 EUR_v_USA    Fst       0.193     
'


#Plotting Dxy

# plot 1 - genome wide plots per comparison
plot_1_dxy <- data %>%
  filter(data_type =='Dxy')%>%
  ggplot(., aes(position*100000, value, col=chromosome)) +
  geom_point(size=1.5) +
  facet_grid(comparison~.) +
  scale_color_npg(guide = guide_legend(override.aes = list(size = 4))) +
  theme_bw() +
  theme(legend.position="top", 
        legend.title = element_blank(),
        axis.title = element_text(size = 32),
        axis.text = element_text(size = 16),
        strip.text = element_text(size = 14),
        legend.text = element_text(size = 20)) +
  labs(x="Genomic Position", y="Dxy")
plot_1_dxy

# plot 2 - density plots of dxy per group
plot_2_dxy <- data %>%
  filter(data_type=='Dxy')%>%
  ggplot(aes(value, chr_type, fill=chr_type), guide="none") +
  geom_density_ridges(quantile_lines=TRUE, quantile_fun=function(x,...)median(x), size=0.5) +
  theme_bw() + theme(legend.position = "top", 
                     legend.title = element_blank(),
                     axis.text.y=element_blank(),
                     axis.title = element_text(size = 32),
                     axis.text = element_text(size = 16),
                     strip.text = element_text(size = 14),
                     legend.text = element_text(size = 20)) +
  facet_grid(comparison~.) +
  scale_fill_npg() +
  labs(x="Dxy", y="Density")
plot_2_dxy

# combine plots
dxy_plot <- plot_1_dxy + plot_2_dxy +  plot_layout(widths = c(5, 1))
ggsave("genomewide_and_density_dxy_region_DOG_ONLY.tif", width=16, height=18, dpi = 300)
ggsave("genomewide_and_density_dxy_region_DOG_ONLY.png", width=16, height=18, dpi = 300)
# higher Dxy = more divergent

#some additional plots
boxplot_dxy <- data %>%
  filter(data_type =='Dxy') %>% 
  ggplot(., aes(comparison, value, col=comparison)) +
  geom_jitter(size = 1, alpha = 0.5) +
  geom_boxplot(fill=NA, col="black") +
  labs(x = "Population" , y = "Dxy", colour = "Population") +
  theme_bw() +
  scale_color_npg() +
  theme(legend.position = "none",
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 18),
        strip.text = element_text(size = 18),
        legend.text = element_text(size = 16)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
boxplot_dxy
ggsave("boxplot_dxy_region_DOG_ONLY.tif", width=14, height = 6, dpi = 300)
ggsave("boxplot_dxy_region_DOG_ONLY.png", width=14, height = 6, dpi = 300)

density_dxy <- data %>%
  filter(data_type =='Dxy') %>% 
  ggplot(., aes(x=value, group = comparison, fill=comparison)) +
  geom_density(adjust=1.5, alpha=.4) +
  labs(x = "Dxy" , y = "Density", colour = "Comparison") +
  theme_bw() +
  scale_fill_npg() +
  theme(legend.position = c(0.8,0.72),
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 18),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
density_dxy
ggsave("density_dxy_region_DOG_ONLY.tif", width = 8, height = 6, dpi = 300)
ggsave("density_dxy_region_DOG_ONLY.png", width = 8, height = 6, dpi = 300)

# AUS_v_ASIA & AUS_v_EUR
scatter_dxy_a <-  data %>%
  filter(data_type == 'Dxy') %>%
  as.data.frame() %>%
  select('position', 'comparison', 'value', 'chromosome')%>%
  as_tibble() %>%
  pivot_wider(names_from = comparison, values_from = value) %>%
  ggplot(., aes(x=AUS_v_ASIA, y = AUS_v_EUR, col=chromosome)) +
  geom_point(alpha=.4) +
  labs(x = "AUS_v_ASIA" , y = "AUS_v_EUR", colour = "Comparison") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
scatter_dxy_a
# can repeat for other pairs

genome_pos_dxy_a <-  data %>%
  filter(data_type == 'Dxy') %>%
  as.data.frame() %>%
  select('position', 'comparison', 'value', 'chromosome')%>%
  as_tibble() %>%
  pivot_wider(names_from = comparison, values_from = value) %>%
  mutate("xx" = AUS_v_ASIA - AUS_v_EUR,
         "yy" = AUS_v_ASIA - EUR_v_ASIA) %>%
  ggplot(., aes(position*100000, xx, col=chromosome)) +
  geom_point(size=1) +
  scale_color_tron() +
  theme_bw() +
  theme(legend.position="top", 
        legend.title = element_blank()) +
  labs(x="Position", y="AUS vs ASIA - AUS vs EUR")
genome_pos_dxy_a
# can repeat for other pairs

genome_pos_dxy_aa <-  data %>%
  filter(data_type == 'Dxy') %>%
  as.data.frame() %>%
  select('position', 'comparison', 'value', 'chromosome')%>%
  as_tibble() %>%
  pivot_wider(names_from = comparison, values_from = value) %>%
  mutate("xx" = AUS_v_ASIA / AUS_v_EUR,
         "yy" = AUS_v_ASIA / EUR_v_ASIA) %>%
  ggplot(., aes(position*100000, xx, col=chromosome)) +
  geom_point(size=1) +
  scale_color_tron() +
  theme_bw() +
  theme(legend.position="top", 
        legend.title = element_blank()) +
  labs(x="Position", y="AUS vs ASIA / AUS vs EUR")
genome_pos_dxy_aa
# can repeat for other pairs

# ggarrange above plots if using them


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
  annotate(geom="text", x=140*100000, y=0.004, label="ChrX", size = 6)+
  annotate(geom="text", x=359.5*100000, y=0.004, label="Chr1", size = 6)+
  annotate(geom="text", x=520*100000, y=0.004, label="Chr2", size = 6)+
  annotate(geom="text", x=665*100000, y=0.004, label="Chr3", size = 6)+
  annotate(geom="text", x=810.5*100000, y=0.004, label="Chr4", size = 6)+
  geom_line() +
  scale_color_npg() +
  theme_bw() +
  theme(legend.position="top", 
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 14)) +
  labs(x="Genomic Position", y="Dxy") +
  guides(col = guide_legend(override.aes = list(linewidth = 2)))
dxy_lineplot
ggsave("lineplot_dxy_region_DOG_ONLY.tif", width=14, height=8, dpi = 300)
ggsave("lineplot_dxy_region_DOG_ONLY.png", width=14, height=8, dpi = 300)


# Dxy by region on a heatmap
# Get median dxy for each comparison
dxy_median <- data %>%
  filter(data_type == "Dxy") %>%
  group_by(comparison) %>%
  summarise(DXY_MEDIAN = median(value, na.rm = TRUE))

# Save as csv file
write.csv(dxy_median, "dxy_median.csv", row.names = TRUE)

dxy_heatmap <- read.csv("dxy_heatmap.csv", header = TRUE)
dxy_heatmap <- as.data.frame(dxy_heatmap)

plot_3_dxy <- ggplot(dxy_heatmap, aes(x = POP1, y=POP2, fill = DXY_MEDIAN)) + 
  geom_tile(color = "white", lwd = 1.5, linetype = 1) +
  coord_fixed() +
  scale_fill_material("indigo", labels = scales::label_number()) +
  guides(fill = guide_colourbar(title = expression(D[XY]), barwidth = 1, barheight = 5)) +
  theme_minimal() +
  labs(x = "Population", y = "Population") +
  theme(axis.text = element_text(size = 18),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title = element_text(size = 24),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank())
plot_3_dxy

ggsave("dxy_heatmap_region_DOG_ONLY.tif", bg = "white", height=5, width = 8, dpi = 300)
ggsave("dxy_heatmap_region_DOG_ONLY.png", bg = "white", height=5, width = 8, dpi = 300)


# Plotting Fst

# plot 1 - genome wide plots per population
plot_1_fst <- data %>%
  filter(data_type=='Fst')%>%
  ggplot(., aes(position*100000, value, col=chromosome)) +
  geom_point(size=1.5) +
  facet_grid(comparison~.) +
  scale_color_npg(guide = guide_legend(override.aes = list(size = 4))) +
  theme_bw() +
  theme(legend.position="top", 
        legend.title = element_blank(),
        axis.title = element_text(size = 32),
        axis.text = element_text(size = 16),
        strip.text = element_text(size = 14),
        legend.text = element_text(size = 20)) +
  labs(x="Genomic Position", y="Fst")
plot_1_fst


# plot 2 - density plots of Fst per group
plot_2_fst <- data %>%
  filter(data_type=='Fst')%>%
  ggplot(aes(value, chr_type, fill=chr_type), guide="none") +
  geom_density_ridges(quantile_lines=TRUE, quantile_fun=function(x,...)median(x), size=0.5) +
  theme_bw() + 
  theme(legend.position = "none", 
        axis.text.y=element_blank(),
        axis.title = element_text(size = 32),
        axis.text = element_text(size = 16),
        strip.text = element_text(size = 14),
        legend.text = element_text(size = 20)) +
  facet_grid(comparison~.) +
  scale_fill_npg() +
  labs(x="Fst", y="Density")
plot_2_fst

# combine plots
plot_1_fst + plot_2_fst +  plot_layout(widths = c(5, 1))
# 1: Removed 35 rows containing missing values or values outside the scale range (`geom_point()`). 
# 2: Removed 35 rows containing non-finite outside the scale range (`stat_density_ridges()`).
ggsave("genomewide_and_density_fst_region_DOG_ONLY.tif", width=16, height=18, dpi = 300)
ggsave("genomewide_and_density_fst_region_DOG_ONLY.png", width=16, height=18, dpi = 300)
# Fst closer to 1 = more genetically distinct
# Fst closer to 0 = more genetically similar


# Region map metadata. Just have one random point representing each region.
location <- read.csv("input/region_metadata.csv", header = TRUE)

## Make world map data
world_map <- map_data("world")

# Set colors for the points
red_palette1 <- brewer.pal(n = 9, name = "YlOrRd")
red_palette2 <- brewer.pal(n=9, name = "YlOrBr")
blue_palette <- brewer.pal(n = 9, name = "Blues")
green_palette1 <- brewer.pal(n = 9, name = "BuGn")
green_palette2 <- brewer.pal(n=9, name = "YlGn")

scale_colour_region <- 
  c('USA' = red_palette2[7],
    'CENAM' = "blueviolet",
    'EUR' = green_palette1[7],
    'ASIA' = "violetred1",
    'AUS' = blue_palette[5])
# chose one random location/point for each region

# Get median fst for each comparison
fst_median <- data %>%
  filter(data_type == "Fst") %>%
  group_by(comparison) %>%
  summarise(FST_MEDIAN = median(value, na.rm = TRUE))
# Save as csv file and then add required metadata such as population longitudes/latitudes
write.csv(fst_median, "fst_median.csv", row.names = TRUE)

# read new csv file with map metadata and median fst between groups
fst_map <- read.csv("fst_map.csv", header = TRUE)

# Fst by region on a map
plot_3_fst <- ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = "grey90") +
  geom_segment(data=fst_map, aes(x=POP1_LONG, y=POP1_LAT, xend=POP2_LONG, yend=POP2_LAT, color = FST_MEDIAN), size = 1) +
  scale_color_gradient(low = "gold", high = "red3", name = expression(F[ST])) +
  geom_point(data = location, aes(x = Longitude, y = Latitude, fill = Region), size = 6, shape = 21) +
  scale_fill_manual(values = scale_colour_region, limits = c("USA", "CENAM", "EUR", "ASIA", "AUS"), name = "Region") +
  geom_text_repel(data = location, aes(x = Longitude, y = Latitude, label = Region), size = 4, fontface = "bold", nudge_y = 7) +
  theme_void() +
  theme(
    legend.position = c(0.1,0.3),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)) +
  guides(fill = "none") +
  coord_fixed()
plot_3_fst

# Save plot
ggsave("fst_map_region_DOG_ONLY.tif", bg = "white", height=5, width=10, dpi = 300)
ggsave("fst_map_region_DOG_ONLY.png", bg = "white", height=5, width=10, dpi = 300)


# Fst by region on a heatmap
fst_heatmap <- read.csv("fst_heatmap.csv", header = TRUE)

plot_4_fst <- ggplot(fst_heatmap, aes(x = POP1, y=POP2, fill = FST_MEDIAN)) + 
  geom_tile(color = "white", lwd = 1.5, linetype = 1) +
  coord_fixed() +
  scale_fill_material("indigo") +
  guides(fill = guide_colourbar(title = expression(F[ST]), barwidth = 1, barheight = 5)) +
  theme_minimal() +
  labs(x = "Population", y = "Population") +
  theme(axis.text = element_text(size = 18),
    axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title = element_text(size = 24),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank())
plot_4_fst

ggsave("fst_heatmap_region_DOG_ONLY.tif", bg = "white", height=5, width = 8, dpi = 300)
ggsave("fst_heatmap_region_DOG_ONLY.png", bg = "white", height=5, width = 8, dpi = 300)
```



### Pixy for Australian dog samples only

```bash
conda activate pixy

module load PaM/environment
module load bsub.py/0.42.1
module load vcftools/0.1.16-c4
module load htslib-1.19/perl-5.38.0

WORKING_DIR=/lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA

cd ${WORKING_DIR}/03_ANALYSIS/04_VARIANTS/FILTER1/NO_OUTGROUPS/FINAL_SETS

# Select Australian dog samples only in vcf
vcftools --gzvcf nuclearSNPssandINVARIANTs_124samples.chrXto4.recode.vcf.gz \
--keep ../../../../05_ANALYSIS/PIXY/AUS_DOG_ONLY_samplelist.keep \
--recode --out AUS_DOG_ONLY_nuclearSNPsandINVARIANTs.chrXto4
# After filtering, kept 44 out of 124 Individuals
# After filtering, kept 84210689 out of a possible 84210689 Sites

# up to here
bgzip AUS_DOG_ONLY_nuclearSNPsandINVARIANTs.chrXto4.recode.vcf;
tabix -p vcf AUS_DOG_ONLY_nuclearSNPsandINVARIANTs.chrXto4.recode.vcf.gz

# run pixy
VCF=${WORKING_DIR}/03_ANALYSIS/04_VARIANTS/FILTER1/NO_OUTGROUPS/FINAL_SETS/AUS_DOG_ONLY_nuclearSNPsandINVARIANTs.chrXto4.recode.vcf.gz

cd ../../../../05_ANALYSIS/PIXY

bsub.py --queue long --threads 10 20 pixy_AUS_DOG_ONLY \
"pixy --stats pi fst dxy \
--vcf ${VCF} \
--populations AUS_DOG_ONLY_pop.list \
--window_size 100000 \
--n_cores 10 \
--output_prefix AUS_DOG_ONLY"
# Successfully completed
```

Plots in R

```R
# Batch 4 - Pixy for Pi, Fst and Dxy for Australian dog cohort only

library(tidyverse)
library(ggsci)
library(ggpubr)
library(patchwork)
library(ggridges)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(purrr)
library(stats)
library(graphics)
library(grDevices)
library(utils)
library(datasets)
library(methods)
library(base)
require(maps)
library(mapdata)
library(readxl)
library(ozmaps) 
library(grid)
library(gridExtra)
library(ggrepel)
library(ggnewscale)
library(reshape)
library(wesanderson)

setwd("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/batch4/FILTER1/NO_OUTGROUPS/pixy/AUS_DOG_ONLY")

##### Nucleotide diversity (Pi) #####

# get nucleotide diversity (pi) data from pixy output
pi_data <- read.table("input/AUS_DOG_ONLY_pi.txt", header=T)
pi_data$chromosome <- str_remove(pi_data$chromosome, 'dirofilaria_immitis_')

# filter pi data to remove small scaffolds not in the linkage groups, and to number the rows per group to help with plotting
pi_data <- pi_data %>%
  group_by(pop) %>%
  mutate(position = 1:n())

#Let's add the chr type variable
pi_data <- pi_data %>%
  mutate(chr_type = ifelse(str_detect(chromosome, "X"), "sexchr", "autosome"))

# calculate the median Pi and check the ratio of sex-to-autosome diversity. Should be about 0.75, as D. immitis is XX/XY
pi_data_sex_median <- pi_data %>%
  group_by(chr_type) %>%
  summarise(median = median(avg_pi, na.rm = TRUE))
pi_data_sex_median

'
# A tibble: 2 × 2
  chr_type   median
  <chr>       <dbl>
1 autosome 0.000479
2 sexchr   0.000128
'
# 0.000128 / 0.000479 = 0.2672234 (way off 0.75)

pi_data_pop_sex_median <- pi_data %>%
  group_by(pop, chr_type) %>%
  summarise(median = median(avg_pi, na.rm = TRUE))
pi_data_pop_sex_median

'
# A tibble: 12 × 3
# Groups:   pop [6]
   pop            chr_type    median
   <chr>          <chr>        <dbl>
 1 Brisbane       autosome 0.000433 
 2 Brisbane       sexchr   0.000106 
 3 Cairns         autosome 0.000451 
 4 Cairns         sexchr   0.0000692
 5 Lockhart_River autosome 0.000812 
 6 Lockhart_River sexchr   0.000294 
 7 Rockhampton    autosome 0.000215 
 8 Rockhampton    sexchr   0.0000312
 9 Sydney         autosome 0.000385 
10 Sydney         sexchr   0.000134 
11 Townsville     autosome 0.000636 
12 Townsville     sexchr   0.000220 
'

# plot 1 - genome wide plots per population
plot_1_pi <- ggplot(pi_data, aes(position*100000, avg_pi, col=chromosome)) +
  geom_point(size=0.8) +
  facet_grid(pop~.) +
  scale_color_npg(guide = guide_legend(override.aes = list(size = 4))) +
  theme_bw() +
  theme(legend.position="top", 
        legend.title = element_blank(),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 10),
        legend.text = element_text(size = 16)) +
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
  theme(legend.position="top",
        legend.title = element_blank(),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 10),
        legend.text = element_text(size = 16)) +
  labs(x="Nucleotide Diversity (Pi)", y="Density")
plot_2_pi

# combine plots
plot_1_pi + plot_2_pi +  plot_layout(widths = c(5, 1))
ggsave("genomewide_and_density_Pi_AUS_DOG_ONLY.tif", width=14, height=8, dpi = 300)
ggsave("genomewide_and_density_Pi_AUS_DOG_ONLY.png", width=14, height=8, dpi = 300)

#Now a boxplot of the pi value per population
#Now a boxplot of the pi value per population
pi_data$pop <- factor(pi_data$pop, 
                      levels = c('Lockhart_River', 'Cairns', 'Townsville', 'Rockhampton', 
                                 'Brisbane', 'Sydney'))

blue_palette <- brewer.pal(n = 9, name = "Blues")

scale_colour_AUS <- function(...){
  ggplot2:::manual_scale(
    'colour', 
    values = setNames(
      c(blue_palette[9],
        blue_palette[8],
        blue_palette[7],
        blue_palette[6],
        blue_palette[5],
        blue_palette[4]),
    c("Lockhart_River", "Cairns", "Townsville", "Rockhampton", "Brisbane", "Sydney")), 
  ...
)
}    


boxplot_pi <- ggplot(pi_data, aes(pop, avg_pi, col=pop)) +
  geom_jitter(size = 1, alpha = 0.5) +
  geom_boxplot(fill=NA, col="black") +
  labs(x = "Population" , y = "Nucleotide diversity (Pi)", colour = "Population") +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 18)
        ) +
  scale_colour_AUS ()

boxplot_pi

ggsave("boxplot_Pi_AUS_DOG_ONLY.tif", width=6, height=6, dpi = 300)
ggsave("boxplot_Pi_AUS_DOG_ONLY.png", width=6, height=6, dpi = 300)
# higher pi value = more diverse population
# lower pi value = less diverse population
## Diversity is pretty similar between cities. Sydney and Brisbane are on par, Townsville is slightly higher.




#### Dxy and Fst ####

# load data
dxy_data <- read.table("input/AUS_DOG_ONLY_dxy.txt", header=T)
fst_data <- read.table("input/AUS_DOG_ONLY_fst.txt", header=T)

# add some columns
dxy_data$data_type <- "Dxy"
dxy_data <- mutate(dxy_data,
                   comparison = paste(pop1, pop2, sep = '_v_'))
fst_data$data_type <- "Fst"
fst_data <- mutate(fst_data,
                   comparison = paste(pop1, pop2, sep = '_v_'))

# subset and merge dataframes
dxy_data_sub <- dxy_data %>% select(comparison,data_type,chromosome,window_pos_1,window_pos_2,avg_dxy)
colnames(dxy_data_sub) <- c("comparison","data_type","chromosome","window_pos_1","window_pos_2","value")

fst_data_sub <- fst_data %>% select(comparison,data_type,chromosome,window_pos_1,window_pos_2,avg_wc_fst)
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
  summarise(median = median(value, na.rm = TRUE)) %>%
  print(n = 30)
'
# A tibble: 30 × 3
# Groups:   comparison [15]
   comparison                   data_type    median
   <chr>                        <chr>         <dbl>
 1 Brisbane_v_Cairns            Dxy        0.000487
 2 Brisbane_v_Cairns            Fst        0.231   
 3 Brisbane_v_Lockhart_River    Dxy        0.000503
 4 Brisbane_v_Lockhart_River    Fst        0.0798  
 5 Brisbane_v_Rockhampton       Dxy        0.000444
 6 Brisbane_v_Rockhampton       Fst        0.104   
 7 Brisbane_v_Sydney            Dxy        0.000471
 8 Brisbane_v_Sydney            Fst        0.219   
 9 Brisbane_v_Townsville        Dxy        0.000468
10 Brisbane_v_Townsville        Fst        0.0626  
11 Cairns_v_Lockhart_River      Dxy        0.000449
12 Cairns_v_Lockhart_River      Fst        0.164   
13 Cairns_v_Rockhampton         Dxy        0.000464
14 Cairns_v_Rockhampton         Fst        0.396   
15 Cairns_v_Sydney              Dxy        0.000443
16 Cairns_v_Sydney              Fst        0.144   
17 Cairns_v_Townsville          Dxy        0.000481
18 Cairns_v_Townsville          Fst        0.0615  
19 Lockhart_River_v_Rockhampton Dxy        0.000492
20 Lockhart_River_v_Rockhampton Fst        0       
21 Lockhart_River_v_Sydney      Dxy        0.000479
22 Lockhart_River_v_Sydney      Fst       -0.0188  
23 Lockhart_River_v_Townsville  Dxy        0.000507
24 Lockhart_River_v_Townsville  Fst       -0.0327  
25 Rockhampton_v_Sydney         Dxy        0.000412
26 Rockhampton_v_Sydney         Fst        0.0468  
27 Rockhampton_v_Townsville     Dxy        0.000488
28 Rockhampton_v_Townsville     Fst        0.0849  
29 Sydney_v_Townsville          Dxy        0.000495
30 Sydney_v_Townsville          Fst        0.130 
'


#Plotting Dxy

# plot 1 - genome wide plots per comparison
plot_1_dxy <- data %>%
  filter(data_type =='Dxy')%>%
  ggplot(., aes(position*100000, value, col=chromosome)) +
  geom_point(size=1.5) +
  facet_grid(comparison~.) +
  scale_color_npg(guide = guide_legend(override.aes = list(size = 4))) +
  theme_bw() +
  theme(legend.position="top", 
        legend.title = element_blank(),
        axis.title = element_text(size = 32),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 20)) +
  labs(x="Genomic Position", y="Dxy")
plot_1_dxy

# plot 2 - density plots of dxy per group
plot_2_dxy <- data %>%
  filter(data_type=='Dxy')%>%
  ggplot(aes(value, chr_type, fill=chr_type), guide="none") +
  geom_density_ridges(quantile_lines=TRUE, quantile_fun=function(x,...)median(x), size=0.5) +
  theme_bw() + theme(legend.position = "top", 
                     legend.title = element_blank(),
                     axis.text.y=element_blank(),
                     axis.title = element_text(size = 32),
                     axis.text = element_text(size = 16),
                     strip.text = element_text(size = 8),
                     legend.text = element_text(size = 20)) +
  facet_grid(comparison~.) +
  scale_fill_npg() +
  labs(x="Dxy", y="Density")
plot_2_dxy

# combine plots
dxy_plot <- plot_1_dxy + plot_2_dxy +  plot_layout(widths = c(5, 1))
ggsave("genomewide_and_density_dxy_AUS_DOG_ONLY.tif", width=20, height=25, dpi = 300)
ggsave("genomewide_and_density_dxy_AUS_DOG_ONLY.png", width=20, height=25, dpi = 300)
# higher Dxy = more divergent

#some additional plots
boxplot_dxy <- data %>%
  filter(data_type =='Dxy') %>% 
  ggplot(., aes(comparison, value, col=comparison)) +
  geom_jitter(size = 1, alpha = 0.5) +
  geom_boxplot(fill=NA, col="black") +
  labs(x = "Population" , y = "Dxy", colour = "Population") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title = element_text(size =24),
        axis.text = element_text(size = 14),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
boxplot_dxy
ggsave("boxplot_dxy_AUS_DOG_ONLY.tif", width=14, height = 6, dpi = 300)
ggsave("boxplot_dxy_AUS_DOG_ONLY.png", width=14, height = 6, dpi = 300)

density_dxy <- data %>%
  filter(data_type =='Dxy') %>% 
  ggplot(., aes(x=value, group = comparison, fill=comparison)) +
  geom_density(adjust=1.5, alpha=.4) +
  labs(x = "Dxy" , y = "Density", colour = "Comparison") +
  theme_bw() +
  theme(legend.position = c(0.78, 0.6),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 18),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
density_dxy
ggsave("density_dxy_AUS_DOG_ONLY.tif", width = 8, height = 6, dpi = 300)
ggsave("density_dxy_AUS_DOG_ONLY.png", width = 8, height = 6, dpi = 300)


#Trying to plot lines
#getting chr positions
chr <- c('X', '1', '2', '3', '4')

for (i in chr) {
  x <- data %>%
    filter(chromosome == paste0('chr', i))
  print(quantile(x$position))
}

data$comparison <- str_replace_all(data$comparison, '_', ' ')

dxy_lineplot <- data %>%
  filter(data_type=='Dxy')%>%
  ggplot(., aes(position*100000, value, col=comparison)) +
  geom_point(size=0.2, alpha = 0.1) +
  geom_vline(xintercept=c(281*100000,
                          438*100000,
                          590*100000,
                          740*100000),
             size=1, linetype="dashed", col='grey41')+
  annotate(geom="text", x=140*100000, y=0.004, label="ChrX", size = 6)+
  annotate(geom="text", x=359.5*100000, y=0.004, label="Chr1", size = 6)+
  annotate(geom="text", x=520*100000, y=0.004, label="Chr2", size = 6)+
  annotate(geom="text", x=665*100000, y=0.004, label="Chr3", size = 6)+
  annotate(geom="text", x=810.5*100000, y=0.004, label="Chr4", size = 6)+
  geom_line() +
  theme_bw() +
  theme(legend.position="top", 
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 18)) +
  labs(x="Genomic Position", y="Dxy") +
  guides(col = guide_legend(override.aes = list(linewidth = 2)))
dxy_lineplot
ggsave("lineplot_dxy_AUS_DOG_ONLY.tif", width=16, height=8, dpi = 300)
ggsave("lineplot_dxy_AUS_DOG_ONLY.png", width=16, height=8, dpi = 300)



# Dxy by city on a heatmap
# Get median dxy for each comparison
dxy_median <- data %>%
  filter(data_type == "Dxy") %>%
  group_by(comparison) %>%
  summarise(DXY_MEDIAN = median(value, na.rm = TRUE))

# Save as csv file
write.csv(dxy_median, "dxy_median.csv", row.names = TRUE)

dxy_heatmap <- read.csv("dxy_heatmap.csv", header = TRUE)
dxy_heatmap <- as.data.frame(dxy_heatmap)

plot_3_dxy <- ggplot(dxy_heatmap, aes(x = POP1, y=POP2, fill = DXY_MEDIAN)) + 
  geom_tile(color = "white", lwd = 1.5, linetype = 1) +
  coord_fixed() +
  scale_fill_material("indigo", labels = scales::label_number()) +
  guides(fill = guide_colourbar(title = expression(D[XY]), barwidth = 1, barheight = 5)) +
  theme_minimal() +
  labs(x = "Population", y = "Population") +
  theme(axis.text = element_text(size = 18),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title = element_text(size = 24),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank())
plot_3_dxy

ggsave("dxy_heatmap_AUS_DOG_ONLY.tif", bg = "white", height=5, width = 8, dpi = 300)
ggsave("dxy_heatmap_AUS_DOG_ONLY.png", bg = "white", height=5, width = 8, dpi = 300)




# Plotting Fst

# plot 1 - genome wide plots per population
plot_1_fst <- data %>%
  filter(data_type=='Fst')%>%
  ggplot(., aes(position*100000, value, col=chromosome)) +
  geom_point(size=1.5) +
  facet_grid(comparison~.) +
  scale_color_npg(guide = guide_legend(override.aes = list(size = 4))) +
  theme_bw() +
  theme(legend.position="top", 
        legend.title = element_blank(),
        axis.title = element_text(size = 32),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 20)) +
  labs(x="Genomic Position", y="Fst")
plot_1_fst


# plot 2 - density plots of Fst per group
plot_2_fst <- data %>%
  filter(data_type=='Fst')%>%
  ggplot(aes(value, chr_type, fill=chr_type), guide="none") +
  geom_density_ridges(quantile_lines=TRUE, quantile_fun=function(x,...)median(x), size=0.5) +
  theme_bw() + 
  theme(legend.position = "none", 
        axis.text.y=element_blank(),
        axis.title = element_text(size = 32),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 20)) +
  facet_grid(comparison~.) +
  scale_fill_npg() +
  labs(x="Fst", y="Density")
plot_2_fst

# combine plots
plot_1_fst + plot_2_fst +  plot_layout(widths = c(5, 1))
# 1: Removed 224 rows containing missing values or values outside the scale range (`geom_point()`). 
# 2: Removed 224 rows containing non-finite outside the scale range (`stat_density_ridges()`). 
## Fst values should only range between 0 and 1. Set anything <0 to 0 and re-do the plots.

# Data cleaning
data_fst_clean <- data %>% filter(data_type == 'Fst') %>% mutate(value = ifelse(value < 0, 0, value)) # this takes my original data, filters it to only include rows where the data_type is Fst. Then mutate function is used to create/modify columns. It modifies the value column, checks if the value is less than 0. If the value is <0, it replaces it with 0. If the value is >0, it keeps the original value.

# Re-do plots with cleaned data
# plot 1 - genome wide plots per population
plot_1_fst_clean <- data_fst_clean %>%
  filter(data_type=='Fst')%>%
  ggplot(., aes(position*100000, value, col=chromosome)) +
  geom_point(size=1) +
  facet_grid(comparison~.) +
  scale_color_tron() +
  theme_bw() +
  theme(legend.position="top", 
        legend.title = element_blank(),
        axis.title = element_text(size = 32),
        axis.text = element_text(size = 16),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 20)) +
  labs(x="Genomic Position", y="Fst")
plot_1_fst_clean


# plot 2 - density plots of Fst per group
plot_2_fst_clean <- data_fst_clean %>%
  filter(data_type=='Fst')%>%
  ggplot(aes(value, chr_type, fill=chr_type), guide="none") +
  geom_density_ridges(quantile_lines=TRUE, quantile_fun=function(x,...)median(x), size=0.5) +
  theme_bw() + theme(legend.position = "none", axis.text.y=element_blank()) +
  facet_grid(comparison~.) +
  scale_fill_npg() +
  labs(x="Fst", y="Density")
plot_2_fst_clean

# combine plots
plot_1_fst_clean + plot_2_fst_clean +  plot_layout(widths = c(5, 1))
# 1: Removed 224 rows containing missing values or values outside the scale range (`geom_point()`). 
# 2: Removed 224 rows containing non-finite outside the scale range (`stat_density_ridges()`).
ggsave("genomewide_and_density_fst_AUS_DOG_ONLY.tif", width=20, height=25, dpi = 300)
ggsave("genomewide_and_density_fst_AUS_DOG_ONLY.png", width=20, height=25, dpi = 300)
# Fst closer to 1 = more genetically distinct
# Fst closer to 0 = more genetically similar
# Some cities only have 1 sample so they look funny.



# AUS map

## Make world map data
world_map <- map_data("world")

# Add map inset to zoom in on Australian dog samples
# Manually specify the coordinates for the area of the world map to show in the inset
aus_xmin <- 141
aus_xmax <- 155
aus_ymin <- -60
aus_ymax <- -9

# Filter the world map data for the inset area
aus_data <- subset(world_map, long >= aus_xmin & long <= aus_xmax & lat >= aus_ymin & lat <= aus_ymax)

# Aus metadata
location_AUS <- read.csv("input/AUS_metadata.csv", header = TRUE)

# Set colors for the points
blue_palette <- brewer.pal(n = 9, name = "Blues")

scale_colour_AUS <- 
  c('Lockhart River' = blue_palette[9],
    'Cairns' = blue_palette[8],
    'Townsville' = blue_palette[7],
    'Rockhampton' = blue_palette[6],
    'Brisbane' = blue_palette[5],
    'Sydney' = blue_palette[4])

# Get median fst for each comparison
fst_median <- data_fst_clean %>%
  filter(data_type == "Fst") %>%
  group_by(comparison) %>%
  summarise(FST_MEDIAN = median(value, na.rm = TRUE))
# Save as csv file and then add required metadata such as population longitudes/latitudes
write.csv(fst_median, "fst_median.csv", row.names = TRUE)

# read new csv file with map metadata and median fst between groups
fst_map <- read.csv("fst_map.csv", header = TRUE)

# Fst on AUS map
plot_3_fst <- ggplot() +
  geom_polygon(data = aus_data, aes(x = long, y = lat, group = group), fill = "grey90") +
  geom_segment(data=fst_map, aes(x=POP1_LONG, y=POP1_LAT, xend=POP2_LONG, yend=POP2_LAT, color = FST_MEDIAN), size = 1) +
  scale_color_gradient(low = "gold", high = "red3", name = expression(F[ST])) +
  geom_point(data = location_AUS, aes(x = Longitude, y = Latitude, fill = City), size = 10, shape = 21) +
  scale_fill_manual(values = scale_colour_AUS, limits = c("Lockhart River", "Cairns", "Townsville", "Rockhampton", "Brisbane", "Sydney"), name = "City") +
  geom_text_repel(data = location_AUS, aes(x = Longitude, y = Latitude, label = City), size = 6, fontface = "bold", nudge_x = 5, nudge_y = 2) +
  theme_void() +
  theme(
    legend.position = c(0.9,0.1),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 14)) +
  guides(fill = "none") +
  coord_fixed(ratio=1)
plot_3_fst

# Save plot
ggsave("fst_map_AUS_DOG_ONLY.tif", bg = "white", height=10, width=5, dpi = 300)
ggsave("fst_map_AUS_DOG_ONLY.png", bg = "white", height=10, width=5, dpi = 300)


# Fst by AUS on a heatmap
fst_heatmap <- read.csv("fst_heatmap.csv", header = TRUE)

plot_4_fst <- ggplot(fst_heatmap, aes(x = POP1, y=POP2, fill = FST_MEDIAN)) + 
  geom_tile(color = "white", lwd = 1.5, linetype = 1) +
  coord_fixed() +
  scale_fill_material("indigo") +
  guides(fill = guide_colourbar(title = expression(F[ST]), barwidth = 1, barheight = 5)) +
  theme_minimal() +
  labs(x = "Population", y = "Population") +
  theme(axis.text = element_text(size = 18),
    axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title = element_text(size = 24),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank())
plot_4_fst

ggsave("fst_heatmap_AUS_DOG_ONLY.tif", bg = "white", height=5, width = 8, dpi = 300)
ggsave("fst_heatmap_AUS_DOG_ONLY.png", bg = "white", height=5, width = 8, dpi = 300)


# Plot comparing genetic vs geographical distance
## Calculated direct distances between cities using: https://www.nhc.noaa.gov/gccalc.shtml
fst_distance <- read.csv("fst_distance.csv", header = TRUE)

plot_5_fst <- ggplot(fst_distance, aes(x = Distance_km, y = FST_MEDIAN)) +
  geom_point(size = 4, col = "cornflowerblue") +
  labs(x = "Distance (km)", y = expression(F[ST])) +
  theme_bw() +
  theme(axis.title = element_text(size = 24),
        axis.text = element_text(size = 18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
plot_5_fst

ggsave("fst_distance_AUS_DOG_ONLY.tif", bg = "white", height = 6, width = 8, dpi = 300)
ggsave("fst_distance_AUS_DOG_ONLY.png", bg = "white", height = 6, width = 8, dpi = 300)
```




### Pixy for USA dog samples only

```bash
conda activate pixy

module load PaM/environment
module load bsub.py/0.42.1
module load vcftools/0.1.16-c4
module load htslib-1.19/perl-5.38.0

WORKING_DIR=/lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA

cd ${WORKING_DIR}/03_ANALYSIS/04_VARIANTS/FILTER1/NO_OUTGROUPS/FINAL_SETS

# Select USA dog samples only in vcf
vcftools --gzvcf nuclearSNPssandINVARIANTs_124samples.chrXto4.recode.vcf.gz \
--keep ../../../../05_ANALYSIS/PIXY/USA_DOG_ONLY_samplelist.keep \
--recode --out USA_DOG_ONLY_nuclearSNPsandINVARIANTs.chrXto4
# After filtering, kept 34 out of 124 Individuals
# After filtering, kept 84210689 out of a possible 84210689 Sites

bgzip USA_DOG_ONLY_nuclearSNPsandINVARIANTs.chrXto4.recode.vcf;
tabix -p vcf USA_DOG_ONLY_nuclearSNPsandINVARIANTs.chrXto4.recode.vcf.gz

# run pixy
VCF=${WORKING_DIR}/03_ANALYSIS/04_VARIANTS/FILTER1/NO_OUTGROUPS/FINAL_SETS/USA_DOG_ONLY_nuclearSNPsandINVARIANTs.chrXto4.recode.vcf.gz

cd ../../../../05_ANALYSIS/PIXY

bsub.py --queue long --threads 10 20 pixy_USA_DOG_ONLY \
"pixy --stats pi fst dxy \
--vcf ${VCF} \
--populations USA_DOG_ONLY_pop.list \
--window_size 100000 \
--n_cores 10 \
--output_prefix USA_DOG_ONLY"
# Successfully completed
```

```R
# Batch 4 - Pixy for Pi, Fst and Dxy for USA dog cohort only

library(tidyverse)
library(ggsci)
library(ggpubr)
library(patchwork)
library(ggridges)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(purrr)
library(stats)
library(graphics)
library(grDevices)
library(utils)
library(datasets)
library(methods)
library(base)
require(maps)
library(mapdata)
library(readxl)
library(ozmaps) 
library(grid)
library(gridExtra)
library(ggrepel)
library(ggnewscale)
library(reshape)
library(wesanderson)

setwd("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/batch4/FILTER1/NO_OUTGROUPS/pixy/USA_DOG_ONLY")

##### Nucleotide diversity (Pi) #####

# get nucleotide diversity (pi) data from pixy output
pi_data <- read.table("input/USA_DOG_ONLY_pi.txt", header=T)
pi_data$chromosome <- str_remove(pi_data$chromosome, 'dirofilaria_immitis_')

# filter pi data to remove small scaffolds not in the linkage groups, and to number the rows per group to help with plotting
pi_data <- pi_data %>%
  group_by(pop) %>%
  mutate(position = 1:n())

#Let's add the chr type variable
pi_data <- pi_data %>%
  mutate(chr_type = ifelse(str_detect(chromosome, "X"), "sexchr", "autosome"))

# calculate the median Pi and check the ratio of sex-to-autosome diversity. Should be about 0.75, as D. immitis is XX/XY
pi_data_sex_median <- pi_data %>%
  group_by(chr_type) %>%
  summarise(median = median(avg_pi, na.rm = TRUE))
pi_data_sex_median

'
# A tibble: 2 × 2
  chr_type   median
  <chr>       <dbl>
1 autosome 0.000652
2 sexchr   0.000257
'
# 0.000257 / 0.000652 = 0.3941718 (way off 0.75)

pi_data_pop_sex_median <- pi_data %>%
  group_by(pop, chr_type) %>%
  summarise(median = median(avg_pi, na.rm = TRUE))
pi_data_pop_sex_median

'
# A tibble: 12 × 3
# Groups:   pop [6]
   pop       chr_type   median
   <chr>     <chr>       <dbl>
 1 Florida   autosome 0.000773
 2 Florida   sexchr   0.000341
 3 Georgia   autosome 0.000732
 4 Georgia   sexchr   0.000330
 5 Illinois  autosome 0.000585
 6 Illinois  sexchr   0.000178
 7 Louisiana autosome 0.000599
 8 Louisiana sexchr   0.000219
 9 Missouri  autosome 0.000516
10 Missouri  sexchr   0.000211
11 Texas     autosome 0.000722
12 Texas     sexchr   0.000305
'

# plot 1 - genome wide plots per population
plot_1_pi <- ggplot(pi_data, aes(position*100000, avg_pi, col=chromosome)) +
  geom_point(size=0.8) +
  facet_grid(pop~.) +
  scale_color_npg(guide = guide_legend(override.aes = list(size = 4))) +
  theme_bw() +
  theme(legend.position="top", 
        legend.title = element_blank(),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 14),
        legend.text = element_text(size = 16)) +
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
  theme(legend.position="top",
        legend.title = element_blank(),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 14),
        legend.text = element_text(size = 16)) +
  labs(x="Nucleotide Diversity (Pi)", y="Density")
plot_2_pi

# combine plots
plot_1_pi + plot_2_pi +  plot_layout(widths = c(5, 1))
ggsave("genomewide_and_density_Pi_USA_DOG_ONLY.tif", width=14, height=8, dpi = 300)
ggsave("genomewide_and_density_Pi_USA_DOG_ONLY.png", width=14, height=8, dpi = 300)

#Now a boxplot of the pi value per population
#Now a boxplot of the pi value per population
pi_data$pop <- factor(pi_data$pop, 
                      levels = c('Florida', 'Georgia', 'Illinois', 'Louisiana', 
                                 'Missouri', 'Texas'))

red_palette1 <- brewer.pal(n = 9, name = "YlOrRd")
red_palette2 <- brewer.pal(n=9, name = "YlOrBr")

scale_colour_USA <- function(...){
  ggplot2:::manual_scale(
    'colour', 
    values = setNames(
      c(red_palette1[5],
        red_palette2[5],
        red_palette1[8],
        red_palette1[6],
        red_palette2[7],
        red_palette1[7]),
    c("Florida", "Georgia", "Illinois", "Louisiana", "Missouri", "Texas")), 
  ...
)
}    


boxplot_pi <- ggplot(pi_data, aes(pop, avg_pi, col=pop)) +
  geom_jitter(size = 1, alpha = 0.5) +
  geom_boxplot(fill=NA, col="black") +
  labs(x = "Population" , y = "Nucleotide diversity (Pi)", colour = "Population") +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 18)
        ) +
  scale_colour_USA ()

boxplot_pi

ggsave("boxplot_Pi_USA_DOG_ONLY.tif", width=6, height=6, dpi = 300)
ggsave("boxplot_Pi_USA_DOG_ONLY.png", width=6, height=6, dpi = 300)
# higher pi value = more diverse population
# lower pi value = less diverse population
## Diversity slightly higher in Florida, Georgia & Texas.




#### Dxy and Fst ####

# load data
dxy_data <- read.table("input/USA_DOG_ONLY_dxy.txt", header=T)
fst_data <- read.table("input/USA_DOG_ONLY_fst.txt", header=T)

# add some columns
dxy_data$data_type <- "Dxy"
dxy_data <- mutate(dxy_data,
                   comparison = paste(pop1, pop2, sep = '_v_'))
fst_data$data_type <- "Fst"
fst_data <- mutate(fst_data,
                   comparison = paste(pop1, pop2, sep = '_v_'))

# subset and merge dataframes
dxy_data_sub <- dxy_data %>% select(comparison,data_type,chromosome,window_pos_1,window_pos_2,avg_dxy)
colnames(dxy_data_sub) <- c("comparison","data_type","chromosome","window_pos_1","window_pos_2","value")

fst_data_sub <- fst_data %>% select(comparison,data_type,chromosome,window_pos_1,window_pos_2,avg_wc_fst)
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
  summarise(median = median(value, na.rm = TRUE)) %>%
  print(n = 30)
'
# A tibble: 30 × 3
# Groups:   comparison [15]
   comparison           data_type   median
   <chr>                <chr>        <dbl>
 1 Florida_v_Georgia    Dxy       0.000569
 2 Florida_v_Georgia    Fst       0.0157  
 3 Florida_v_Illinois   Dxy       0.000610
 4 Florida_v_Illinois   Fst       0.0784  
 5 Florida_v_Louisiana  Dxy       0.000607
 6 Florida_v_Louisiana  Fst       0.0601  
 7 Florida_v_Missouri   Dxy       0.000615
 8 Florida_v_Missouri   Fst       0.0910  
 9 Florida_v_Texas      Dxy       0.000627
10 Florida_v_Texas      Fst       0.0181  
11 Georgia_v_Illinois   Dxy       0.000513
12 Georgia_v_Illinois   Fst       0.169   
13 Georgia_v_Louisiana  Dxy       0.000533
14 Georgia_v_Louisiana  Fst       0.118   
15 Georgia_v_Missouri   Dxy       0.000532
16 Georgia_v_Missouri   Fst       0.225   
17 Georgia_v_Texas      Dxy       0.000540
18 Georgia_v_Texas      Fst       0.000749
19 Illinois_v_Louisiana Dxy       0.000520
20 Illinois_v_Louisiana Fst       0.149   
21 Illinois_v_Missouri  Dxy       0.000560
22 Illinois_v_Missouri  Fst       0.237   
23 Illinois_v_Texas     Dxy       0.000567
24 Illinois_v_Texas     Fst       0.0474  
25 Louisiana_v_Missouri Dxy       0.000514
26 Louisiana_v_Missouri Fst       0.221   
27 Louisiana_v_Texas    Dxy       0.000551
28 Louisiana_v_Texas    Fst       0.0296  
29 Missouri_v_Texas     Dxy       0.000574
30 Missouri_v_Texas     Fst       0.0856
'


#Plotting Dxy

# plot 1 - genome wide plots per comparison
plot_1_dxy <- data %>%
  filter(data_type =='Dxy')%>%
  ggplot(., aes(position*100000, value, col=chromosome)) +
  geom_point(size=1.5) +
  facet_grid(comparison~.) +
  scale_color_npg(guide = guide_legend(override.aes = list(size = 4))) +
  theme_bw() +
  theme(legend.position="top", 
        legend.title = element_blank(),
        axis.title = element_text(size = 32),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 10),
        legend.text = element_text(size = 20)) +
  labs(x="Genomic Position", y="Dxy")
plot_1_dxy

# plot 2 - density plots of dxy per group
plot_2_dxy <- data %>%
  filter(data_type=='Dxy')%>%
  ggplot(aes(value, chr_type, fill=chr_type), guide="none") +
  geom_density_ridges(quantile_lines=TRUE, quantile_fun=function(x,...)median(x), size=0.5) +
  theme_bw() + theme(legend.position = "top", 
                     legend.title = element_blank(),
                     axis.text.y=element_blank(),
                     axis.title = element_text(size = 32),
                     axis.text = element_text(size = 16),
                     strip.text = element_text(size = 10),
                     legend.text = element_text(size = 20)) +
  facet_grid(comparison~.) +
  scale_fill_npg() +
  labs(x="Dxy", y="Density")
plot_2_dxy

# combine plots
dxy_plot <- plot_1_dxy + plot_2_dxy +  plot_layout(widths = c(5, 1))
ggsave("genomewide_and_density_dxy_USA_DOG_ONLY.tif", width=20, height=25, dpi = 300)
ggsave("genomewide_and_density_dxy_USA_DOG_ONLY.png", width=20, height=25, dpi = 300)
# higher Dxy = more divergent

#some additional plots
boxplot_dxy <- data %>%
  filter(data_type =='Dxy') %>% 
  ggplot(., aes(comparison, value, col=comparison)) +
  geom_jitter(size = 1, alpha = 0.5) +
  geom_boxplot(fill=NA, col="black") +
  labs(x = "Population" , y = "Dxy", colour = "Population") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title = element_text(size =24),
        axis.text = element_text(size = 14),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
boxplot_dxy
ggsave("boxplot_dxy_USA_DOG_ONLY.tif", width=14, height = 6, dpi = 300)
ggsave("boxplot_dxy_USA_DOG_ONLY.png", width=14, height = 6, dpi = 300)

density_dxy <- data %>%
  filter(data_type =='Dxy') %>% 
  ggplot(., aes(x=value, group = comparison, fill=comparison)) +
  geom_density(adjust=1.5, alpha=.4) +
  labs(x = "Dxy" , y = "Density", colour = "Comparison") +
  theme_bw() +
  theme(legend.position = c(0.78, 0.7),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 18),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
density_dxy
ggsave("density_dxy_USA_DOG_ONLY.tif", width = 8, height = 8, dpi = 300)
ggsave("density_dxy_USA_DOG_ONLY.png", width = 8, height = 8, dpi = 300)


#Trying to plot lines
#getting chr positions
chr <- c('X', '1', '2', '3', '4')

for (i in chr) {
  x <- data %>%
    filter(chromosome == paste0('chr', i))
  print(quantile(x$position))
}

data$comparison <- str_replace_all(data$comparison, '_', ' ')

dxy_lineplot <- data %>%
  filter(data_type=='Dxy')%>%
  ggplot(., aes(position*100000, value, col=comparison)) +
  geom_point(size=0.2, alpha = 0.1) +
  geom_vline(xintercept=c(281*100000,
                          438*100000,
                          590*100000,
                          740*100000),
             size=1, linetype="dashed", col='grey41')+
  annotate(geom="text", x=140*100000, y=0.004, label="ChrX", size = 6)+
  annotate(geom="text", x=359.5*100000, y=0.004, label="Chr1", size = 6)+
  annotate(geom="text", x=520*100000, y=0.004, label="Chr2", size = 6)+
  annotate(geom="text", x=665*100000, y=0.004, label="Chr3", size = 6)+
  annotate(geom="text", x=810.5*100000, y=0.004, label="Chr4", size = 6)+
  geom_line() +
  theme_bw() +
  theme(legend.position="top", 
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 18)) +
  labs(x="Genomic Position", y="Dxy") +
  guides(col = guide_legend(override.aes = list(linewidth = 2)))
dxy_lineplot
ggsave("lineplot_dxy_USA_DOG_ONLY.tif", width=16, height=8, dpi = 300)
ggsave("lineplot_dxy_USA_DOG_ONLY.png", width=16, height=8, dpi = 300)



# Dxy by city on a heatmap
# Get median dxy for each comparison
dxy_median <- data %>%
  filter(data_type == "Dxy") %>%
  group_by(comparison) %>%
  summarise(DXY_MEDIAN = median(value, na.rm = TRUE))

# Save as csv file
write.csv(dxy_median, "dxy_median.csv", row.names = TRUE)

dxy_heatmap <- read.csv("dxy_heatmap.csv", header = TRUE)
dxy_heatmap <- as.data.frame(dxy_heatmap)

plot_3_dxy <- ggplot(dxy_heatmap, aes(x = POP1, y=POP2, fill = DXY_MEDIAN)) + 
  geom_tile(color = "white", lwd = 1.5, linetype = 1) +
  coord_fixed() +
  scale_fill_material("indigo", labels = scales::label_number()) +
  guides(fill = guide_colourbar(title = expression(D[XY]), barwidth = 1, barheight = 5)) +
  theme_minimal() +
  labs(x = "Population", y = "Population") +
  theme(axis.text = element_text(size = 18),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title = element_text(size = 24),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank())
plot_3_dxy

ggsave("dxy_heatmap_USA_DOG_ONLY.tif", bg = "white", height=5, width = 8, dpi = 300)
ggsave("dxy_heatmap_USA_DOG_ONLY.png", bg = "white", height=5, width = 8, dpi = 300)




# Plotting Fst

# plot 1 - genome wide plots per population
plot_1_fst <- data %>%
  filter(data_type=='Fst')%>%
  ggplot(., aes(position*100000, value, col=chromosome)) +
  geom_point(size=1.5) +
  facet_grid(comparison~.) +
  scale_color_npg(guide = guide_legend(override.aes = list(size = 4))) +
  theme_bw() +
  theme(legend.position="top", 
        legend.title = element_blank(),
        axis.title = element_text(size = 32),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 20)) +
  labs(x="Genomic Position", y="Fst")
plot_1_fst


# plot 2 - density plots of Fst per group
plot_2_fst <- data %>%
  filter(data_type=='Fst')%>%
  ggplot(aes(value, chr_type, fill=chr_type), guide="none") +
  geom_density_ridges(quantile_lines=TRUE, quantile_fun=function(x,...)median(x), size=0.5) +
  theme_bw() + 
  theme(legend.position = "none", 
        axis.text.y=element_blank(),
        axis.title = element_text(size = 32),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 20)) +
  facet_grid(comparison~.) +
  scale_fill_npg() +
  labs(x="Fst", y="Density")
plot_2_fst

# combine plots
plot_1_fst + plot_2_fst +  plot_layout(widths = c(5, 1))
# 1: Removed 22 rows containing missing values or values outside the scale range (`geom_point()`). 
# 2: Removed 22 rows containing non-finite outside the scale range (`stat_density_ridges()`). 
## Fst values should only range between 0 and 1. Set anything <0 to 0 and re-do the plots.

# Data cleaning
data_fst_clean <- data %>% filter(data_type == 'Fst') %>% mutate(value = ifelse(value < 0, 0, value)) # this takes my original data, filters it to only include rows where the data_type is Fst. Then mutate function is used to create/modify columns. It modifies the value column, checks if the value is less than 0. If the value is <0, it replaces it with 0. If the value is >0, it keeps the original value.

# Re-do plots with cleaned data
# plot 1 - genome wide plots per population
plot_1_fst_clean <- data_fst_clean %>%
  filter(data_type=='Fst')%>%
  ggplot(., aes(position*100000, value, col=chromosome)) +
  geom_point(size=1) +
  facet_grid(comparison~.) +
  scale_color_tron() +
  theme_bw() +
  theme(legend.position="top", 
        legend.title = element_blank(),
        axis.title = element_text(size = 32),
        axis.text = element_text(size = 16),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 20)) +
  labs(x="Genomic Position", y="Fst")
plot_1_fst_clean


# plot 2 - density plots of Fst per group
plot_2_fst_clean <- data_fst_clean %>%
  filter(data_type=='Fst')%>%
  ggplot(aes(value, chr_type, fill=chr_type), guide="none") +
  geom_density_ridges(quantile_lines=TRUE, quantile_fun=function(x,...)median(x), size=0.5) +
  theme_bw() + theme(legend.position = "none", axis.text.y=element_blank()) +
  facet_grid(comparison~.) +
  scale_fill_npg() +
  labs(x="Fst", y="Density")
plot_2_fst_clean

# combine plots
plot_1_fst_clean + plot_2_fst_clean +  plot_layout(widths = c(5, 1))
# 1: Removed 22 rows containing missing values or values outside the scale range (`geom_point()`). 
# 2: Removed 22 rows containing non-finite outside the scale range (`stat_density_ridges()`).
ggsave("genomewide_and_density_fst_USA_DOG_ONLY.tif", width=20, height=25, dpi = 300)
ggsave("genomewide_and_density_fst_USA_DOG_ONLY.png", width=20, height=25, dpi = 300)
# Fst closer to 1 = more genetically distinct
# Fst closer to 0 = more genetically similar
# Some cities only have 1 sample so they look funny.



# USA map

## Make world map data
world_map <- map_data("world")

# Add map inset to zoom in on USA dog samples
# Manually specify the coordinates for the area of the world map to show in the inset
USA_xmin <- -105
USA_xmax <- -70
USA_ymin <- 20
USA_ymax <- 50

# Filter the world map data for the inset area
USA_data <- subset(world_map, long >= USA_xmin & long <= USA_xmax & lat >= USA_ymin & lat <= USA_ymax)

# USA metadata
location_USA <- read.csv("input/USA_metadata.csv", header = TRUE)

# Set colors for the points
red_palette1 <- brewer.pal(n = 9, name = "YlOrRd")
red_palette2 <- brewer.pal(n=9, name = "YlOrBr")

scale_colour_USA <- 
  c('Florida' = red_palette1[5],
    'Georgia' = red_palette2[5],
    'Illinois' = red_palette1[8],
    'Louisiana' = red_palette1[6],
    'Missouri' = red_palette2[7],
    'Texas' = red_palette1[7])

# Get median fst for each comparison
fst_median <- data_fst_clean %>%
  filter(data_type == "Fst") %>%
  group_by(comparison) %>%
  summarise(FST_MEDIAN = median(value, na.rm = TRUE))
# Save as csv file and then add required metadata such as population longitudes/latitudes
write.csv(fst_median, "fst_median.csv", row.names = TRUE)

# read new csv file with map metadata and median fst between groups
fst_map <- read.csv("fst_map.csv", header = TRUE)

# Fst on USA map
plot_3_fst <- ggplot() +
  geom_polygon(data = USA_data, aes(x = long, y = lat, group = group), fill = "grey90") +
  geom_segment(data=fst_map, aes(x=POP1_LONG, y=POP1_LAT, xend=POP2_LONG, yend=POP2_LAT, color = FST_MEDIAN), size = 1) +
  scale_color_gradient(low = "gold", high = "red3", name = expression(F[ST])) +
  geom_point(data = location_USA, aes(x = Longitude, y = Latitude, fill = State), size = 10, shape = 21) +
  scale_fill_manual(values = scale_colour_USA, limits = c("Florida", "Georgia", "Illinois", "Louisiana", "Missouri", "Texas"), name = "State") +
  geom_text_repel(data = location_USA, aes(x = Longitude, y = Latitude, label = State), size = 6, fontface = "bold", nudge_x = 2, nudge_y = 1) +
  theme_void() +
  theme(
    legend.position = c(0.9,0.1),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 14)) +
  guides(fill = "none") +
  coord_fixed(ratio=1)
plot_3_fst

# Save plot
ggsave("fst_map_USA_DOG_ONLY.tif", bg = "white", height=8, width=10, dpi = 300)
ggsave("fst_map_USA_DOG_ONLY.png", bg = "white", height=8, width=10, dpi = 300)


# Fst by USA on a heatmap
fst_heatmap <- read.csv("fst_heatmap.csv", header = TRUE)

plot_4_fst <- ggplot(fst_heatmap, aes(x = POP1, y=POP2, fill = FST_MEDIAN)) + 
  geom_tile(color = "white", lwd = 1.5, linetype = 1) +
  coord_fixed() +
  scale_fill_material("indigo") +
  guides(fill = guide_colourbar(title = expression(F[ST]), barwidth = 1, barheight = 5)) +
  theme_minimal() +
  labs(x = "Population", y = "Population") +
  theme(axis.text = element_text(size = 18),
    axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title = element_text(size = 24),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank())
plot_4_fst

ggsave("fst_heatmap_USA_DOG_ONLY.tif", bg = "white", height=5, width = 8, dpi = 300)
ggsave("fst_heatmap_USA_DOG_ONLY.png", bg = "white", height=5, width = 8, dpi = 300)
```