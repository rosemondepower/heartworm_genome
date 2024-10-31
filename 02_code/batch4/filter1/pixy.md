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

# Run pixy per host
bsub.py --queue long --threads 10 40 pixy_host \
"pixy --stats pi fst dxy \
--vcf ${VCF} \
--populations host_pop.list \
--window_size 100000 \
--n_cores 10 \
--output_prefix host"
## Keep in mind that the sample sizes for non-dog hosts are quite low, but will still be interesting to explore
```

Download the file and use it as input in R to estimate the values of pi, fst and dxy and generate plots.

```R
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
```

### Pi, Dxy & Fst by region

```R
############ By region #############

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

#Let's add the chr type variable
pi_data <- pi_data %>%
  mutate(chr_type = ifelse(str_detect(chromosome, "X"), "sexchr", "autosome"))

# calculate the median Pi and check the ratio of sex-to-autosome diversity. Should be about 0.75, as D. immitis is XX/XY
pi_data_sex_median <- pi_data %>%
  group_by(chr_type) %>%
  summarise(median = median(avg_pi, na.rm = TRUE))

pi_data_pop_sex_median <- pi_data %>%
  group_by(pop, chr_type) %>%
  summarise(median = median(avg_pi, na.rm = TRUE))

# plot 1 - genome wide plots per population
plot_1 <- ggplot(pi_data, aes(position*100000, avg_pi, col=chromosome)) +
  geom_point(size=1) +
  facet_grid(pop~.) +
  scale_color_simpsons() +
  theme_bw() +
  theme(legend.position="top", 
        legend.title = element_blank()) +
  labs(x="Genomic Position", y="Nucleotide Diversity (Pi)")


# plot 2 - density plots of pi per group
plot_2 <- ggplot(pi_data, aes(avg_pi, chr_type, fill=chr_type), guide="none") +
  geom_density_ridges(quantile_lines=TRUE, quantile_fun=function(x,...)median(x), size=0.5) +
  theme_bw() + theme(legend.position = "none", axis.text.y=element_blank()) +
  facet_grid(pop~.) +
  xlim(0, 0.005) +
  scale_x_continuous(breaks = c(0, 0.0025, 0.005)) +
  scale_fill_npg() +
  labs(x="Nucleotide Diversity (Pi)", y="Density")

# combine plots
plot_1 + plot_2 +  plot_layout(widths = c(5, 1))
ggsave("genomewide_and_density_Pi_region.tif", width=9, height=6, dpi = 300)

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
        red_palette2[7]), 
        "blueviolet",
        green_palette1[7]),
     c('AUS', 'ASIA', 'USA', 'EUR', 'CENAM')), 
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
  scale_colour_region ()+
  ylim(0, 0.003)

boxplot_pi

ggsave("boxplot_Pi_region.tif", width=4, height=4, dpi = 300)
# higher pi value = more diverse population
# lower pi value = less diverse population


#let's check the statistical significance of the differences between pi values
#Let's explore the normality
pi_data %>% filter(chromosome != 'chrX') %>%
ggplot(., aes(x = avg_pi, colour = pop)) +
  geom_histogram() +
  theme_bw() +
  facet_grid(. ~ pop)

#Confirm the absence of normality by shapiro
pi_data_AUS <- pi_data_AUS %>% filter(chromosome != 'chrX')
AUS_shapiro <- shapiro.test(pi_data_AUS$avg_pi)
print(AUS_shapiro)

pi_data_ASIA <- pi_data_ASIA %>% filter(chromosome != 'chrX')
ASIA_shapiro <- shapiro.test(pi_data_ASIA$avg_pi)
print(ASIA_shapiro)

pi_data_USA <- pi_data_USA %>% filter(chromosome != 'chrX')
USA_shapiro <- shapiro.test(pi_data_USA$avg_pi)
print(USA_shapiro)

pi_data_EUR <- pi_data_EUR %>% filter(chromosome != 'chrX')
EUR_shapiro <- shapiro.test(pi_data_EUR$avg_pi)
print(EUR_shapiro)

pi_data_CENAM <- pi_data_CENAM %>% filter(chromosome != 'chrX')
CENAM_shapiro <- shapiro.test(pi_data_CENAM$avg_pi)
print(CENAM_shapiro)

#Willcoxson test for all pairs of pops
wilcox.test(pi_data_AUS$avg_pi, pi_data_ASIA$avg_pi)
wilcox.test(pi_data_AUS$avg_pi, pi_data_USA$avg_pi)
wilcox.test(pi_data_AUS$avg_pi, pi_data_EUR$avg_pi) 
wilcox.test(pi_data_AUS$avg_pi, pi_data_CENAM$avg_pi)
wilcox.test(pi_data_ASIA$avg_pi, pi_data_USA$avg_pi)
wilcox.test(pi_data_ASIA$avg_pi, pi_data_EUR$avg_pi)
wilcox.test(pi_data_ASIA$avg_pi, pi_data_CENAM$avg_pi)
wilcox.test(pi_data_USA$avg_pi, pi_data_EUR$avg_pi)
wilcox.test(pi_data_USA$avg_pi, pi_data_CENAM$avg_pi)
wilcox.test(pi_data_EUR$avg_pi, pi_data_CENAM$avg_pi)


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



#Plotting Dxy

# plot 1 - genome wide plots per comparison
plot_1_dxy <- data %>%
  filter(data_type =='Dxy')%>%
  ggplot(., aes(position*100000, value, col=chromosome)) +
  geom_point(size=1) +
  facet_grid(comparison~.) +
  scale_color_simpsons() +
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

# combine plots
dxy_plot <- plot_1_dxy + plot_2_dxy +  plot_layout(widths = c(5, 1))
ggsave("genomewide_and_density_dxy_region.tif", width=9, height=6, dpi = 300)

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

# AUS_v_ASIA & AUS_v_EUR
scatter_dxy_a <-  data %>%
  filter(data_type == 'Dxy') %>%
  as.data.frame() %>%
  select('position', 'comparison', 'value', 'chromosome')%>%
  as_tibble() %>%
  pivot_wider(names_from = comparison, values_from = value) %>%
  ggplot(., aes(x=AUS_v_ITL, y = AUS_v_USA, col=chromosome)) +
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
         "yy" = AUS_v_ASIA - ASIA_v_EUR) %>%
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
         "yy" = AUS_v_ASIA / ASIA_v_EUR) %>%
  ggplot(., aes(position*100000, xx, col=chromosome)) +
  geom_point(size=1) +
  scale_color_tron() +
  theme_bw() +
  theme(legend.position="top", 
        legend.title = element_blank()) +
  labs(x="Position", y="AUS vs ASIA / AUS vs EUR")
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
ggsave("dxy_lineplot_region.tif", width=9, height=6, dpi = 300)





# Plotting Fst

# plot 1 - genome wide plots per population
plot_1_fst <- data %>%
  filter(data_type=='Fst')%>%
  ggplot(., aes(position*100000, value, col=chromosome)) +
  geom_point(size=1) +
  facet_grid(comparison~.) +
  scale_color_simpsons() +
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

# combine plots
plot_1_fst + plot_2_fst +  plot_layout(widths = c(5, 1))
ggsave("genomewide_and_density_fst_region.tif", width=9, height=6, dpi = 300)
# Fst closer to 1 = more genetically distinct
# Fst closer to 0 = more genetically similar


# Region map metadata. Just have one random point representing each region.
location <- read.csv("region_metadata.csv", header = TRUE)

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

# Fst by region on a map
plot_3_fst <- ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = "grey90") +
  geom_segment(data=fst, aes(x=POP1_LONG, y=POP1_LAT, xend=POP2_LONG, yend=POP2_LAT, color = cut(Fst, breaks = c(0, 0.1, 0.5, 1), labels = c("low", "mid", "high")), size = 0.5))
  geom_point(data = location, aes(x = Longitude, y = Latitude, color = Region), size = 1.5) +
  geom_text_repel(data = location, aes(x = Longitude, y = Latitude, label = Region), size = 8, fontface = "bold", max.overlaps = 20) +  # Add 
  scale_color_manual(values = c("low" = "lightsalmon", "mid" = "tomato" , "orangered3" = "red")) +
  new_scale_color() +
  scale_color_manual(values = scale_colour_region, limits = c("USA", "CENAM", "EUR", "ASIA", "AUS")) +
  theme_void() +
  theme(
    legend.position = "none"
  )
plot_3_fst

# Save plot
ggsave("fst_map_region.tif", height=5, width=10, dpi = 300)


# Fst by region on a heatmap
plot_4_fst <- ggplot(fst, aes(x = POP1, y=POP2, fill = Fst)) + 
  geom_tile(color = "white", lwd = 1.5, linetype = 1) +
  coord_fixed() +
  scale_fill_gradient(low = "lightsalmon", high = "orangered3") +
  guides(fill = guide_colourbar(title = expression(F[ST]), barwidth = 0.5, barheight = 20)) +
  theme_minimal() +
  labs(x = "Population", y = "Population") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot_4_fst

ggsave(fst_heatmap_region.tif, height=5, width = 8, dpi = 300)
```


Next...

Run pixy by host

Run pixy by city for AUS only

Do the same for USA only