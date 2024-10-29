# Dirofilaria immitis WGS Lab Book - Pixy

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

cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/04_VARIANTS/ALLSITES

# make list of vcfs
ls -1 *.list.vcf.gz | sort -n vcf_files.list
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
#How many did it keep?

#select nuclear invariants
bsub.py 1 select_nuclearINVARIANTs "gatk SelectVariants \
--reference ${REFERENCE} \
--variant ${VCF} \
--select-type-to-include NO_VARIATION \
--exclude-intervals ${MIT_CONTIG} \
--exclude-intervals ${WB_CONTIG} \
--output ${VCF%.vcf.gz}.nuclearINVARIANTs.vcf"

vcftools --vcf DIMMITIS_POPGEN_ALLSITES.nuclearINVARIANTs.vcf --remove-indels
#How many did it keep?
```

### Merge these nuclear invariants with the nuclear SNPs I previously filtered

```bash
#Keep n=128 samples in the nuclear dataset (nuclear_samplelist.keep) and merge these with the SNPs
bsub.py 1 filter_nuclear_INVARIANTvcftools "vcftools --vcf ${VCF%.vcf.gz}.nuclearINVARIANTs.vcf \
--keep /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/04_VARIANTS/FILTER1/NO_OUTGROUPS/nuclear_samplelist.keep \
--recode --out nuclearINVARIANTS_128samples"
#How many did it keep?

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
#How many did it keep?

bgzip nuclearSNPssandINVARIANTs_124samples.chrXto4.recode.vcf;
tabix -p vcf nuclearSNPssandINVARIANTs_124samples.chrXto4.recode.vcf.gz
mv nuclearSNPssandINVARIANTs_124samples.chrXto4.recode.vcf.* ../FILTER1/NO_OUTGROUPS/FINAL_SETS/
```

### Run pixy on nuclear variants & invariants

```bash
# Create pixy environment
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

# Population files per country and city
# country_pop.list
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
CRI_SJO_AD_001	CRI
GRC_XAN_AD_001	GRC
GRC_XAN_AD_002	GRC
GRC_XAN_AD_003	GRC
GRC_XAN_AD_004	GRC
GRC_XAN_AD_005	GRC
GRC_XAN_AD_006	GRC
GRC_XAN_AD_007	GRC
GRC_XAN_AD_008	GRC
GRC_XAN_AD_009	GRC
GRC_XAN_AD_010	GRC
GRC_XAN_AD_011	GRC
GRC_XAN_AD_012	GRC
ITA_NEA_AD_001	ITA
ITA_NEA_AD_002	ITA
ITA_NEA_AD_003	ITA
ITA_PAV_AD_001	ITA
MYS_SEL_AD_001	MYS
PAN_BOC_AD_001	PAN
PAN_BOC_AD_002	PAN
PAN_BOC_AD_003	PAN
PAN_BOC_AD_004	PAN
PAN_BOC_AD_005	PAN
PAN_PUE_AD_001	PAN
PAN_PUE_AD_002	PAN
PAN_PUE_AD_003	PAN
PAN_PUE_AD_004	PAN
PAN_PUE_AD_005	PAN
PAN_PUE_AD_006	PAN
PAN_SLO_AD_001	PAN
PAN_SLO_AD_002	PAN
PAN_SLO_AD_003	PAN
ROU_BUC_AD_001	ROU
ROU_COM_AD_001	ROU
ROU_GIU_AD_001	ROU
THA_BKK_AD_001	THA
THA_BKK_AD_002	THA
THA_BKK_AD_003	THA
THA_BKK_AD_004	THA
THA_BKK_AD_005	THA
THA_BKK_AD_006	THA
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

# city_pop.list
'
AUS_BNE_AD_001	Brisbane
AUS_BNE_AD_002	Brisbane
AUS_BNE_AD_003	Brisbane
AUS_BNE_AD_004	Brisbane
AUS_BNE_AD_006	Brisbane
AUS_BNE_AD_008	Brisbane
AUS_BNE_AD_009	Brisbane
AUS_CNS_AD_001	Cairns
AUS_CNS_AD_002	Cairns
AUS_LHR_AD_001	Lockhart River
AUS_ROK_AD_001	Rockhampton
AUS_SYD_AD_001	Sydney
AUS_SYD_AD_002	Sydney
AUS_SYD_AD_003	Sydney
AUS_SYD_AD_004	Sydney
AUS_SYD_AD_005	Sydney
AUS_SYD_AD_006	Sydney
AUS_SYD_AD_007	Sydney
AUS_SYD_AD_008	Sydney
AUS_SYD_AD_009	Sydney
AUS_SYD_AD_010	Sydney
AUS_SYD_AD_011	Sydney
AUS_SYD_AD_012	Sydney
AUS_SYD_AD_013	Sydney
AUS_SYD_AD_014	Sydney
AUS_SYD_AD_015	Sydney
AUS_SYD_AD_017	Sydney
AUS_TVS_AD_001	Townsville
AUS_TVS_AD_002	Townsville
AUS_TVS_AD_003	Townsville
AUS_TVS_AD_004	Townsville
AUS_TVS_AD_005	Townsville
AUS_TVS_AD_006	Townsville
AUS_TVS_AD_007	Townsville
AUS_TVS_AD_008	Townsville
AUS_TVS_AD_010	Townsville
AUS_TVS_AD_011	Townsville
AUS_TVS_AD_012	Townsville
AUS_TVS_AD_013	Townsville
AUS_TVS_AD_014	Townsville
AUS_TVS_AD_015	Townsville
AUS_TVS_AD_016	Townsville
AUS_TVS_AD_017	Townsville
AUS_TVS_AD_018	Townsville
AUS_TVS_AD_019	Townsville
CRI_SJO_AD_001	San Jose
GRC_XAN_AD_001	Thessaloniki/Xanthi
GRC_XAN_AD_002	Thessaloniki/Xanthi
GRC_XAN_AD_003	Thessaloniki/Xanthi
GRC_XAN_AD_004	Thessaloniki/Xanthi
GRC_XAN_AD_005	Thessaloniki/Xanthi
GRC_XAN_AD_006	Thessaloniki/Xanthi
GRC_XAN_AD_007	Thessaloniki/Xanthi
GRC_XAN_AD_008	Thessaloniki/Xanthi
GRC_XAN_AD_009	Thessaloniki/Xanthi
GRC_XAN_AD_010	Thessaloniki/Xanthi
GRC_XAN_AD_011	Thessaloniki/Xanthi
GRC_XAN_AD_012	Thessaloniki/Xanthi
ITA_NEA_AD_001	Northeast
ITA_NEA_AD_002	Northeast
ITA_NEA_AD_003	Northeast
ITA_PAV_AD_001	Pavia
MYS_SEL_AD_001	Selangor
PAN_BOC_AD_001	Boca Chica
PAN_BOC_AD_002	Boca Chica
PAN_BOC_AD_003	Boca Chica
PAN_BOC_AD_004	Boca Chica
PAN_BOC_AD_005	Boca Chica
PAN_PUE_AD_001	Puerto Armuelles
PAN_PUE_AD_002	Puerto Armuelles
PAN_PUE_AD_003	Puerto Armuelles
PAN_PUE_AD_004	Puerto Armuelles
PAN_PUE_AD_005	Puerto Armuelles
PAN_PUE_AD_006	Puerto Armuelles
PAN_SLO_AD_001	San Lorenzo
PAN_SLO_AD_002	San Lorenzo
PAN_SLO_AD_003	San Lorenzo
ROU_BUC_AD_001	Bucharest
ROU_COM_AD_001	Comana
ROU_GIU_AD_001	Giurgiu
THA_BKK_AD_001	Bangkok
THA_BKK_AD_002	Bangkok
THA_BKK_AD_003	Bangkok
THA_BKK_AD_004	Bangkok
THA_BKK_AD_005	Bangkok
THA_BKK_AD_006	Bangkok
USA_FLO_AD_001	Florida
USA_FLO_AD_002	Florida
USA_FLO_AD_003	Florida
USA_FLO_AD_004	Florida
USA_FLO_AD_005	Florida
USA_FLO_AD_006	Florida
USA_FLO_AD_007	Florida
USA_FLO_AD_008	Florida
USA_FLO_AD_009	Florida
USA_FLO_AD_010	Florida
USA_FLO_AD_011	Florida
USA_FLO_AD_012	Florida
USA_FLO_AD_013	Florida
USA_FLO_AD_014	Florida
USA_FLO_AD_015	Florida
USA_FLO_AD_016	Florida
USA_GEO_AD_001	Georgia
USA_ILL_AD_001	Illinois
USA_ILL_AD_002	Illinois
USA_LOU_AD_001	Louisiana
USA_LOU_AD_002	Louisiana
USA_MIS_AD_001	Missouri
USA_MIS_AD_002	Missouri
USA_TEX_AD_001	Texas
USA_TEX_AD_002	Texas
USA_TEX_AD_003	Texas
USA_TEX_AD_004	Texas
USA_TEX_AD_005	Texas
USA_TEX_AD_006	Texas
USA_TEX_AD_007	Texas
USA_TEX_AD_008	Texas
USA_TEX_AD_009	Texas
USA_TEX_AD_010	Texas
USA_TEX_AD_011	Texas
USA_TEX_AD_012	Texas
USA_TEX_AD_013	Texas
USA_TEX_AD_014	Texas
USA_TEX_AD_015	Texas
'


# Run pixy per country
bsub.py --queue long --threads 20 20 pixy_country \
"pixy --stats pi fst dxy \
--vcf ${VCF} \
--populations country_pop.list \
--window_size 100000 \
--n_cores 20 \
--output_prefix country"

# Run pixy per city
bsub.py --queue long --threads 20 20 pixy_city \
"pixy --stats pi fst dxy \
--vcf ${VCF} \
--populations city_pop.list \
--window_size 100000 \
--n_cores 20 \
--output_prefix city"
```

Download the file and use it as input in R to estimate the values of pi, fst and dxy and generate plots.

```R
library(tidyverse)
library(ggsci)
library(ggpubr)
library(patchwork)
library(ggridges)

# Nucleotide diversity (pi)

##### PI #####

# get nucleotide diversity (pi) data from pixy output
pi_data <- read.table("input/city_pi.txt", header=T)
pi_data$chromosome <- str_remove(pi_data$chromosome, 'dirofilaria_immitis_')

# subset
pi_data_AUS <- pi_data %>%
  filter(pop=="AUS")
pi_data_CRI <- pi_data %>%
  filter(pop=="CRI")
pi_data_GRC <- pi_data %>%
  filter(pop=="GRC")
pi_data_ITA <- pi_data %>%
  filter(pop=="ITA")
pi_data_MYS <- pi_data %>%
  filter(pop=="MYS")
pi_data_PAN <- pi_data %>%
  filter(pop=="PAN")
pi_data_ROU <- pi_data %>%
  filter(pop=="ROU")
pi_data_THA <- pi_data %>%
  filter(pop=="THA")
pi_data_USA <- pi_data %>%
  filter(pop=="USA")

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

# calculate the median Pi and checking the ratio of sex-to-autosome diversity. Should be about 0.75, as D. immitis is XX/XY
pi_data_sex_median <- pi_data %>%
  group_by(chr_type) %>%
  summarise(median = median(avg_pi, na.rm = TRUE)) ############# UP TO HERE
'
# A tibble: 2 × 2
chr_type   median
<chr>       <dbl>
1 autosome 0.000269
2 sexchr   0.000133
'
# 0.000133 / 0.000269 = 0.4944238 (far off 0.75 expected of diversity on sex chromosome relative to autosome)

pi_data_pop_sex_median <- pi_data %>%
  group_by(pop, chr_type) %>%
  summarise(median = median(avg_pi, na.rm = TRUE))
'
# A tibble: 6 × 3
# Groups:   pop [3]
pop   chr_type     median
<chr> <chr>         <dbl>
1 AUS   autosome 0.0000391 
2 AUS   sexchr   0.00000376
3 ITL   autosome 0.0000864 
4 ITL   sexchr   0.0000322 
5 USA   autosome 0.000736  
6 USA   sexchr   0.000426  
'
# plot 1 - genome wide plots per population
plot_1 <- ggplot(pi_data, aes(position*100000, avg_pi, col=chromosome)) +
  geom_point(size=1) +
  facet_grid(pop~.) +
  scale_color_tron() +
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

# bring it together
plot_1 + plot_2 +  plot_layout(widths = c(5, 1))
ggsave("plots_genomewide_and_density_Pi.png", width=9, height=6)

#Now a boxplot of the pi value per population

scale_colour_pop <- function(...){
  ggplot2:::manual_scale(
    'colour', 
    values = setNames(
      c('purple2',
        'turquoise3',
        'green4',
        'orange3',
        'red3'), 
      c('Lockhart River Cooktown', 'Cairns', 
        'Townsville',
        'Brisbane',
        'Sydney')), 
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

ggsave("plots_boxplot_pop_Pi.png", width=4, height=4)

#let's check the statistical significance of the differences between pi values
#Let's explore the normality
pi_data %>% filter(chromosome != 'chrX') %>%
ggplot(., aes(x = avg_pi, colour = pop)) +
  geom_histogram() +
  theme_bw() +
  facet_grid(. ~ pop)

#Confirm the absence of normality by shapiro
pi_data_Brisbane <- pi_data_Brisbane %>% filter(chromosome != 'chrX')
Brisbane_shapiro <- shapiro.test(pi_data_Brisbane$avg_pi)
print(Brisbane_shapiro)#W = 0.69062, p-value < 2.2e-16

pi_data_Cairns <- pi_data_Cairns %>% filter(chromosome != 'chrX')
Cairns_shapiro <- shapiro.test(pi_data_Cairns$avg_pi)
print(Cairns_shapiro)#W = 0.69994, p-value < 2.2e-16

pi_data_Townsville <- pi_data_Townsville %>% filter(chromosome != 'chrX')
Townsville_shapiro <- shapiro.test(pi_data_Townsville$avg_pi)
print(Townsville_shapiro)#W = 0.90598, p-value < 2.2e-16

pi_data_Sydney <- pi_data_Sydney %>% filter(chromosome != 'chrX')
Sydney_shapiro <- shapiro.test(pi_data_Sydney$avg_pi)
print(Sydney_shapiro)

pi_data_Cooktown <- pi_data_Cooktown %>% filter(chromosome != 'chrX')
Cooktown_shapiro <- shapiro.test(pi_data_Cooktown$avg_pi)
print(Cooktown_shapiro)

#Willcoxson test for everyone
wilcox.test(pi_data_Brisbane$avg_pi, pi_data_Cairns$avg_pi)
#W = 152048, p-value = 1.16e-06
wilcox.test(pi_data_Brisbane$avg_pi, pi_data_Townsville$avg_pi)
#W = 41288, p-value < 2.2e-16
wilcox.test(pi_data_Brisbane$avg_pi, pi_data_Sydney$avg_pi)
#W = 50720, p-value < 2.2e-16
wilcox.test(pi_data_Brisbane$avg_pi, pi_data_Cooktown$avg_pi)
wilcox.test(pi_data_Cairns$avg_pi, pi_data_Townsville$avg_pi)
wilcox.test(pi_data_Cairns$avg_pi, pi_data_Sydneyavg_pi)
wilcox.test(pi_data_Cairns$avg_pi, pi_data_Cooktown$avg_pi)
wilcox.test(pi_data_Townsville$avg_pi, pi_data_Sydney$avg_pi)
wilcox.test(pi_data_Townsville$avg_pi, pi_data_Cooktown$avg_pi)
wilcox.test(pi_data_Sydney$avg_pi, pi_data_Cooktown$avg_pi)


#let's generate some dataframe for the estatistics
# subset
pi_data_Brisbane <- pi_data %>%
  filter(pop=="Brisbane") %>%
  filter(chromosome != 'chrX')

pi_data_Cairns <- pi_data %>%
  filter(pop=="Cairns") %>%
  filter(chromosome != 'chrX') 

pi_data_Townsville <- pi_data %>%
  filter(pop=="Townsville") %>%
  filter(chromosome != 'chrX')

pi_data_Sydney <- pi_data %>%
  filter(pop=="Sydney") %>%
  filter(chromosome != 'chrX')

pi_data_Cooktown <- pi_data %>%
  filter(pop=="Cooktown") %>%
  filter(chromosome != 'chrX')
```