# Dirofilaria immitis WGS Lab Book - Pixy NOTES

## Generate an ALL SITES variant set for running pixy

```bash
module load bsub.py/0.42.1
module load gatk/4.1.4.1

cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/04_VARIANTS/ALLSITES

# Create cohort gvcf

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

# try running another job with threads to be faster (still took a while, disregard this)
for file in $SPLIT; do
 base=$(basename $file .txt)
  bsub.py --queue basement --threads 8 20 gatk_combinegvcfs_${base}_NEW \
  "gatk CombineGVCFs -R /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/01_REF/dimmitis_WSI_2.2.fa --intervals ${file} ${args} -O /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/04_VARIANTS/ALLSITES_NEW/DIMMITIS_POPGEN_ALLSITES_${base}.g.vcf.gz" ;
done


# Run genotyping on each chromosome separately
for file in $SPLIT; do
 base=$(basename $file .txt)
  bsub.py --queue long 20 gatk_genotypevcfs_${base} \
  "gatk GenotypeGVCFs -R /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/01_REF/dimmitis_WSI_2.2.fa -V /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/04_VARIANTS/ALLSITES/DIMMITIS_POPGEN_ALLSITES_${base}.g.vcf.gz --intervals ${file} --all-sites -O DIMMITIS_POPGEN_ALLSITES_${base}.vcf.gz -G StandardAnnotation -G AS_StandardAnnotation" ;
done






#!/bin/bash

# Export constants
export PREFIX=DIMMITIS_POPGEN_ALLSITES  # prefix for output files
export REFERENCE=/lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/01_REF/dimmitis_WSI_2.2.fa  # path to reference genome
export GVCF_LIST=/lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/04_VARIANTS/ALLSITES/gvcf.list  # path to list of GVCF files

# Load GATK module
module load gatk/4.1.4.1
module load fastaq/3.17.0-docker3
module load samtools/1.14--hb421002_0

# Define directories
export LOG_FILES="${PWD}/gatk_hc_${PREFIX}/LOG_FILES_RP"  # directory for log files
export REFERENCE_FILES="${PWD}/gatk_hc_${PREFIX}/REFERENCE_FILES_RP"  # directory for reference files
export GATK_HC_MERGED="${PWD}/gatk_hc_${PREFIX}/GATK_HC_MERGED_RP"  # directory for merged haplotype caller files

# Create directories if they don't exist
mkdir -p ${LOG_FILES}
mkdir -p ${REFERENCE_FILES}
mkdir -p ${GATK_HC_MERGED}

# Save current script in run folder to reproduce the exact output
cp ${PWD}/run_combinegvcfs.sh ${PWD}/gatk_hc_${PREFIX}/commands.$(date -Iminutes).txt

#-------------------------------------------------------------------------------
### 01. Prepare reference files
#-------------------------------------------------------------------------------

func_build_reference() {
    # Check if the reference genome file already exists
    if [ -f "${REFERENCE_FILES}/REF.fa" ]; then
        echo "Reference is already setup. Moving on."
        exit 0
    else
        # Copy the reference genome file to the REFERENCE_FILES directory
        cp "${REFERENCE}" "${REFERENCE_FILES}/REF.fa"
        # Create an index file for the reference genome
        samtools faidx "${REFERENCE_FILES}/REF.fa"
        # Create a dictionary file for the reference genome

        # Split the reference genome file into chunks of approximately 10 Mb in size
        fastaq split_by_base_count "${REFERENCE_FILES}/REF.fa" "${REFERENCE_FILES}/REFsplit" 10000000

        # For each chunked genome section, create a list of contig/scaffold names
        for i in "${REFERENCE_FILES}/REFsplit"* ; do
            # Extract the name of the chunked genome section
            NAME=$( echo "${i}" | awk -F '/' '{print $NF}' )
            # Extract the contig/scaffold names from the chunked genome section
            grep ">" "${i}" | sed 's/>//g' > "${REFERENCE_FILES}/${NAME}.list"
        done
    fi
}

export -f func_build_reference




#-------------------------------------------------------------------------------
### 03. Merge GVCFs
#-------------------------------------------------------------------------------
func_merge_gvcf() {
  
  mkdir -p ${GATK_HC_MERGED}/LOGFILES

  n=1
  for SEQUENCE in ${REFERENCE_FILES}/REFsplit*list; do
    SEQUENCE=$( echo ${SEQUENCE} | awk -F '/' '{print $NF}' )
    echo -e "gatk CombineGVCFs -R ${REFERENCE_FILES}/REF.fa --intervals ${REFERENCE_FILES}/${SEQUENCE} \\" > ${GATK_HC_MERGED}/run_merge_gvcfs.tmp.job_${n}
    while read SAMPLE; do
        echo -e "--variant ${SAMPLE} \\" >> ${GATK_HC_MERGED}/run_merge_gvcfs.tmp.job_${n}
    done < ${GVCF_LIST}
    echo -e "--output ${GATK_HC_MERGED}/${n}.cohort.tmp.g.vcf.gz" >> ${GATK_HC_MERGED}/run_merge_gvcfs.tmp.job_${n}
    let "n+=1"
  done

  chmod a+x ${GATK_HC_MERGED}/run_merge_gvcfs.tmp.job_*

  JOBS=$( ls -1 ${GATK_HC_MERGED}/run_merge_gvcfs.tmp.job_* | wc -l )
  ID="U$(date +%s)"

  bsub -q long -R'span[hosts=1] select[mem>30000] rusage[mem=30000]' -n 10 -M30000 -J "gatk_merge_gvcf_[1-$JOBS]%100" -e "${GATK_HC_MERGED}/LOGFILES/gatk_merge_gvcf_[1-$JOBS].e" -o "${GATK_HC_MERGED}/LOGFILES/gatk_merge_gvcf_[1-$JOBS].o" "${GATK_HC_MERGED}/run_merge_gvcfs.tmp.job_\$LSB_JOBINDEX"

  rm ${GATK_HC_MERGED}/MERGE_ARRAY_FINISHED
  bsub -w "done(gatk_merge_gvcf_)" -q normal -R'span[hosts=1] select[mem>100] rusage[mem=100]' -n 1 -M100 -J "gatk_merge_gvcf_finish" -e "${GATK_HC_MERGED}/LOGFILES/gatk_merge_gvcf_finish.e" -o "${GATK_HC_MERGED}/LOGFILES/gatk_merge_gvcf_finish.o" "touch ${GATK_HC_MERGED}/MERGE_ARRAY_FINISHED"

  until [ -f "${GATK_HC_MERGED}/MERGE_ARRAY_FINISHED" ]; do
    sleep 10
  done
}

export -f func_merge_gvcf


#-------------------------------------------------------------------------------
# running the pipeline
#-------------------------------------------------------------------------------

# func_build_reference
bsub -E 'test -e /nfs/users/nfs_r/rp24' -R "select[mem>1000] rusage[mem=1000]" -M1000 -o ${LOG_FILES}/gatk_01_build_reference.o -e ${LOG_FILES}/gatk_01_build_reference.e -J gatk_01_build_reference_${PREFIX} func_build_reference

# func_merge_gvcf
bsub -w "done(gatk_01_build_reference_${PREFIX})" -E 'test -e /nfs/users/nfs_r/rp24' -R "select[mem>50000] rusage[mem=50000]" -q long -M50000 -n20 -o ${LOG_FILES}/gatk_03_merge_gvcf.o -e ${LOG_FILES}/gatk_03_merge_gvcf.e -J gatk_03_merge_gvcf_${PREFIX} func_merge_gvcf
















































## Run pixy in Artemis


```bash
#Generate a population file
#pop.list
'JS6277	Sydney
JS6278	Sydney
JS6279	Sydney
JS6280	Sydney
JS6281	Lockhart River Cooktown
JS6342	Brisbane
JS6343	Brisbane
JS6344	Brisbane
JS6345	Brisbane
JS6346	Brisbane
JS6347	Brisbane
JS6349	Cairns
JS6350	Cairns
JS6351	Townsville
JS6352	Townsville
JS6353	Townsville
JS6354	Townsville
JS6355	Townsville
JS6356	Townsville
JS6357	Townsville
JS6358	Townsville
JS6359	Townsville
JS6360	Townsville
JS6368	Townsville
JS6369	Townsville
JS6370	Townsville
SRR13154013	Sydney
SRR13154014	Sydney
SRR13154015	Sydney
SRR13154016	Sydney
SRR13154017	Sydney'
```


```bash
#!/bin/bash

# PBS directives 
#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N pixy
#PBS -l select=1:ncpus=20:mem=50GB
#PBS -l walltime=10:00:00
#PBS -m e
#PBS -q defaultQ
#PBS -o pixy.txt
#PBS -M rosemonde.power@sydney.edu.au

# qsub ../pixy.pbs

# Load modules
module load anaconda3

# Activate myenv
conda activate myenv
module load htslib/1.14

# Environmental variables
WORKING_DIR=/scratch/RDS-FSC-Heartworm_MLR-RW/mapping
#VCF=${WORKING_DIR}/pixy/nuclearSNPssandINVARIANTs.chrxto4.recode.vcf.gz

cd ${WORKING_DIR}/pixy

# Diversity per city
pixy --stats pi fst dxy \
--vcf ${VCF} \
--populations /project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_ALL/data/analysis/pop.list \
--window_size 100000 \
--n_cores 20 \
--output_prefix city
```

Download the file and use it as input in R to estimate the values and generate plots.


## Nucleotide diversity in R

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
pi_data_Brisbane <- pi_data %>%
  filter(pop=="Brisbane")
pi_data_Cairns <- pi_data %>%
  filter(pop=="Cairns")
pi_data_Townsville <- pi_data %>%
  filter(pop=="Townsville")
  pi_data_Sydney <- pi_data %>%
  filter(pop=="Sydney")
  pi_data_Cooktown <- pi_data %>%
  filter(pop=="Lockhart River Cooktown")

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

# calculate the median Pi and checking the ratio of sex-to-autosome diversity. Should be about 0.75, as Trichuris is XX/XY
pi_data_sex_median <- pi_data %>%
  group_by(chr_type) %>%
  summarise(median = median(avg_pi, na.rm = TRUE))
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

## Dxy and Fst

```R

#### Dxy and Fst ####

# load data
dxy_data <- read.table("input/city_dxy.txt", header=T)
fst_data <- read.table("input/city_fst.txt", header=T)

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
# A tibble: 6 × 3
# Groups:   comparison [3]
comparison data_type   median
<chr>      <chr>        <dbl>
1 AUS_v_ITL  Dxy       0.000608
2 AUS_v_ITL  Fst       0.876   
3 AUS_v_USA  Dxy       0.000643
4 AUS_v_USA  Fst       0.327   
5 ITL_v_USA  Dxy       0.000666
6 ITL_v_USA  Fst       0.161
'

#Plotting Dxy

# plot 1 - genome wide plots per comparison
plot_1 <- data %>%
  filter(data_type =='Dxy')%>%
  ggplot(., aes(position*100000, value, col=chromosome)) +
  geom_point(size=1) +
  facet_grid(comparison~.) +
  scale_color_tron() +
  theme_bw() +
  theme(legend.position="top", 
        legend.title = element_blank()) +
  labs(x="Genomic Position", y="Dxy")


# plot 2 - density plots of dxy per group
plot_2 <- data %>%
  filter(data_type=='Dxy')%>%
  ggplot(aes(value, chr_type, fill=chr_type), guide="none") +
  geom_density_ridges(quantile_lines=TRUE, quantile_fun=function(x,...)median(x), size=0.5) +
  theme_bw() + theme(legend.position = "none", axis.text.y=element_blank()) +
  facet_grid(comparison~.) +
  scale_fill_npg() +
  labs(x="Dxy", y="Density")

# bring it together
dxy_plot <- plot_1 + plot_2 +  plot_layout(widths = c(5, 1))
ggsave("plots_genomewide_and_density_dxy.png", width=9, height=6)

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

density_dxy <- data %>%
  filter(data_type =='Dxy') %>% 
  ggplot(., aes(x=value, group = comparison, fill=comparison)) +
  geom_density(adjust=1.5, alpha=.4) +
  labs(x = "Dxy" , y = "density", colour = "Comparison") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


# Edit this
scatter_dxy_a <-  data %>%
  filter(data_type == 'Dxy') %>%
  as.data.frame() %>%
  select('position', 'comparison', 'value', 'chromosome')%>%
  as_tibble() %>%
  pivot_wider(names_from = comparison, values_from = value) %>%
  ggplot(., aes(x=AUS_v_ITL, y = AUS_v_USA, col=chromosome)) +
  geom_point(alpha=.4) +
  labs(x = "AUS_v_ITL" , y = "AUS_v_USA", colour = "Comparison") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
scatter_dxy_b <-  data %>%
  filter(data_type == 'Dxy') %>%
  as.data.frame() %>%
  select('position', 'comparison', 'value', 'chromosome')%>%
  as_tibble() %>%
  pivot_wider(names_from = comparison, values_from = value) %>%
  ggplot(., aes(x=AUS_v_ITL, y = ITL_v_USA, col=chromosome)) +
  geom_point(alpha=.4) +
  labs(x = "AUS_v_ITL" , y = "ITL_v_USA", colour = "Comparison") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

genome_pos_dxy_a <-  data %>%
  filter(data_type == 'Dxy') %>%
  as.data.frame() %>%
  select('position', 'comparison', 'value', 'chromosome')%>%
  as_tibble() %>%
  pivot_wider(names_from = comparison, values_from = value) %>%
  mutate("xx" = AUS_v_ITL - AUS_v_USA,
         "yy" = AUS_v_ITL - ITL_v_USA) %>%
  ggplot(., aes(position*100000, xx, col=chromosome)) +
  geom_point(size=1) +
  scale_color_tron() +
  theme_bw() +
  theme(legend.position="top", 
        legend.title = element_blank()) +
  labs(x="Position", y="AUS vs ITL - AUS vs USA")

genome_pos_dxy_b <-  data %>%
  filter(data_type == 'Dxy') %>%
  as.data.frame() %>%
  select('position', 'comparison', 'value', 'chromosome')%>%
  as_tibble() %>%
  pivot_wider(names_from = comparison, values_from = value) %>%
  mutate("xx" = AUS_v_ITL - AUS_v_USA,
         "yy" = AUS_v_ITL - ITL_v_USA) %>%
  ggplot(., aes(position*100000, yy, col=chromosome)) +
  geom_point(size=1) +
  scale_color_tron() +
  theme_bw() +
  theme(legend.position="top", 
        legend.title = element_blank()) +
  labs(x="Position", y="AUS vs ITL - ITL vs USA")

genome_pos_dxy_aa <-  data %>%
  filter(data_type == 'Dxy') %>%
  as.data.frame() %>%
  select('position', 'comparison', 'value', 'chromosome')%>%
  as_tibble() %>%
  pivot_wider(names_from = comparison, values_from = value) %>%
  mutate("xx" = AUS_v_ITL / AUS_v_USA,
         "yy" = AUS_v_ITL / ITL_v_USA) %>%
  ggplot(., aes(position*100000, xx, col=chromosome)) +
  geom_point(size=1) +
  scale_color_tron() +
  theme_bw() +
  theme(legend.position="top", 
        legend.title = element_blank()) +
  labs(x="Position", y="AUS vs ITL / AUS vs USA")

genome_pos_dxy_bb <-  data %>%
  filter(data_type == 'Dxy') %>%
  as.data.frame() %>%
  select('position', 'comparison', 'value', 'chromosome')%>%
  as_tibble() %>%
  pivot_wider(names_from = comparison, values_from = value) %>%
  mutate("xx" = AUS_v_ITL / AUS_v_USA,
         "yy" = AUS_v_ITL / ITL_v_USA) %>%
  ggplot(., aes(position*100000, yy, col=chromosome)) +
  geom_point(size=1) +
  scale_color_tron() +
  theme_bw() +
  theme(legend.position="top", 
        legend.title = element_blank()) +
  labs(x="Position", y="AUS vs ITL / ITL vs USA")

genome_pos_dxy_cc <-  data %>%
  filter(data_type == 'Dxy') %>%
  as.data.frame() %>%
  select('position', 'comparison', 'value', 'chromosome')%>%
  as_tibble() %>%
  pivot_wider(names_from = comparison, values_from = value) %>%
  mutate("xx" = AUS_v_ITL / AUS_v_USA,
         "yy" = AUS_v_ITL / ITL_v_USA,
         "zz" = AUS_v_USA / ITL_v_USA) %>%
  ggplot(., aes(position*100000, zz, col=chromosome)) +
  geom_point(size=1) +
  scale_color_tron() +
  theme_bw() +
  theme(legend.position="top", 
        legend.title = element_blank()) +
  labs(x="Position", y="AUS_v_USA / ITL vs USA")
  
ggarrange(scatter_dxy_a, scatter_dxy_b, common.legend = T)
ggarrange(genome_pos_dxy_a, genome_pos_dxy_b, common.legend = T, ncol = 1)
ggarrange(genome_pos_dxy_aa, genome_pos_dxy_bb, genome_pos_dxy_cc,
          common.legend = T, ncol = 1)

ggsave("nosequecosa.png", width=9, height=6)
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
# Plotting Fst

# plot 1 - genome wide plots per population
plot_1 <- data %>%
  filter(data_type=='Fst')%>%
  ggplot(., aes(position*100000, value, col=chromosome)) +
  geom_point(size=1) +
  facet_grid(comparison~.) +
  scale_color_tron() +
  theme_bw() +
  theme(legend.position="top", 
        legend.title = element_blank()) +
  labs(x="Genomic Position", y="Fst")


# plot 2 - density plots of Fst per group
plot_2 <- data %>%
  filter(data_type=='Fst')%>%
  ggplot(aes(value, chr_type, fill=chr_type), guide="none") +
  geom_density_ridges(quantile_lines=TRUE, quantile_fun=function(x,...)median(x), size=0.5) +
  theme_bw() + theme(legend.position = "none", axis.text.y=element_blank()) +
  facet_grid(comparison~.) +
  scale_fill_npg() +
  labs(x="Fst", y="Density")

# bring it together
plot_1 + plot_2 +  plot_layout(widths = c(5, 1))
ggsave("plots_genomewide_and_density_fst.png", width=9, height=6)

```