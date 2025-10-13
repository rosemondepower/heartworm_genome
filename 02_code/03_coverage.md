# Genome coverage

## Coverage on bam files

Adopted the code from Javier's paper.

```bash
module load PaM/environment
module load bsub.py/0.42.1
module load bamtools/2.5.1--he860b03_5
module load bedtools/2.31.0--hf5e1c6e_3
module load samtools/1.14--hb421002_0
bsub.py -q long --threads 4 50 coverage "./coverage.sh"
```

```bash
#!/bin/bash

WORKING_DIR=/lustre/scratch127/pam/teams/team333/rp24/dirofilaria_immitis/R2/coverage
cd ${WORKING_DIR}

WINDOW='100000'



for i in /lustre/scratch127/pam/teams/team333/sd21/dirofilaria_immitis/project_global-diversity/04_analysis/01_mapping/results_dirofilaria_only/*.bam; do

base_name=$(basename ${i%.bam})

bamtools header -in ${i} | grep "^@SQ" | awk -F'[:\t]' '{printf $3"\t"1"\t"$5"\n"}' OFS="\t" > ${base_name}.chr.bed
bamtools header -in ${i} | grep "^@SQ" | awk -F'[:\t]' '{printf $3"\t"$5"\n"}' OFS="\t" > ${base_name}.chr.genome

bedtools makewindows -g ${base_name}.chr.genome -w ${WINDOW} > ${base_name}.${WINDOW}_window.bed

samtools bedcov -Q 20 ${base_name}.chr.bed ${i} | awk -F'\t' '{printf $1"\t"$2"\t"$3"\t"$4"\t"$4/($3-$2)"\n"}' OFS="\t" > ${base_name}.chr.cov
samtools bedcov -Q 20 ${base_name}.${WINDOW}_window.bed ${i} | awk -F'\t' '{printf $1"\t"$2"\t"$3"\t"$4"\t"$4/($3-$2)"\n"}' OFS="\t" > ${base_name}.${WINDOW}_window.cov

rm ${base_name}.chr.bed ${base_name}.${WINDOW}_window.bed ${base_name}.chr.genome;

done



for i in *.chr.cov; do 

printf "${i}\n" > ${i}.tmp | awk '{print $5}' OFS="\t" ${i} >> ${i}.tmp;

done
```
Successfully completed

```bash
paste *.tmp > coverage_stats.summary
#rm *.tmp

mkdir COV_STATS
mv *.chr.cov *_window.cov *.cov coverage_stats.summary COV_STATS
cd COV_STATS
```


## Generate quantitative stats on coverage for supplementary tables etc
Extract mtDNA, Wb and nuclear (mean & stddev) data

For nuclear, we will select only the defined Chr (chrX and chr1 to chr4)

```bash
# extract mtDNA and nuclear (mean & stddev) data
for i in *.chr.cov; do
	name=${i%.chr.cov};
	nuc=$(grep -v "scaffold\|Wb\|Mt\|NC\|NW" ${i%.chr.cov}.100000_window.cov | datamash mean 5 sstdev 5 );
	mtDNA=$(grep "chrMtDNA" ${i} | cut -f5 );
	Wb=$(grep 'chrWb' ${i} | cut -f5 ); 
	echo -e "${name}\t${nuc}\t${mtDNA}\t${Wb}";
done > 'mito_wolb_cov.stats'
```
Now we'll generate some plots and stats.


## Plot coverage in R

- Coverage ratios
- Coverage across nuclear genome
- Sex of samples

```R
# HW WGS Coverage

# load libraries
library(tidyverse)
library(ggpubr)
library(ggsci)
library(stringr)

setwd("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/R2_Sep25/coverage")

#first, I have to read the nuclear cov stat and estimate the mean and sd
#then, to add it to 'mito_wolb_cov.stats

nuc_mito_wb_cov <- read.table('input/mito_wolb_cov.stats', header = F) %>% as_tibble()

colnames(nuc_mito_wb_cov) <- c('ID', 'nuc_cov', 'sd_nuc_cov', 'mito_cov', 'wb_cov')

write_csv(nuc_mito_wb_cov, 'nuc_mit_wb_cov.csv')

# nuclear, mitochondrial and Wb DNA coverage ratio

n_m <- ggplot(nuc_mito_wb_cov, aes(x=ID, y=mito_cov/nuc_cov)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 12),
        axis.text.y = element_text(size = 18),
        axis.title.x = element_text(size = 28),
        axis.title.y = element_text(size = 28),
        plot.title = element_text(size = 32)
        ) +
  labs(title = "Nuclear to mitochondrial genome coverage ratio", y = "Coverage Ratio")
n_m

n_wb <- ggplot(nuc_mito_wb_cov, aes(x=ID, y=wb_cov/nuc_cov)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 12),
        axis.text.y = element_text(size = 18),
        axis.title.x = element_text(size = 28),
        axis.title.y = element_text(size = 28),
        plot.title = element_text(size = 32)
        ) +
  labs(title = "Nuclear to Wolbachia genome coverage ratio", y = "Coverage Ratio")
n_wb

m_wb <- ggplot(nuc_mito_wb_cov, aes(x=ID, y=mito_cov/wb_cov)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 12),
        axis.text.y = element_text(size = 18),
        axis.title.x = element_text(size = 28),
        axis.title.y = element_text(size = 28),
        plot.title = element_text(size = 32)
        ) +
  labs(title = "Mitochondrial to Wolbachia genome coverage ratio", y = "Coverage Ratio")
m_wb

ggarrange(n_m, n_wb, m_wb, nrow = 3)
ggsave("cov_ratios.png", height=15, width=20)
ggsave("cov_ratios.tif", height=15, width=20, dpi = 300)




# list file names
file_names.window <- list.files(path = "C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/R2_Sep25/coverage/input",pattern = ".100000_window.cov")

# load data using file names, and make a formatted data frame
setwd("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/R2_Sep25/coverage/input")

data <- purrr::map_df(file_names.window, function(x) {
  data <- read.delim(x, header = F, sep="\t")
  data$V1 <- str_replace(data$V1, 'dirofilaria_immitis_', '')
  data <- tibble::rowid_to_column(data, "NUM")
  cbind(sample_name = gsub(".100000_window.cov","",x), data)
})
colnames(data) <- c("ID", "NUM", "CHR", "START", "END", 
                    "RAW_COVERAGE", "PROPORTION_COVERAGE")







# D. immitis coverage:

# remove scaffolds, mitochondrial and wolbachia genome
data_nuc <- dplyr::filter(data, !grepl("scaffold|MtDNA|Wb|NC|NW",CHR)) %>%
  group_by(ID) %>%
  mutate(RowNo = row_number()) %>%
  ungroup()


# data$SEX <- str_detect(data$SCAF,"Trichuris_trichiura_1_")


# plot the general cov for each sample
ggplot(data_nuc, aes(RowNo, PROPORTION_COVERAGE/(median(PROPORTION_COVERAGE)), group = ID, col = CHR)) +
  geom_point(size=0.5) +
  labs( x = "Genome position" , y = "Relative coverage per 100kb window") +
  theme_bw() + theme(strip.text.x = element_text(size = 6)) +
  facet_wrap(~ID, scales = "free_y")
# lots of samples, can't fit into a single plot. Break up into 2 separate plots.

# part1
ggplot(subset(data_nuc, ID %in% unique(ID)[1:72]), aes(RowNo, PROPORTION_COVERAGE/(median(PROPORTION_COVERAGE)), group = ID, col = CHR)) +
  geom_point(size=0.5) +
  labs( x = "Window number" , y = "Relative coverage per 100kb window") +
  theme_bw() + theme(strip.text.x = element_text(size = 8),
                     axis.text.x = element_text(size = 10),
                     axis.text.y = element_text(size = 10),
                     axis.title.x = element_text(size = 20),
                     axis.title.y = element_text(size = 20),
                     plot.title = element_text(size = 24),
                     legend.title = element_text(size = 20),
                     legend.text = element_text(size = 18)) +
  facet_wrap(~ID, scales = "free_y", nrow = 12, ncol = 6)

ggsave("../ALL_genomewide_RELATIVEcoverage_allsamples_part1.png", height=16, width=12)
ggsave("../ALL_genomewide_RELATIVEcoverage_allsamples_part1.tif", height=16, width=12, dpi = 300)

# part2
ggplot(subset(data_nuc, ID %in% unique(ID)[73:138]), aes(RowNo, PROPORTION_COVERAGE/(median(PROPORTION_COVERAGE)), group = ID, col = CHR)) +
  geom_point(size=0.5) +
  labs( x = "Window number" , y = "Relative coverage per 100kb window") +
  theme_bw() + theme(strip.text.x = element_text(size = 8),
                     axis.text.x = element_text(size = 10),
                     axis.text.y = element_text(size = 10),
                     axis.title.x = element_text(size = 20),
                     axis.title.y = element_text(size = 20),
                     plot.title = element_text(size = 24),
                     legend.title = element_text(size = 20),
                     legend.text = element_text(size = 18)) +
  facet_wrap(~ID, scales = "free_y", nrow = 12, ncol = 6)

ggsave("../ALL_genomewide_RELATIVEcoverage_allsamples_part2.png", height=16, width=12)
ggsave("../ALL_genomewide_RELATIVEcoverage_allsamples_part2.tif", height=16, width=12, dpi = 300)





# the above plots show the proportion of coverage relative to the median coverage of all samples. Let's also just look at the coverage by itself (not in relation to anything else). Set the same y-axes for comparison between samples (some data points will be cut off but that's ok, it's just to visualise).

# part1
ggplot(subset(data_nuc, ID %in% unique(ID)[1:72]), aes(RowNo, PROPORTION_COVERAGE, group = ID, col = CHR)) +
  geom_point(size=0.5) +
  labs( x = "Window number" , y = "Coverage per 100kb window") +
  theme_bw() + theme(strip.text.x = element_text(size = 8),
                     axis.text.x = element_text(size = 10),
                     axis.text.y = element_text(size = 10),
                     axis.title.x = element_text(size = 20),
                     axis.title.y = element_text(size = 20),
                     plot.title = element_text(size = 24),
                     legend.title = element_text(size = 20),
                     legend.text = element_text(size = 18)) +
  ylim(0, 200) +
  facet_wrap(~ID, scales = "fixed", nrow = 12, ncol = 6)

ggsave("../ALL_genomewide_PERBPcoverage_allsamples_part1.png", height=16, width=12)
ggsave("../ALL_genomewide_PERBPcoverage_allsamples_part1.tif", height=16, width=12, dpi = 300)

# part2
ggplot(subset(data_nuc, ID %in% unique(ID)[73:138]), aes(RowNo, PROPORTION_COVERAGE, group = ID, col = CHR)) +
  geom_point(size=0.5) +
  labs( x = "Window number" , y = "Coverage per 100kb window") +
  theme_bw() + theme(strip.text.x = element_text(size = 8),
                     axis.text.x = element_text(size = 10),
                     axis.text.y = element_text(size = 10),
                     axis.title.x = element_text(size = 20),
                     axis.title.y = element_text(size = 20),
                     plot.title = element_text(size = 24),
                     legend.title = element_text(size = 20),
                     legend.text = element_text(size = 18)) +
  ylim(0, 200) +
  facet_wrap(~ID, scales = "fixed", nrow = 12, ncol = 6)

ggsave("../ALL_genomewide_PERBPcoverage_allsamples_part2.png", height=16, width=12)
ggsave("../ALL_genomewide_PERBPcoverage_allsamples_part2.tif", height=16, width=12, dpi = 300)







# Let's see only the chrX to explore the sex of the sample
#Plotting with the chr1 helps to see differences. Males would have ~ half the coverage of the X chromosome compared to the autosome.

# part 1
data_nuc %>%
  filter(., grepl("chrX|chr1",CHR) & ID %in% unique(ID)[1:72]) %>%
  ggplot(aes(RowNo, PROPORTION_COVERAGE/(median(PROPORTION_COVERAGE)), group = ID, col = CHR)) +
  geom_point(size=0.2) +
  labs( x = "Window number" , y = "Relative coverage per 100kb window") +
  theme_bw() + theme(strip.text.x = element_text(size = 8),
                     axis.text.x = element_text(size = 10),
                     axis.text.y = element_text(size = 10),
                     axis.title.x = element_text(size = 20),
                     axis.title.y = element_text(size = 20),
                     plot.title = element_text(size = 24),
                     legend.title = element_text(size = 20),
                     legend.text = element_text(size = 18)) +
  facet_wrap(~ID, scales = "free_y", nrow = 12, ncol = 6)

ggsave("../chrXtochr1_genomewide_coverage_part1.png", height=16, width=12)
ggsave("../chrXtochr1_genomewide_coverage_part1.tif", height=16, width=12, dpi = 300)

# part2
data_nuc %>%
  filter(., grepl("chrX|chr1",CHR) & ID %in% unique(ID)[73:138]) %>%
  ggplot(aes(RowNo, PROPORTION_COVERAGE/(median(PROPORTION_COVERAGE)), group = ID, col = CHR)) +
  geom_point(size=0.2) +
  labs( x = "Window number" , y = "Relative coverage per 100kb window") +
  theme_bw() + theme(strip.text.x = element_text(size = 8),
                     axis.text.x = element_text(size = 10),
                     axis.text.y = element_text(size = 10),
                     axis.title.x = element_text(size = 20),
                     axis.title.y = element_text(size = 20),
                     plot.title = element_text(size = 24),
                     legend.title = element_text(size = 20),
                     legend.text = element_text(size = 18)) +
  facet_wrap(~ID, scales = "free_y", nrow = 12, ncol = 6)

ggsave("../chrXtochr1_genomewide_coverage_part2.png", height=16, width=12)
ggsave("../chrXtochr1_genomewide_coverage_part2.tif", height=16, width=12, dpi = 300)
```
