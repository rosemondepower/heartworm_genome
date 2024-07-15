# ADMIXTOOLS via admixr

## Convert vcf file to eigenstrat format

Following these tutorials:
- https://speciationgenomics.github.io/ADMIXTOOLS_admixr/
- https://bodkan.net/admixr/articles/01-tutorial.html

ADMIXTOOLS requires eigenstrat input files. To convert a vcf file to eigenstrat, we will use a conversion script written by Joana. That script just requires the name of the vcf file (NB - the script can use uncompressed vcfs too).
- downloaded convertVCFtoEigenstrat.sh from: https://github.com/joanam/scripts/blob/e8c6aa4b919b58d69abba01e7b7e38a892587111/convertVCFtoEigenstrat.sh

Also looked at Steve's ancient trichuris code for help getting the eigenstrat format files.

I am doing this on the Sanger farm.

```bash
# create conda environment
conda create --name admixtools
conda init --all
conda activate admixtools

# install admixtools and eigensoft
conda install bioconda::admixtools
conda install bioconda::eigensoft

# load other modules
module load bsub.py/0.42.1
#module load compilers/gcc/12.2
module load vcftools/0.1.16-c4
#module load perl-velvetoptimiser/2.2.6
module load plink/1.90b6.18--h516909a_0
module load bcftools/1.14--h88f3f91_0

cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/04_VARIANTS/SCRATCH/FILTER2/WITH_OUTGROUPS/ADMIXTOOLS

# create link
ln -s /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/04_VARIANTS/SCRATCH/FILTER2/WITH_OUTGROUPS/FINAL_SETS/nuclear_samples3x_missing0.8.chr1to4.recode.vcf

# generate "chrom-map.txt" for chromosome names
bcftools view -H nuclear_samples3x_missing0.8.chr1to4.recode.vcf | cut -f 1 | uniq | awk '{print $0"\t"$0}' > chrom-map.txt
## This command appends each chromosome name ($0) with a tab character ("\t") and then repeats the chromosome name again ($0) before redirecting the output to "chrom-map.txt".


# Convert vcf to eigenstrat format
chmod a+x convertVCFtoEigenstrat_hw.sh
bsub.py 4 eigenstrat "./convertVCFtoEigenstrat_hw.sh nuclear_samples3x_missing0.8.chr1to4.recode"
## modify to allow non-standard chromosomes because eigenstrat format expects chromosomes to be called "1, 2, 3 etc". Indicate a "chom-map.txt" file in the vcftools command.
# my file wasn't working, had a look at Steve's to see where I'm going wrong

# this produced 3 important files I need for admixtools:
# .snp
# .eigenstratgeno # i changed this to .geno so that eigenstrat would recognise it
# .ind - modified col 3 of .ind file to have country groups. Make sure it is all 1 word e.g. needed to change "Costa Rica" to Costa_Rica".
```

## f3 stats

```bash
# make a new populations file
admixtools_pops.txt

# URSI outgroup
# set the outgroup
OUTGROUP=URSI

# loop throguh the populations to generate the pop file as input to admixtools
for i in URSI AUS ASIA EUR CENAM USA; do
     for j in  URSI AUS ASIA EUR CENAM USA; do
          if [[ "$i" == "$j" ]] || [[ "$i" == "$OUTGROUP" ]] || [[ "$j" == "$OUTGROUP" ]]; then
               :
               else
               echo -e "${i}\t${j}\t${OUTGROUP}" >> admixtools_pops.txt;
          fi;
     done;
done

# remove any duplicates



# REPENS outgroup
# set the outgroup
OUTGROUP=REPENS

# loop throguh the populations to generate the pop file as input to admixtools
for i in REPENS AUS ASIA EUR CENAM USA; do
     for j in  REPENS AUS ASIA EUR CENAM USA; do
          if [[ "$i" == "$j" ]] || [[ "$i" == "$OUTGROUP" ]] || [[ "$j" == "$OUTGROUP" ]]; then
               :
               else
               echo -e "${i}\t${j}\t${OUTGROUP}" >> admixtools_pops.txt;
          fi;
     done;
done



# THAI outgroup
# set the outgroup
OUTGROUP=THAI

# loop throguh the populations to generate the pop file as input to admixtools
for i in THAI AUS ASIA EUR CENAM USA; do
     for j in  THAI AUS ASIA EUR CENAM USA; do
          if [[ "$i" == "$j" ]] || [[ "$i" == "$OUTGROUP" ]] || [[ "$j" == "$OUTGROUP" ]]; then
               :
               else
               echo -e "${i}\t${j}\t${OUTGROUP}" >> admixtools_pops.txt;
          fi;
     done;
done




# run admixtools to generate f3 stats
qp3Pop -p PARAMETER_FILE > qp3Pop.out
```

PARAMETER_FILE:

```bash
genotypename:   nuclear_samples3x_missing0.8.chr1to4.recode.eigenstratgeno (in eigenstrat format)
snpname:        nuclear_samples3x_missing0.8.chr1to4.recode.snp      (in eigenstrat format)
indivname:      nuclear_samples3x_missing0.8.chr1to4.recode.ind    (in eigenstrat format)
popfilename:    admixtools_pops.txt
inbreed: YES
```

```bash
grep "result" qp3Pop.out | awk '{print $2,$3,$4,$5,$6,$7,$8}' OFS="\t" > qp3Pop.clean.out
# didnt work using these other species as outgroups. Maybe just compare within the D. immitis.
```


```bash
# AUS as outgroup
# set the outgroup
OUTGROUP=AUS

# loop throguh the populations to generate the pop file as input to admixtools
for i in AUS ASIA EUR CENAM USA; do
     for j in  AUS ASIA EUR CENAM USA; do
          if [[ "$i" == "$j" ]] || [[ "$i" == "$OUTGROUP" ]] || [[ "$j" == "$OUTGROUP" ]]; then
               :
               else
               echo -e "${i}\t${j}\t${OUTGROUP}" >> admixtools_pops.txt;
          fi;
     done;
done

# ASIA as outgroup
# set the outgroup
OUTGROUP=ASIA

# loop throguh the populations to generate the pop file as input to admixtools
for i in AUS ASIA EUR CENAM USA; do
     for j in  AUS ASIA EUR CENAM USA; do
          if [[ "$i" == "$j" ]] || [[ "$i" == "$OUTGROUP" ]] || [[ "$j" == "$OUTGROUP" ]]; then
               :
               else
               echo -e "${i}\t${j}\t${OUTGROUP}" >> admixtools_pops.txt;
          fi;
     done;
done

# EUR as outgroup
# set the outgroup
OUTGROUP=EUR

# loop throguh the populations to generate the pop file as input to admixtools
for i in AUS ASIA EUR CENAM USA; do
     for j in  AUS ASIA EUR CENAM USA; do
          if [[ "$i" == "$j" ]] || [[ "$i" == "$OUTGROUP" ]] || [[ "$j" == "$OUTGROUP" ]]; then
               :
               else
               echo -e "${i}\t${j}\t${OUTGROUP}" >> admixtools_pops.txt;
          fi;
     done;
done

# CENAM as outgroup
# set the outgroup
OUTGROUP=CENAM

# loop throguh the populations to generate the pop file as input to admixtools
for i in AUS ASIA EUR CENAM USA; do
     for j in  AUS ASIA EUR CENAM USA; do
          if [[ "$i" == "$j" ]] || [[ "$i" == "$OUTGROUP" ]] || [[ "$j" == "$OUTGROUP" ]]; then
               :
               else
               echo -e "${i}\t${j}\t${OUTGROUP}" >> admixtools_pops.txt;
          fi;
     done;
done

# USA as outgroup
# set the outgroup
OUTGROUP=USA

# loop through the populations to generate the pop file as input to admixtools
for i in AUS ASIA EUR CENAM USA; do
     for j in  AUS ASIA EUR CENAM USA; do
          if [[ "$i" == "$j" ]] || [[ "$i" == "$OUTGROUP" ]] || [[ "$j" == "$OUTGROUP" ]]; then
               :
               else
               echo -e "${i}\t${j}\t${OUTGROUP}" >> admixtools_pops.txt;
          fi;
     done;
done

# remove duplicates

# run admixtools to generate f3 stats
qp3Pop -p PARAMETER_FILE > qp3Pop.out

```



## Plot f3 stats

```R
# load libraries
library(tidyverse)

setwd("")

# read data
data <- read.delim("qp3Pop.clean.out", header=F, sep="\t")

# fix headings
colnames(data) <- c("Source_1", "Source_2", "Outgroup", "f_3", "std_err", "Z_score", "SNPs")

# make a plot
ggplot(data,aes(f_3, reorder(paste0(Source_1,",",Source_2), -f_3))) +
     geom_point(size = 2) +
     geom_segment(aes(x = f_3-std_err, y = paste0(Source_1,",",Source_2), xend = f_3+std_err, yend = paste0(Source_1,",",Source_2))) +
     theme_bw() + xlim(0,1) +
     labs(x = "f3(Source1,Source2;Outgroup)" , y = "") +
     facet_grid(Outgroup~., scale="free_y", space = "free_y")

# save it
ggsave("plot_admixtools_f3_statistics.png")
ggsave("plot_admixtools_f3_statistics.pdf", height = 4, width = 5, useDingbats = FALSE)
```

## D stats

```bash
# loop through the populations to generate the pop file as input to admixtools
OUTGROUP=URSI
for i in AUS ASIA EUR CENAM USA; do
     for j in  AUS ASIA EUR CENAM USA; do
          for k in AUS ASIA EUR CENAM USA; do
               if [[ "$i" == "$j" ]] || [[ "$i" == "$OUTGROUP" ]] || [[ "$j" == "$OUTGROUP" ]] || [[ "$i" == "$k" ]] || [[ "$j" == "$k" ]] || [[ "$k" == "$OUTGROUP" ]]; then
                    :
                    else
                    echo -e "${i}\t${j}\t${k}\t${OUTGROUP}" >> admixtools_pops_4.txt;
               fi;
          done;
     done;
done

# run admixtools to generate D-stats
qpDstat -p PARAMETER_FILE_dstat > qpDstat.out
```


## f4 ratio

```bash
# run admixtools to generate f4 ratios
qpF4ratio -p PARAMETER_FILE_dstat >qpF4ratio.out
```























# Sanger isn't connecting, try in Artemis for now

```bash
module load anaconda3

# create conda environment
conda create --name admixr
conda init --all
```
```bash
#!/bin/bash

# PBS directives 
#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N admixr
#PBS -l select=1:ncpus=1:mem=10GB
#PBS -l walltime=01:00:00
#PBS -m abe
#PBS -q defaultQ
#PBS -o admixr.txt

module load anaconda3
module load vcftools
module load plink
module load bcftools

conda activate admixr

cd /scratch/RDS-FSC-Heartworm_MLR-RW/batch4

# generate "chrom-map.txt" for chromosome names
bcftools view -H nuclear_samples3x_missing0.8.chr1to4.recode.vcf | cut -f 1 | uniq | awk '{print $0"\t"$0}' > chrom-map.txt
## This command appends each chromosome name ($0) with a tab character ("\t") and then repeats the chromosome name again ($0) before redirecting the output to "chrom-map.txt".


# Convert vcf to eigenstrat format
chmod a+x convertVCFtoEigenstrat_hw.sh
./convertVCFtoEigenstrat_hw.sh nuclear_samples3x_missing0.8.chr1to4.recode
## modify to allow non-standard chromosomes because eigenstrat format expects chromosomes to be called "1, 2, 3 etc". Indicate a "chrom-map.txt" file in the vcftools command.
# my file wasn't working, had a look at Steve's to see where I'm going wrong

# this produced 3 important files I need for admixtools:
# .snp
# .eigenstratgeno # i changed this to .geno so that eigenstrat would recognise it
# .ind - modified col 3 of .ind file to have country groups. Make sure it is all 1 word e.g. needed to change "Costa Rica" to Costa_Rica".
```



## admixr in R

Admixr is an R package that provides a more straightforward way to use ADMIXTOOLS.

```bash
# open R on command line
R
##R version 4.1.0 (2021-05-18) -- "Camp Pontanezen"
```

```R
# install
install.packages("admixr")
# made personal library '~/R/x86_64-pc-linux-gnu-library/4.1’ and downloaded the admixtools package there
# selected CRAN mirror for use in this session:  2: Australia (Canberra) [https]

# make sure R can find admixtools binaries on the $PATH

export PATH="/nfs/users/nfs_r/rp24/miniconda3/envs/admixtools/bin"

# load packages
library(tidyverse)
library(admixr)
library(ggplot2)

# set data prefix
data_prefix <- "/lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/BATCH3/02_VARIANTS/ADMIXTOOLS/nuclear_samples3x_missing0.9.chr1to4.recode.RENAMED"

# read in data
snps <- eigenstrat(data_prefix)
snps

#count SNPs
##The count_snps function can be useful for quality control, weighting of admixture statistics (D, f4, etc.) in regression analyses etc.
snp_info <- count_snps(snps)
snp_info
'# A tibble: 81 × 4
   id        sex   label     present
   <chr>     <chr> <chr>       <int>
 1 ERR034940 U     USA        155325
 2 ERR034941 U     Italy      171468
 3 Dog2.1    U     Australia  174564
 4 Dog2.2    U     Australia  174348
 5 Dog2.3    U     Australia  174557
 6 Dog2.4    U     Australia  174555
 7 JS6281    U     Australia  174585
 8 JS6342    U     Australia  174374
 9 JS6343    U     Australia  174540
10 JS6344    U     Australia  174533
# … with 71 more rows
'

# Now we can calculate some statistics...




# D statistic (ABBA-BABA test)

## "Significant departure of D from zero indicates an excess of allele sharing between the first and the third population (positive D), or an excess of allele sharing between the second and the third population (negative D). If we get D that is not significantly different from 0, this suggests that the first and second populations form a clade, and don’t differ in the rate of allele sharing with the third population (this is the null hypothesis that the data is compared against)."


# Q1: Do Australian heartworms and Thailand heartworms share more alleles with one another than expected by chance? Use D. repens as outgroup.
##Using the admixr package we can then calculate our D statistic simply by running:

# D(non-Europe pops W,African,Neanderthal,Chimp).
result <- d(W = "Australia", X = "USA", Y = "Thailand", Z = "D.repens", data = snps)
# Z is the outgroup


# Q2: Do Central American heartworms show evidence of closer genetic affinity with pops from Europe?
result <- d(W = "Panama", X = "Thailand", Y = "Romania", Z = "D.repens", data = snps)
# Z is the outgroup

head(result)

# plot
plot_dstat <- ggplot(result, aes(fct_reorder(W, D), D, color = abs(Zscore) > 2)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_errorbar(aes(ymin = D - 2 * stderr, ymax = D + 2 * stderr))

ggsave("plot_dstat.pdf", plot = plot_dstat)

# flip axes
plot_dstat2 <- ggplot(result, aes(fct_reorder(W, D), D, color = abs(Zscore) > 2)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_errorbar(aes(ymin = D - 2 * stderr, ymax = D + 2 * stderr)) +
  coord_flip()

ggsave("plot_dstat2.pdf", plot = plot_dstat2)

## We can see that the D values for Australia are significantly different from 0, meaning that the data rejects the null hypothesis of no Thailand ancestry in Australians. So this suggests that Thailand admixed with the ancestors of present-day Australian worms.






# f4 statistic 

## Very similar to D statistic
## "Significant departure of f4 from 0 can be interpreted as evidence of gene flow."
## "You might be wondering why we have both f4 and D if they are so similar. The truth is that f4 is, among other things, directly informative about the amount of shared genetic drift (“branch length”) between pairs of populations, which is a very useful theoretical property. Other than that, it’s often a matter of personal preference and so admixr provides functions for calculating both."

result <- f4(W = pops, X = "Yoruba", Y = "Vindija", Z = "Chimp", data = snps)

head(result)

'

'

# f4-ratio statistic

## If the D and f4 stats told us that a certain population carries ancestry from another population, what if we want to know how MUCH of this ancestry they have? What proportion of their genomes is of this origin?

result <- f4ratio(X = pops, A = "Altai", B = "Vindija", C = "Yoruba", O = "Chimp", data = snps)

head(result)
# The ancestry proportion (a number between 0 and 1) is given in the alpha column

ggplot(result, aes(fct_reorder(X, alpha), alpha, color = abs(Zscore) > 2)) +
  geom_point() +
  geom_errorbar(aes(ymin = alpha - 2 * stderr, ymax = alpha + 2 * stderr)) +
  geom_hline(yintercept = 0, linetype = 2) +
  labs(y = "Neandertal ancestry proportion", x = "present-day individual")




# f3 statistic

## Let's say we are interested in relative divergence times between pairs of HW populations, and we want to know in which approximate order they split off from each other. To address this question, we can use f3 statistic by fixing the C outgroup as D. repens, and calculating pairwise f3 stats between all D. immitis populations.

pops <- c("Australia", "USA", "Thailand", "Panama", "Costa_Rica", "Italy", "Romania")

result <- f3(A = pops, B = pops, C = "D_repens", data = snps)
#or
result <- f3(A = "Europe", B = pops, C = "D_repens", data = snps)

head(result)

# sort the population labels according to an increasing f3 value relative to French/lowest one
ordered <- filter(result, A == "Mbuti", B != "Mbuti") %>% arrange(f3) %>% .[["B"]] %>% c("Mbuti")

# plot heatmap of pairwise f3 values
result %>%
  filter(A != B) %>%
  mutate(A = factor(A, levels = ordered),
         B = factor(B, levels = ordered)) %>%
  ggplot(aes(A, B)) + geom_tile(aes(fill = f3))



```






















# Convert PLINK format to .geno, .ind, and .snp files
plink --bfile hw --recodeAD --out hw
## This command will produce hw.raw file containing genotype information.

./plink2geno.py hw.raw > hw.geno
./plink2ind.py hw.raw > hw.ind
./plink2snp.py hw.raw > hw.snp


# Convert plink format to Eigenstrat format
convertf -p par.PED.EIGENSTRAT
#par.PED.EIGENSTRAT:
'
genotypename:	data_hw.bed
snpname:	data_hw.bim
indivname:	data_hw.fam
outputformat:	EIGENSTRAT
genotypeoutname:	data_hw.geno
snpoutname:	data_hw.snp
indivoutname:	data_hw.ind
'

# some samples are being ignored because eigenstrat format requires sex information. Add sex info into the .fam file column 6.



# this produces 3 important files I need for admixtools:
# .snp
# .eigenstratgeno
# .ind
```

chmod a+x nuclear_samples3x_missing0.9.chr1to4.recode.RENAMED.ind
chmod a+x nuclear_samples3x_missing0.9.chr1to4.recode.RENAMED.snp
chmod a+x nuclear_samples3x_missing0.9.chr1to4.recode.RENAMED.eigenstratgeno

