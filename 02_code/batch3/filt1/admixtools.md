# ADMIXTOOLS and admixr

## Convert vcf file to eigenstrat format

Following these tutorials:
- https://speciationgenomics.github.io/ADMIXTOOLS_admixr/
- https://bodkan.net/admixr/articles/tutorial.html

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

cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/BATCH3/02_VARIANTS/ADMIXTOOLS

# create link
ln -s ../FINAL_SETS/nuclear_samples3x_missing0.9.chr1to4.recode.RENAMED.vcf

# generate "chrom-map.txt" for chromosome names
bcftools view -H nuclear_samples3x_missing0.9.chr1to4.recode.RENAMED.vcf | cut -f 1 | uniq | awk '{print $0"\t"$0}' > chrom-map.txt
## This command appends each chromosome name ($0) with a tab character ("\t") and then repeats the chromosome name again ($0) before redirecting the output to "chrom-map.txt".


# Convert vcf to eigenstrat format
chmod a+x convertVCFtoEigenstrat_hw.sh
./convertVCFtoEigenstrat_hw.sh nuclear_samples3x_missing0.9.chr1to4.recode.RENAMED
## modify to allow non-standard chromosomes because eigenstrat format expects chromosomes to be called "1, 2, 3 etc". Indicate a "chom-map.txt" file in the vcftools command.
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

# Q1: Do pops from Central America show evidence of closer genetic affinity with pops from Europe compared to elsewhere?
##Using the admixr package we can then calculate our D statistic simply by running:
result <- d(W = "", X = "", Y = "", Z = "", data = snps)
# Z is the outgroup


# Q2: Do Australian heartworms show evidence of closer genetic affinity with pops from Asia?
result <- d(W = "Australia", X = "Romania", Y = "Thailand", Z = "Panama", data = snps)
# Z is the outgroup

head(result)

#plot
plot_dstat <- ggplot(result, aes(fct_reorder(W, D), D, color = abs(Zscore) > 2)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_errorbar(aes(ymin = D - 2 * stderr, ymax = D + 2 * stderr))

ggsave("plot_dstat.pdf", plot = plot_dstat)

#We can see that the D values for Australia are significantly different from 0, meaning that the data rejects the null hypothesis of no Thailand ancestry in Australians. So this suggests that Thailand admixed with the ancestors of present-day Australian worms.






# f4 statistic 





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

