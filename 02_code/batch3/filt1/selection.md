# Selection analysis

Finding loci that are under divergent selection - genetic outlier analysis.

## Prep

```bash
DIR=/lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/BATCH3/02_VARIANTS/SELECTION/PCADAPT
VCF=/lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/BATCH3/02_VARIANTS/FINAL_SETS/nuclear_samples3x_missing0.9.chr1to4.recode.RENAMED.vcf
METADATA=/lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/BATCH3/02_VARIANTS/SELECTION/METADATA.txt
NAME=nuclear_samples3x_missing0.9.chr1to4.recode.RENAMED

cd $DIR

# we will need the genetic data to be in several different formats. Let's prepare that now. First we convert the VCF to PLINK, and then to BED.

# Load modules
module load vcftools/0.1.16-c4
module load plink/1.90b6.18--h516909a_0
module load  bcftools/1.14--h88f3f91_0

# Convert VCF to PLINK
bcftools view -H $VCF | cut -f 1 | uniq | awk '{print $0"\t"$0}' > hw.chrom-map.txt
vcftools --vcf $VCF --plink --chrom-map hw.chrom-map.txt --out $NAME.plink

# Convert PLINK to BED
plink --file $NAME.plink --make-bed --noweb --out $NAME --allow-extra-chr
```



## PCAdapt

PCAdapt can detect genetic marker outliers without having population designation using a PCA approach. It uses an ordination approach to find sites in a data set that are outliers with respect to background population structure.

```R
R
'R version 4.1.0'

install.packages("pcadapt")
library("pcadapt")

setwd("/lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/BATCH3/02_VARIANTS/SELECTION/PCADAPT")

# Load data
hw_bed <- "../nuclear_samples3x_missing0.9.chr1to4.recode.RENAMED.bed"
hw_pcadapt <- read.pcadapt(hw_bed, type = "bed")

# Produce k-plot
hw_pcadapt_kplot <- pcadapt(input = hw_pcadapt, K = 20)
pdf("pcadapt_hw_kplot.pdf")
plot(hw_pcadapt_kplot, option = "screeplot")
dev.off()
# a k value of 4 is most appropriate, as this is the value of k after which the curve starts to flatten out more.

hw_pcadapt_pca <- pcadapt(hw_pcadapt, K = 4)
summary(hw_pcadapt_pca)

# Investigate axis projections
poplist.names <- readLines("poplist.names.txt")
print(poplist.names)

pdf("pcadapt_hw_projection1v2.pdf")
plot(hw_pcadapt_kplot, option = "scores", i = 1, j = 2, pop = poplist.names)
dev.off()

pdf("pcadapt_hw_projection5v7.pdf")
plot(hw_pcadapt_kplot, option = "scores", i = 5, j = 7, pop = poplist.names)

# # axes 5 vs 7 - There really isn't much population structure -maybe k = 6 is definitely too much, so we made the correct decision by making k = 4.

# Investigate Manhattan and Q-Qplot
## Manhattan plots are a way to visualize the GWAS (genome-wide association study) p-values (or other statistical values) at each SNP locus along the genome
## Q-Qplots plots are a quick way to check if your residuals are normally distributed.
pdf("pcadapt_hw_manhattan.pdf")
plot(hw_pcadapt_pca, option = "manhattan")
dev.off()
# on the x-axes we have position in the genome, on y-axis we have log10 of the p-value. y= interestingness, the higher it is, the more interesting it is. We will be grabbing a few of the high ones and looking at those as outliers. But we haven't looked at the disribution of p-values themselves.

pdf("pcadapt_hw_qqplot.pdf")
plot(hw_pcadapt_pca, option = "qqplot")
dev.off()
# looks at expected p-values against observed values. P-value starts deviating away pretty quickly but steadily, at the more extreme end are becoming more small (more interesting) compared to what we would expect. The shape of this plot will be very important. It's more extreme than the starling example, but not super extreme.

# Plotting and adjusting the p-values
hw_pcadapt_pvalues <- as.data.frame(hw_pcadapt_pca$pvalues)

library("ggplot2")

pdf("pcadapt_hw_pvalues.pdf")
hist(hw_pcadapt_pca$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "green") # putting p-values into 50 bins from values of 0.0 to 1.0 and seeing how many values exist within that bin. Fairly flat distribution, but spike towards the left and right ends. Will pull out outliers at the left end, small values.
dev.off()

#p-value of 0.05 is so that only 1 time out of 20, will you have accidentally a false positive. That works fairly well when you're working with ecology datasets, but when you're working with genetics data, you're running 5007 diff statistical tests, each SNP is running through its wown test to see whether it's an outlier. Having 5% false positive threshold seems like it's not really high enough, you might end up having 10s or 100s of false positive SNPs and that'll drown out the true signal in your data. So that's why you have p-value corrections and multiple testing. Here we're doing that...



hw_pcadapt_padj <- p.adjust(hw_pcadapt_pca$pvalues,method="bonferroni") # most common corrrection is bonferroni, it's pretty strict but not too strict.
alpha <- 0.1 # this is basically p-value on top of p-value. We will make it a bit higher so that we're not too strict with our correction. It reduces false positives. It's traditional - if you go higher to 0.2 that's less strict, if you go down to 0.05 that's more strict. alpha of 0.1 is likely the one that reviewers don't question, if you choose higher or lower you may need to justify why.
outliers <- which(hw_pcadapt_padj < alpha)
length(outliers) 
'[1] 144'
# I have 144. They have allele freq that are really wild, they're weird compared to the rest of the genome, we're not associating them with any phenotype. Could get where the individuals are sampled on a map, get the allele freq (more red, more blue), and perhaps see a gradient - see if outliers based on distance.

write.table(outliers, file="hw_pcadapt_outliers.txt") #save the position of these outlier SNPs
```

```bash
# Mapping outliers: PCAdapt

#The first thing we will do is create a list of SNPs in VCF, and then assign line numbers that can be used to find matching line numbers in outliers (SNP IDs are lost in PCadapt & Bayescan, line numbers are used as signifiers).

cd ..
grep -v "^#" $VCF | cut -f1-3 | awk '{print $0"\t"NR}' > hw_populations_SNPs.txt # -v is the opposite (don't grab)
# ^ lines that start with
# so grab any line that doesn't start with #. We're basically chopping off the header lines. prints the existing output and introduced a tab/new column, then prints the row number. So we want list of SNP IDs, position/chromosome, and what order they occur in in our files.
# BUT my HW data doesn't seem to have specific SNP IDs assigned, it just has chromosome and position. So I'll extract field 1 and 2 for chr and pos.

cd PCADAPT
awk '{print $2}' hw_pcadapt_outliers.txt > hw_pcadapt_outliers_numbers.txt 

awk 'FNR==NR{a[$1];next} (($4) in a)' hw_pcadapt_outliers_numbers.txt ../hw_populations_SNPs.txt   | cut -f1,2 > pcadapt_outlierSNPIDs.txt # we want first column to match the 4th. FNR refers to the record number (typically the line number) in the current file. NR refers to the total record number. The operator == is a comparison operator, which returns true when the two surrounding operands are equal. Col 1&2 in our file -obtained from looking at vcf file header.
head pcadapt_outlierSNPIDs.txt
```

## VCFtools windowed Fst

```bash
# Enter data directory
DIR=/lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/BATCH3/02_VARIANTS/SELECTION
cd $DIR
mkdir VCFTOOLS_FST

# Subset metadata
grep "USA" $METADATA | awk '{print $1}' > individuals_USA.txt
grep "Europe" $METADATA | awk '{print $1}' > individuals_Europe.txt
grep "Australia" $METADATA | awk '{print $1}' > individuals_Australia.txt
grep "Asia" $METADATA | awk '{print $1}' > individuals_Asia.txt
grep "CenAmerica" $METADATA | awk '{print $1}' > individuals_CenAmerica.txt

cd VCFTOOLS_FST

# USA vs Australia
vcftools --vcf $VCF --weir-fst-pop $DIR/individuals_USA.txt --weir-fst-pop $DIR/individuals_Australia.txt --out USA_Australia # pairwise analysis, will do 1 population at a time. The fst numbers are higher than the starling example. If it's a low fst value, then maybe the populations haven't been separated very long. If the fst is higher (closer to 1) then there is high population structure.
# Fst is measure of gen differentiation between diff populations, Comparing allele freq between USA and Australia and saying how genetically different are these?

head -n 5 USA_Australia.weir.fst 
wc -l USA_Australia.weir.fst # extra line for the header
'175199' # lost some SNPs from VCFtools?



awk '{print $0"\t"NR}' ./USA_Australia.weir.fst  > USA_Australia.weir.fst.edit # introduced new column numbering each line - but I am not plotting the windowed fst, I am plotting the actual SNPs so I will use the actual SNP position.
# these row numbers don't line up...

R

```

```R
library("ggplot2")

setwd("/lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/BATCH3/02_VARIANTS/SELECTION/VCFTOOLS_FST")

fst <- read.table("USA_Australia.weir.fst.edit", sep="\t", header=TRUE)
str(fst)

quantile(fst$WEIR_AND_COCKERHAM_FST, probs = c(.95, .99, .999), na.rm = TRUE) 
'  95%      99%    99.9%
0.649214 0.808868 0.913254'

# exploring these arbitrary thresholds. If we were to grab top 5% fst we would grab any SNP with fst of 0.649 etc. If you just want to explore patterns, PCA, admixture analysis etc then slightly lower threshold is probs better.

# Choose the quantile threshold above which SNPs will be classified as outliers. In the code block below, we chose 99% (i.e., the top 1% of SNPs).

#plot
pdf("USA_Australia_fst_hw.pdf", width=10, height=10)
ggplot(fst, aes(x=POS, y=WEIR_AND_COCKERHAM_FST)) + 
  geom_point(color = "cornflowerblue", size = 0.2) + 
  geom_hline(yintercept=0.81, linetype="dashed", color = "red")+
  labs(x = "SNP position") +
  ggtitle("USA & Australia") +
  theme(axis.text = element_text(size = 10),
                axis.title = element_text(size = 12, face = "bold"),
                plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
                legend.position = "none",
                panel.background = element_rect(fill = "white"),
                plot.background = element_rect(fill = "white")) +
  facet_wrap(~ CHROM, scales = "free") +
  scale_x_continuous(labels = scales::number_format()) # don't want x axis labels to be in scientific notation

dev.off()

q()
# Anything above that line is an outlier SNP.Each dot is fst
```

```bash
# Select all outlier SNPs
cat USA_Australia.weir.fst.edit | awk '$3>0.808868 && $3 != "-nan"' > USA_Australia.outliers
wc -l USA_Australia.outliers # how many outlier windows we have in total


#Identify the regions of the genome that were found to be outliers and subset the VCF file 
awk 'NR > 1 {print $1 "\t" $2 "\t" $2}' USA_Australia.outliers > USA_Australia.outliers.bed #created a bed file

vcftools --vcf $VCF --bed USA_Australia.outliers.bed --out USA_Australia_fst_outliers --recode # expects a bed file, subset vcf file to only have snps that are outliers
'  Read 1580 BED file entries.
After filtering, kept 1579 out of a possible 175415 Sites
'

#Create a list of outlier SNP IDs
grep -v "#" USA_Australia_fst_outliers.recode.vcf | awk '{print $1 "\t" $2}' > USA_Australia_fst_outlier_pos.txt # from this list of SNP files, just grab the SNP positions. Cutting off all the header. There's 1579 SNPs.
wc -l USA_Australia_fst_outlier_pos.txt # yep we've got 61, which is what we were expecting.
```



# Bayescan

BayeScan identifies outlier SNPs based on allele frequencies.

We need to convert into these old Bayescan filetypes, but we still use them very frequently. To do this, we will use PGDspider.

```bash
cd ..
mkdir BAYESCAN
cut -f1,2 $METADATA > hw_metadata_INDPOP.txt
cd BAYESCAN
nano VCF_PGD.spid

'# VCF Parser questions 
PARSER_FORMAT=VCF
# Only output SNPs with a phred-scaled quality of at least: 
VCF_PARSER_QUAL_QUESTION=
# Select population definition file:
VCF_PARSER_POP_FILE_QUESTION=../hw_metadata_INDPOP.txt
# What is the ploidy of the data?
VCF_PARSER_PLOIDY_QUESTION=DIPLOID
# Do you want to include a file with population definitions?
VCF_PARSER_POP_QUESTION=true
# Output genotypes as missing if the phred-scale genotype quality is below:
VCF_PARSER_GTQUAL_QUESTION=
# Do you want to include non-polymorphic SNPs?
VCF_PARSER_MONOMORPHIC_QUESTION=false
# Only output following individuals (ind1, ind2, ind4, ...):
VCF_PARSER_IND_QUESTION=
# Only input following regions (refSeqName:start:end, multiple regions: whitespace separated):
VCF_PARSER_REGION_QUESTION=
# Output genotypes as missing if the read depth of a position for the sample is below:
VCF_PARSER_READ_QUESTION=
# Take most likely genotype if "PL" or "GL" is given in the genotype field?
VCF_PARSER_PL_QUESTION=false
# Do you want to exclude loci with only missing data?
VCF_PARSER_EXC_MISSING_LOCI_QUESTION=false

# PGD Writer questions
WRITER_FORMAT=PGD'

# this is pretty standard code, but you will need to update the file locations. Press ctrl+x, y, enter


#load module
module load quay.io/biocontainers/pgdspider/2.1.1.5--hdfd78af_1/module

#run the step 1 of the conversion
PGDSpider2-cli -inputfile $VCF -inputformat VCF -outputfile starling_3populations.pgd -outputformat  PGD -spid VCF_PGD.spid #first it converts to universal PHD spdier format, then have long list of about 30 diff types you can convert into. This one covnerts to intermediate input file. We provide a spid so it knows which individuals belong to which populations. "PGDSpider configuration file not found! Loading default configuration." --> that's ok, we gave it some information.

#run the step 2 of the conversion
PGDSpider2-cli -inputfile starling_3populations.pgd -inputformat PGD -outputfile starling_3populations.bs -outputformat GESTE_BAYE_SCAN #there's a big list of what you need to convert it to for bayescan format. They have automatically made a new spid file for us with default setting which is good enough for us, we don't need to tell it that it's diploidy etc. They made a template for us.

head starling_3populations.bs
#Looks like this:
'[loci]=5007

[populations]=3

[pop]=1
 1      12      2       9 3
 2      20      2       11 9
 3      18      2       15 3
 4      20      2       0 20
 5      22      2       2 20
'
#kind of doesn't matter in the context of this program which is defined ancestor allele which one is ref/alt etc. But usually 1st = ref, 2nd = alternate. 


#load bayescan
module load quay.io/biocontainers/bayescan/2.0.1--h9f5acd7_4

#run bayescan. 
bayescan2 ./starling_3populations.bs -od ./ -threads 2 -n 5000 -thin 10 -nbp 20 -pilot 5000 -burn 50000 -pr_odds 10
# Bayesian works cyclically, it looks at experimental evidence and feedbacks to improve its prediction - convergence - how we know if it works, we feed into it and it gets smarter and smarter. Prediciton of which sites in your genome are likely to be under selection and which ones are not. If it doesn't reach covnergece, it just can't reach an answer.
# Typically ignore the first 50,000 iterations because it takes a while for the program to find tis feet and start converging onto something that makes sense.

# If you're having trouble converging, change the chain. Only fiddle with these if you know what you're doing and have a reason.

# takes a while to run, so we're cutting it off and using pre-cooked data
cp $DIR/workshop_material/Ev1_SelectionMetaAnalysis/workshop_files/backup_files/starling_3population_AccRte.txt .
cp $DIR/workshop_material/Ev1_SelectionMetaAnalysis/workshop_files/backup_files/starling_3population_Verif.txt .
cp $DIR/workshop_material/Ev1_SelectionMetaAnalysis/workshop_files/backup_files/starling_3population_fst.txt .
cp $DIR/workshop_material/Ev1_SelectionMetaAnalysis/workshop_files/backup_files/starling_3population.sel .

wc -l starling_3population.sel
'5000'
# representing fsts across all diff cyles

less starling_3population_AccRte.txt
# things are starting to stabilise, oscillated up and down. If it doesn't stagnate, sign that something has gone wrong with the covnergece.

wc -l starling_3population_fst.txt
'5008' each line must correspond to a SNP

R
```