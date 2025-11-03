# ADMIXTOOLS

## Convert vcf file to eigenstrat format

Useful tutorials:
- https://speciationgenomics.github.io/ADMIXTOOLS_admixr/
- https://bodkan.net/admixr/articles/01-tutorial.html

ADMIXTOOLS requires eigenstrat input files. To convert a vcf file to eigenstrat, we will use a conversion script written by Joana. That script just requires the name of the vcf file (NB - the script can use uncompressed vcfs too).
- downloaded convertVCFtoEigenstrat.sh from: https://github.com/joanam/scripts/blob/e8c6aa4b919b58d69abba01e7b7e38a892587111/convertVCFtoEigenstrat.sh


```bash
# create conda environment
module load ISG/conda
conda create --name admixtools
conda init --all
conda activate admixtools

# install admixtools and eigensoft
conda install bioconda::admixtools
conda install bioconda::eigensoft

# load other modules
module load PaM/environment
module load bsub.py/0.42.1
#module load compilers/gcc/12.2
module load vcftools/0.1.16-c4
#module load perl-velvetoptimiser/2.2.6
module load plink/1.90b6.18--h516909a_0
module load bcftools/1.14--h88f3f91_0

cd /lustre/scratch127/pam/teams/team333/rp24/dirofilaria_immitis/R2/dstats

# create link
zcat ../treemix/treemix.vcf.gz > admixtools.vcf

# generate "chrom-map.txt" for chromosome names
bcftools view -H admixtools.vcf | cut -f 1 | uniq | awk '{print $0"\t"$0}' > chrom-map.txt
## This command appends each chromosome name ($0) with a tab character ("\t") and then repeats the chromosome name again ($0) before redirecting the output to "chrom-map.txt".


# Convert vcf to eigenstrat format
chmod a+x convertVCFtoEigenstrat_hw.sh
bsub.py 4 eigenstrat "./convertVCFtoEigenstrat_hw.sh admixtools"
## modify to allow non-standard chromosomes because eigenstrat format expects chromosomes to be called "1, 2, 3 etc". Indicate a "chom-map.txt" file in the vcftools command.

# this produced 3 important files I need for admixtools:
# .snp
# .eigenstratgeno
# .ind - modified col 3 of .ind file to have continent groups. Make sure it is all 1 word, otherwise things won't work.
```


## D statistic

The 4-population test, implemented here as D-statistics, is also a formal test for admixture based on a four taxon 4 statistic, which can provide some information about the direction of gene ﬂow.
For any 4 populations (W, X, Y, Z), qpDstat computes the D-statistics as - 
num = (w − x)(y − z )
den = (w + x − 2wx)(y + z − 2yz )

D = num/ den

The output of qpDstat is informative about the direction of gene flow. So for 4 populations (W, X, Y, Z) as follows - 
If the Z-score is +ve, then the gene flow occured either between W and Y or X and Z 
If the Z-score is -ve, then the gene flow occured either between W and Z or X and Y. 

Run using different outgroups:

- D. ursi
- D. repens
- D. 'Thai'

And run using all the D. immitis pops from different continents as source populations.

```bash
# loop through the populations to generate the pop file as input
OUTGROUP=URSI
for i in AUS ASIA EUR CENAM USA; do
     for j in  AUS ASIA EUR CENAM USA; do
          for k in AUS ASIA EUR CENAM USA; do
               if [[ "$i" == "$j" ]] || [[ "$i" == "$OUTGROUP" ]] || [[ "$j" == "$OUTGROUP" ]] || [[ "$i" == "$k" ]] || [[ "$j" == "$k" ]] || [[ "$k" == "$OUTGROUP" ]]; then
                    :
                    else
                    echo -e "${k}\t${j}\t${i}\t${OUTGROUP}" >> admixtools_pops_dstat_ursi.txt;
               fi;
          done;
     done;
done

# run admixtools to generate D-stats
qpDstat -p PARAMETER_FILE_dstat_ursi > qpDstat_ursi.out

# PARAMETER_FILE_dstat:
genotypename:   admixtools.eigenstratgeno (in eigenstrat format)
snpname:        admixtools.snp      (in eigenstrat format)
indivname:      admixtools.ind    (in eigenstrat format)
popfilename:    admixtools_pops_dstat_ursi.txt
inbreed: YES
printsd: YES

# extract the relevant output lines
grep "result" qpDstat_ursi.out | awk '{print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11}' OFS="\t" > qpDstat_ursi.clean.out
```

Plot in R:

```R
# Plot D statistics

# load libraries
library(tidyverse)
library(RColorBrewer)
library(ggstance)

setwd("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/R2_Sep25/admixtools")

set.seed(3)
############################### URSI AS OUTGROUP
# read data
data <- read.delim("qpDstat_ursi.clean.out", header=F, sep="\t")

# headings
colnames(data) <- c("Pop1", "Pop2", "Pop3", "Outgroup", "D_stat", "SE", "Z", "BABA", "ABBA", "SNPs")

# specify whether absolute value of Z is > 3
data <- mutate(data, SIG = ifelse(abs(Z) > 3, "YES", "NO"))

# colours
red_palette1 <- brewer.pal(n = 9, name = "YlOrRd")
red_palette2 <- brewer.pal(n=9, name = "YlOrBr")
blue_palette <- brewer.pal(n = 9, name = "Blues")
green_palette1 <- brewer.pal(n = 9, name = "BuGn")
green_palette2 <- brewer.pal(n=9, name = "YlGn")

colour_map <- c("AUS" = "cornflowerblue",
             "ASIA" = "hotpink",
             "EUR" = green_palette1[6],
             "CENAM" = "darkviolet",
             "USA" = red_palette1[6])


# Replace population names with corresponding color names
data$Pop2Color <- colour_map[data$Pop2]
# Adjust colors so that if the abs(Z) < 3, the points appear duller
data$Color <- mapply(function(sig, colour) {
  if (sig == "YES") {
    return(colour)
  } else {
    return(adjustcolor(colour, alpha.f = 0.2))
  }
}, data$SIG, data$Pop2Color)



## POP1: AUS
# subset data
pop3_AUS <- "AUS"
data_AUS <- data[data$Pop3 == pop3_AUS, ]

qpDstat_ursi_plot_AUS <- ggplot(data_AUS, aes(x = D_stat, y = Pop1, color = Color)) +
  geom_point(position = position_dodge(width = 0.3), size = 4) + # Add points
  geom_errorbar(aes(xmin = D_stat - SE, xmax = D_stat + SE, y = Pop1), 
                width = 0.0,
                position = position_dodge(width = 0.3)) + # Add error bars
  geom_vline(xintercept = 0, color = "grey", linetype = "dashed", linewidth = 0.5) + # Add dashed line at x = 0
  theme_bw() +
  labs(x = "D-stat", y = "Pop1") +
  facet_wrap(~Pop3, nrow = 1) +
  scale_color_identity(guide = "legend", name = "Pop2", breaks = colour_map, labels = names(colour_map)) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
qpDstat_ursi_plot_AUS
ggsave("qpDstat_ursi_plot_AUS.png", qpDstat_ursi_plot_AUS, height=8, width=5)




## POP1: ASIA
# subset data
pop3_ASIA <- "ASIA"
data_ASIA <- data[data$Pop3 == pop3_ASIA, ]

qpDstat_ursi_plot_ASIA <- ggplot(data_ASIA, aes(x = D_stat, y = Pop1, color = Color)) +
  geom_point(position = position_dodge(width = 0.3), size = 4) + # Add points
  geom_errorbar(aes(xmin = D_stat - SE, xmax = D_stat + SE, y = Pop1), 
                width = 0.0,
                position = position_dodge(width = 0.3)) + # Add error bars
  geom_vline(xintercept = 0, color = "grey", linetype = "dashed", linewidth = 0.5) + # Add dashed line at x = 0
  theme_bw() +
  labs(x = "D-stat", y = "Pop1") +
  facet_wrap(~Pop3, nrow = 1) +
  scale_color_identity(guide = "legend", name = "Pop2", breaks = colour_map, labels = names(colour_map)) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
qpDstat_ursi_plot_ASIA
ggsave("qpDstat_ursi_plot_ASIA.png", qpDstat_ursi_plot_ASIA, height=8, width=5)


## POP1: EUR
# subset data
pop3_EUR <- "EUR"
data_EUR <- data[data$Pop3 == pop3_EUR, ]

qpDstat_ursi_plot_EUR <- ggplot(data_EUR, aes(x = D_stat, y = Pop1, color = Color)) +
  geom_point(position = position_dodge(width = 0.3), size = 4) + # Add points
  geom_errorbar(aes(xmin = D_stat - SE, xmax = D_stat + SE, y = Pop1), 
                width = 0.0,
                position = position_dodge(width = 0.3)) + # Add error bars
  geom_vline(xintercept = 0, color = "grey", linetype = "dashed", linewidth = 0.5) + # Add dashed line at x = 0
  theme_bw() +
  labs(x = "D-stat", y = "Pop1") +
  facet_wrap(~Pop3, nrow = 1) +
  scale_color_identity(guide = "legend", name = "Pop2", breaks = colour_map, labels = names(colour_map)) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
qpDstat_ursi_plot_EUR
ggsave("qpDstat_ursi_plot_EUR.png", qpDstat_ursi_plot_EUR, height=8, width=5)


## POP1: CENAM
# subset data
pop3_CENAM <- "CENAM"
data_CENAM <- data[data$Pop3 == pop3_CENAM, ]

qpDstat_ursi_plot_CENAM <- ggplot(data_CENAM, aes(x = D_stat, y = Pop1, color = Color)) +
  geom_point(position = position_dodge(width = 0.3), size = 4) + # Add points
  geom_errorbar(aes(xmin = D_stat - SE, xmax = D_stat + SE, y = Pop1), 
                width = 0.0,
                position = position_dodge(width = 0.3)) + # Add error bars
  geom_vline(xintercept = 0, color = "grey", linetype = "dashed", linewidth = 0.5) + # Add dashed line at x = 0
  theme_bw() +
  labs(x = "D-stat", y = "Pop1") +
  facet_wrap(~Pop3, nrow = 1) +
  scale_color_identity(guide = "legend", name = "Pop2", breaks = colour_map, labels = names(colour_map)) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
qpDstat_ursi_plot_CENAM
ggsave("qpDstat_ursi_plot_CENAM.png", qpDstat_ursi_plot_CENAM, height=8, width=5)


## POP1: USA
# subset data
pop3_USA <- "USA"
data_USA <- data[data$Pop3 == pop3_USA, ]

qpDstat_ursi_plot_USA <- ggplot(data_USA, aes(x = D_stat, y = Pop1, color = Color)) +
  geom_point(position = position_dodge(width = 0.3), size = 4) + # Add points
  geom_errorbar(aes(xmin = D_stat - SE, xmax = D_stat + SE, y = Pop1), 
                width = 0.0,
                position = position_dodge(width = 0.3)) + # Add error bars
  geom_vline(xintercept = 0, color = "grey", linetype = "dashed", linewidth = 0.5) + # Add dashed line at x = 0
  theme_bw() +
  labs(x = "D-stat", y = "Pop1") +
  facet_wrap(~Pop3, nrow = 1) +
  scale_color_identity(guide = "legend", name = "Pop2", breaks = colour_map, labels = names(colour_map)) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
qpDstat_ursi_plot_USA
ggsave("qpDstat_ursi_plot_USA.png", qpDstat_ursi_plot_USA, height=8, width=5)

```

![](images/qpDstat_ursi_plot.png)

Excess allele sharing for AUS:
- ASIA
- CENAM
- USA

Excess allele sharing for ASIA:
- AUS
- USA

Excess allele sharing for CENAM:
- AUS
- EUR

Excess allele sharing for EUR:
- CENAM

Excess allele sharing for USA:
- AUS
- ASIA
- CENAM

This shows that the D. immitis dispersal is not super super recent, because if that were the case, then there would be allele sharing between all populations because they would basically be identical. But it's also not super super ancient, because there is still some strong shared ancestry.

# Other tests (not included in manuscript)

The below code was for testing purposes only, and were not included in the final manuscript.

## outgroup f3 stats

From this tutorial: https://compvar-workshop.readthedocs.io/en/latest/contents/03_f3stats/f3stats.html#f4-statistics

F3 stats can measure admixture, as well as how closely two populations (A and B) are related to one another using a distant outgroup (C).

```bash
# make a new populations file

# URSI outgroup
# set the outgroup
OUTGROUP=URSI

# loop throguh the populations to generate the pop file as input to admixtools
for i in URSI AUS ASIA EUR CENAM USA; do
     for j in  URSI AUS ASIA EUR CENAM USA; do
          if [[ "$i" == "$j" ]] || [[ "$i" == "$OUTGROUP" ]] || [[ "$j" == "$OUTGROUP" ]]; then
               :
               else
               echo -e "${i}\t${j}\t${OUTGROUP}" >> admixtools_pops_f3.txt;
          fi;
     done;
done

# remove any duplicates in pops A and B


# REPENS outgroup
# set the outgroup
OUTGROUP=REPENS

# loop throguh the populations to generate the pop file as input to admixtools
for i in REPENS AUS ASIA EUR CENAM USA; do
     for j in  REPENS AUS ASIA EUR CENAM USA; do
          if [[ "$i" == "$j" ]] || [[ "$i" == "$OUTGROUP" ]] || [[ "$j" == "$OUTGROUP" ]]; then
               :
               else
               echo -e "${i}\t${j}\t${OUTGROUP}" >> admixtools_pops_f3.txt;
          fi;
     done;
done



# THAI outgroup
# set the outgroup
OUTGROUP=DTHAI

# loop throguh the populations to generate the pop file as input to admixtools
for i in DTHAI AUS ASIA EUR CENAM USA; do
     for j in  DTHAI AUS ASIA EUR CENAM USA; do
          if [[ "$i" == "$j" ]] || [[ "$i" == "$OUTGROUP" ]] || [[ "$j" == "$OUTGROUP" ]]; then
               :
               else
               echo -e "${i}\t${j}\t${OUTGROUP}" >> admixtools_pops_f3.txt;
          fi;
     done;
done




# run admixtools to generate f3 stats
qp3Pop -p PARAMETER_FILE_f3 > qp3Pop.out
```

PARAMETER_FILE:

```bash
genotypename:   nuclear_samples3x_missing0.9.chr1to4.recode.eigenstratgeno (in eigenstrat format)
snpname:        nuclear_samples3x_missing0.9.chr1to4.recode.snp      (in eigenstrat format)
indivname:      nuclear_samples3x_missing0.9.chr1to4.recode.ind    (in eigenstrat format)
popfilename:    admixtools_pops_f3.txt
inbreed: YES
printsd: YES
```

```bash
grep "result" qp3Pop.out | awk '{print $2,$3,$4,$5,$6,$7,$8}' OFS="\t" > qp3Pop.clean.out
# didnt work using these other species as outgroups.
```






## qpGraph

```bash
# install qpBrute
conda env create --name qpbrute --file https://raw.githubusercontent.com/ekirving/qpbrute/master/environment.yaml
# activate env
conda activate qpbrute

# parameter file
DIR:   ./
S1:              sim1           
indivname:       ../nuclear_samples3x_missing0.9.chr1to4.recode.ind  
snpname:         ../nuclear_samples3x_missing0.9.chr1to4.recode.snp       
genotypename:    ../nuclear_samples3x_missing0.9.chr1to4.recode.eigenstratgeno  
outpop:         NULL     
blgsize:        0.05
lsqmode:       YES
diag:          .0001
hires:         YES
initmix:      1000 
precision:    .0001  
zthresh:      3.0
terse:        NO
useallsnps:   YES

# graph file - start very simple first. we know CENAM is from EUR, AUS is from ASIA.
# qpgraph1:
root     R
label    AUS AUS 
label    ASIA   ASIA     
label    EUR   EUR    
label    CENAM   CENAM   
label    USA   USA
edge root_usa R USA
edge root_a1 R A1
edge A1_eur A1 EUR
edge A1_asia A1 ASIA
edge eur_cenam EUR CENAM
edge asia_aus ASIA AUS

qpGraph -p PARAMETER_FILE_qpgraph -g qpgraph1 -o qpgraph1.ggg -d qpgraph1.dot > qpgraph1.out

# qpgraph2:
root     R
label    AUS AUS 
label    ASIA   ASIA     
label    EUR   EUR    
label    CENAM   CENAM   
label    USA   USA
edge root_usa R USA
edge root_a1 R A1
edge A1_eur A1 EUR
edge A1_asia A1 ASIA
edge eur_cenam EUR CENAM
edge asia_aus ASIA AUS
admix mix_usa_asia USA ASIA
admix mix_cenam_aus CENAM AUS
admix mix_usa_cenam USA CENAM
admix mix_usa_eur USA EUR

qpGraph -p PARAMETER_FILE_qpgraph -g qpgraph2 -o qpgraph2.ggg -d qpgraph2.dot > qpgraph2.out

# go simpler

#qpgraph3:
root     R
label    ASIA   ASIA     
label    EUR   EUR    
edge root_A1 R A1
edge A1_asia A1 ASIA
edge A1_usa A1 USA
edge root_eur R EUR

qpGraph -p PARAMETER_FILE_qpgraph -g qpgraph3 -o qpgraph3.ggg -d qpgraph3.dot > qpgraph3.out

#qpgraph4: Now add EUR -> CENAM
root     R
label    ASIA   ASIA     
label    EUR   EUR    
edge root_A1 R A1
edge A1_asia A1 ASIA
edge A1_usa A1 USA
edge root_eur R EUR
edge eur_cenam EUR CENAM

qpGraph -p PARAMETER_FILE_qpgraph -g qpgraph4 -o qpgraph4.ggg -d qpgraph4.dot > qpgraph4.out

#qpgraph5: Now add ASIA -> AUS
root     R
label    ASIA   ASIA     
label    EUR   EUR    
edge root_A1 R A1
edge A1_asia A1 ASIA
edge A1_usa A1 USA
edge root_eur R EUR
edge eur_cenam EUR CENAM
edge asia_aus ASIA AUS

qpGraph -p PARAMETER_FILE_qpgraph -g qpgraph5 -o qpgraph5.ggg -d qpgraph5.dot > qpgraph5.out

# no outliers!
```

## qpBrute

Now try qpbrute - it searches for the best model automatically, so it will be a better option. Use URSI as outgroup.

```bash
# exhaustive search - will search all possible models
bsub.py --threads 10 10 qpBrute "qpBrute --threads 10 --par qpbrute_run1.par --prefix qpbrute_run1 --pops AUS ASIA EUR CENAM USA --out URSI"
# THIS WORKED
# FINISHED: Found 49 unique solution(s) from a total of 5,144 unique graphs!


# heuristic search - focuses on the most promising models
bsub.py --threads 10 10 qpBrute_heur "qpBrute --threads 10 --par qpbrute_run2.par --prefix qpbrute_run2 --pops AUS ASIA EUR CENAM USA --out URSI --heuristic"
# ERROR: Cannot resolve the graph from any permutation of the given nodes.
#FINISHED: Found 0 unique solution(s) from a total of 107 unique graphs!
# Ok so heuristic search doesn't work for this dataset...
```

# qpBayes

To calculate Bayes factors

```bash
module load cellgen/R/4.3.1
R

install.packages("optimbase_1.0-10.tar.gz", repos = NULL, type = "source") # worked
install.packages("optimsimplex_1.0-8.tar.gz", repos = NULL, type = "source") # worked
install.packages("neldermead_1.0-12.tar.gz", repos = NULL, type = "source") # worked
install.packages('foreach') # worked
install.packages('pracma') # worked
install.packages("admixturegraph_1.0.2.tar.gz", repos = NULL, type = "source") # worked?


bsub.py --threads 10 20 qpBayes "qpBayes --threads 10 \
    --geno ../nuclear_samples3x_missing0.9.chr1to4.recode.eigenstratgeno \
    --ind ../nuclear_samples3x_missing0.9.chr1to4.recode.ind   \
    --snp ../nuclear_samples3x_missing0.9.chr1to4.recode.snp  \
    --prefix qpbrute_run1 \
    --pops AUS ASIA EUR CENAM USA \
    --out URSI"

    ea6ffc108913
```
   



## Qpwave

```bash
# parameter file
genotypename:   nuclear_samples3x_missing0.9.chr1to4.recode.eigenstratgeno
snpname:       nuclear_samples3x_missing0.9.chr1to4.recode.snp
indivname:     nuclear_samples3x_missing0.9.chr1to4.recode.ind
popleft:       popleft_AUS.txt
popright:      popright_AUS.txt
details:       YES 
allsnps: YES
# did this for each pop

qpWave -p PARAMETER_FILE_qpwave > qpwave
```



## F4 ratio test

F4 ratio estimation allows inference of the mixing proportions of an admixture event, even without access to accurate surrogates for the ancestral populations.

For any 5 populations that are related to each other, we can compute the admixture proportion \alpha as - 

\alpha = f4(A,O; X,C)/ f4(A,O; B, C) 

```bash
qpF4ratio -p PARAMETER_FILE_F4ratio_ >logfile

A:
B:
X:
C:
O:

# loop through the populations to generate the pop file as input
OUTGROUP=URSI
for i in AUS ASIA EUR CENAM USA; do
     for j in  AUS ASIA EUR CENAM USA; do
          for k in AUS ASIA EUR CENAM USA; do
               if [[ "$i" == "$j" ]] || [[ "$i" == "$OUTGROUP" ]] || [[ "$j" == "$OUTGROUP" ]] || [[ "$i" == "$k" ]] || [[ "$j" == "$k" ]] || [[ "$k" == "$OUTGROUP" ]]; then
                    :
                    else
                    echo -e "${k}\t${j}\t${i}\t${OUTGROUP}" >> admixtools_pops_dstat_ursi.txt;
               fi;
          done;
     done;
done

# run admixtools to generate D-stats
qpDstat -p PARAMETER_FILE_dstat_ursi > qpDstat_ursi.out

# PARAMETER_FILE_dstat:
genotypename:   nuclear_samples3x_missing0.9.chr1to4.recode.eigenstratgeno (in eigenstrat format)
snpname:        nuclear_samples3x_missing0.9.chr1to4.recode.snp      (in eigenstrat format)
indivname:      nuclear_samples3x_missing0.9.chr1to4.recode.ind    (in eigenstrat format)
popfilename:    admixtools_pops_dstat_ursi.txt
inbreed: YES

# extract the relevant output lines
grep "result" qpDstat_ursi.out | awk '{print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11}' OFS="\t" > qpDstat_ursi.clean.out
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
printsd: YES
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




## f4 ratio

```bash
# run admixtools to generate f4 ratios
qpF4ratio -p PARAMETER_FILE_dstat >qpF4ratio.out
```


