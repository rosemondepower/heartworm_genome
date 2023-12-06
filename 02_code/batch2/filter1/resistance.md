


bsub.py 1 resistance "../resistance.sh"
```bash
# Load modules
module load vcftools/0.1.16-c4
module load common-apps/htslib/1.9.229

cd /lustre/scratch125/pam/teams/team333/rp24/Diro/data/VARIANTS/RESISTANCE

vcftools --gzvcf ../ALLSITES/nuclearVARIANTsandINVARIANTs_allsamples.recode.RENAMED.COPY.vcf.gz --positions resistance_snps.txt --recode --out nuclearVARIANTsandINVARIANTs_allsamples.recode.RENAMED.RESISTANCE

# Summary stats
## reports p-value for each site from a hardy-weinberg equilibrium test
vcftools --gzvcf ../ALLSITES/nuclearVARIANTsandINVARIANTs_allsamples.recode.RENAMED.COPY.vcf.gz --positions resistance_snps.txt --hardy --out nuclearVARIANTsandINVARIANTs_allsamples.recode.RENAMED.RESISTANCE
## outputs allele frequency for each site
vcftools --gzvcf ../ALLSITES/nuclearVARIANTsandINVARIANTs_allsamples.recode.RENAMED.COPY.vcf.gz --positions resistance_snps.txt --freq --out nuclearVARIANTsandINVARIANTs_allsamples.recode.RENAMED.RESISTANCE
## measures nucleotide divergency on a per-site basis
vcftools --gzvcf ../ALLSITES/nuclearVARIANTsandINVARIANTs_allsamples.recode.RENAMED.COPY.vcf.gz --positions resistance_snps.txt --site-pi --out nuclearVARIANTsandINVARIANTs_allsamples.recode.RENAMED.RESISTANCE

```

bsub.py 1 resistance2 "../resistance2.sh"
```bash
# Load modules
module load vcftools/0.1.16-c4
module load common-apps/htslib/1.9.229

cd /lustre/scratch125/pam/teams/team333/rp24/Diro/data/VARIANTS/RESISTANCE

# genotype matrix
## .012 file contains genotypes of each individual on a separate line. Gnotypes are represented as 0, 1, 2 (number represents that number of non-reference alleles. Missing genotypes are -1)
## .012.indv file contains the individuals included in the main file
## .012.pos file contains site locations included in the main file
vcftools --gzvcf ../ALLSITES/nuclearVARIANTsandINVARIANTs_allsamples.recode.RENAMED.COPY.vcf.gz --positions resistance_snps.txt --012 --out nuclearVARIANTsandINVARIANTs_allsamples.recode.RENAMED.RESISTANCE
```


```bash
# Transpose positions file
 awk '{for (i=1; i<=NF; i++) a[i,NR]=$i} END {for (i=1; i<=NF; i++) {for (j=1; j<=NR; j++) {printf a[i,j]; if (j < NR) printf "\t"; else printf "\n"}}}' /lustre/scratch125/pam/teams/team333/rp24/Diro/data/VARIANTS/RESISTANCE/nuclearVARIANTsandINVARIANTs_allsamples.recode.RENAMED.RESISTANCE.012.pos >  /lustre/scratch125/pam/teams/team333/rp24/Diro/data/VARIANTS/RESISTANCE/transposed.pos

# Cut ID column out
cut -f 2- < nuclearVARIANTsandINVARIANTs_allsamples.recode.RENAMED.RESISTANCE.012 > 0.12.cut

# Combine
cat transposed.pos 0.12.cut > combined.txt

# insert column at start
{ printf "\n\n"; cat nuclearVARIANTsandINVARIANTs_allsamples.recode.RENAMED.RESISTANCE.012.indv; } | paste - combined.txt > genotype.txt

# Make csv file
sed 's/\t/,/g' genotype.txt > genotype.csv
```

Discard everything above, I will try again using bcftools below to get the genotypes of the 42 SNPs.

bsub.py 1 resistance3 "../resistance3.sh"

```bash
module load bcftools/1.14--h88f3f91_0

cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/02_VARIANTS/RESISTANCE
bcftools query -R resistance_snps.txt -H -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' ../ALLSITES/nuclearVARIANTsandINVARIANTs_allsamples.recode.RENAMED.COPY.vcf.gz > bcftools_output.tsv

awk '{gsub("0/0", $3"/"$3); gsub("0/1", $3"/"$4); gsub("1/1", $4"/"$4); gsub("0|0", $3"/"$3); gsub("0|1", $3"/"$4); gsub("1|1", $4"/"$4)}1' bcftools_output.tsv > bcftools_alleles.tsv
# replace 0/0 with ref/ref
# replace 0/1 with ref/alt
# replace 1/1 with alt/alt
# sometimes it can give | e.g. 0|0 so replace those too
```