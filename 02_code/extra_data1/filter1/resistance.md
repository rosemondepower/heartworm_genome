


bsub.py 1 resistance "../resistance.sh"
```bash
# Load modules
module load vcftools/0.1.16-c4
module load common-apps/htslib/1.9.229

cd /lustre/scratch125/pam/teams/team333/rp24/Diro/data/VARIANTS/RESISTANCE

vcftools --gzvcf ../ALLSITES/nuclearVARIANTsandINVARIANTs_allsamples.recode.RENAMED.COPY.vcf.gz --positions resistance_snps.txt --recode --out nuclearVARIANTsandINVARIANTs_allsamples.recode.RENAMED.RESISTANCE

vcftools --gzvcf ../ALLSITES/nuclearVARIANTsandINVARIANTs_allsamples.recode.RENAMED.COPY.vcf.gz --positions resistance_snps.txt --hardy --out nuclearVARIANTsandINVARIANTs_allsamples.recode.RENAMED.RESISTANCE

vcftools --gzvcf ../ALLSITES/nuclearVARIANTsandINVARIANTs_allsamples.recode.RENAMED.COPY.vcf.gz --positions resistance_snps.txt --freq --out nuclearVARIANTsandINVARIANTs_allsamples.recode.RENAMED.RESISTANCE

vcftools --gzvcf ../ALLSITES/nuclearVARIANTsandINVARIANTs_allsamples.recode.RENAMED.COPY.vcf.gz --positions resistance_snps.txt --site-pi --out nuclearVARIANTsandINVARIANTs_allsamples.recode.RENAMED.RESISTANCE
```
