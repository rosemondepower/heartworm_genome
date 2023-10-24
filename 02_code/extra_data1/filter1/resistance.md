

bsub.py 1 resistance "../resistance.sh"

```bash
# load modules
module load bcftools/1.14--h88f3f91_0
module load common-apps/htslib/1.9.229

cd /lustre/scratch125/pam/teams/team333/rp24/Diro/data/VARIANTS/RESISTANCE

cp ../ALLSITES/nuclearVARIANTsandINVARIANTs_allsamples.recode.RENAMED.vcf ../ALLSITES/nuclearVARIANTsandINVARIANTs_allsamples.recode.RENAMED.COPY.vcf
bgzip ../ALLSITES/nuclearVARIANTsandINVARIANTs_allsamples.recode.RENAMED.COPY.vcf

bcftools query -R resistance.bed -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' ../ALLSITES/nuclearVARIANTsandINVARIANTs_allsamples.recode.RENAMED.COPY.vcf.gz > extracted_alleles.txt
```