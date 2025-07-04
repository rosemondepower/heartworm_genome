# Re-run BEAST2 analysis for mtDNA

Plan:
- make an alignment with the published Diro/Oncho mitogenomes
- clean up a bit
- get in NEXUS format
- Run BEAST2





## BEAST2

- Sequences: Published mitogenomes of D. immitis & Oncho
- Mutation rate: 1.6E-7 (upper 1.91E-7, lower 1.29E-7) (Denver 2000) -----> SLOW IT DOWN TO GENERATION = 4 YEARS
- C. elegans generation time is let's say ~ 1 week so D. immitis generation time of 4 years (208 weeks) is 208 times slower. So mutation rate would be 7.692307692E-10 (upper 9.182692308E-10, lower 6.201923077E-10)
- Model: Yule
- So slightly longer chain length and less frequent logging - Chain length 100000000, log every 500000






# SCRATCH

Plan:
- Extract whole mtDNA from the mapped genomic data, as using Get organelle may have introduced errors
- Make an alignment
- Run BEAST2

## Extract whole mitogenomes

```bash
# Load modules
module load PaM/environment
module load bsub.py/0.42.1
module load vcftools/0.1.16-c4
module load bcftools/1.17--h3cc50cf_1
module load gatk/4.1.4.1
module load samtools/1.14--hb421002_0

cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/05_ANALYSIS/MITO/GATK_Consensus

ln -s /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/04_VARIANTS/FILTER1/NO_OUTGROUPS/FINAL_SETS/mito_samples3x_missing0.9.recode.vcf

cp /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/04_VARIANTS/FILTER1/NO_OUTGROUPS/mito_samplelist.keep .
cd Sample_vcf

# Get a single sample per VCF
while read SAMPLE; do
vcftools --vcf ../mito_samples3x_missing0.9.recode.vcf \
--indv ${SAMPLE} \
--recode --out ${SAMPLE}_mito_samples3x_missing0.9
done < ../mito_samplelist.keep
cat *.log > combined.log

# gzip
while read SAMPLE; do
bgzip -f ${SAMPLE}_mito_samples3x_missing0.9.recode.vcf
done < ../mito_samplelist.keep

# index
while read SAMPLE; do
tabix -p vcf ${SAMPLE}_mito_samples3x_missing0.9.recode.vcf.gz
done < ../mito_samplelist.keep

# Filter out any rows that are 0/0 (only keep genotypes with the alternate allele)
while read SAMPLE; do
bcftools view -i 'GT!~"^0/0"' ${SAMPLE}_mito_samples3x_missing0.9.recode.vcf.gz -Oz -o ${SAMPLE}_mito_samples3x_missing0.9.filtered.vcf.gz
done < ../mito_samplelist.keep

# index
while read SAMPLE; do
tabix -p vcf ${SAMPLE}_mito_samples3x_missing0.9.filtered.vcf.gz
done < ../mito_samplelist.keep



# prep mtDNA reference
REFERENCE=/lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/01_REF/dimmitis_WSI_2.2_chrMtDNA.fa
samtools dict ${REFERENCE} > /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/01_REF/dimmitis_WSI_2.2_chrMtDNA.dict

cd ../Sample_fasta
while read SAMPLE; do
bsub.py 2 ${SAMPLE}_fasta \
"gatk FastaAlternateReferenceMaker \
-R ${REFERENCE} \
-O ${SAMPLE}_mito.fasta \
-L dirofilaria_immitis_chrMtDNA \
-V ../Sample_vcf/${SAMPLE}_mito_samples3x_missing0.9.filtered.vcf.gz"
done < ../mito_samplelist.keep

# incorporate sample ID into header
while read SAMPLE; do
sed "s/^>1 dirofilaria_immitis_chrMtDNA:1-13814/>${SAMPLE}/" ${SAMPLE}_mito.fasta > ${SAMPLE}_mito_ID.fasta
done < ../mito_samplelist.keep

Moved these fasta sequences into R_analysis>batch4>mito>GATK_Dec24>consensus_fasta

```

Made an alignment with other published Oncho mitogenomes so it's ready for BEAST.