# Re-run BEAST2 analysis for mtDNA

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

cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/05_ANALYSIS/MITO/Bcftools_Consensus

ln -s /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/04_VARIANTS/FILTER1/NO_OUTGROUPS/FINAL_SETS/mito_samples3x_missing0.9.recode.vcf

cd Sample_fasta

# Get a single sample per VCF
while read SAMPLE; do
vcftools --vcf ../mito_samples3x_missing0.9.recode.vcf \
--indv ${SAMPLE} \
--recode --out ${SAMPLE}_mito_samples3x_missing0.9
done < mito sample list


gatk FastaAlternateReferenceMaker \
-R reference.fasta \
-O output.fasta \
--intervals input.intervals \
-V input.vcf



``