

bsub.py 10 hetero "../hetero.sh"

```bash
#vcftools
module load vcftools/0.1.16-c4

# working dir
WORKING_DIR=/lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA
cd ${WORKING_DIR}/02_VARIANTS/HETERO

# set environment variables
# VCF file
vcf_file=${WORKING_DIR}/02_VARIANTS/FINAL_SETS/nuclear_samples3x_missing0.9.chr1to4.recode.RENAMED.vcf
# Sample location file (format: sample_name location)
sample_location=sample_locations.txt
# Output file
output=observed_heterozygosity_by_location.txt

# Calculate observed heterozygosity by location
vcftools --vcf ${vcf_file} --het --out observed_heterozygosity

# Join the observed heterozygosity file with sample locations
join -1 1 -2 1 <(sort -k1,1 observed_heterozygosity.het) <(sort -k1,1 sample_locations.txt) > ${output}
```
