# SMC++


```bash
module load smcpp/1.15.3-c1
module load common-apps/htslib/1.9.229
module load bcftools/1.14--h88f3f91_0

# prep VCF file
bgzip /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/OTHER/BATCH3/02_VARIANTS/FINAL_SETS/nuclear_samples3x_missing0.9.chr1to4.recode.RENAMED.vcf
tabix -p vcf /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/OTHER/BATCH3/02_VARIANTS/FINAL_SETS/nuclear_samples3x_missing0.9.chr1to4.recode.RENAMED.vcf.gz
VCF=/lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/OTHER/BATCH3/02_VARIANTS/FINAL_SETS/nuclear_samples3x_missing0.9.chr1to4.recode.RENAMED.vcf.gz

cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/OTHER/BATCH3/02_VARIANTS/SMCPP

# Get sample names for each population
ASIA=$(awk -F'\t' '$2 == "Thailand" {print $1}' location.txt | paste -sd ',')
AUS=$(awk -F'\t' '$2 == "Australia" {print $1}' location.txt | paste -sd ',')
CAM=$(awk -F'\t' '$2 == "Panama" || $2 == "Costa Rica" {print $1}' location.txt | paste -sd ',')
EUR=$(awk -F'\t' '$2 == "Italy" || $2 == "Romania" {print $1}' location.txt | paste -sd ',')
USA=$(awk -F'\t' '$2 == "USA" {print $1}' location.txt | paste -sd ',')

## There may be missing genotypes designated as '.', whereas smc++ expects it to be './.' - so see if there are any lines like this:
zcat ${VCF} | \
  bcftools view | \
  grep -v '^#' | \
  awk '{
    found = 0;
    for (i=10; i<=90; i++) {
      split($i, gt, ":");
      if (gt[1] == "." && length(gt[1]) == 1 && found == 0) {
        print "Found \".\" in line", NR, "in column", i;
        found = 1;
      }
    }
  }'
# There are.

# view the actual lines
zcat ${VCF} | \
  bcftools view | \
  grep -v '^#' | \
  awk '{
    for (i=10; i<=90; i++) {
      split($i, gt, ":");
      if (gt[1] == "." && length(gt[1]) == 1) {
        print "Found \".\" in line", NR ":", $0;
        break;
      }
    }
  }'

## Replace them as './.'
zcat ${VCF} | \
  awk '
  BEGIN {OFS="\t"}
  /^#/ {print; next}
  {
    for (i=10; i<=NF; i++) {
      split($i, gt, ":");
      if (gt[1] == "." && length(gt[1]) == 1) {
        gt[1] = "./."
      }
      $i = gt[1]
      for (j=2; j<=length(gt); j++) {
        $i = $i ":" gt[j]
      }
    }
    print
  }' | bgzip > modified.vcf.gz

# check if it worked
zcat modified.vcf.gz | \
  bcftools view | \
  grep -v '^#' | \
  awk '{
    found = 0;
    for (i=10; i<=90; i++) {
      split($i, gt, ":");
      if (gt[1] == "." && length(gt[1]) == 1 && found == 0) {
        print "Found \".\" in line", NR, "in column", i;
        found = 1;
      }
    }
  }'

zcat modified.vcf.gz | \
  bcftools view | \
  grep -v '^#' | \
  awk '{
    for (i=10; i<=90; i++) {
      split($i, gt, ":");
      if (gt[1] == "." && length(gt[1]) == 1) {
        print "Found \".\" in line", NR ":", $0;
        break;
      }
    }
  }'
# Yep, none anymore! Continue with smc++ analysis.
tabix -p vcf modified.vcf.gz
VCF_MOD=modified.vcf.gz

# Convert VCF to the SMC++ input format with vcf2smc
smc++ vcf2smc ${VCF_MOD} chr1_ASIA.smc.gz dirofilaria_immitis_chr1 ASIA:${ASIA}
smc++ vcf2smc ${VCF_MOD} chr1_AUS.smc.gz dirofilaria_immitis_chr1 AUS:${AUS}
smc++ vcf2smc ${VCF_MOD} chr1_CAM.smc.gz dirofilaria_immitis_chr1 CAM:${CAM}
smc++ vcf2smc ${VCF_MOD} chr1_EUR.smc.gz dirofilaria_immitis_chr1 EUR:${EUR}
smc++ vcf2smc ${VCF_MOD} chr1_USA.smc.gz dirofilaria_immitis_chr1 USA:${USA}

smc++ vcf2smc ${VCF} chr2_ASIA.smc.gz chr2 ASIA:${ASIA}
smc++ vcf2smc ${VCF} chr2_AUS.smc.gz chr2 AUS:${AUS}
smc++ vcf2smc ${VCF} chr2_CAM.smc.gz chr2 CAM:${CAM}
smc++ vcf2smc ${VCF} chr2_EUR.smc.gz chr2 EUR:${EUR}
smc++ vcf2smc ${VCF} chr2_USA.smc.gz chr2 USA:${USA}

smc++ vcf2smc ${VCF} chr3_ASIA.smc.gz chr3 ASIA:${ASIA}
smc++ vcf2smc ${VCF} chr3_AUS.smc.gz chr3 AUS:${AUS}
smc++ vcf2smc ${VCF} chr3_CAM.smc.gz chr3 CAM:${CAM}
smc++ vcf2smc ${VCF} chr3_EUR.smc.gz chr3 EUR:${EUR}
smc++ vcf2smc ${VCF} chr3_USA.smc.gz chr3 USA:${USA}

smc++ vcf2smc ${VCF} chr4_ASIA.smc.gz chr4 ASIA:${ASIA}
smc++ vcf2smc ${VCF} chr4_AUS.smc.gz chr4 AUS:${AUS}
smc++ vcf2smc ${VCF} chr4_CAM.smc.gz chr4 CAM:${CAM}
smc++ vcf2smc ${VCF} chr4_EUR.smc.gz chr4 EUR:${EUR}
smc++ vcf2smc ${VCF} chr4_USA.smc.gz chr4 USA:${USA}
# each call to vcf2smc processes a single contig. VCFs containing multiple contigs should be processed via multiple separate runs.

# Fit the model using Estimate and plot
# chr 1
smc++ estimate -o /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/OTHER/BATCH3/02_VARIANTS/SMCPP 6.66e-9 chr1_ASIA.smc.gz
smc++ plot plot.png model.final.json




smc++ estimate -o /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/OTHER/BATCH3/02_VARIANTS/SMC 6.66e-9 chr1_AUS.smc.gz
smc++ estimate -o /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/OTHER/BATCH3/02_VARIANTS/SMC 6.66e-9 chr1_CAM.smc.gz 
smc++ estimate -o /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/OTHER/BATCH3/02_VARIANTS/SMC 6.66e-9 chr1_EUR.smc.gz 
smc++ estimate -o /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/OTHER/BATCH3/02_VARIANTS/SMC 6.66e-9 chr1_USA.smc.gz

# chr 2
smc++ estimate -o /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/OTHER/BATCH3/02_VARIANTS/SMC 6.66e-9 chr2_ASIA.smc.gz chr2_AUS.smc.gz chr2_CAM.smc.gz chr2_EUR.smc.gz chr2_USA.smc.gz
# chr 3
smc++ estimate -o /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/OTHER/BATCH3/02_VARIANTS/SMC 6.66e-9 chr3_ASIA.smc.gz chr3_AUS.smc.gz chr3_CAM.smc.gz chr3_EUR.smc.gz chr3_USA.smc.gz
# chr 4
smc++ estimate -o /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/OTHER/BATCH3/02_VARIANTS/SMC 6.66e-9 chr4_ASIA.smc.gz chr4_AUS.smc.gz chr4_CAM.smc.gz chr4_EUR.smc.gz chr1_USA.smc.gz

```