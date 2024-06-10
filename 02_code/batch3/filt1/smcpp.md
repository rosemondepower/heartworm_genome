# SMC++


```bash
module load smcpp/1.15.3-c1
module load common-apps/htslib/1.9.229
module load bcftools/1.14--h88f3f91_0

cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/OTHER/BATCH3/02_VARIANTS/SMCPP

# prep VCF file
bgzip /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/OTHER/BATCH3/02_VARIANTS/FINAL_SETS/nuclear_samples3x_missing0.9.chr1to4.recode.RENAMED.vcf
tabix -p vcf /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/OTHER/BATCH3/02_VARIANTS/FINAL_SETS/nuclear_samples3x_missing0.9.chr1to4.recode.RENAMED.vcf.gz
VCF=/lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/OTHER/BATCH3/02_VARIANTS/FINAL_SETS/nuclear_samples3x_missing0.9.chr1to4.recode.RENAMED.vcf.gz

# Get sample names for each population
ASIA=$(awk -F'\t' '$2 == "Thailand" {print $1}' location.txt | paste -sd ',')
AUS=$(awk -F'\t' '$2 == "Australia" {print $1}' location.txt | paste -sd ',')
CAM=$(awk -F'\t' '$2 == "Panama" || $2 == "Costa Rica" {print $1}' location.txt | paste -sd ',')
EUR=$(awk -F'\t' '$2 == "Italy" || $2 == "Romania" {print $1}' location.txt | paste -sd ',')
USA=$(awk -F'\t' '$2 == "USA" {print $1}' location.txt | paste -sd ',')
declare -A populations
populations=(
  [ASIA]=$ASIA
  [AUS]=$AUS
  [CAM]=$CAM
  [EUR]=$EUR
  [USA]=$USA
)

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
for chr in {1..4}; do
  smc++ vcf2smc ${VCF_MOD} chr${chr}_ASIA.smc.gz dirofilaria_immitis_chr${chr} ASIA:${ASIA}
  smc++ vcf2smc ${VCF_MOD} chr${chr}_AUS.smc.gz dirofilaria_immitis_chr${chr} AUS:${AUS}
  smc++ vcf2smc ${VCF_MOD} chr${chr}_CAM.smc.gz dirofilaria_immitis_chr${chr} CAM:${CAM}
  smc++ vcf2smc ${VCF_MOD} chr${chr}_EUR.smc.gz dirofilaria_immitis_chr${chr} EUR:${EUR}
  smc++ vcf2smc ${VCF_MOD} chr${chr}_USA.smc.gz dirofilaria_immitis_chr${chr} USA:${USA};
done
# each call to vcf2smc processes a single contig. VCFs containing multiple contigs should be processed via multiple separate runs.

# Fit the model using Estimate
## Mutation rate for C. elegans is 2.1x10^-8 mutations/site/generation (https://doi.org/10.1038/nature02697).
for chr in {1..4}; do
  smc++ estimate --timepoints 0 100000 -o /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/OTHER/BATCH3/02_VARIANTS/SMCPP/chr${chr}_ASIA_model 2.1e-8 chr${chr}_ASIA.smc.gz
  smc++ estimate --timepoints 0 100000 -o /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/OTHER/BATCH3/02_VARIANTS/SMCPP/chr${chr}_AUS_model 2.1e-8 chr${chr}_AUS.smc.gz
  smc++ estimate --timepoints 0 100000 -o /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/OTHER/BATCH3/02_VARIANTS/SMCPP/chr${chr}_CAM_model 2.1e-8 chr${chr}_CAM.smc.gz
  smc++ estimate --timepoints 0 100000 -o /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/OTHER/BATCH3/02_VARIANTS/SMCPP/chr${chr}_EUR_model 2.1e-8 chr${chr}_EUR.smc.gz
  smc++ estimate --timepoints 0 100000 -o /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/OTHER/BATCH3/02_VARIANTS/SMCPP/chr${chr}_USA_model 2.1e-8 chr${chr}_USA.smc.gz;
done

# plot models
## Use generation time of D. immitis ~ 1
# separately
for chr in {1..4}; do
  smc++ plot -g 1 chr${chr}_ASIA_model/chr${chr}_ASIA_plot.png chr${chr}_ASIA_model/model.final.json
  smc++ plot -g 1 chr${chr}_AUS_model/chr${chr}_AUS_plot.png chr${chr}_AUS_model/model.final.json
  smc++ plot -g 1 chr${chr}_CAM_model/chr${chr}_CAM_plot.png chr${chr}_CAM_model/model.final.json
  smc++ plot -g 1 chr${chr}_EUR_model/chr${chr}_EUR_plot.png chr${chr}_EUR_model/model.final.json
  smc++ plot -g 1 chr${chr}_USA_model/chr${chr}_USA_plot.png chr${chr}_USA_model/model.final.json;
done

# combined per chromosome
for chr in {1..4}; do
  smc++ plot -g 1 chr${chr}_plot.png chr${chr}_ASIA_model/model.final.json chr${chr}_AUS_model/model.final.json chr${chr}_CAM_model/model.final.json chr${chr}_EUR_model/model.final.json chr${chr}_USA_model/model.final.json;
done
```


# Clean split models

The split command fits two-population clean split models (assumes no ongoing gene flow between the two populations after they diverged)

```bash
# Create datasets containing the joint frequency spectrum for both populations. Run split.
## ASIA & AUS
for chr in {1..4}; do
mkdir split/ASIA_AUS/chr${chr}
smc++ vcf2smc ${VCF_MOD} split/ASIA_AUS/chr${chr}/chr${chr}_ASIA_AUS.smc.gz dirofilaria_immitis_chr${chr} ASIA:${ASIA} AUS:${AUS}
smc++ vcf2smc ${VCF_MOD} split/ASIA_AUS/chr${chr}/chr${chr}_AUS_ASIA.smc.gz dirofilaria_immitis_chr${chr} AUS:${AUS} ASIA:${ASIA}
smc++ split --timepoints 0 100000 -o split/ASIA_AUS/chr${chr} chr${chr}_ASIA_model/model.final.json chr${chr}_AUS_model.json/model.final.json split/ASIA_AUS/chr${chr}/*.smc.gz
smc++ plot split/ASIA_AUS/chr${chr}/chr${chr}_ASIA_AUS.pdf split/ASIA_AUS/chr${chr}/model.final.json;
done

## ASIA & CAM
for chr in {1..4}; do
mkdir split/ASIA_CAM/chr${chr}
smc++ vcf2smc ${VCF_MOD} split/ASIA_CAM/chr${chr}/chr${chr}_ASIA_CAM.smc.gz dirofilaria_immitis_chr${chr} ASIA:${ASIA} CAM:${CAM}
smc++ vcf2smc ${VCF_MOD} split/ASIA_CAM/chr${chr}/chr${chr}_CAM_ASIA.smc.gz dirofilaria_immitis_chr${chr} CAM:${CAM} ASIA:${ASIA}
smc++ split -o split/ASIA_CAM/chr${chr} chr${chr}_ASIA_model.json/model.final.json chr${chr}_CAM_model.json/model.final.json split/ASIA_CAM/chr${chr}/*.smc.gz
smc++ plot split/ASIA_CAM/chr${chr}/chr${chr}_ASIA_CAM.pdf split/ASIA_CAM/chr${chr}/model.final.json;
done

## ASIA & EUR

## ASIA & USA

## AUS & CAM

## AUS & EUR

## AUS & USA

## CAM & EUR

## CAM & USA

## EUR & USA





```