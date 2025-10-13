# Filter 1 - Moderate stringency - OUTGROUPS

Filter the VCF file with moderate stringency, but this time include the outgroups.

## SNPs QC

```bash
# Load modules
module load vcftools/0.1.16-c4
module load bsub.py/0.42.1
module load common-apps/htslib/1.9.229

VCF=/lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/04_VARIANTS/FILTER1/DIMMITIS_POPGEN.cohort.2024-06-27.vcf.gz

cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/04_VARIANTS/FILTER1/OUTGROUPS
```

### Querying SNP and INDEL QC profiles to determine thresholds for filters

Adopted from Javier's paper.

```bash
bsub.py 10 run_snps_qc "run_snps_qc.sh"
```

```bash
# Load modules
module load gatk/4.1.4.1

WORKING_DIR=/lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/04_VARIANTS/FILTER1/OUTGROUPS

# set reference, vcf, and mitochondrial and Wb contig
REFERENCE=/lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/01_REF/dimmitis_WSI_2.2.fa
VCF=/lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/04_VARIANTS/FILTER1/DIMMITIS_POPGEN.cohort.2024-06-27.vcf.gz
MIT_CONTIG=dirofilaria_immitis_chrMtDNA
WB_CONTIG=dirofilaria_immitis_chrWb
cd ${WORKING_DIR}

# select nuclear SNPs
gatk SelectVariants \
--reference ${REFERENCE} \
--variant ${VCF} \
--select-type-to-include SNP \
--exclude-intervals ${MIT_CONTIG} \
--exclude-intervals ${WB_CONTIG} \
--output ${VCF%.vcf.gz}.nuclearSNPs.vcf

# select nuclear INDELs
gatk SelectVariants \
--reference ${REFERENCE} \
--variant ${VCF} \
--select-type-to-include INDEL \
--exclude-intervals ${MIT_CONTIG} \
--exclude-intervals ${WB_CONTIG} \
--output ${VCF%.vcf.gz}.nuclearINDELs.vcf

# select mitochondrial SNPs
gatk SelectVariants \
--reference ${REFERENCE} \
--variant ${VCF} \
--select-type-to-include SNP \
--intervals ${MIT_CONTIG} \
--output ${VCF%.vcf.gz}.mitoSNPs.vcf

# select mitochondrial INDELs
gatk SelectVariants \
--reference ${REFERENCE} \
--variant ${VCF} \
--select-type-to-include INDEL \
--intervals ${MIT_CONTIG} \
--output ${VCF%.vcf.gz}.mitoINDELs.vcf

# select WB SNPs
gatk SelectVariants \
--reference ${REFERENCE} \
--variant ${VCF} \
--select-type-to-include SNP \
--intervals ${WB_CONTIG} \
--output ${VCF%.vcf.gz}.WbSNPs.vcf

# select WB INDELs
gatk SelectVariants \
--reference ${REFERENCE} \
--variant ${VCF} \
--select-type-to-include INDEL \
--intervals ${WB_CONTIG} \
--output ${VCF%.vcf.gz}.WbINDELs.vcf

# won't make tables for the various SNP/indel data because I'm just going to stick with the same filtering parameters used in 'filter1'.

mv *SNPs* *INDELs* OUTGROUPS
```

Won't be making density plots this time, just use the ones obtained from filter1.



### Attending to the quantiles, thresholds for specific parameters are established

```bash
# load modules
module load bsub.py/0.42.1
module load gatk/4.1.4.1

WORKING_DIR=/lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/04_VARIANTS/FILTER1/OUTGROUPS
cd ${WORKING_DIR}

# set reference
REFERENCE=/lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/01_REF/dimmitis_WSI_2.2.fa

#Nuclear
bsub.py 1 filter_nuclearSNPs "gatk VariantFiltration \
--reference ${REFERENCE} \
--variant DIMMITIS_POPGEN.cohort.2024-06-27.nuclearSNPs.vcf \
--filter-expression 'QUAL < 40 || DP < 1837 || DP > 11972 || MQ < 21.52 || SOR > 7.112 || QD < 0.780 || FS > 8.275 || MQRankSum < -5.870 || ReadPosRankSum < -2.792 || ReadPosRankSum > 2.240' \
--filter-name "SNP_filtered" \
--output DIMMITIS_POPGEN.cohort.2024-06-27.nuclearSNPs.filtered.vcf"

bsub.py 1 filter_nuclearINDELS "gatk VariantFiltration \
--reference ${REFERENCE} \
--variant DIMMITIS_POPGEN.cohort.2024-06-27.nuclearINDELs.vcf \
--filter-expression 'QUAL < 46 || DP < 787 || DP > 12600 || MQ < 22.83 || SOR > 6.736 || QD < 0.990 || FS > 7.324 || MQRankSum < -4.797 || ReadPosRankSum < -4.178 || ReadPosRankSum > 1.960' \
--filter-name "INDEL_filtered" \
--output DIMMITIS_POPGEN.cohort.2024-06-27.nuclearINDELs.filtered.vcf"

#Mitochondrial
bsub.py 1 filter_mitoSNPs "gatk VariantFiltration \
--reference ${REFERENCE} \
--variant DIMMITIS_POPGEN.cohort.2024-06-27.mitoSNPs.vcf \
--filter-expression 'QUAL < 77 || DP < 199695 || DP > 659990 || MQ < 20.92 || SOR > 12.402 || QD < 3.320 || FS > 41.869 || MQRankSum < -5.622 || ReadPosRankSum < -3.003 || ReadPosRankSum > 5.452' \
--filter-name "SNP_filtered" \
--output DIMMITIS_POPGEN.cohort.2024-06-27.mitoSNPs.filtered.vcf"

bsub.py 1 filter_mitoINDELS "gatk VariantFiltration \
--reference ${REFERENCE} \
--variant DIMMITIS_POPGEN.cohort.2024-06-27.mitoINDELs.vcf \
--filter-expression 'QUAL < 61 || DP < 200280 || DP > 655589 || MQ < 21.36 || SOR > 14.646 || QD < 0.288 || FS > 82.622 || MQRankSum < -7.928 || ReadPosRankSum < -7.416 || ReadPosRankSum > 6.170' \
--filter-name "INDEL_filtered" \
--output DIMMITIS_POPGEN.cohort.2024-06-27.mitoINDELs.filtered.vcf"

#Wolbachia
bsub.py 1 filter_WbSNPs "gatk VariantFiltration \
--reference ${REFERENCE} \
--variant DIMMITIS_POPGEN.cohort.2024-06-27.WbSNPs.vcf \
--filter-expression 'QUAL < 51 || DP < 93526 || DP > 140940 || MQ < 22.00 || SOR > 11.661 || QD < 1.790 || FS > 46.426 || MQRankSum < -15.506 || ReadPosRankSum < -6.165 || ReadPosRankSum > 13.108' \
--filter-name "SNP_filtered" \
--output DIMMITIS_POPGEN.cohort.2024-06-27.WbSNPs.filtered.vcf"

bsub.py 1 filter_WbINDELS "gatk VariantFiltration \
--reference ${REFERENCE} \
--variant DIMMITIS_POPGEN.cohort.2024-06-27.WbINDELs.vcf \
--filter-expression 'QUAL < 42 || DP < 93456 || DP > 145612 || MQ < 22.017 || SOR > 11.268 || QD < 1.210 || FS > 59.220 || MQRankSum < -6.716 || ReadPosRankSum < -6.413 || ReadPosRankSum > 6.014' \
--filter-name "INDEL_filtered" \
--output DIMMITIS_POPGEN.cohort.2024-06-27.WbINDELs.filtered.vcf"

# once done, count the filtered sites - funny use of "|" allows direct markdown table format
echo -e "| Filtered_VCF | Variants_PASS | Variants_FILTERED |\n| -- | -- | -- | " > filter.stats

for i in *filtered.vcf; do
     name=${i}; pass=$( grep -E 'PASS' ${i} | wc -l ); filter=$( grep -E 'filter' ${i} | wc -l );
     echo -e "| ${name} | ${pass} | ${filter} |" >> filter.stats
done

mv *.o *.e LOGS
```
This is the summary of the filtered variants ('filter.stats'):

| Filtered_VCF | Variants_PASS | Variants_FILTERED |
| -- | -- | -- | 
| DIMMITIS_POPGEN.cohort.2024-06-27.mitoINDELs.filtered.vcf | 240 | 30 |
| DIMMITIS_POPGEN.cohort.2024-06-27.mitoSNPs.filtered.vcf | 1511 | 141 |
| DIMMITIS_POPGEN.cohort.2024-06-27.nuclearINDELs.filtered.vcf | 748192 | 85599 |
| DIMMITIS_POPGEN.cohort.2024-06-27.OUTGROUPS.nuclearSNPs.filtered.vcf | 1775848 | 189001 |
| DIMMITIS_POPGEN.cohort.2024-06-27.WbINDELs.filtered.vcf | 8897 | 731 |
| DIMMITIS_POPGEN.cohort.2024-06-27.WbSNPs.filtered.vcf | 69952 | 5110 |


Filtering looks ok.

Usually we would merge the SNP and INDEL files together to make a joined VCF. However, I want to disregard the indels moving forward and only focus on the SNPs. Some downstream tools don't like having indels in there.


### Filter genotypes based on x3 depth per genotype

```bash
#Nuclear 
bsub.py 1 filter_nuclear_GT \
"gatk VariantFiltration \
--reference ${REFERENCE} \
--variant DIMMITIS_POPGEN.cohort.2024-06-27.nuclearSNPs.filtered.vcf \
--genotype-filter-expression ' DP < 3 '  \
--genotype-filter-name "DP_lt3" \
--output DIMMITIS_POPGEN.cohort.2024-06-27.nuclearSNPs.DPfiltered.vcf"

bsub.py --done "filter_nuclear_GT" 1 filter_nuclear_GT2 \
"gatk SelectVariants \
--reference ${REFERENCE} \
--variant DIMMITIS_POPGEN.cohort.2024-06-27.nuclearSNPs.DPfiltered.vcf \
--set-filtered-gt-to-nocall \
--output DIMMITIS_POPGEN.cohort.2024-06-27.nuclearSNPs.DPfilterNoCall.vcf"

#Mito
bsub.py 1 filter_mito_GT \
"gatk VariantFiltration \
--reference ${REFERENCE} \
--variant DIMMITIS_POPGEN.cohort.2024-06-27.mitoSNPs.filtered.vcf \
--genotype-filter-expression ' DP < 3 '  \
--genotype-filter-name "DP_lt3" \
--output DIMMITIS_POPGEN.cohort.2024-06-27.mitoSNPs.DPfiltered.vcf"

bsub.py --done "filter_mito_GT" 1 filter_mito_GT2 \
"gatk SelectVariants \
--reference ${REFERENCE} \
--variant DIMMITIS_POPGEN.cohort.2024-06-27.mitoSNPs.DPfiltered.vcf \
--set-filtered-gt-to-nocall \
--output DIMMITIS_POPGEN.cohort.2024-06-27.mitoSNPs.DPfilterNoCall.vcf"

#wolbachia
bsub.py 1 filter_Wb_GT \
"gatk VariantFiltration \
--reference ${REFERENCE} \
--variant DIMMITIS_POPGEN.cohort.2024-06-27.WbSNPs.filtered.vcf \
--genotype-filter-expression ' DP < 3 '  \
--genotype-filter-name "DP_lt3" \
--output DIMMITIS_POPGEN.cohort.2024-06-27.WbSNPs.DPfiltered.vcf"

bsub.py --done "filter_Wb_GT" 1 filter_Wb_GT2 \
"gatk SelectVariants \
--reference ${REFERENCE} \
--variant DIMMITIS_POPGEN.cohort.2024-06-27.WbSNPs.DPfiltered.vcf \
--set-filtered-gt-to-nocall \
--output DIMMITIS_POPGEN.cohort.2024-06-27.WbSNPs.DPfilterNoCall.vcf"
```


### Now we apply a set of standard filters for population genomics

bsub.py 10 run_standard_filt "run_standard_filt.sh"

```bash
# load modules
module load vcftools/0.1.16-c4

cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/04_VARIANTS/FILTER1/OUTGROUPS
VCF=/lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/04_VARIANTS/FILTER1/OUTGROUPS/DIMMITIS_POPGEN.cohort.2024-06-27.vcf.gz

#Nuclear variants
vcftools \
--vcf DIMMITIS_POPGEN.cohort.2024-06-27.nuclearSNPs.DPfilterNoCall.vcf \
--remove-filtered-geno-all \
--remove-filtered-all \
--min-alleles 2 \
--max-alleles 2 \
--hwe 1e-06 \
--maf 0.02 \
--recode \
--recode-INFO-all \
--out DIMMITIS_POPGEN.cohort.2024-06-27.nuclear_SNPs.final
After filtering, kept 143 out of 143 Individuals
Outputting VCF file...
After filtering, kept 268135 out of a possible 1964841 Sites

#--- nuclear SNPs
vcftools --vcf DIMMITIS_POPGEN.cohort.2024-06-27.nuclear_SNPs.final.recode.vcf --remove-indels

#--- nuclear  INDELs
vcftools --vcf DIMMITIS_POPGEN.cohort.2024-06-27.nuclear_SNPs.final.recode.vcf --keep-only-indels


#Mitochondrial variants
vcftools \
--vcf DIMMITIS_POPGEN.cohort.2024-06-27.mitoSNPs.DPfilterNoCall.vcf \
--remove-filtered-geno-all \
--remove-filtered-all \
--min-alleles 2 \
--max-alleles 2 \
--maf 0.02 \
--recode \
--recode-INFO-all \
--out DIMMITIS_POPGEN.cohort.2024-06-27.mito_SNPs.final
#After filtering, kept 143 out of 143 Individuals
Outputting VCF file...
After filtering, kept 500 out of a possible 1644 Sites

#--- mito SNPs
vcftools --vcf DIMMITIS_POPGEN.cohort.2024-06-27.mito_SNPs.final.recode.vcf --remove-indels

#--- mito INDELs
vcftools --vcf DIMMITIS_POPGEN.cohort.2024-06-27.mito_SNPs.final.recode.vcf --keep-only-indels


#wolbachia variants
vcftools \
--vcf DIMMITIS_POPGEN.cohort.2024-06-27.WbSNPs.DPfilterNoCall.vcf \
--remove-filtered-geno-all \
--remove-filtered-all \
--min-alleles 2 \
--max-alleles 2 \
--maf 0.02 \
--recode \
--recode-INFO-all \
--out DIMMITIS_POPGEN.cohort.2024-06-27.Wb_SNPs.final
#After filtering, kept 143 out of 143 Individuals
Outputting VCF file...
After filtering, kept 25563 out of a possible 75054 Sites

#--- Wb SNPs
vcftools --vcf DIMMITIS_POPGEN.cohort.2024-06-27.Wb_SNPs.final.recode.vcf --remove-indels

#--- Wb INDELs
vcftools --vcf DIMMITIS_POPGEN.cohort.2024-06-27.Wb_SNPs.final.recode.vcf --keep-only-indels
```

### Now, we are filtering by missingness

```bash
#determine missingness per individual
vcftools --vcf DIMMITIS_POPGEN.cohort.2024-06-27.nuclear_SNPs.final.recode.vcf --out nuclear --missing-indv
vcftools --vcf DIMMITIS_POPGEN.cohort.2024-06-27.mito_SNPs.final.recode.vcf --out mito --missing-indv
vcftools --vcf DIMMITIS_POPGEN.cohort.2024-06-27.Wb_SNPs.final.recode.vcf --out wb --missing-indv
```

### Check the missingess in R

```R
 # Check missingness
  data_nuclear <- read.delim("nuclear.imiss", header=T)
  data_mito <- read.delim("mito.imiss", header=T)
  data_wb <- read.delim("wb.imiss", header=T)
  
  #creating the function - per sample
  fun_plot_missingness <- function(data,title) {
    
    plot <- ggplot(data, aes(INDV, 1-F_MISS)) +
      geom_boxplot(color = 'brown') +
      geom_point(size = 1, color = 'brown4') +
      theme_bw() +
      labs(x="Sample ID", y="Proportion of total variants present (1-missingness)")+
      ggtitle(title) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    print(plot)
    ggsave(paste0("plot_missingness_figure_",title,".png"))
  }
  
  # plotting for each dataset
  fun_plot_missingness(data_nuclear, "nuclear_variants")

  fun_plot_missingness(data_mito,"mitochondrial_variants")

  fun_plot_missingness(data_wb, "wb_variants")
```

```bash
# For nuclear (n=133) - nuclear_samplelist.keep
## remove samples with missingness less than ~0.5 (AUS_BNE_AD_005, AUS_SYD_AD_016, AUS_TVS_AD_009, MYS_SEL_AD_001_R, ROU_GIU_AD_001_R, USA_GEO_MF_001, USA_GEO_MF_002, USA_ILL_MF_001, USA_LOU_MF_001, USA_MIP_MF_001)


# For mithochondiral (n=137) - mito_samplelist.keep
# removed: (ROU_GIU_AD_001_R, USA_GEO_MF_001, USA_GEO_MF_002, USA_ILL_MF_001, USA_LOU_MF_001, USA_MIP_MF_001)


# For wb (n=130) - wb_samplelist.keep
# remove: (AUS_SYD_AD_002, GRC_XAN_AD_011, MYS_SEL_AD_001, MYS_SEL_AD_001_R, PAN_BOC_AD_001, PAN_BOC_AD_003, ROU_COM_AD_001, ROU_GIU_AD_001_R, USA_GEO_MF_001, USA_GEO_MF_002, USA_ILL_MF_001, USA_LOU_MF_001, USA_MIP_MF_001)

```


### Let's check different thresholds for each dataset

bsub.py 2 run_check_thresh "run_check_thresh.sh"

```bash
# load modules
module load vcftools/0.1.16-c4

cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/04_VARIANTS/FILTER1/OUTGROUPS

# For nuclear variants
for i in 0.7 0.8 0.9 1; do
     vcftools --vcf DIMMITIS_POPGEN.cohort.2024-06-27.nuclear_SNPs.final.recode.vcf --keep nuclear_samplelist.keep --max-missing ${i} ;
done

# max-missing = 0.7
#After filtering, kept 133 out of 143 Individuals
#After filtering, kept 267405 out of a possible 268135 Sites


# max-missing = 0.8
#After filtering, kept 133 out of 143 Individuals
#After filtering, kept 266404 out of a possible 268135 Sites


# max-missing = 0.9
#After filtering, kept 133 out of 143 Individuals
#After filtering, kept 263414 out of a possible 268135 Sites


# max-missing = 1
#After filtering, kept 133 out of 143 Individuals
#After filtering, kept 4 out of a possible 268135 Sites


# For mito variants
for i in 0.7 0.8 0.9 1; do
     vcftools --vcf DIMMITIS_POPGEN.cohort.2024-06-27.mito_SNPs.final.recode.vcf --keep mito_samplelist.keep --max-missing ${i} ;
done

# max-missing = 0.7
#After filtering, kept 137 out of 143 Individuals
#After filtering, kept 500 out of a possible 500 Sites


# max-missing = 0.8
#After filtering, kept 137 out of 143 Individuals
#After filtering, kept 500 out of a possible 500 Sites


# max-missing = 0.9
#After filtering, kept 137 out of 143 Individuals
#After filtering, kept 500 out of a possible 500 Sites


# max-missing = 1
#After filtering, kept 137 out of 143 Individuals
#After filtering, kept 207 out of a possible 500 Sites


# For Wb variants
for i in 0.7 0.8 0.9 1; do
     vcftools --vcf DIMMITIS_POPGEN.cohort.2024-06-27.Wb_SNPs.final.recode.vcf --keep wb_samplelist.keep --max-missing ${i} ;
done

# max-missing = 0.7
#After filtering, kept 130 out of 143 Individuals
#After filtering, kept 25563 out of a possible 25563 Sites

# max-missing = 0.8
#After filtering, kept 130 out of 143 Individuals
#After filtering, kept 25563 out of a possible 25563 Sites

# max-missing = 0.9
#After filtering, kept 130 out of 143 Individuals
#After filtering, kept 25561 out of a possible 25563 Sites

# max-missing = 1
#After filtering, kept 130 out of 143 Individuals
#After filtering, kept 1378 out of a possible 25563 Sites
```

Selecting a max missingness of 0.9 for nuclear, 0.9 for mito and 0.9 for Wb is sensible.

bsub.py 1 run_thresh "run_thresh.sh"

```bash
# load modules
module load vcftools/0.1.16-c4

cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/04_VARIANTS/FILTER1/OUTGROUPS

mkdir FINAL_SETS

# For nuclear
vcftools --vcf DIMMITIS_POPGEN.cohort.2024-06-27.nuclear_SNPs.final.recode.vcf \
     --keep nuclear_samplelist.keep \
     --max-missing 0.9 \
     --recode --recode-INFO-all \
     --out FINAL_SETS/nuclear_samples3x_missing0.9

# For mito
vcftools --vcf DIMMITIS_POPGEN.cohort.2024-06-27.mito_SNPs.final.recode.vcf \
     --keep mito_samplelist.keep \
     --max-missing 0.9 \
     --recode --recode-INFO-all \
     --out FINAL_SETS/mito_samples3x_missing0.9

# For wb
vcftools --vcf DIMMITIS_POPGEN.cohort.2024-06-27.Wb_SNPs.final.recode.vcf \
     --keep wb_samplelist.keep \
     --max-missing 0.9 \
     --recode --recode-INFO-all \
     --out FINAL_SETS/wb_samples3x_missing0.9
```


### Select only the variants in the chr 1 to chr4, avoiding the chrX and the scaffolds. Select chr1-4 and chrX separately.

bsub.py 4 run_select_chr "run_select_chr.sh"

```bash
# load modules
module load vcftools/0.1.16-c4

cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/04_VARIANTS/FILTER1/OUTGROUPS/FINAL_SETS

# chr1-4
vcftools --vcf nuclear_samples3x_missing0.9.recode.vcf \
--chr dirofilaria_immitis_chr1 \
--chr dirofilaria_immitis_chr2 \
--chr dirofilaria_immitis_chr3 \
--chr dirofilaria_immitis_chr4 \
--recode --out nuclear_samples3x_missing0.9.chr1to4
#After filtering, kept 133 out of 133 Individuals
#After filtering, kept 195900 out of a possible 263414 Sites


vcftools --vcf nuclear_samples3x_missing0.9.chr1to4.recode.vcf --remove-indels
#After filtering, kept 133 out of 133 Individuals
#After filtering, kept 195900 out of a possible 195900 Sites

vcftools --vcf nuclear_samples3x_missing0.9.chr1to4.recode.vcf --keep-only-indels
#After filtering, kept 133 out of 133 Individuals
#After filtering, kept 0 out of a possible 195900 Sites- this makes sense because I removed the indels earlier and only focused on the SNPs.
```

Now I can run downstream analyses with these outgroup samples, e.g. admixtools.


## Outgroups SNP stats

```bash
module load bcftools/1.17--h3cc50cf_1

cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/06_ANALYSIS/OUTGROUPS

vcftools --vcf /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/04_VARIANTS/FILTER1/OUTGROUPS/FINAL_SETS/nuclear_samples3x_missing0.9.chr1to4.recode.vcf --out outgroup_stats --depth --site-mean-depth --geno --het --counts

bcftools stats \
--samples 'ITA_NEA_AD_004,ITA_NEA_AD_005,USA_WIS_AD_001,VNM_HAI_AD_001,VNM_HAN_AD_001' \
/lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/04_VARIANTS/FILTER1/OUTGROUPS/FINAL_SETS/nuclear_samples3x_missing0.9.chr1to4.recode.vcf > outgroup.stats

grep "PSC" vcf_outgroups.stats > vcf_outgroups.clean.stats

```




