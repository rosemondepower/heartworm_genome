# Dirofilaria immitis WGS Lab Book - Pixy

Ideas for running pixy:
- Populations = region (AUS, ASIA, USA, EUR, CENAM) (code below - adopted from Javi's paper)
- Populations = host (dog, cat, fox, ferret etc) (will adapt below code for this if everything looks ok)
- Populations = city (all samples worldwide) or perhaps just focus on specific countries, e.g. run pixy for Australia only and USA only.                                            

## Generate an ALL SITES variant set for running pixy

Steve has already done some. Continue making an -all-sites vcf.


```bash
#!/bin/bash

#-------------------------------------------------------------------------------
# run_gatk_hc.sh
#-------------------------------------------------------------------------------

# stephen doyle
# Jan 2023

# Export environment variables
export PREFIX=dirofilaria_global  # prefix for output files
#export REFERENCE=/nfs/users/nfs_s/sd21/lustre_link/dirofilaria_immitis/project_global-diversity/01_references/dimmitis_WSI_2.2.fa  # path to reference genome
#export BAM_LIST=/nfs/users/nfs_s/sd21/lustre_link/dirofilaria_immitis/project_global-diversity/04_analysis/01_mapping/results_dirofilaria_only/bam.list  # path to list of BAM files

# Load modules
module load PaM/environment
module load bsub.py/0.42.1
module load gatk/4.1.4.1
module load htslib-1.19/perl-5.38.0
module load vcftools/0.1.16-c4

cd /lustre/scratch127/pam/teams/team333/rp24/dirofilaria_immitis/R2/variants_ALLSITES

# Define file locations
export LOG_FILES="${PWD}/gatk_hc_${PREFIX}/LOG_FILES"  # directory for log files
export REFERENCE_FILES="${PWD}/gatk_hc_${PREFIX}/REFERENCE_FILES"  # directory for reference files
export GATK_HC_GVCFs="${PWD}/gatk_hc_${PREFIX}/GATK_HC_GVCFs"  # directory for GATK HC GVCF files
export GATK_HC_MERGED="${PWD}/gatk_hc_${PREFIX}/GATK_HC_MERGED"  # directory for merged haplotype caller files

# Create directories if they don't exist
[ -d ${LOG_FILES} ] || mkdir -p ${LOG_FILES}
[ -d ${REFERENCE_FILES} ] || mkdir -p ${REFERENCE_FILES}
[ -d ${GATK_HC_GVCFs} ] || mkdir -p ${GATK_HC_GVCFs}
[ -d ${GATK_HC_MERGED} ] || mkdir -p ${GATK_HC_MERGED}



# Save current script in run folder to reproduce the exact output
cp ${PWD}/run_gatk_hc.sh ${PWD}/gatk_hc_${PREFIX}/commands.$(date -Iminutes).txt




#-------------------------------------------------------------------------------
### 03. Merge GVCFs
#-------------------------------------------------------------------------------
func_merge_gvcf() { 

#ls -1 ${GATK_HC_GVCFs}/*complete/*gz > ${GATK_HC_MERGED}/gvcf.list

[ -d ${GATK_HC_MERGED}/LOGFILES ] || mkdir -p ${GATK_HC_MERGED}/LOGFILES


n=1
for SEQUENCE in /lustre/scratch127/pam/teams/team333/sd21/dirofilaria_immitis/project_global-diversity/04_analysis/02_variants/gatk_hc_dirofilaria_global/REFERENCE_FILES/REFsplit*list; do
    SEQUENCE=$( echo ${SEQUENCE} | awk -F '/' '{print $NF}' )
    echo -e "gatk CombineGVCFs -R /lustre/scratch127/pam/teams/team333/sd21/dirofilaria_immitis/project_global-diversity/04_analysis/02_variants/gatk_hc_dirofilaria_global/REFERENCE_FILES/REF.fa --intervals /lustre/scratch127/pam/teams/team333/sd21/dirofilaria_immitis/project_global-diversity/04_analysis/02_variants/gatk_hc_dirofilaria_global/REFERENCE_FILES/${SEQUENCE} \\" > ${GATK_HC_MERGED}/run_merge_gvcfs.tmp.job_${n}
    while read SAMPLE; do
        echo -e "--variant ${SAMPLE} \\" >> ${GATK_HC_MERGED}/run_merge_gvcfs.tmp.job_${n};
   done < /lustre/scratch127/pam/teams/team333/sd21/dirofilaria_immitis/project_global-diversity/04_analysis/02_variants/gatk_hc_dirofilaria_global/GATK_HC_MERGED/gvcf.list
   echo -e "--output ${GATK_HC_MERGED}/${n}.cohort.tmp.g.vcf.gz" >> ${GATK_HC_MERGED}/run_merge_gvcfs.tmp.job_${n};
   let "n+=1"; 
done

chmod a+x ${GATK_HC_MERGED}/run_merge_gvcfs.tmp.job_*

# setup job conditions
JOBS=$( ls -1 ${GATK_HC_MERGED}/run_merge_gvcfs.tmp.job_* | wc -l )
ID="U$(date +%s)"

#submit job array to call variants put scaffold / contig
bsub -q long -R'span[hosts=1] select[mem>10000] rusage[mem=10000]' -n 10 -M10000 -J "gatk_merge_gvcf_[1-$JOBS]%100" -e "${GATK_HC_MERGED}/LOGFILES/gatk_merge_gvcf_[1-$JOBS].e" -o "${GATK_HC_MERGED}/LOGFILES/gatk_merge_gvcf_[1-$JOBS].o" "${GATK_HC_MERGED}/run_merge_gvcfs.tmp.job_\$LSB_JOBINDEX"

rm ${GATK_HC_MERGED}/MERGE_ARRAY_FINISHED
bsub -w "done(gatk_merge_gvcf_)" -q normal -R'span[hosts=1] select[mem>100] rusage[mem=100]' -n 1 -M100 -J "gatk_merge_gvcf_finish" -e "${GATK_HC_MERGED}/LOGFILES/gatk_merge_gvcf_finish.e" -o "${GATK_HC_MERGED}/LOGFILES/gatk_merge_gvcf_finish.o" "touch ${GATK_HC_MERGED}/MERGE_ARRAY_FINISHED"

until [ -f "${GATK_HC_MERGED}/MERGE_ARRAY_FINISHED" ]
do
     sleep 10
done

}

export -f func_merge_gvcf





#-------------------------------------------------------------------------------
### 04. Genotype GVCFs
#-------------------------------------------------------------------------------

func_genotype_gvcfs() { 

# split each chromosome up into separate jobs, and run genotyping on each individually.   
n=1
for SEQUENCE in /lustre/scratch127/pam/teams/team333/sd21/dirofilaria_immitis/project_global-diversity/04_analysis/02_variants/gatk_hc_dirofilaria_global/REFERENCE_FILES/REFsplit*list; do
    SEQUENCE=$( echo ${SEQUENCE} | awk -F '/' '{print $NF}' )
    echo -e "gatk GenotypeGVCFs \
    -all-sites \
    -R /lustre/scratch127/pam/teams/team333/sd21/dirofilaria_immitis/project_global-diversity/04_analysis/02_variants/gatk_hc_dirofilaria_global/REFERENCE_FILES/REF.fa \
    -V ${GATK_HC_MERGED}/${n}.cohort.tmp.g.vcf.gz \
    --intervals /lustre/scratch127/pam/teams/team333/sd21/dirofilaria_immitis/project_global-diversity/04_analysis/02_variants/gatk_hc_dirofilaria_global/REFERENCE_FILES/${SEQUENCE} \
    -O ${GATK_HC_MERGED}/${n}.cohort.tmp.vcf.gz -G StandardAnnotation -G AS_StandardAnnotation" > ${GATK_HC_MERGED}/run_hc_genotype.tmp.job_${n};
    let "n+=1"; 
done

chmod a+x ${GATK_HC_MERGED}/run_hc_genotype*

# setup job conditions
JOBS=$( ls -1 ${GATK_HC_MERGED}/run_hc_genotype* | wc -l )
ID="U$(date +%s)"

bsub -q long -R'span[hosts=1] select[mem>100000] rusage[mem=100000]' -n 6 -M100000 -J "gatk_genotype_cohort_gvcf_[1-$JOBS]" -e "${GATK_HC_MERGED}/LOGFILES/gatk_genotype_cohort_gvcf_[1-$JOBS].e" -o "${GATK_HC_MERGED}/LOGFILES/gatk_genotype_cohort_gvcf_[1-$JOBS].o" "${GATK_HC_MERGED}/run_hc_genotype.tmp.job_*\$LSB_JOBINDEX"

rm ${GATK_HC_MERGED}/GENOTYPE_ARRAY_FINISHED
bsub -w "done(gatk_genotype_cohort_gvcf_)" -q normal -R'span[hosts=1] select[mem>100] rusage[mem=100]' -n 1 -M100 -J "gatk_genotype_cohort_gvcf_finish" -e "${GATK_HC_MERGED}/LOGFILES/gatk_genotype_cohort_gvcf_finish.e" -o "${GATK_HC_MERGED}/LOGFILES/gatk_genotype_cohort_gvcf_finish.o" "touch ${GATK_HC_MERGED}/GENOTYPE_ARRAY_FINISHED"

until [ -f "${GATK_HC_MERGED}/GENOTYPE_ARRAY_FINISHED" ]
do
     sleep 10
done

}

export -f func_genotype_gvcfs



#-------------------------------------------------------------------------------
### 05. Finish making VCF and cleanup
#-------------------------------------------------------------------------------


func_finish_vcf() {

    # 
    ls ${GATK_HC_MERGED}/*.cohort.tmp.vcf.gz > ${GATK_HC_MERGED}/cohort.vcf.list

    # concatenate the vcf files in the list
    vcf-concat --files ${GATK_HC_MERGED}/cohort.vcf.list > ${GATK_HC_MERGED}/${PREFIX}.cohort.$(date -I).vcf

    # Compress the combined VCF file with bgzip
    bgzip -f ${GATK_HC_MERGED}/${PREFIX}.cohort.$(date -I).vcf

    # Create a tabix index for the compressed combined VCF file
    tabix -f ${GATK_HC_MERGED}/${PREFIX}.cohort.$(date -I).vcf.gz
    
    # Remove all files in the directory specified by GATK_HC_MERGED that match the pattern *tmp*
    #rm ${GATK_HC_MERGED}/*tmp*

}

export -f func_finish_vcf





#-------------------------------------------------------------------------------
# running the pipeline
#-------------------------------------------------------------------------------

# func_build_reference
#bsub -E 'test -e /nfs/users/nfs_s/sd21' -R "select[mem>1000] rusage[mem=1000]" -M1000 -o ${LOG_FILES}/gatk_01_build_reference.o -e ${LOG_FILES}/gatk_01_build_reference.e -J gatk_01_build_reference_${PREFIX} func_build_reference

# func_make_gvcf
#bsub -w "done(gatk_01_build_reference_${PREFIX})" -E 'test -e /nfs/users/nfs_s/sd21' -R "select[mem>50000] rusage[mem=50000]" -q long -M50000 -n20 -o ${LOG_FILES}/gatk_02_make_gvcf.o -e ${LOG_FILES}/gatk_02_make_gvcf.e -J gatk_02_make_gvcf_${PREFIX} func_make_gvcf

# func_merge_gvcf
bsub -E 'test -e /nfs/users/nfs_r/rp24' -R "select[mem>50000] rusage[mem=50000]" -q long -M50000 -n20 -o ${LOG_FILES}/gatk_03_merge_gvcf.o -e ${LOG_FILES}/gatk_03_merge_gvcf.e -J gatk_03_merge_gvcf_${PREFIX} func_merge_gvcf
# done

# func_genotype_gvcfs
bsub -E 'test -e /nfs/users/nfs_r/rp24' -R "select[mem>100000] rusage[mem=100000]" -q long -M100000 -n20 -o ${LOG_FILES}/gatk_04_genotype_gvcfs.o -e ${LOG_FILES}/gatk_04_genotype_gvcfs.e -J gatk_04_genotype_gvcfs_${PREFIX} func_genotype_gvcfs
# done

# func_finish_vcf
bsub -E 'test -e /nfs/users/nfs_r/rp24' -R "select[mem>50000] rusage[mem=50000]" -q long -M50000 -n1 -o ${LOG_FILES}/gatk_05_finish_vcf.o -e ${LOG_FILES}/gatk_05_finish_vcf.e -J gatk_05_finish_vcf_${PREFIX} func_finish_vcf
# done

# renamed the all sites VCF file to include _ALLSITES.vcf.gz
```






## Filter nuclear variants and invariants

```bash
VCF=/lustre/scratch127/pam/teams/team333/rp24/dirofilaria_immitis/R2/variants_ALLSITES/gatk_hc_dirofilaria_global/GATK_HC_MERGED/dirofilaria_global.cohort.2025-10-10_ALLSITES.vcf.gz
MIT_CONTIG=dirofilaria_immitis_chrMtDNA
WB_CONTIG=dirofilaria_immitis_chrWb
REFERENCE=/lustre/scratch127/pam/teams/team333/sd21/dirofilaria_immitis/project_global-diversity/01_references/dimmitis_WSI_2.2.fa

cd /lustre/scratch127/pam/teams/team333/rp24/dirofilaria_immitis/R2/variants_ALLSITES/gatk_hc_dirofilaria_global/GATK_HC_MERGED

vcftools --gzvcf dirofilaria_global.cohort.2025-10-10_ALLSITES.vcf.gz --remove-indels
#After filtering, kept 138 out of 138 Individuals
#After filtering, kept 87124750 out of a possible 88135118 Sites


#select nuclear invariants
bsub.py 1 select_nuclearINVARIANTs "gatk SelectVariants \
--reference ${REFERENCE} \
--variant ${VCF} \
--select-type-to-include NO_VARIATION \
--exclude-intervals ${MIT_CONTIG} \
--exclude-intervals ${WB_CONTIG} \
--output ${VCF%.vcf.gz}.nuclearINVARIANTs.vcf"
# done

vcftools --vcf dirofilaria_global.cohort.2025-10-10_ALLSITES.nuclearINVARIANTs.vcf --remove-indels
#After filtering, kept 138 out of 138 Individuals
#After filtering, kept 84268830 out of a possible 84352254 Sites
```

### Merge these nuclear invariants with the nuclear SNPs I previously filtered

```bash
#Keep n=129 samples in the nuclear dataset (keep_samples.dimmitis-only.nuclear_snps.list) and merge these with the SNPs
bsub.py 1 filter_nuclear_INVARIANTvcftools "vcftools --vcf ${VCF%.vcf.gz}.nuclearINVARIANTs.vcf \
--keep /lustre/scratch127/pam/teams/team333/sd21/dirofilaria_immitis/project_global-diversity/03_metadata/keep_samples.dimmitis-only.nuclear_snps.list \
--recode --out nuclearINVARIANTS_129samples"
# done
#After filtering, kept 129 out of 138 Individuals
#After filtering, kept 84352254 out of a possible 84352254 Sites

# Merge with the already filtered variants
bsub.py --done "filter_nuclear_INVARIANTvcftools" 1 merge_nuclear_VARIANTsandINVARIANTs "gatk MergeVcfs \
--INPUT /lustre/scratch127/pam/teams/team333/rp24/dirofilaria_immitis/R2/vcf/dirofilaria_global.cohort.2025-06-18.nuclearSNPs..keep_samples.dimmitis-only.recode.vcf \
--INPUT nuclearINVARIANTS_129samples.recode.vcf \
--OUTPUT nuclearVARIANTsandINVARIANTs_129samples.recode.vcf"
# done



## BY HOST (ALL SAMPLES)
# Select only the chrX to ch4 (avoid the scaffolds) and remove replicate samples
cd /lustre/scratch127/pam/teams/team333/rp24/dirofilaria_immitis/R2/pixy/host
bsub.py --done "merge_nuclear_VARIANTsandINVARIANTs" 1 filter_CHR "vcftools \
--vcf /lustre/scratch127/pam/teams/team333/rp24/dirofilaria_immitis/R2/variants_ALLSITES/gatk_hc_dirofilaria_global/GATK_HC_MERGED/nuclearVARIANTsandINVARIANTs_129samples.recode.vcf \
--chr dirofilaria_immitis_chrX \
--chr dirofilaria_immitis_chr1 \
--chr dirofilaria_immitis_chr2 \
--chr dirofilaria_immitis_chr3 \
--chr dirofilaria_immitis_chr4 \
--keep /lustre/scratch127/pam/teams/team333/sd21/dirofilaria_immitis/project_global-diversity/03_metadata/keep_samples.dimmitis-only.nuclear_snps.no-reps.list \
--remove-indels \
--recode --out nuclearSNPssandINVARIANTs_124samples.chrXto4_HOST"
# done
# After filtering, kept 124 out of 129 Individuals
# After filtering, kept 84270097 out of a possible 84653564 Sites
bgzip nuclearSNPssandINVARIANTs_124samples.chrXto4_HOST.recode.vcf;
tabix -p vcf nuclearSNPssandINVARIANTs_124samples.chrXto4_HOST.recode.vcf.gz


# BY REGION (DOG ONLY)
# Select only the chrX to ch4 (avoid the scaffolds) and remove replicate samples
cd /lustre/scratch127/pam/teams/team333/rp24/dirofilaria_immitis/R2/pixy/region_dog
bsub.py 1 filter_CHR_REGION_DOG "vcftools \
--vcf /lustre/scratch127/pam/teams/team333/rp24/dirofilaria_immitis/R2/variants_ALLSITES/gatk_hc_dirofilaria_global/GATK_HC_MERGED/nuclearVARIANTsandINVARIANTs_129samples.recode.vcf \
--chr dirofilaria_immitis_chrX \
--chr dirofilaria_immitis_chr1 \
--chr dirofilaria_immitis_chr2 \
--chr dirofilaria_immitis_chr3 \
--chr dirofilaria_immitis_chr4 \
--keep keep_samples.dimmitis-only.nuclear_snps.no-reps_DOG.list \
--remove-indels \
--recode --out nuclearSNPssandINVARIANTs_113samples.chrXto4_REGION_DOG"
# submitted
#After filtering, kept 113 out of 129 Individuals
#After filtering, kept 84270097 out of a possible 84653564 Sites

bgzip nuclearSNPssandINVARIANTs_113samples.chrXto4_REGION_DOG.recode.vcf;
tabix -p vcf nuclearSNPssandINVARIANTs_113samples.chrXto4_REGION_DOG.recode.vcf.gz


# AUSTRALIAN SAMPLES (DOG)
# Select only the chrX to ch4 (avoid the scaffolds) and remove replicate samples
cd /lustre/scratch127/pam/teams/team333/rp24/dirofilaria_immitis/R2/pixy/aus_dog
bsub.py 1 filter_CHR_AUS_DOG "vcftools \
--vcf /lustre/scratch127/pam/teams/team333/rp24/dirofilaria_immitis/R2/variants_ALLSITES/gatk_hc_dirofilaria_global/GATK_HC_MERGED/nuclearVARIANTsandINVARIANTs_129samples.recode.vcf \
--chr dirofilaria_immitis_chrX \
--chr dirofilaria_immitis_chr1 \
--chr dirofilaria_immitis_chr2 \
--chr dirofilaria_immitis_chr3 \
--chr dirofilaria_immitis_chr4 \
--keep keep_samples.dimmitis-only.nuclear_snps.no-reps_AUS_DOG.list \
--remove-indels \
--recode --out nuclearSNPssandINVARIANTs_44samples.chrXto4_AUS_DOG"
# done
#After filtering, kept 44 out of 129 Individuals
#After filtering, kept 84270097 out of a possible 84653564 Sites

bgzip nuclearSNPssandINVARIANTs_44samples.chrXto4_AUS_DOG.recode.vcf;
tabix -p vcf nuclearSNPssandINVARIANTs_44samples.chrXto4_AUS_DOG.recode.vcf.gz

```

### Run pixy on nuclear variants & invariants

```bash
# Create pixy environment
module load ISG/conda
conda create --name pixy
conda init --all
conda activate pixy

# Install pixy & htslib
conda install -c conda-forge pixy
conda install -c bioconda htslib

# Environmental variables
WORKING_DIR=/lustre/scratch127/pam/teams/team333/rp24/dirofilaria_immitis/R2

cd ${WORKING_DIR}/pixy

# Population files per region and host
# region_dog_pop.list
'
AUS_BNE_AD_001	AUS
AUS_BNE_AD_002	AUS
AUS_BNE_AD_003	AUS
AUS_BNE_AD_004	AUS
AUS_BNE_AD_008	AUS
AUS_BNE_AD_009	AUS
AUS_CNS_AD_001	AUS
AUS_CNS_AD_002	AUS
AUS_LHR_AD_001	AUS
AUS_ROK_AD_001	AUS
AUS_SYD_AD_001	AUS
AUS_SYD_AD_002	AUS
AUS_SYD_AD_003	AUS
AUS_SYD_AD_004	AUS
AUS_SYD_AD_005	AUS
AUS_SYD_AD_006	AUS
AUS_SYD_AD_007	AUS
AUS_SYD_AD_008	AUS
AUS_SYD_AD_009	AUS
AUS_SYD_AD_010	AUS
AUS_SYD_AD_011	AUS
AUS_SYD_AD_012	AUS
AUS_SYD_AD_013	AUS
AUS_SYD_AD_014	AUS
AUS_SYD_AD_015	AUS
AUS_SYD_AD_017	AUS
AUS_TVS_AD_001	AUS
AUS_TVS_AD_002	AUS
AUS_TVS_AD_003	AUS
AUS_TVS_AD_004	AUS
AUS_TVS_AD_005	AUS
AUS_TVS_AD_006	AUS
AUS_TVS_AD_007	AUS
AUS_TVS_AD_008	AUS
AUS_TVS_AD_010	AUS
AUS_TVS_AD_011	AUS
AUS_TVS_AD_012	AUS
AUS_TVS_AD_013	AUS
AUS_TVS_AD_014	AUS
AUS_TVS_AD_015	AUS
AUS_TVS_AD_016	AUS
AUS_TVS_AD_017	AUS
AUS_TVS_AD_018	AUS
AUS_TVS_AD_019	AUS
GRC_XAN_AD_001	EUR
GRC_XAN_AD_002	EUR
GRC_XAN_AD_003	EUR
GRC_XAN_AD_004	EUR
GRC_XAN_AD_005	EUR
GRC_XAN_AD_006	EUR
GRC_XAN_AD_007	EUR
GRC_XAN_AD_008	EUR
GRC_XAN_AD_009	EUR
GRC_XAN_AD_010	EUR
GRC_XAN_AD_011	EUR
GRC_XAN_AD_012	EUR
ITA_NEA_AD_001	EUR
ITA_PAV_AD_001	EUR
MYS_SEL_AD_001	ASIA
PAN_BOC_AD_001	CENAM
PAN_BOC_AD_002	CENAM
PAN_BOC_AD_003	CENAM
PAN_BOC_AD_004	CENAM
PAN_BOC_AD_005	CENAM
PAN_PUE_AD_001	CENAM
PAN_PUE_AD_002	CENAM
PAN_PUE_AD_003	CENAM
PAN_PUE_AD_004	CENAM
PAN_PUE_AD_005	CENAM
PAN_PUE_AD_006	CENAM
PAN_SLO_AD_001	CENAM
PAN_SLO_AD_002	CENAM
PAN_SLO_AD_003	CENAM
THA_BKK_AD_001	ASIA
THA_BKK_AD_002	ASIA
THA_BKK_AD_003	ASIA
THA_BKK_AD_004	ASIA
THA_BKK_AD_005	ASIA
THA_BKK_AD_006	ASIA
USA_FLO_AD_001	USA
USA_FLO_AD_002	USA
USA_FLO_AD_003	USA
USA_FLO_AD_004	USA
USA_FLO_AD_005	USA
USA_FLO_AD_006	USA
USA_FLO_AD_007	USA
USA_FLO_AD_008	USA
USA_FLO_AD_009	USA
USA_FLO_AD_010	USA
USA_FLO_AD_011	USA
USA_FLO_AD_012	USA
USA_FLO_AD_013	USA
USA_FLO_AD_014	USA
USA_FLO_AD_015	USA
USA_FLO_AD_016	USA
USA_GEO_AD_001	USA
USA_ILL_AD_001	USA
USA_ILL_AD_002	USA
USA_LOU_AD_001	USA
USA_LOU_AD_002	USA
USA_MIS_AD_001	USA
USA_MIS_AD_002	USA
USA_TEX_AD_002	USA
USA_TEX_AD_003	USA
USA_TEX_AD_004	USA
USA_TEX_AD_005	USA
USA_TEX_AD_006	USA
USA_TEX_AD_007	USA
USA_TEX_AD_008	USA
USA_TEX_AD_009	USA
USA_TEX_AD_010	USA
USA_TEX_AD_013	USA
USA_TEX_AD_014	USA
'

# host_pop.list
'
AUS_BNE_AD_001	Dog
AUS_BNE_AD_002	Dog
AUS_BNE_AD_003	Dog
AUS_BNE_AD_004	Dog
AUS_BNE_AD_006	Fox
AUS_BNE_AD_008	Dog
AUS_BNE_AD_009	Dog
AUS_CNS_AD_001	Dog
AUS_CNS_AD_002	Dog
AUS_LHR_AD_001	Dog
AUS_ROK_AD_001	Dog
AUS_SYD_AD_001	Dog
AUS_SYD_AD_002	Dog
AUS_SYD_AD_003	Dog
AUS_SYD_AD_004	Dog
AUS_SYD_AD_005	Dog
AUS_SYD_AD_006	Dog
AUS_SYD_AD_007	Dog
AUS_SYD_AD_008	Dog
AUS_SYD_AD_009	Dog
AUS_SYD_AD_010	Dog
AUS_SYD_AD_011	Dog
AUS_SYD_AD_012	Dog
AUS_SYD_AD_013	Dog
AUS_SYD_AD_014	Dog
AUS_SYD_AD_015	Dog
AUS_SYD_AD_017	Dog
AUS_TVS_AD_001	Dog
AUS_TVS_AD_002	Dog
AUS_TVS_AD_003	Dog
AUS_TVS_AD_004	Dog
AUS_TVS_AD_005	Dog
AUS_TVS_AD_006	Dog
AUS_TVS_AD_007	Dog
AUS_TVS_AD_008	Dog
AUS_TVS_AD_010	Dog
AUS_TVS_AD_011	Dog
AUS_TVS_AD_012	Dog
AUS_TVS_AD_013	Dog
AUS_TVS_AD_014	Dog
AUS_TVS_AD_015	Dog
AUS_TVS_AD_016	Dog
AUS_TVS_AD_017	Dog
AUS_TVS_AD_018	Dog
AUS_TVS_AD_019	Dog
CRI_SJO_AD_001	Cat
GRC_XAN_AD_001	Dog
GRC_XAN_AD_002	Dog
GRC_XAN_AD_003	Dog
GRC_XAN_AD_004	Dog
GRC_XAN_AD_005	Dog
GRC_XAN_AD_006	Dog
GRC_XAN_AD_007	Dog
GRC_XAN_AD_008	Dog
GRC_XAN_AD_009	Dog
GRC_XAN_AD_010	Dog
GRC_XAN_AD_011	Dog
GRC_XAN_AD_012	Dog
ITA_NEA_AD_001	Dog
ITA_NEA_AD_002	Fox
ITA_NEA_AD_003	Fox
ITA_PAV_AD_001	Dog
MYS_SEL_AD_001	Dog
PAN_BOC_AD_001	Dog
PAN_BOC_AD_002	Dog
PAN_BOC_AD_003	Dog
PAN_BOC_AD_004	Dog
PAN_BOC_AD_005	Dog
PAN_PUE_AD_001	Dog
PAN_PUE_AD_002	Dog
PAN_PUE_AD_003	Dog
PAN_PUE_AD_004	Dog
PAN_PUE_AD_005	Dog
PAN_PUE_AD_006	Dog
PAN_SLO_AD_001	Dog
PAN_SLO_AD_002	Dog
PAN_SLO_AD_003	Dog
ROU_BUC_AD_001	Leopard
ROU_COM_AD_001	Wildcat
ROU_GIU_AD_001	Golden_jackal
THA_BKK_AD_001	Dog
THA_BKK_AD_002	Dog
THA_BKK_AD_003	Dog
THA_BKK_AD_004	Dog
THA_BKK_AD_005	Dog
THA_BKK_AD_006	Dog
USA_FLO_AD_001	Dog
USA_FLO_AD_002	Dog
USA_FLO_AD_003	Dog
USA_FLO_AD_004	Dog
USA_FLO_AD_005	Dog
USA_FLO_AD_006	Dog
USA_FLO_AD_007	Dog
USA_FLO_AD_008	Dog
USA_FLO_AD_009	Dog
USA_FLO_AD_010	Dog
USA_FLO_AD_011	Dog
USA_FLO_AD_012	Dog
USA_FLO_AD_013	Dog
USA_FLO_AD_014	Dog
USA_FLO_AD_015	Dog
USA_FLO_AD_016	Dog
USA_GEO_AD_001	Dog
USA_ILL_AD_001	Dog
USA_ILL_AD_002	Dog
USA_LOU_AD_001	Dog
USA_LOU_AD_002	Dog
USA_MIS_AD_001	Dog
USA_MIS_AD_002	Dog
USA_TEX_AD_001	Ferret
USA_TEX_AD_002	Dog
USA_TEX_AD_003	Dog
USA_TEX_AD_004	Dog
USA_TEX_AD_005	Dog
USA_TEX_AD_006	Dog
USA_TEX_AD_007	Dog
USA_TEX_AD_008	Dog
USA_TEX_AD_009	Dog
USA_TEX_AD_010	Dog
USA_TEX_AD_011	Cat
USA_TEX_AD_012	Cat
USA_TEX_AD_013	Dog
USA_TEX_AD_014	Dog
USA_TEX_AD_015	Cat
'

# aus_dog_pop.list:
'
AUS_BNE_AD_001	Brisbane
AUS_BNE_AD_002	Brisbane
AUS_BNE_AD_003	Brisbane
AUS_BNE_AD_004	Brisbane
AUS_BNE_AD_008	Brisbane
AUS_BNE_AD_009	Brisbane
AUS_CNS_AD_001	Cairns
AUS_CNS_AD_002	Cairns
AUS_LHR_AD_001	Lockhart_River
AUS_ROK_AD_001	Rockhampton
AUS_SYD_AD_001	Sydney
AUS_SYD_AD_002	Sydney
AUS_SYD_AD_003	Sydney
AUS_SYD_AD_004	Sydney
AUS_SYD_AD_005	Sydney
AUS_SYD_AD_006	Sydney
AUS_SYD_AD_007	Sydney
AUS_SYD_AD_008	Sydney
AUS_SYD_AD_009	Sydney
AUS_SYD_AD_010	Sydney
AUS_SYD_AD_011	Sydney
AUS_SYD_AD_012	Sydney
AUS_SYD_AD_013	Sydney
AUS_SYD_AD_014	Sydney
AUS_SYD_AD_015	Sydney
AUS_SYD_AD_017	Sydney
AUS_TVS_AD_001	Townsville
AUS_TVS_AD_002	Townsville
AUS_TVS_AD_003	Townsville
AUS_TVS_AD_004	Townsville
AUS_TVS_AD_005	Townsville
AUS_TVS_AD_006	Townsville
AUS_TVS_AD_007	Townsville
AUS_TVS_AD_008	Townsville
AUS_TVS_AD_010	Townsville
AUS_TVS_AD_011	Townsville
AUS_TVS_AD_012	Townsville
AUS_TVS_AD_013	Townsville
AUS_TVS_AD_014	Townsville
AUS_TVS_AD_015	Townsville
AUS_TVS_AD_016	Townsville
AUS_TVS_AD_017	Townsville
AUS_TVS_AD_018	Townsville
AUS_TVS_AD_019	Townsville
'

# Run pixy per region (dog only)
VCF=region_dog/nuclearSNPssandINVARIANTs_113samples.chrXto4_REGION_DOG.recode.vcf.gz
bsub.py --queue long --threads 10 40 pixy_region_dog \
"pixy --stats pi fst dxy \
--vcf ${VCF} \
--populations region_dog/region_dog_pop.list \
--window_size 100000 \
--n_cores 10 \
--output_folder region_dog \
--output_prefix region_dog"
# Submitted

# Run pixy per host (all samples)
VCF=host/nuclearSNPssandINVARIANTs_124samples.chrXto4.recode.vcf.gz
bsub.py --queue long --threads 10 40 pixy_host \
"pixy --stats pi fst dxy \
--vcf ${VCF} \
--populations host/host_pop.list \
--window_size 100000 \
--n_cores 10 \
--output_folder host \
--output_prefix host"
# Done
## Keep in mind that the sample sizes for non-dog hosts are quite low, but will still be interesting to explore

# Run pixy AUS (dog only)
VCF=aus_dog/nuclearSNPssandINVARIANTs_44samples.chrXto4_AUS_DOG.recode.vcf.gz
bsub.py --queue long --threads 10 40 pixy_aus_dog \
"pixy --stats pi fst dxy \
--vcf ${VCF} \
--populations aus_dog/aus_dog_pop.list \
--window_size 100000 \
--n_cores 10 \
--output_folder aus_dog \
--output_prefix aus_dog"
# Submitted

```

Download the files and use it as input in R to estimate the values of pi, fst and dxy and generate plots.

### Pixy by region (dog)

```R
# Batch 4 - Pixy for Pi, Fst and Dxy by REGION for dog samples only

library(tidyverse)
library(ggsci)
library(ggpubr)
library(patchwork)
library(ggridges)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(purrr)
library(stats)
library(graphics)
library(grDevices)
library(utils)
library(datasets)
library(methods)
library(base)
require(maps)
library(mapdata)
library(readxl)
library(ozmaps) 
library(grid)
library(gridExtra)
library(ggrepel)
library(ggnewscale)
library(reshape)
library(wesanderson)

setwd("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/R2_Sep25/pixy/region_dog")

##### Nucleotide diversity (Pi) #####

# get nucleotide diversity (pi) data from pixy output
pi_data <- read.table("input/region_dog_pi.txt", header=T)
pi_data$chromosome <- str_remove(pi_data$chromosome, 'dirofilaria_immitis_')

# subset
pi_data_AUS <- pi_data %>%
  filter(pop=="AUS")
pi_data_ASIA <- pi_data %>%
  filter(pop=="ASIA")
pi_data_USA <- pi_data %>%
  filter(pop=="USA")
pi_data_EUR <- pi_data %>%
  filter(pop=="EUR")
pi_data_CENAM <- pi_data %>%
  filter(pop=="CENAM")

# filter pi data to remove small scaffolds not in the linkage groups, and to number the rows per group to help with plotting
pi_data <- pi_data %>%
  group_by(pop) %>%
  mutate(position = 1:n())

# get position for vertical lines used in plot to delineate the linkage groups
pi_data %>%
  group_by(chromosome) %>%
  summarise(max = max(position, na.rm = TRUE))

'
# A tibble: 5 × 2
  chromosome   max
  <chr>      <int>
1 chr1         440
2 chr2         592
3 chr3         744
4 chr4         885
5 chrX         283
'

#Let's add the chr type variable
pi_data <- pi_data %>%
  mutate(chr_type = ifelse(str_detect(chromosome, "X"), "sexchr", "autosome"))

# calculate the median Pi and check the ratio of sex-to-autosome diversity. Should be about 0.75, as D. immitis is XX/XY
pi_data_sex_median <- pi_data %>%
  group_by(chr_type) %>%
  summarise(median = median(avg_pi, na.rm = TRUE))
pi_data_sex_median

'
# A tibble: 2 × 2
  chr_type   median
  <chr>       <dbl>
1 autosome 0.000492
2 sexchr   0.000254
'
# 0.000254 / 0.000492 = 0.5162602 (bit off 0.75)

pi_data_pop_sex_median <- pi_data %>%
  group_by(pop, chr_type) %>%
  summarise(median = median(avg_pi, na.rm = TRUE))
pi_data_pop_sex_median

'
# A tibble: 10 × 3
# Groups:   pop [5]
   pop   chr_type   median
   <chr> <chr>       <dbl>
 1 ASIA  autosome 0.000312
 2 ASIA  sexchr   0.000144
 3 AUS   autosome 0.000629
 4 AUS   sexchr   0.000275
 5 CENAM autosome 0.000385
 6 CENAM sexchr   0.000222
 7 EUR   autosome 0.000461
 8 EUR   sexchr   0.000254
 9 USA   autosome 0.000777
10 USA   sexchr   0.000440
'

# plot 1 - genome wide plots per population
plot_1_pi <- ggplot(pi_data, aes(position*100000, avg_pi, col=chromosome)) +
  geom_point(size=0.8) +
  facet_grid(pop~.) +
  scale_color_npg(guide = guide_legend(override.aes = list(size = 4))) +
  theme_bw() +
  theme(legend.position="top", 
        legend.title = element_blank(),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 12, color = "black"),
        strip.text = element_text(size = 18),
        legend.text = element_text(size = 16),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 2)) +
  labs(x="Genomic Position", y="Nucleotide Diversity (Pi)")
plot_1_pi


# plot 2 - density plots of pi per group
plot_2_pi <- ggplot(pi_data, aes(avg_pi, chr_type, fill=chr_type), guide="none") +
  geom_density_ridges(quantile_lines=TRUE, quantile_fun=function(x,...)median(x), size=0.5) +
  theme_bw() + theme(legend.position = "none", axis.text.y=element_blank()) +
  facet_grid(pop~.) +
  xlim(0, 0.005) +
  scale_x_continuous(breaks = c(0, 0.0025, 0.005)) +
  scale_fill_npg() +
  theme(legend.position="top",
        legend.title = element_blank(),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 12, color = "black"),
        strip.text = element_text(size = 18),
        legend.text = element_text(size = 16),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 2)) +
  labs(x="Nucleotide Diversity (Pi)", y="Density")
plot_2_pi

# combine plots
plot_1_pi + plot_2_pi +  plot_layout(widths = c(5, 1))
ggsave("genomewide_and_density_Pi_region_dog.tif", width=14, height=8, dpi = 300)
ggsave("genomewide_and_density_Pi_region_dog.png", width=14, height=8, dpi = 300)
ggsave("genomewide_and_density_Pi_region_dog.pdf", width=14, height=8, dpi = 300)

#Now a boxplot of the pi value per population

red_palette1 <- brewer.pal(n = 9, name = "YlOrRd")
red_palette2 <- brewer.pal(n=9, name = "YlOrBr")
blue_palette <- brewer.pal(n = 9, name = "Blues")
green_palette1 <- brewer.pal(n = 9, name = "BuGn")
green_palette2 <- brewer.pal(n=9, name = "YlGn")

scale_colour_region <- function(...){
  ggplot2:::manual_scale(
    'colour', 
    values = setNames(
      c(blue_palette[5],
        "violetred1",
        red_palette1[7], 
        green_palette1[7],
        "blueviolet"),
    c("AUS", "ASIA", "USA", "EUR", "CENAM")), 
  ...
)
}    

boxplot_pi <- ggplot(pi_data, aes(pop, avg_pi, col=pop)) +
  geom_jitter(size = 1, alpha = 0.5) +
  geom_boxplot(fill=NA, col="black") +
  labs(x = "Population" , y = "Nucleotide diversity (Pi)", colour = "Population") +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 18, color = "black"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 2)
        ) +
  scale_colour_region ()

boxplot_pi

ggsave("boxplot_Pi_region_dog_grid.tif", width=4, height=6, dpi = 300)
ggsave("boxplot_Pi_region_dog_grid.png", width=4, height=6, dpi = 300)
ggsave("boxplot_Pi_region_dog_grid.pdf", width=4, height=6, dpi = 300)
# higher pi value = more diverse population
# lower pi value = less diverse population
## Looks like USA and AUS have higher diversity.


#let's check the statistical significance of the differences between pi values
#Let's explore the normality
pi_data %>% filter(chromosome != 'chrX') %>%
  ggplot(., aes(x = avg_pi, colour = pop)) +
  geom_histogram() +
  theme_bw() +
  facet_grid(. ~ pop)

#Confirm the absence of normality by shapiro
## If p > 0.05 the data is normal, if p < 0.05 the data is NOT normal
pi_data_AUS <- pi_data_AUS %>% filter(chromosome != 'chrX')
AUS_shapiro <- shapiro.test(pi_data_AUS$avg_pi)
print(AUS_shapiro)
# W = 0.9118, p-value < 2.2e-16

pi_data_ASIA <- pi_data_ASIA %>% filter(chromosome != 'chrX')
ASIA_shapiro <- shapiro.test(pi_data_ASIA$avg_pi)
print(ASIA_shapiro)
# W = 0.86847, p-value < 2.2e-16

pi_data_USA <- pi_data_USA %>% filter(chromosome != 'chrX')
USA_shapiro <- shapiro.test(pi_data_USA$avg_pi)
print(USA_shapiro)
# W = 0.91081, p-value < 2.2e-16

pi_data_EUR <- pi_data_EUR %>% filter(chromosome != 'chrX')
EUR_shapiro <- shapiro.test(pi_data_EUR$avg_pi)
print(EUR_shapiro)
# W = 0.91439, p-value < 2.2e-16

pi_data_CENAM <- pi_data_CENAM %>% filter(chromosome != 'chrX')
CENAM_shapiro <- shapiro.test(pi_data_CENAM$avg_pi)
print(CENAM_shapiro)
# W = 0.88403, p-value < 2.2e-16

#Wilcoxson test for all pairs of pops
## This tests whether the mean values of 2 dependent groups differ significantly from each other
## If p < 0.05 they are sig different
wilcox.test(pi_data_AUS$avg_pi, pi_data_ASIA$avg_pi)
# W = 254902, p-value < 2.2e-16
wilcox.test(pi_data_AUS$avg_pi, pi_data_USA$avg_pi)
# W = 149019, p-value = 9.559e-08
wilcox.test(pi_data_AUS$avg_pi, pi_data_EUR$avg_pi)
# W = 221303, p-value = 2.985e-11
wilcox.test(pi_data_AUS$avg_pi, pi_data_CENAM$avg_pi)
# W = 235273, p-value < 2.2e-16
wilcox.test(pi_data_ASIA$avg_pi, pi_data_USA$avg_pi)
# W = 86465, p-value < 2.2e-16
wilcox.test(pi_data_ASIA$avg_pi, pi_data_EUR$avg_pi)
# W = 144675, p-value = 1.404e-09
wilcox.test(pi_data_ASIA$avg_pi, pi_data_CENAM$avg_pi)
# W = 156031, p-value = 3.012e-05
wilcox.test(pi_data_USA$avg_pi, pi_data_EUR$avg_pi)
# W = 247376, p-value < 2.2e-16
wilcox.test(pi_data_USA$avg_pi, pi_data_CENAM$avg_pi)
# W = 259934, p-value < 2.2e-16
wilcox.test(pi_data_EUR$avg_pi, pi_data_CENAM$avg_pi)
# W = 195316, p-value = 0.01931

## They're all sig different


#let's generate some dataframes for the estatistics
# subset
pi_data_AUS <- pi_data %>%
  filter(pop=="AUS") %>%
  filter(chromosome != 'chrX')

pi_data_ASIA <- pi_data %>%
  filter(pop=="ASIA") %>%
  filter(chromosome != 'chrX') 

pi_data_USA <- pi_data %>%
  filter(pop=="USA") %>%
  filter(chromosome != 'chrX')

pi_data_EUR <- pi_data %>%
  filter(pop=="EUR") %>%
  filter(chromosome != 'chrX')

pi_data_CENAM <- pi_data %>%
  filter(pop=="CENAM") %>%
  filter(chromosome != 'chrX')




#### Dxy and Fst ####

# load data
dxy_data <- read.table("input/region_dog_dxy.txt", header=T)
fst_data <- read.table("input/region_dog_fst.txt", header=T)

# add some columns
dxy_data$data_type <- "Dxy"
dxy_data <- mutate(dxy_data,
                   comparison = paste(pop1, pop2, sep = '_v_'))
fst_data$data_type <- "Fst"
fst_data <- mutate(fst_data,
                   comparison = paste(pop1, pop2, sep = '_v_'))

# subset and merge dataframes
dxy_data_sub <- dxy_data %>% select(comparison,data_type,chromosome,window_pos_1,window_pos_2,avg_dxy)
colnames(dxy_data_sub) <- c("comparison","data_type","chromosome","window_pos_1","window_pos_2","value")

fst_data_sub <- fst_data %>% select(comparison,data_type,chromosome,window_pos_1,window_pos_2,avg_wc_fst)
colnames(fst_data_sub) <- c("comparison","data_type","chromosome","window_pos_1","window_pos_2","value")

# cheeky fix to get matching rows from both datasets
tmp_dxy <- semi_join(dxy_data_sub, fst_data_sub, by=c("comparison", "chromosome", "window_pos_1", "window_pos_2"))
tmp_fst <- semi_join(fst_data_sub, dxy_data_sub,  by=c("comparison", "chromosome", "window_pos_1", "window_pos_2"))

# join the datasets together to create a single dataframe
data <- full_join(tmp_dxy, tmp_fst)
data$chromosome <- str_remove(data$chromosome, 'dirofilaria_immitis_')

# add numbering to help with plotting
data <- data %>%
  group_by(comparison, data_type) %>%
  mutate(position = 1:n())

# add sex chromosome information
data <- data %>%
  mutate(chr_type = ifelse(str_detect(chromosome, "X"), "sexchr", "autosome"))


# summarise median data for Fst and Dxy
data %>%
  group_by(comparison,data_type) %>%
  summarise(median = median(value, na.rm = TRUE))
'
# A tibble: 20 × 3
# Groups:   comparison [10]
   comparison   data_type   median
   <chr>        <chr>        <dbl>
 1 ASIA_v_CENAM Dxy       0.000747
 2 ASIA_v_CENAM Fst       0.512   
 3 ASIA_v_USA   Dxy       0.000661
 4 ASIA_v_USA   Fst       0.204   
 5 AUS_v_ASIA   Dxy       0.000581
 6 AUS_v_ASIA   Fst       0.224   
 7 AUS_v_CENAM  Dxy       0.000702
 8 AUS_v_CENAM  Fst       0.328   
 9 AUS_v_EUR    Dxy       0.000705
10 AUS_v_EUR    Fst       0.315   
11 AUS_v_USA    Dxy       0.000713
12 AUS_v_USA    Fst       0.179   
13 CENAM_v_USA  Dxy       0.000724
14 CENAM_v_USA  Fst       0.240   
15 EUR_v_ASIA   Dxy       0.000737
16 EUR_v_ASIA   Fst       0.490   
17 EUR_v_CENAM  Dxy       0.000593
18 EUR_v_CENAM  Fst       0.289   
19 EUR_v_USA    Dxy       0.000719
20 EUR_v_USA    Fst       0.225   
'


#Plotting Dxy

# plot 1 - genome wide plots per comparison
plot_1_dxy <- data %>%
  filter(data_type =='Dxy')%>%
  ggplot(., aes(position*100000, value, col=chromosome)) +
  geom_point(size=1.5) +
  facet_grid(comparison~.) +
  scale_color_npg(guide = guide_legend(override.aes = list(size = 4))) +
  theme_bw() +
  theme(legend.position="top", 
        legend.title = element_blank(),
        axis.title = element_text(size = 32),
        axis.text = element_text(size = 16, color = "black"),
        strip.text = element_text(size = 14),
        legend.text = element_text(size = 20),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 2)) +
  labs(x="Genomic Position", y="Dxy")
plot_1_dxy

# plot 2 - density plots of dxy per group
plot_2_dxy <- data %>%
  filter(data_type=='Dxy')%>%
  ggplot(aes(value, chr_type, fill=chr_type), guide="none") +
  geom_density_ridges(quantile_lines=TRUE, quantile_fun=function(x,...)median(x), size=0.5) +
  theme_bw() + theme(legend.position = "top", 
                     legend.title = element_blank(),
                     axis.text.y=element_blank(),
                     axis.title = element_text(size = 32),
                     axis.text = element_text(size = 16, color = "black"),
                     strip.text = element_text(size = 14),
                     legend.text = element_text(size = 20),
                     panel.grid = element_blank(),
                     panel.border = element_rect(color = "black", fill = NA, linewidth = 2)) +
  facet_grid(comparison~.) +
  scale_fill_npg() +
  labs(x="Dxy", y="Density")
plot_2_dxy

# combine plots
dxy_plot <- plot_1_dxy + plot_2_dxy +  plot_layout(widths = c(5, 1))
ggsave("genomewide_and_density_dxy_region_dog.tif", width=16, height=18, dpi = 300)
ggsave("genomewide_and_density_dxy_region_dog.png", width=16, height=18, dpi = 300)
ggsave("genomewide_and_density_dxy_region_dog.pdf", width=16, height=18, dpi = 300)
# higher Dxy = more divergent

#some additional plots
boxplot_dxy <- data %>%
  filter(data_type =='Dxy') %>% 
  ggplot(., aes(comparison, value, col=comparison)) +
  geom_jitter(size = 1, alpha = 0.5) +
  geom_boxplot(fill=NA, col="black") +
  labs(x = "Population" , y = "Dxy", colour = "Population") +
  theme_bw() +
  scale_color_npg() +
  theme(legend.position = "none",
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 18, color = "black"),
        strip.text = element_text(size = 18),
        legend.text = element_text(size = 16),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 2))
boxplot_dxy
ggsave("boxplot_dxy_region_dog.tif", width=14, height = 6, dpi = 300)
ggsave("boxplot_dxy_region_dog.png", width=14, height = 6, dpi = 300)
ggsave("boxplot_dxy_region_dog.pdf", width=14, height = 6, dpi = 300)

density_dxy <- data %>%
  filter(data_type =='Dxy') %>% 
  ggplot(., aes(x=value, group = comparison, fill=comparison)) +
  geom_density(adjust=1.5, alpha=.4) +
  labs(x = "Dxy" , y = "Density", colour = "Comparison") +
  theme_bw() +
  scale_fill_npg() +
  theme(legend.position = c(0.8,0.72),
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 18, color = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 2))
density_dxy
ggsave("density_dxy_region_dog.tif", width = 8, height = 6, dpi = 300)
ggsave("density_dxy_region_dog.png", width = 8, height = 6, dpi = 300)
ggsave("density_dxy_region_dog.pdf", width = 8, height = 6, dpi = 300)

# AUS_v_ASIA & AUS_v_EUR
scatter_dxy_a <-  data %>%
  filter(data_type == 'Dxy') %>%
  as.data.frame() %>%
  select('position', 'comparison', 'value', 'chromosome')%>%
  as_tibble() %>%
  pivot_wider(names_from = comparison, values_from = value) %>%
  ggplot(., aes(x=AUS_v_ASIA, y = AUS_v_EUR, col=chromosome)) +
  geom_point(alpha=.4) +
  labs(x = "AUS_v_ASIA" , y = "AUS_v_EUR", colour = "Comparison") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
scatter_dxy_a
# can repeat for other pairs

genome_pos_dxy_a <-  data %>%
  filter(data_type == 'Dxy') %>%
  as.data.frame() %>%
  select('position', 'comparison', 'value', 'chromosome')%>%
  as_tibble() %>%
  pivot_wider(names_from = comparison, values_from = value) %>%
  mutate("xx" = AUS_v_ASIA - AUS_v_EUR,
         "yy" = AUS_v_ASIA - EUR_v_ASIA) %>%
  ggplot(., aes(position*100000, xx, col=chromosome)) +
  geom_point(size=1) +
  scale_color_tron() +
  theme_bw() +
  theme(legend.position="top", 
        legend.title = element_blank()) +
  labs(x="Position", y="AUS vs ASIA - AUS vs EUR")
genome_pos_dxy_a
# can repeat for other pairs

genome_pos_dxy_aa <-  data %>%
  filter(data_type == 'Dxy') %>%
  as.data.frame() %>%
  select('position', 'comparison', 'value', 'chromosome')%>%
  as_tibble() %>%
  pivot_wider(names_from = comparison, values_from = value) %>%
  mutate("xx" = AUS_v_ASIA / AUS_v_EUR,
         "yy" = AUS_v_ASIA / EUR_v_ASIA) %>%
  ggplot(., aes(position*100000, xx, col=chromosome)) +
  geom_point(size=1) +
  scale_color_tron() +
  theme_bw() +
  theme(legend.position="top", 
        legend.title = element_blank()) +
  labs(x="Position", y="AUS vs ASIA / AUS vs EUR")
genome_pos_dxy_aa
# can repeat for other pairs

# ggarrange above plots if using them


#Trying to plot lines
#getting chr positions
chr <- c('X', '1', '2', '3', '4')

for (i in chr) {
  x <- data %>%
    filter(chromosome == paste0('chr', i))
  print(quantile(x$position))
}

data$comparison <- str_replace_all(data$comparison, '_', ' ')
data$comparison <- str_replace_all(data$comparison, 'v', 'vs')

dxy_lineplot <- data %>%
  filter(data_type=='Dxy')%>%
  ggplot(., aes(position*100000, value, col=comparison)) +
  geom_point(size=0.2, alpha = 0.1) +
  geom_vline(xintercept=c(281*100000,
                          438*100000,
                          590*100000,
                          740*100000),
             size=1, linetype="dashed", col='grey41')+
  annotate(geom="text", x=140*100000, y=0.004, label="ChrX", size = 6)+
  annotate(geom="text", x=359.5*100000, y=0.004, label="Chr1", size = 6)+
  annotate(geom="text", x=520*100000, y=0.004, label="Chr2", size = 6)+
  annotate(geom="text", x=665*100000, y=0.004, label="Chr3", size = 6)+
  annotate(geom="text", x=810.5*100000, y=0.004, label="Chr4", size = 6)+
  geom_line() +
  scale_color_npg() +
  theme_bw() +
  theme(legend.position="top", 
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 14)) +
  labs(x="Genomic Position", y="Dxy") +
  guides(col = guide_legend(override.aes = list(linewidth = 2)))
dxy_lineplot
ggsave("lineplot_dxy_region_dog.tif", width=14, height=8, dpi = 300)
ggsave("lineplot_dxy_region_dog.png", width=14, height=8, dpi = 300)
ggsave("lineplot_dxy_region_dog.pdf", width=14, height=8, dpi = 300)


# Dxy by region on a heatmap
# Get median dxy for each comparison
dxy_median <- data %>%
  filter(data_type == "Dxy") %>%
  group_by(comparison) %>%
  summarise(DXY_MEDIAN = median(value, na.rm = TRUE))

# Save as csv file
write.csv(dxy_median, "dxy_median.csv", row.names = TRUE)

dxy_heatmap <- read.csv("dxy_heatmap.csv", header = TRUE)
dxy_heatmap <- as.data.frame(dxy_heatmap)

plot_3_dxy <- ggplot(dxy_heatmap, aes(x = POP1, y=POP2, fill = DXY_MEDIAN)) + 
  geom_tile(color = "white", lwd = 1.5, linetype = 1) +
  coord_fixed() +
  scale_fill_material("indigo", labels = scales::label_number()) +
  guides(fill = guide_colourbar(title = expression(D[XY]), barwidth = 1, barheight = 5)) +
  theme_minimal() +
  labs(x = "Population", y = "Population") +
  theme(axis.text = element_text(size = 18),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title = element_text(size = 24),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank())
plot_3_dxy

ggsave("dxy_heatmap_region_dog.tif", bg = "white", height=5, width = 8, dpi = 300)
ggsave("dxy_heatmap_region_dog.png", bg = "white", height=5, width = 8, dpi = 300)
ggsave("dxy_heatmap_region_dog.pdf", bg = "white", height=5, width = 8, dpi = 300)


# Plotting Fst

# plot 1 - genome wide plots per population
plot_1_fst <- data %>%
  filter(data_type=='Fst')%>%
  ggplot(., aes(position*100000, value, col=chromosome)) +
  geom_point(size=1.5) +
  facet_grid(comparison~.) +
  scale_color_npg(guide = guide_legend(override.aes = list(size = 4))) +
  theme_bw() +
  theme(legend.position="top", 
        legend.title = element_blank(),
        axis.title = element_text(size = 32),
        axis.text = element_text(size = 16, color = "black"),
        strip.text = element_text(size = 14),
        legend.text = element_text(size = 20),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 2)) +
  labs(x="Genomic Position", y="Fst")
plot_1_fst


# plot 2 - density plots of Fst per group
plot_2_fst <- data %>%
  filter(data_type=='Fst')%>%
  ggplot(aes(value, chr_type, fill=chr_type), guide="none") +
  geom_density_ridges(quantile_lines=TRUE, quantile_fun=function(x,...)median(x), size=0.5) +
  theme_bw() + 
  theme(legend.position = "none", 
        axis.text.y=element_blank(),
        axis.title = element_text(size = 32),
        axis.text = element_text(size = 16, color = "black"),
        strip.text = element_text(size = 14),
        legend.text = element_text(size = 20),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 2)) +
  facet_grid(comparison~.) +
  scale_fill_npg() +
  labs(x="Fst", y="Density")
plot_2_fst

# combine plots
plot_1_fst + plot_2_fst +  plot_layout(widths = c(5, 1))
# 1: Removed 35 rows containing missing values or values outside the scale range (`geom_point()`). 
# 2: Removed 35 rows containing non-finite outside the scale range (`stat_density_ridges()`).
ggsave("genomewide_and_density_fst_region_dog.tif", width=16, height=18, dpi = 300)
ggsave("genomewide_and_density_fst_region_dog.png", width=16, height=18, dpi = 300)
ggsave("genomewide_and_density_fst_region_dog.pdf", width=16, height=18, dpi = 300)
# Fst closer to 1 = more genetically distinct
# Fst closer to 0 = more genetically similar


# Region map metadata. Just have one random point representing each region.
location <- read.csv("input/region_metadata.csv", header = TRUE)

## Make world map data
world_map <- map_data("world")

# Set colors for the points
red_palette1 <- brewer.pal(n = 9, name = "YlOrRd")
red_palette2 <- brewer.pal(n=9, name = "YlOrBr")
blue_palette <- brewer.pal(n = 9, name = "Blues")
green_palette1 <- brewer.pal(n = 9, name = "BuGn")
green_palette2 <- brewer.pal(n=9, name = "YlGn")

scale_colour_region <- 
  c('USA' = red_palette2[7],
    'CENAM' = "blueviolet",
    'EUR' = green_palette1[7],
    'ASIA' = "violetred1",
    'AUS' = blue_palette[5])
# chose one random location/point for each region

# Get median fst for each comparison
fst_median <- data %>%
  filter(data_type == "Fst") %>%
  group_by(comparison) %>%
  summarise(FST_MEDIAN = median(value, na.rm = TRUE))
# Save as csv file and then add required metadata such as population longitudes/latitudes
write.csv(fst_median, "fst_median.csv", row.names = TRUE)

# read new csv file with map metadata and median fst between groups
fst_map <- read.csv("fst_map.csv", header = TRUE)

# Fst by region on a map
plot_3_fst <- ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = "grey90") +
  geom_segment(data=fst_map, aes(x=POP1_LONG, y=POP1_LAT, xend=POP2_LONG, yend=POP2_LAT, color = FST_MEDIAN), size = 1) +
  scale_color_gradient(low = "gold", high = "red3", name = expression(F[ST])) +
  geom_point(data = location, aes(x = Longitude, y = Latitude, fill = Region), size = 6, shape = 21) +
  scale_fill_manual(values = scale_colour_region, limits = c("USA", "CENAM", "EUR", "ASIA", "AUS"), name = "Region") +
  geom_text_repel(data = location, aes(x = Longitude, y = Latitude, label = Region), size = 4, fontface = "bold", nudge_y = 7) +
  theme_void() +
  theme(
    legend.position = c(0.1,0.3),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)) +
  guides(fill = "none") +
  coord_fixed()
plot_3_fst

# Save plot
ggsave("fst_map_region_dog.tif", bg = "white", height=5, width=10, dpi = 300)
ggsave("fst_map_region_dog.png", bg = "white", height=5, width=10, dpi = 300)
ggsave("fst_map_region_dog.pdf", bg = "white", height=5, width=10, dpi = 300)


# Fst by region on a heatmap
fst_heatmap <- read.csv("fst_heatmap.csv", header = TRUE)

plot_4_fst <- ggplot(fst_heatmap, aes(x = POP1, y=POP2, fill = FST_MEDIAN)) + 
  geom_tile(color = "white", lwd = 1.5, linetype = 1) +
  coord_fixed() +
  scale_fill_material("indigo") +
  guides(fill = guide_colourbar(title = expression(F[ST]), barwidth = 1, barheight = 5)) +
  theme_minimal() +
  labs(x = "Population", y = "Population") +
  theme(axis.text = element_text(size = 18),
    axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title = element_text(size = 24),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank())
plot_4_fst

ggsave("fst_heatmap_region_dog.tif", bg = "white", height=5, width = 8, dpi = 300)
ggsave("fst_heatmap_region_dog.png", bg = "white", height=5, width = 8, dpi = 300)
ggsave("fst_heatmap_region_dog.pdf", bg = "white", height=5, width = 8, dpi = 300)
```



### Pixy by host (all samples)

```R
# Pixy for Pi, Fst and Dxy by HOST

library(tidyverse)
library(ggsci)
library(ggpubr)
library(patchwork)
library(ggridges)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(purrr)
library(stats)
library(graphics)
library(grDevices)
library(utils)
library(datasets)
library(methods)
library(base)
require(maps)
library(mapdata)
library(readxl)
library(ozmaps) 
library(grid)
library(gridExtra)
library(ggrepel)
library(ggnewscale)
library(reshape)
library(wesanderson)
library(scales)

setwd("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/R2_Sep25/pixy/host")

##### Nucleotide diversity (Pi) #####

# get nucleotide diversity (pi) data from pixy output
pi_data <- read.table("input/host_pi.txt", header=T)
pi_data$chromosome <- str_remove(pi_data$chromosome, 'dirofilaria_immitis_')

# subset
pi_data_Ferret <- pi_data %>%
  filter(pop=="Ferret")
pi_data_Golden_jackal <- pi_data %>%
  filter(pop=="Golden_jackal")
pi_data_Wildcat <- pi_data %>%
  filter(pop=="Wildcat")
pi_data_Leopard <- pi_data %>%
  filter(pop=="Leopard")
pi_data_Cat <- pi_data %>%
  filter(pop=="Cat")
pi_data_Fox <- pi_data %>%
  filter(pop=="Fox")
pi_data_Dog <- pi_data %>%
  filter(pop=="Dog")

# filter pi data to remove small scaffolds not in the linkage groups, and to number the rows per group to help with plotting
pi_data <- pi_data %>%
  group_by(pop) %>%
  mutate(position = 1:n())

# get position for vertical lines used in plot to delineate the linkage groups
pi_data %>%
  group_by(chromosome) %>%
  summarise(max = max(position, na.rm = TRUE))

'
# A tibble: 5 × 2
  chromosome   max
  <chr>      <int>
1 chr1         440
2 chr2         592
3 chr3         744
4 chr4         885
5 chrX         283
'

#Let's add the chr type variable
pi_data <- pi_data %>%
  mutate(chr_type = ifelse(str_detect(chromosome, "X"), "sexchr", "autosome"))

# calculate the median Pi and check the ratio of sex-to-autosome diversity. Should be about 0.75, as D. immitis is XX/XY
pi_data_sex_median <- pi_data %>%
  group_by(chr_type) %>%
  summarise(median = median(avg_pi, na.rm = TRUE))
pi_data_sex_median

'
# A tibble: 2 × 2
  chr_type   median
  <chr>       <dbl>
1 autosome 0.000524
2 sexchr   0.000183
'
# 0.000267 / 0.000561 = 0.4759358 (way off 0.75)

pi_data_pop_sex_median <- pi_data %>%
  group_by(pop, chr_type) %>%
  summarise(median = median(avg_pi, na.rm = TRUE))
pi_data_pop_sex_median

'
# A tibble: 14 × 3
# Groups:   pop [7]
   pop           chr_type    median
   <chr>         <chr>        <dbl>
 1 Cat           autosome 0.000721 
 2 Cat           sexchr   0.000411 
 3 Dog           autosome 0.000761 
 4 Dog           sexchr   0.000441 
 5 Ferret        autosome 0.00127  
 6 Ferret        sexchr   0.000611 
 7 Fox           autosome 0.000676 
 8 Fox           sexchr   0.000373 
 9 Golden_jackal autosome 0.000204 
10 Golden_jackal sexchr   0.0000214
11 Leopard       autosome 0.000146 
12 Leopard       sexchr   0.0000414
13 Wildcat       autosome 0.000209 
14 Wildcat       sexchr   0.0000432
'

# plot 1 - genome wide plots per population
plot_1_pi <- ggplot(pi_data, aes(position*100000, avg_pi, col=chromosome)) +
  geom_point(size=0.8) +
  facet_grid(pop~.) +
  scale_color_npg(guide = guide_legend(override.aes = list(size = 4))) +
  theme_bw() +
  theme(legend.position="top", 
        legend.title = element_blank(),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 12, color = "black"),
        strip.text = element_text(size = 10),
        legend.text = element_text(size = 16),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 2)) +
  labs(x="Genomic Position", y="Nucleotide Diversity (Pi)")
plot_1_pi


# plot 2 - density plots of pi per group
plot_2_pi <- ggplot(pi_data, aes(avg_pi, chr_type, fill=chr_type), guide="none") +
  geom_density_ridges(quantile_lines=TRUE, quantile_fun=function(x,...)median(x), size=0.5) +
  theme_bw() + theme(legend.position = "none", axis.text.y=element_blank()) +
  facet_grid(pop~.) +
  xlim(0, 0.005) +
  scale_x_continuous(breaks = c(0, 0.0025, 0.005)) +
  scale_fill_npg() +
  theme(legend.position="top",
        legend.title = element_blank(),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 12, color = "black"),
        strip.text = element_text(size = 10),
        legend.text = element_text(size = 16),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 2)) +
  labs(x="Nucleotide Diversity (Pi)", y="Density")
plot_2_pi

# combine plots
plot_1_pi + plot_2_pi +  plot_layout(widths = c(5, 1))
ggsave("genomewide_and_density_Pi_host.tif", width=14, height=8, dpi = 300)
ggsave("genomewide_and_density_Pi_host.png", width=14, height=8, dpi = 300)
ggsave("genomewide_and_density_Pi_host.pdf", width=14, height=8, dpi = 300)

#Now a boxplot of the pi value per population
pi_data_renamed <- pi_data
pi_data_renamed$pop <- gsub("Golden_jackal", "Jackal", pi_data_renamed$pop)

boxplot_pi <- ggplot(pi_data_renamed, aes(pop, avg_pi, col=pop)) +
  geom_jitter(size = 1, alpha = 0.5) +
  geom_boxplot(fill=NA, col="black") +
  labs(x = "Population" , y = "Nucleotide diversity (Pi)", colour = "Population") +
  scale_color_npg() +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 18, color = "black"),
        #panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 2)
        )

boxplot_pi

ggsave("boxplot_Pi_host.tif", width=6, height=6, dpi = 300)
ggsave("boxplot_Pi_host.png", width=6, height=6, dpi = 300)
ggsave("boxplot_Pi_host.pdf", width=6, height=6, dpi = 300)
# higher pi value = more diverse population
# lower pi value = less diverse population
## Looks like Ferret has highest diversity. Jackal, leopard and wildcat have the lowest (but keep in mind very small sample size).


#let's check the statistical significance of the differences between pi values
#Let's explore the normality
pi_data %>% filter(chromosome != 'chrX') %>%
  ggplot(., aes(x = avg_pi, colour = pop)) +
  geom_histogram() +
  theme_bw() +
  facet_grid(. ~ pop)

#Confirm the absence of normality by shapiro
## If p > 0.05 the data is normal, if p < 0.05 the data is NOT normal
pi_data_Cat <- pi_data_Cat %>% filter(chromosome != 'chrX')
Cat_shapiro <- shapiro.test(pi_data_Cat$avg_pi)
print(Cat_shapiro)
# W = 0.90509, p-value < 2.2e-16

pi_data_Dog <- pi_data_Dog %>% filter(chromosome != 'chrX')
Dog_shapiro <- shapiro.test(pi_data_Dog$avg_pi)
print(Dog_shapiro)
# W = 0.91281, p-value < 2.2e-16

pi_data_Ferret <- pi_data_Ferret %>% filter(chromosome != 'chrX')
Ferret_shapiro <- shapiro.test(pi_data_Ferret$avg_pi)
print(Ferret_shapiro)
# W = 0.92422, p-value < 2.2e-16

pi_data_Fox <- pi_data_Fox %>% filter(chromosome != 'chrX')
Fox_shapiro <- shapiro.test(pi_data_Fox$avg_pi)
print(Fox_shapiro)
# W = 0.92589, p-value < 2.2e-16

pi_data_Golden_jackal <- pi_data_Golden_jackal %>% filter(chromosome != 'chrX')
Golden_jackal_shapiro <- shapiro.test(pi_data_Golden_jackal$avg_pi)
print(Golden_jackal_shapiro)
# W = 0.79487, p-value < 2.2e-16

pi_data_Leopard <- pi_data_Leopard %>% filter(chromosome != 'chrX')
Leopard_shapiro <- shapiro.test(pi_data_Leopard$avg_pi)
print(Leopard_shapiro)
# W = 0.74904, p-value < 2.2e-16

pi_data_Wildcat <- pi_data_Wildcat %>% filter(chromosome != 'chrX')
Wildcat_shapiro <- shapiro.test(pi_data_Wildcat$avg_pi)
print(Wildcat_shapiro)
# W = 0.7881, p-value < 2.2e-16

#Wilcoxson test for all pairs of pops
## This tests whether the mean values of 2 dependent groups differ significantly from each other
## If p < 0.05 they are sig different
wilcox.test(pi_data_Cat$avg_pi, pi_data_Dog$avg_pi)
# W = 170525, p-value = 0.07676
wilcox.test(pi_data_Cat$avg_pi, pi_data_Ferret$avg_pi)
# W = 119004, p-value < 2.2e-16
wilcox.test(pi_data_Cat$avg_pi, pi_data_Fox$avg_pi)
# W = 192171, p-value = 0.06903
wilcox.test(pi_data_Cat$avg_pi, pi_data_Golden_jackal$avg_pi)
# W = 269116, p-value < 2.2e-16
wilcox.test(pi_data_Cat$avg_pi, pi_data_Leopard$avg_pi)
# W = 282078, p-value < 2.2e-16
wilcox.test(pi_data_Cat$avg_pi, pi_data_Wildcat$avg_pi)
# W = 267850, p-value < 2.2e-16
wilcox.test(pi_data_Dog$avg_pi, pi_data_Ferret$avg_pi)
# W = 125792, p-value < 2.2e-16
wilcox.test(pi_data_Dog$avg_pi, pi_data_Fox$avg_pi)
# W = 202856, p-value = 0.0003314
wilcox.test(pi_data_Dog$avg_pi, pi_data_Golden_jackal$avg_pi)
# W = 274823, p-value < 2.2e-16
wilcox.test(pi_data_Dog$avg_pi, pi_data_Leopard$avg_pi)
# W = 288228, p-value < 2.2e-16
wilcox.test(pi_data_Dog$avg_pi, pi_data_Wildcat$avg_pi)
# W = 273142, p-value < 2.2e-16
wilcox.test(pi_data_Ferret$avg_pi, pi_data_Fox$avg_pi)
# W = 252126, p-value < 2.2e-16
wilcox.test(pi_data_Ferret$avg_pi, pi_data_Golden_jackal$avg_pi)
# W = 302341, p-value < 2.2e-16
wilcox.test(pi_data_Ferret$avg_pi, pi_data_Leopard$avg_pi)
# W = 311633, p-value < 2.2e-16
wilcox.test(pi_data_Ferret$avg_pi, pi_data_Wildcat$avg_pi)
# W = 301603, p-value < 2.2e-16
wilcox.test(pi_data_Fox$avg_pi, pi_data_Golden_jackal$avg_pi)
# W = 262138, p-value < 2.2e-16
wilcox.test(pi_data_Fox$avg_pi, pi_data_Leopard$avg_pi)
# W = 275757, p-value < 2.2e-16
wilcox.test(pi_data_Fox$avg_pi, pi_data_Wildcat$avg_pi)
# W = 261208, p-value < 2.2e-16
wilcox.test(pi_data_Golden_jackal$avg_pi, pi_data_Leopard$avg_pi)
# W = 193552, p-value = 0.03998
wilcox.test(pi_data_Golden_jackal$avg_pi, pi_data_Wildcat$avg_pi)
# W = 186628, p-value = 0.3665
wilcox.test(pi_data_Leopard$avg_pi, pi_data_Wildcat$avg_pi)
# W = 174546, p-value = 0.2672


#let's generate some dataframes for the estatistics
# subset
pi_data_Cat <- pi_data %>%
  filter(pop=="Cat") %>%
  filter(chromosome != 'chrX')

pi_data_Dog <- pi_data %>%
  filter(pop=="Dog") %>%
  filter(chromosome != 'chrX') 

pi_data_Ferret <- pi_data %>%
  filter(pop=="Ferret") %>%
  filter(chromosome != 'chrX')

pi_data_Fox <- pi_data %>%
  filter(pop=="Fox") %>%
  filter(chromosome != 'chrX')

pi_data_Golden_jackal <- pi_data %>%
  filter(pop=="Golden_jackal") %>%
  filter(chromosome != 'chrX')

pi_data_Leopard <- pi_data %>%
  filter(pop=="Leopard") %>%
  filter(chromosome != 'chrX')

pi_data_Wildcat <- pi_data %>%
  filter(pop=="Wildcat") %>%
  filter(chromosome != 'chrX')


#### Dxy and Fst ####

# load data
dxy_data <- read.table("input/host_dxy.txt", header=T)
fst_data <- read.table("input/host_fst.txt", header=T)

# add some columns
dxy_data$data_type <- "Dxy"
dxy_data <- mutate(dxy_data,
                   comparison = paste(pop1, pop2, sep = '_v_'))
fst_data$data_type <- "Fst"
fst_data <- mutate(fst_data,
                   comparison = paste(pop1, pop2, sep = '_v_'))

# subset and merge dataframes
dxy_data_sub <- dxy_data %>% select(comparison,data_type,chromosome,window_pos_1,window_pos_2,avg_dxy)
colnames(dxy_data_sub) <- c("comparison","data_type","chromosome","window_pos_1","window_pos_2","value")

fst_data_sub <- fst_data %>% select(comparison,data_type,chromosome,window_pos_1,window_pos_2,avg_wc_fst)
colnames(fst_data_sub) <- c("comparison","data_type","chromosome","window_pos_1","window_pos_2","value")

# cheeky fix to get matching rows from both datasets
tmp_dxy <- semi_join(dxy_data_sub, fst_data_sub, by=c("comparison", "chromosome", "window_pos_1", "window_pos_2"))
tmp_fst <- semi_join(fst_data_sub, dxy_data_sub,  by=c("comparison", "chromosome", "window_pos_1", "window_pos_2"))

# join the datasets together to create a single dataframe
data <- full_join(tmp_dxy, tmp_fst)
data$chromosome <- str_remove(data$chromosome, 'dirofilaria_immitis_')

# add numbering to help with plotting
data <- data %>%
  group_by(comparison, data_type) %>%
  mutate(position = 1:n())

# add sex chromosome information
data <- data %>%
  mutate(chr_type = ifelse(str_detect(chromosome, "X"), "sexchr", "autosome"))


# summarise median data for Fst and Dxy
data %>%
  group_by(comparison,data_type) %>%
  summarise(median = median(value, na.rm = TRUE)) %>%
  print(n = 42)
'
# A tibble: 42 × 3
# Groups:   comparison [21]
   comparison              data_type    median
   <chr>                   <chr>         <dbl>
 1 Cat_v_Ferret            Dxy        0.000675
 2 Cat_v_Ferret            Fst       -0.114   
 3 Cat_v_Golden_jackal     Dxy        0.000691
 4 Cat_v_Golden_jackal     Fst        0.0791  
 5 Cat_v_Leopard           Dxy        0.000671
 6 Cat_v_Leopard           Fst        0.0763  
 7 Cat_v_Wildcat           Dxy        0.000684
 8 Cat_v_Wildcat           Fst        0.0758  
 9 Dog_v_Cat               Dxy        0.000645
10 Dog_v_Cat               Fst        0.0213  
11 Dog_v_Ferret            Dxy        0.000725
12 Dog_v_Ferret            Fst       -0.0616  
13 Dog_v_Fox               Dxy        0.000628
14 Dog_v_Fox               Fst       -0.00231 
15 Dog_v_Golden_jackal     Dxy        0.000645
16 Dog_v_Golden_jackal     Fst        0.0362  
17 Dog_v_Leopard           Dxy        0.000649
18 Dog_v_Leopard           Fst        0.0437  
19 Dog_v_Wildcat           Dxy        0.000643
20 Dog_v_Wildcat           Fst        0.0371  
21 Fox_v_Cat               Dxy        0.000642
22 Fox_v_Cat               Fst        0.0504  
23 Fox_v_Ferret            Dxy        0.000727
24 Fox_v_Ferret            Fst        0.0431  
25 Fox_v_Golden_jackal     Dxy        0.000593
26 Fox_v_Golden_jackal     Fst        0.104   
27 Fox_v_Leopard           Dxy        0.000607
28 Fox_v_Leopard           Fst        0.110   
29 Fox_v_Wildcat           Dxy        0.000577
30 Fox_v_Wildcat           Fst        0.100   
31 Golden_jackal_v_Ferret  Dxy        0.000755
32 Golden_jackal_v_Ferret  Fst        0       
33 Leopard_v_Ferret        Dxy        0.000754
34 Leopard_v_Ferret        Fst        0       
35 Leopard_v_Golden_jackal Dxy        0.000292
36 Leopard_v_Golden_jackal Fst        0       
37 Leopard_v_Wildcat       Dxy        0.000285
38 Leopard_v_Wildcat       Fst        0       
39 Wildcat_v_Ferret        Dxy        0.000744
40 Wildcat_v_Ferret        Fst        0       
41 Wildcat_v_Golden_jackal Dxy        0.000293
42 Wildcat_v_Golden_jackal Fst        0    
'


#Plotting Dxy

# plot 1 - genome wide plots per comparison
plot_1_dxy <- data %>%
  filter(data_type =='Dxy')%>%
  ggplot(., aes(position*100000, value, col=chromosome)) +
  geom_point(size=1.5) +
  facet_grid(comparison~.) +
  scale_color_npg(guide = guide_legend(override.aes = list(size = 4))) +
  theme_bw() +
  theme(legend.position="top", 
        legend.title = element_blank(),
        axis.title = element_text(size = 32),
        axis.text = element_text(size = 16, color = "black"),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 20),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 2)) +
  labs(x="Genomic Position", y="Dxy")
plot_1_dxy

# plot 2 - density plots of dxy per group
plot_2_dxy <- data %>%
  filter(data_type=='Dxy')%>%
  ggplot(aes(value, chr_type, fill=chr_type), guide="none") +
  geom_density_ridges(quantile_lines=TRUE, quantile_fun=function(x,...)median(x), size=0.5) +
  theme_bw() + theme(legend.position = "top", 
                     legend.title = element_blank(),
                     axis.text.y=element_blank(),
                     axis.title = element_text(size = 32),
                     axis.text = element_text(size = 16, color = "black"),
                     strip.text = element_text(size = 8),
                     legend.text = element_text(size = 20),
                     panel.grid = element_blank(),
                     panel.border = element_rect(color = "black", fill = NA, linewidth = 2)) +
  facet_grid(comparison~.) +
  scale_fill_npg() +
  labs(x="Dxy", y="Density")
plot_2_dxy

# combine plots
dxy_plot <- plot_1_dxy + plot_2_dxy +  plot_layout(widths = c(5, 1))
ggsave("genomewide_and_density_dxy_host.tif", width=20, height=30, dpi = 300)
ggsave("genomewide_and_density_dxy_host.png", width=20, height=30, dpi = 300)
ggsave("genomewide_and_density_dxy_host.pdf", width=20, height=30, dpi = 300)
# higher Dxy = more divergent

#some additional plots
boxplot_dxy <- data %>%
  filter(data_type =='Dxy') %>% 
  ggplot(., aes(comparison, value, col=comparison)) +
  geom_jitter(size = 1, alpha = 0.5) +
  geom_boxplot(fill=NA, col="black") +
  labs(x = "Host" , y = "Dxy", colour = "Population") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 12, color = "black"),
        strip.text = element_text(size = 18),
        legend.text = element_text(size = 16),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 2)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
boxplot_dxy
ggsave("boxplot_dxy_host.tif", width=20, height = 6, dpi = 300)
ggsave("boxplot_dxy_host.png", width=20, height = 6, dpi = 300)
ggsave("boxplot_dxy_host.pdf", width=20, height = 6, dpi = 300)

density_dxy <- data %>%
  filter(data_type =='Dxy') %>% 
  ggplot(., aes(x=value, group = comparison, fill=comparison)) +
  geom_density(adjust=1.5, alpha=.4) +
  labs(x = "Dxy" , y = "Density", colour = "Comparison") +
  theme_bw() +
  theme(legend.position = c(0.73,0.77),
        legend.title = element_blank(),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 12, color = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 2))
density_dxy
ggsave("density_dxy_host.tif", width = 8, height = 8, dpi = 300)
ggsave("density_dxy_host.png", width = 8, height = 8, dpi = 300)
ggsave("density_dxy_host.pdf", width = 8, height = 8, dpi = 300)


#Trying to plot lines
#getting chr positions
chr <- c('X', '1', '2', '3', '4')

for (i in chr) {
  x <- data %>%
    filter(chromosome == paste0('chr', i))
  print(quantile(x$position))
}

data$comparison <- str_replace_all(data$comparison, '_', ' ')
data$comparison <- str_replace_all(data$comparison, 'v', 'vs')

dxy_lineplot <- data %>%
  filter(data_type=='Dxy')%>%
  ggplot(., aes(position*100000, value, col=comparison)) +
  geom_point(size=0.2, alpha = 0.1) +
  geom_vline(xintercept=c(281*100000,
                          438*100000,
                          590*100000,
                          740*100000),
             size=1, linetype="dashed", col='grey41')+
  annotate(geom="text", x=140*100000, y=0.004, label="ChrX", size = 6)+
  annotate(geom="text", x=359.5*100000, y=0.004, label="Chr1", size = 6)+
  annotate(geom="text", x=520*100000, y=0.004, label="Chr2", size = 6)+
  annotate(geom="text", x=665*100000, y=0.004, label="Chr3", size = 6)+
  annotate(geom="text", x=810.5*100000, y=0.004, label="Chr4", size = 6)+
  geom_line() +
  theme_bw() +
  theme(legend.position="top", 
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 14)) +
  labs(x="Genomic Position", y="Dxy") +
  guides(col = guide_legend(override.aes = list(linewidth = 2)))
dxy_lineplot
ggsave("lineplot_dxy_host.tif", width=14, height=8, dpi = 300)
ggsave("lineplot_dxy_host.png", width=14, height=8, dpi = 300)
ggsave("lineplot_dxy_host.pdf", width=14, height=8, dpi = 300)

# Dxy by host on a heatmap
# Get median dxy for each comparison
dxy_median <- data %>%
  filter(data_type == "Dxy") %>%
  group_by(comparison) %>%
  summarise(DXY_MEDIAN = median(value, na.rm = TRUE))

# Save as csv file
write.csv(dxy_median, "dxy_median.csv", row.names = TRUE)

dxy_heatmap <- read.csv("dxy_heatmap.csv", header = TRUE)
dxy_heatmap <- as.data.frame(dxy_heatmap)

plot_3_dxy <- ggplot(dxy_heatmap, aes(x = POP1, y=POP2, fill = DXY_MEDIAN)) + 
  geom_tile(color = "white", lwd = 1.5, linetype = 1) +
  coord_fixed() +
  scale_fill_material("indigo", labels = scales::label_number()) +
  guides(fill = guide_colourbar(title = expression(D[XY]), barwidth = 1, barheight = 5)) +
  theme_minimal() +
  labs(x = "Host", y = "Host") +
  theme(axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 10),
        panel.grid.major = element_blank())
plot_3_dxy

ggsave("dxy_heatmap_host.tif", bg = "white", height=5, width = 8, dpi = 300)
ggsave("dxy_heatmap_host.png", bg = "white", height=5, width = 8, dpi = 300)
ggsave("dxy_heatmap_host.pdf", bg = "white", height=5, width = 8, dpi = 300)

# Plotting Fst

# plot 1 - genome wide plots per population
plot_1_fst <- data %>%
  filter(data_type=='Fst')%>%
  ggplot(., aes(position*100000, value, col=chromosome)) +
  geom_point(size=1) +
  facet_grid(comparison~.) +
  scale_color_npg(guide = guide_legend(override.aes = list(size = 4))) +
  theme_bw() +
  theme(legend.position="top", 
        legend.title = element_blank(),
        axis.title = element_text(size = 32),
        axis.text = element_text(size = 16, color = "black"),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 20),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 2)) +
  labs(x="Genomic Position", y="Fst")
plot_1_fst
# Removed 338 rows containing missing values or values outside the scale range (`geom_point()`). 

# plot 2 - density plots of Fst per group
plot_2_fst <- data %>%
  filter(data_type=='Fst')%>%
  ggplot(aes(value, chr_type, fill=chr_type), guide="none") +
  geom_density_ridges(quantile_lines=TRUE, quantile_fun=function(x,...)median(x), size=0.5) +
  theme_bw() + 
  theme(legend.position = "none", 
        axis.text.y=element_blank(),
        axis.title = element_text(size = 32),
        axis.text = element_text(size = 16, color = "black"),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 20),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 2)) +
  facet_grid(comparison~.) +
  scale_fill_npg() +
  labs(x="Fst", y="Density")
plot_2_fst
# Removed 338 rows containing non-finite outside the scale range (`stat_density_ridges()`). 

# combine plots
plot_1_fst + plot_2_fst +  plot_layout(widths = c(5, 1))
## Fst values should only range between 0 and 1. Set anything <0 to 0 and re-do the plots.

# Data cleaning
data_fst_clean <- data %>% filter(data_type == 'Fst') %>% mutate(value = ifelse(value < 0, 0, value)) # this takes my original data, filters it to only include rows where the data_type is Fst. Then mutate function is used to create/modify columns. It modifies the value column, checks if the value is less than 0. If the value is <0, it replaces it with 0. If the value is >0, it keeps the original value.

# Re-do plots with cleaned data
# plot 1 - genome wide plots per population
plot_1_fst_clean <- data_fst_clean %>%
  filter(data_type=='Fst')%>%
  ggplot(., aes(position*100000, value, col=chromosome)) +
  geom_point(size=1) +
  facet_grid(comparison~.) +
  scale_color_tron() +
  theme_bw() +
  theme(legend.position="top", 
        legend.title = element_blank(),
        axis.title = element_text(size = 32),
        axis.text = element_text(size = 16, color = "black"),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 20),
        panel.grid = element_blank(),
panel.border = element_rect(color = "black", fill = NA, linewidth = 2)) +
  labs(x="Genomic Position", y="Fst")
plot_1_fst_clean


# plot 2 - density plots of Fst per group
plot_2_fst_clean <- data_fst_clean %>%
  filter(data_type=='Fst')%>%
  ggplot(aes(value, chr_type, fill=chr_type), guide="none") +
  geom_density_ridges(quantile_lines=TRUE, quantile_fun=function(x,...)median(x), size=0.5) +
  theme_bw() + theme(legend.position = "none", axis.text.y=element_blank(), panel.grid = element_blank(), panel.border = element_rect(color = "black", fill = NA, linewidth = 2)) +
  facet_grid(comparison~.) +
  scale_fill_npg() +
  labs(x="Fst", y="Density")
plot_2_fst_clean

# combine plots
plot_1_fst_clean + plot_2_fst_clean +  plot_layout(widths = c(5, 1))
ggsave("genomewide_and_density_fst_host.tif", width=20, height=30, dpi = 300)
ggsave("genomewide_and_density_fst_host.png", width=20, height=30, dpi = 300)
ggsave("genomewide_and_density_fst_host.pdf", width=20, height=30, dpi = 300)
# Fst closer to 1 = more genetically distinct
# Fst closer to 0 = more genetically similar
# Doesn't look too great because some hosts only have 1 sample, and there's a bunch of NA data and negative values I had to clean.

# Fst by host on a heatmap
# Get median fst for each comparison
fst_median <- data_fst_clean %>%
  filter(data_type == "Fst") %>%
  group_by(comparison) %>%
  summarise(FST_MEDIAN = median(value, na.rm = TRUE))
# Save as csv file
write.csv(fst_median, "fst_median.csv", row.names = TRUE)

fst_heatmap <- read.csv("fst_heatmap.csv", header = TRUE)

plot_4_fst <- ggplot(fst_heatmap, aes(x = POP1, y=POP2, fill = FST_MEDIAN)) + 
  geom_tile(color = "white", lwd = 1.5, linetype = 1) +
  coord_fixed() +
  scale_fill_material("indigo") +
  guides(fill = guide_colourbar(title = expression(F[ST]), barwidth = 1, barheight = 5)) +
  theme_minimal() +
  labs(x = "Host", y = "Host") +
  theme(axis.text = element_text(size = 12),
    axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 10),
        panel.grid.major = element_blank())
plot_4_fst

ggsave("fst_heatmap_host.tif", bg = "white", height=5, width = 8, dpi = 300)
ggsave("fst_heatmap_host.png", bg = "white", height=5, width = 8, dpi = 300)
ggsave("fst_heatmap_host.pdf", bg = "white", height=5, width = 8, dpi = 300)
```



### Pixy for AUS samples (dog)

```R
# Batch 4 - Pixy for Pi, Fst and Dxy for Australian dog cohort only

library(tidyverse)
library(ggsci)
library(ggpubr)
library(patchwork)
library(ggridges)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(purrr)
library(stats)
library(graphics)
library(grDevices)
library(utils)
library(datasets)
library(methods)
library(base)
require(maps)
library(mapdata)
library(readxl)
library(ozmaps) 
library(grid)
library(gridExtra)
library(ggrepel)
library(ggnewscale)
library(reshape)
library(wesanderson)

setwd("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/R2_Sep25/pixy/aus_dog")

##### Nucleotide diversity (Pi) #####

# get nucleotide diversity (pi) data from pixy output
pi_data <- read.table("input/aus_dog_pi.txt", header=T)
pi_data$chromosome <- str_remove(pi_data$chromosome, 'dirofilaria_immitis_')

# subset
pi_data_LHR <- pi_data %>%
  filter(pop=="Lockhart_River")
pi_data_CNS <- pi_data %>%
  filter(pop=="Cairns")
pi_data_TSV <- pi_data %>%
  filter(pop=="Townsville")
pi_data_ROK <- pi_data %>%
  filter(pop=="Rockhampton")
pi_data_BNE <- pi_data %>%
  filter(pop=="Brisbane")
pi_data_SYD <- pi_data %>%
  filter(pop=="Sydney")

# filter pi data to remove small scaffolds not in the linkage groups, and to number the rows per group to help with plotting
pi_data <- pi_data %>%
  group_by(pop) %>%
  mutate(position = 1:n())

#Let's add the chr type variable
pi_data <- pi_data %>%
  mutate(chr_type = ifelse(str_detect(chromosome, "X"), "sexchr", "autosome"))

# calculate the median Pi and check the ratio of sex-to-autosome diversity. Should be about 0.75, as D. immitis is XX/XY
pi_data_sex_median <- pi_data %>%
  group_by(chr_type) %>%
  summarise(median = median(avg_pi, na.rm = TRUE))
pi_data_sex_median

'
# A tibble: 2 × 2
  chr_type   median
  <chr>       <dbl>
1 autosome 0.000517
2 sexchr   0.000209
'
# 0.000209 / 0.000517 = 0.4042553 (way off 0.75)

pi_data_pop_sex_median <- pi_data %>%
  group_by(pop, chr_type) %>%
  summarise(median = median(avg_pi, na.rm = TRUE))
pi_data_pop_sex_median

'
# A tibble: 12 × 3
# Groups:   pop [6]
   pop            chr_type    median
   <chr>          <chr>        <dbl>
 1 Brisbane       autosome 0.000480 
 2 Brisbane       sexchr   0.000207 
 3 Cairns         autosome 0.000474 
 4 Cairns         sexchr   0.0000985
 5 Lockhart_River autosome 0.000887 
 6 Lockhart_River sexchr   0.000391 
 7 Rockhampton    autosome 0.000260 
 8 Rockhampton    sexchr   0.0000331
 9 Sydney         autosome 0.000419 
10 Sydney         sexchr   0.000229 
11 Townsville     autosome 0.000689 
12 Townsville     sexchr   0.000306 
'

# plot 1 - genome wide plots per population
plot_1_pi <- ggplot(pi_data, aes(position*100000, avg_pi, col=chromosome)) +
  geom_point(size=0.8) +
  facet_grid(pop~.) +
  scale_color_npg(guide = guide_legend(override.aes = list(size = 4))) +
  theme_bw() +
  theme(legend.position="top", 
        legend.title = element_blank(),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 12, color = "black"),
        strip.text = element_text(size = 10),
        legend.text = element_text(size = 16),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 2)) +
  labs(x="Genomic Position", y="Nucleotide Diversity (Pi)")
plot_1_pi


# plot 2 - density plots of pi per group
plot_2_pi <- ggplot(pi_data, aes(avg_pi, chr_type, fill=chr_type), guide="none") +
  geom_density_ridges(quantile_lines=TRUE, quantile_fun=function(x,...)median(x), size=0.5) +
  theme_bw() + theme(legend.position = "none", axis.text.y=element_blank()) +
  facet_grid(pop~.) +
  xlim(0, 0.005) +
  scale_x_continuous(breaks = c(0, 0.0025, 0.005)) +
  scale_fill_npg() +
  theme(legend.position="top",
        legend.title = element_blank(),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 12, color = "black"),
        strip.text = element_text(size = 10),
        legend.text = element_text(size = 16),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 2)) +
  labs(x="Nucleotide Diversity (Pi)", y="Density")
plot_2_pi

# combine plots
plot_1_pi + plot_2_pi +  plot_layout(widths = c(5, 1))
ggsave("genomewide_and_density_Pi_aus_dog.tif", width=14, height=8, dpi = 300)
ggsave("genomewide_and_density_Pi_aus_dog.png", width=14, height=8, dpi = 300)
ggsave("genomewide_and_density_Pi_aus_dog.pdf", width=14, height=8, dpi = 300)

#Now a boxplot of the pi value per population
pi_data$pop <- factor(pi_data$pop, 
                      levels = c('Lockhart_River', 'Cairns', 'Townsville', 'Rockhampton', 
                                 'Brisbane', 'Sydney'))

blue_palette <- brewer.pal(n = 9, name = "Blues")

scale_colour_AUS <- function(...){
  ggplot2:::manual_scale(
    'colour', 
    values = setNames(
      c(blue_palette[9],
        blue_palette[8],
        blue_palette[7],
        blue_palette[6],
        blue_palette[5],
        blue_palette[4]),
    c("Lockhart_River", "Cairns", "Townsville", "Rockhampton", "Brisbane", "Sydney")), 
  ...
)
}    


boxplot_pi <- ggplot(pi_data, aes(pop, avg_pi, col=pop)) +
  geom_jitter(size = 1, alpha = 0.5) +
  geom_boxplot(fill=NA, col="black") +
  labs(x = "Population" , y = "Nucleotide diversity (Pi)", colour = "Population") +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 18, color = "black"),
        #panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 2)
        ) +
  scale_colour_AUS ()

boxplot_pi

ggsave("boxplot_Pi_aus_dog.tif", width=4, height=6, dpi = 300)
ggsave("boxplot_Pi_aus_dog.png", width=4, height=6, dpi = 300)
ggsave("boxplot_Pi_aus_dog.pdf", width=4, height=6, dpi = 300)
# higher pi value = more diverse population
# lower pi value = less diverse population
## Diversity is pretty similar between cities. Sydney and Brisbane are on par, Townsville is slightly higher.


#let's check the statistical significance of the differences between pi values
#Let's explore the normality
pi_data %>% filter(chromosome != 'chrX') %>%
  ggplot(., aes(x = avg_pi, colour = pop)) +
  geom_histogram() +
  theme_bw() +
  facet_grid(. ~ pop)

#Confirm the absence of normality by shapiro
## If p > 0.05 the data is normal, if p < 0.05 the data is NOT normal
pi_data_LHR <- pi_data_LHR %>% filter(chromosome != 'chrX')
LHR_shapiro <- shapiro.test(pi_data_LHR$avg_pi)
print(LHR_shapiro)
# W = 0.89779, p-value < 2.2e-16

pi_data_CNS <- pi_data_CNS %>% filter(chromosome != 'chrX')
CNS_shapiro <- shapiro.test(pi_data_CNS$avg_pi)
print(CNS_shapiro)
# W = 0.88412, p-value < 2.2e-16

pi_data_TSV <- pi_data_TSV %>% filter(chromosome != 'chrX')
TSV_shapiro <- shapiro.test(pi_data_TSV$avg_pi)
print(TSV_shapiro)
# W = 0.91773, p-value < 2.2e-16

pi_data_ROK <- pi_data_ROK %>% filter(chromosome != 'chrX')
ROK_shapiro <- shapiro.test(pi_data_ROK$avg_pi)
print(ROK_shapiro)
# W = 0.79065, p-value < 2.2e-16

pi_data_BNE <- pi_data_BNE %>% filter(chromosome != 'chrX')
BNE_shapiro <- shapiro.test(pi_data_BNE$avg_pi)
print(BNE_shapiro)
# W = 0.91712, p-value < 2.2e-16

pi_data_SYD <- pi_data_SYD %>% filter(chromosome != 'chrX')
SYD_shapiro <- shapiro.test(pi_data_SYD$avg_pi)
print(SYD_shapiro)
# W = 0.89084, p-value < 2.2e-16

#Wilcoxson test for all pairs of pops
## This tests whether the mean values of 2 dependent groups differ significantly from each other
## If p < 0.05 they are sig different
wilcox.test(pi_data_LHR$avg_pi, pi_data_CNS$avg_pi)
# W = 244000, p-value < 2.2e-16
wilcox.test(pi_data_LHR$avg_pi, pi_data_TSV$avg_pi)
# W = 211096, p-value = 7.222e-07
wilcox.test(pi_data_LHR$avg_pi, pi_data_ROK$avg_pi)
# W = 273059, p-value < 2.2e-16
wilcox.test(pi_data_LHR$avg_pi, pi_data_BNE$avg_pi)
# W = 250541, p-value < 2.2e-16
wilcox.test(pi_data_LHR$avg_pi, pi_data_SYD$avg_pi)
# W = 257315, p-value < 2.2e-16
wilcox.test(pi_data_CNS$avg_pi, pi_data_TSV$avg_pi)
# W = 137951, p-value = 7.52e-13
wilcox.test(pi_data_CNS$avg_pi, pi_data_ROK$avg_pi)
# W = 228283, p-value = 5.587e-15
wilcox.test(pi_data_CNS$avg_pi, pi_data_BNE$avg_pi)
# W = 186599, p-value = 0.3711
wilcox.test(pi_data_CNS$avg_pi, pi_data_SYD$avg_pi)
# W = 196079, p-value = 0.01366
wilcox.test(pi_data_TSV$avg_pi, pi_data_ROK$avg_pi)
# W = 255917, p-value < 2.2e-16
wilcox.test(pi_data_TSV$avg_pi, pi_data_BNE$avg_pi)
# W = 232458, p-value < 2.2e-16
wilcox.test(pi_data_TSV$avg_pi, pi_data_SYD$avg_pi)
# W = 241180, p-value < 2.2e-16
wilcox.test(pi_data_ROK$avg_pi, pi_data_BNE$avg_pi)
# W = 136737, p-value = 1.61e-13
wilcox.test(pi_data_ROK$avg_pi, pi_data_SYD$avg_pi)
# W = 143136, p-value = 2.68e-10
wilcox.test(pi_data_BNE$avg_pi, pi_data_SYD$avg_pi)
# W = 191409, p-value = 0.09068

## Quite a few are not sig different


#let's generate some dataframes for the estatistics
# subset
pi_data_LHR <- pi_data %>%
  filter(pop=="Lockhart_River") %>%
  filter(chromosome != 'chrX')

pi_data_CNS <- pi_data %>%
  filter(pop=="Cairns") %>%
  filter(chromosome != 'chrX') 

pi_data_TSV <- pi_data %>%
  filter(pop=="Townsville") %>%
  filter(chromosome != 'chrX')

pi_data_ROK <- pi_data %>%
  filter(pop=="Rockhampton") %>%
  filter(chromosome != 'chrX')

pi_data_BNE <- pi_data %>%
  filter(pop=="Brisbane") %>%
  filter(chromosome != 'chrX')

pi_data_SYD <- pi_data %>%
  filter(pop=="Sydney") %>%
  filter(chromosome != 'chrX')



#### Dxy and Fst ####

# load data
dxy_data <- read.table("input/aus_dog_dxy.txt", header=T)
fst_data <- read.table("input/aus_dog_fst.txt", header=T)

# add some columns
dxy_data$data_type <- "Dxy"
dxy_data <- mutate(dxy_data,
                   comparison = paste(pop1, pop2, sep = '_v_'))
fst_data$data_type <- "Fst"
fst_data <- mutate(fst_data,
                   comparison = paste(pop1, pop2, sep = '_v_'))

# subset and merge dataframes
dxy_data_sub <- dxy_data %>% select(comparison,data_type,chromosome,window_pos_1,window_pos_2,avg_dxy)
colnames(dxy_data_sub) <- c("comparison","data_type","chromosome","window_pos_1","window_pos_2","value")

fst_data_sub <- fst_data %>% select(comparison,data_type,chromosome,window_pos_1,window_pos_2,avg_wc_fst)
colnames(fst_data_sub) <- c("comparison","data_type","chromosome","window_pos_1","window_pos_2","value")

# cheeky fix to get matching rows from both datasets
tmp_dxy <- semi_join(dxy_data_sub, fst_data_sub, by=c("comparison", "chromosome", "window_pos_1", "window_pos_2"))
tmp_fst <- semi_join(fst_data_sub, dxy_data_sub,  by=c("comparison", "chromosome", "window_pos_1", "window_pos_2"))

# join the datasets together to create a single dataframe
data <- full_join(tmp_dxy, tmp_fst)
data$chromosome <- str_remove(data$chromosome, 'dirofilaria_immitis_')

# add numbering to help with plotting
data <- data %>%
  group_by(comparison, data_type) %>%
  mutate(position = 1:n())

# add sex chromosome information
data <- data %>%
  mutate(chr_type = ifelse(str_detect(chromosome, "X"), "sexchr", "autosome"))


# summarise median data for Fst and Dxy
data %>%
  group_by(comparison,data_type) %>%
  summarise(median = median(value, na.rm = TRUE)) %>%
  print(n = 30)
'
# A tibble: 30 × 3
# Groups:   comparison [15]
   comparison                   data_type    median
   <chr>                        <chr>         <dbl>
 1 Brisbane_v_Cairns            Dxy        0.000530
 2 Brisbane_v_Cairns            Fst        0.235   
 3 Brisbane_v_Lockhart_River    Dxy        0.000559
 4 Brisbane_v_Lockhart_River    Fst        0.0937  
 5 Brisbane_v_Rockhampton       Dxy        0.000493
 6 Brisbane_v_Rockhampton       Fst        0.127   
 7 Brisbane_v_Sydney            Dxy        0.000547
 8 Brisbane_v_Sydney            Fst        0.225   
 9 Brisbane_v_Townsville        Dxy        0.000543
10 Brisbane_v_Townsville        Fst        0.0672  
11 Cairns_v_Lockhart_River      Dxy        0.000493
12 Cairns_v_Lockhart_River      Fst        0.178   
13 Cairns_v_Rockhampton         Dxy        0.000521
14 Cairns_v_Rockhampton         Fst        0.392   
15 Cairns_v_Sydney              Dxy        0.000498
16 Cairns_v_Sydney              Fst        0.138   
17 Cairns_v_Townsville          Dxy        0.000548
18 Cairns_v_Townsville          Fst        0.0633  
19 Lockhart_River_v_Rockhampton Dxy        0.000536
20 Lockhart_River_v_Rockhampton Fst        0       
21 Lockhart_River_v_Sydney      Dxy        0.000550
22 Lockhart_River_v_Sydney      Fst       -0.0217  
23 Lockhart_River_v_Townsville  Dxy        0.000568
24 Lockhart_River_v_Townsville  Fst       -0.0328  
25 Rockhampton_v_Sydney         Dxy        0.000477
26 Rockhampton_v_Sydney         Fst        0.0395  
27 Rockhampton_v_Townsville     Dxy        0.000565
28 Rockhampton_v_Townsville     Fst        0.0847  
29 Sydney_v_Townsville          Dxy        0.000570
30 Sydney_v_Townsville          Fst        0.131   
'


#Plotting Dxy

# plot 1 - genome wide plots per comparison
plot_1_dxy <- data %>%
  filter(data_type =='Dxy')%>%
  ggplot(., aes(position*100000, value, col=chromosome)) +
  geom_point(size=1.5) +
  facet_grid(comparison~.) +
  scale_color_npg(guide = guide_legend(override.aes = list(size = 4))) +
  theme_bw() +
  theme(legend.position="top", 
        legend.title = element_blank(),
        axis.title = element_text(size = 32),
        axis.text = element_text(size = 10, color = "black"),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 20),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 2)) +
  labs(x="Genomic Position", y="Dxy")
plot_1_dxy

# plot 2 - density plots of dxy per group
plot_2_dxy <- data %>%
  filter(data_type=='Dxy')%>%
  ggplot(aes(value, chr_type, fill=chr_type), guide="none") +
  geom_density_ridges(quantile_lines=TRUE, quantile_fun=function(x,...)median(x), size=0.5) +
  theme_bw() + theme(legend.position = "top", 
                     legend.title = element_blank(),
                     axis.text.y=element_blank(),
                     axis.title = element_text(size = 32),
                     axis.text = element_text(size = 16, color = "black"),
                     strip.text = element_text(size = 8),
                     legend.text = element_text(size = 20),
                     panel.grid = element_blank(),
                     panel.border = element_rect(color = "black", fill = NA, linewidth = 2)) +
  facet_grid(comparison~.) +
  scale_fill_npg() +
  labs(x="Dxy", y="Density")
plot_2_dxy

# combine plots
dxy_plot <- plot_1_dxy + plot_2_dxy +  plot_layout(widths = c(5, 1))
ggsave("genomewide_and_density_dxy_aus_dog.tif", width=20, height=25, dpi = 300)
ggsave("genomewide_and_density_dxy_aus_dog.png", width=20, height=25, dpi = 300)
ggsave("genomewide_and_density_dxy_aus_dog.pdf", width=20, height=25, dpi = 300)
# higher Dxy = more divergent

#some additional plots
boxplot_dxy <- data %>%
  filter(data_type =='Dxy') %>% 
  ggplot(., aes(comparison, value, col=comparison)) +
  geom_jitter(size = 1, alpha = 0.5) +
  geom_boxplot(fill=NA, col="black") +
  labs(x = "Population" , y = "Dxy", colour = "Population") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title = element_text(size =24),
        axis.text = element_text(size = 14, color = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 2))
boxplot_dxy
ggsave("boxplot_dxy_aus_dog.tif", width=14, height = 6, dpi = 300)
ggsave("boxplot_dxy_aus_dog.png", width=14, height = 6, dpi = 300)
ggsave("boxplot_dxy_aus_dog.pdf", width=14, height = 6, dpi = 300)

density_dxy <- data %>%
  filter(data_type =='Dxy') %>% 
  ggplot(., aes(x=value, group = comparison, fill=comparison)) +
  geom_density(adjust=1.5, alpha=.4) +
  labs(x = "Dxy" , y = "Density", colour = "Comparison") +
  theme_bw() +
  theme(legend.position = c(0.78, 0.6),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 18, color = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 2))
density_dxy
ggsave("density_dxy_aus_dog.tif", width = 8, height = 6, dpi = 300)
ggsave("density_dxy_aus_dog.png", width = 8, height = 6, dpi = 300)
ggsave("density_dxy_aus_dog.pdf", width = 8, height = 6, dpi = 300)


#Trying to plot lines
#getting chr positions
chr <- c('X', '1', '2', '3', '4')

for (i in chr) {
  x <- data %>%
    filter(chromosome == paste0('chr', i))
  print(quantile(x$position))
}

data$comparison <- str_replace_all(data$comparison, '_', ' ')

dxy_lineplot <- data %>%
  filter(data_type=='Dxy')%>%
  ggplot(., aes(position*100000, value, col=comparison)) +
  geom_point(size=0.2, alpha = 0.1) +
  geom_vline(xintercept=c(281*100000,
                          438*100000,
                          590*100000,
                          740*100000),
             size=1, linetype="dashed", col='grey41')+
  annotate(geom="text", x=140*100000, y=0.004, label="ChrX", size = 6)+
  annotate(geom="text", x=359.5*100000, y=0.004, label="Chr1", size = 6)+
  annotate(geom="text", x=520*100000, y=0.004, label="Chr2", size = 6)+
  annotate(geom="text", x=665*100000, y=0.004, label="Chr3", size = 6)+
  annotate(geom="text", x=810.5*100000, y=0.004, label="Chr4", size = 6)+
  geom_line() +
  theme_bw() +
  theme(legend.position="top", 
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 18)) +
  labs(x="Genomic Position", y="Dxy") +
  guides(col = guide_legend(override.aes = list(linewidth = 2)))
dxy_lineplot
ggsave("lineplot_dxy_aus_dog.tif", width=16, height=8, dpi = 300)
ggsave("lineplot_dxy_aus_dog.png", width=16, height=8, dpi = 300)
ggsave("lineplot_dxy_aus_dog.pdf", width=16, height=8, dpi = 300)



# Dxy by city on a heatmap
# Get median dxy for each comparison
dxy_median <- data %>%
  filter(data_type == "Dxy") %>%
  group_by(comparison) %>%
  summarise(DXY_MEDIAN = median(value, na.rm = TRUE))

# Save as csv file
write.csv(dxy_median, "dxy_median.csv", row.names = TRUE)

dxy_heatmap <- read.csv("dxy_heatmap.csv", header = TRUE)
dxy_heatmap <- as.data.frame(dxy_heatmap)

plot_3_dxy <- ggplot(dxy_heatmap, aes(x = POP1, y=POP2, fill = DXY_MEDIAN)) + 
  geom_tile(color = "white", lwd = 1.5, linetype = 1) +
  coord_fixed() +
  scale_fill_material("indigo", labels = scales::label_number()) +
  guides(fill = guide_colourbar(title = expression(D[XY]), barwidth = 1, barheight = 5)) +
  theme_minimal() +
  labs(x = "Population", y = "Population") +
  theme(axis.text = element_text(size = 18),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title = element_text(size = 24),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank())
plot_3_dxy

ggsave("dxy_heatmap_aus_dog.tif", bg = "white", height=5, width = 8, dpi = 300)
ggsave("dxy_heatmap_aus_dog.png", bg = "white", height=5, width = 8, dpi = 300)
ggsave("dxy_heatmap_aus_dog.pdf", bg = "white", height=5, width = 8, dpi = 300)




# Plotting Fst

# plot 1 - genome wide plots per population
plot_1_fst <- data %>%
  filter(data_type=='Fst')%>%
  ggplot(., aes(position*100000, value, col=chromosome)) +
  geom_point(size=1.5) +
  facet_grid(comparison~.) +
  scale_color_npg(guide = guide_legend(override.aes = list(size = 4))) +
  theme_bw() +
  theme(legend.position="top", 
        legend.title = element_blank(),
        axis.title = element_text(size = 32),
        axis.text = element_text(size = 10, color = "black"),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 20),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 2)) +
  labs(x="Genomic Position", y="Fst")
plot_1_fst


# plot 2 - density plots of Fst per group
plot_2_fst <- data %>%
  filter(data_type=='Fst')%>%
  ggplot(aes(value, chr_type, fill=chr_type), guide="none") +
  geom_density_ridges(quantile_lines=TRUE, quantile_fun=function(x,...)median(x), size=0.5) +
  theme_bw() + 
  theme(legend.position = "none", 
        axis.text.y=element_blank(),
        axis.title = element_text(size = 32),
        axis.text = element_text(size = 10, color = "black"),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 20),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 2)) +
  facet_grid(comparison~.) +
  scale_fill_npg() +
  labs(x="Fst", y="Density")
plot_2_fst

# combine plots
plot_1_fst + plot_2_fst +  plot_layout(widths = c(5, 1))
# 1: Removed 114 rows containing missing values or values outside the scale range (`geom_point()`). 
# 2: Removed 114 rows containing non-finite outside the scale range (`stat_density_ridges()`). 
## Fst values should only range between 0 and 1. Set anything <0 to 0 and re-do the plots.

# Data cleaning
data_fst_clean <- data %>% filter(data_type == 'Fst') %>% mutate(value = ifelse(value < 0, 0, value)) # this takes my original data, filters it to only include rows where the data_type is Fst. Then mutate function is used to create/modify columns. It modifies the value column, checks if the value is less than 0. If the value is <0, it replaces it with 0. If the value is >0, it keeps the original value.

# Re-do plots with cleaned data
# plot 1 - genome wide plots per population
plot_1_fst_clean <- data_fst_clean %>%
  filter(data_type=='Fst')%>%
  ggplot(., aes(position*100000, value, col=chromosome)) +
  geom_point(size=1) +
  facet_grid(comparison~.) +
  scale_color_tron() +
  theme_bw() +
  theme(legend.position="top", 
        legend.title = element_blank(),
        axis.title = element_text(size = 32),
        axis.text = element_text(size = 16, color = "black"),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 20),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 2)) +
  labs(x="Genomic Position", y="Fst")
plot_1_fst_clean


# plot 2 - density plots of Fst per group
plot_2_fst_clean <- data_fst_clean %>%
  filter(data_type=='Fst')%>%
  ggplot(aes(value, chr_type, fill=chr_type), guide="none") +
  geom_density_ridges(quantile_lines=TRUE, quantile_fun=function(x,...)median(x), size=0.5) +
  theme_bw() + theme(legend.position = "none", axis.text.y=element_blank(), panel.grid = element_blank(), panel.border = element_rect(color = "black", fill = NA, linewidth = 2)) +
  facet_grid(comparison~.) +
  scale_fill_npg() +
  labs(x="Fst", y="Density")
plot_2_fst_clean

# combine plots
plot_1_fst_clean + plot_2_fst_clean +  plot_layout(widths = c(5, 1))
# 1: Removed 114 rows containing missing values or values outside the scale range (`geom_point()`). 
# 2: Removed 114 rows containing non-finite outside the scale range (`stat_density_ridges()`). 
ggsave("genomewide_and_density_fst_aus_dog.tif", width=20, height=25, dpi = 300)
ggsave("genomewide_and_density_fst_aus_dog.png", width=20, height=25, dpi = 300)
ggsave("genomewide_and_density_fst_aus_dog.pdf", width=20, height=25, dpi = 300)
# Fst closer to 1 = more genetically distinct
# Fst closer to 0 = more genetically similar
# Some cities only have 1 sample so they look funny.



# AUS map

## Make world map data
world_map <- map_data("world")

# Add map inset to zoom in on Australian dog samples
# Manually specify the coordinates for the area of the world map to show in the inset
aus_xmin <- 141
aus_xmax <- 155
aus_ymin <- -60
aus_ymax <- -9

# Filter the world map data for the inset area
aus_data <- subset(world_map, long >= aus_xmin & long <= aus_xmax & lat >= aus_ymin & lat <= aus_ymax)

# Aus metadata
location_AUS <- read.csv("input/AUS_metadata.csv", header = TRUE)

# Set colors for the points
blue_palette <- brewer.pal(n = 9, name = "Blues")

scale_colour_AUS <- 
  c('Lockhart River' = blue_palette[9],
    'Cairns' = blue_palette[8],
    'Townsville' = blue_palette[7],
    'Rockhampton' = blue_palette[6],
    'Brisbane' = blue_palette[5],
    'Sydney' = blue_palette[4])

# Get median fst for each comparison
fst_median <- data_fst_clean %>%
  filter(data_type == "Fst") %>%
  group_by(comparison) %>%
  summarise(FST_MEDIAN = median(value, na.rm = TRUE))
# Save as csv file and then add required metadata such as population longitudes/latitudes
write.csv(fst_median, "fst_median.csv", row.names = TRUE)

# read new csv file with map metadata and median fst between groups
fst_map <- read.csv("fst_map.csv", header = TRUE)

# Fst on AUS map
plot_3_fst <- ggplot() +
  geom_polygon(data = aus_data, aes(x = long, y = lat, group = group), fill = "grey90") +
  geom_segment(data=fst_map, aes(x=POP1_LONG, y=POP1_LAT, xend=POP2_LONG, yend=POP2_LAT, color = FST_MEDIAN), size = 1) +
  scale_color_gradient(low = "gold", high = "red3", name = expression(F[ST])) +
  geom_point(data = location_AUS, aes(x = Longitude, y = Latitude, fill = City), size = 10, shape = 21) +
  scale_fill_manual(values = scale_colour_AUS, limits = c("Lockhart River", "Cairns", "Townsville", "Rockhampton", "Brisbane", "Sydney"), name = "City") +
  geom_text_repel(data = location_AUS, aes(x = Longitude, y = Latitude, label = City), size = 6, fontface = "bold", nudge_x = 5, nudge_y = 2) +
  theme_void() +
  theme(
    legend.position = c(0.9,0.1),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 14)) +
  guides(fill = "none") +
  coord_fixed(ratio=1)
plot_3_fst

# Save plot
ggsave("fst_map_aus_dog.tif", bg = "white", height=10, width=5, dpi = 300)
ggsave("fst_map_aus_dog.png", bg = "white", height=10, width=5, dpi = 300)
ggsave("fst_map_aus_dog.pdf", bg = "white", height=10, width=5, dpi = 300)


# Fst by AUS on a heatmap
fst_heatmap <- read.csv("fst_heatmap.csv", header = TRUE)

plot_4_fst <- ggplot(fst_heatmap, aes(x = POP1, y=POP2, fill = FST_MEDIAN)) + 
  geom_tile(color = "white", lwd = 1.5, linetype = 1) +
  coord_fixed() +
  scale_fill_material("indigo") +
  guides(fill = guide_colourbar(title = expression(F[ST]), barwidth = 1, barheight = 5)) +
  theme_minimal() +
  labs(x = "Population", y = "Population") +
  theme(axis.text = element_text(size = 18),
    axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title = element_text(size = 24),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank())
plot_4_fst

ggsave("fst_heatmap_aus_dog.tif", bg = "white", height=5, width = 8, dpi = 300)
ggsave("fst_heatmap_aus_dog.png", bg = "white", height=5, width = 8, dpi = 300)
ggsave("fst_heatmap_aus_dog.pdf", bg = "white", height=5, width = 8, dpi = 300)
```