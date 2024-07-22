# Troubleshooting outgroups

The outgroup data looked a bit funny when attempting to run admixtools F3. Let's investigate.

![](images/Mito_outgroups_Differences_vs_distance.png)

There's 5-14% difference in the mtDNA, so the nuclear is unlikely to be as high. We would expect there to have been more SNPs, so perhaps it is a technical reason it hasn't worked well.

![](images/Outgroups_fastqc.png)

The duplicated sequences are very high and most would've been thrown away in deduplication

![](images/Outgroups_mapping.png)

The mapping % weren't that great either. Low mapping rates could be due to diversity, or there could be little worm DNA in there. Check for contamination.

## kraken to test for contamination in my outgroup samples

```bash
# load modules
module load bsub.py/0.42.1
module load kraken2/2.1.2

# build custom database
WORKING_DIR=/lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA
cd ${WORKING_DIR}/03_ANALYSIS/01_PREP/KRAKEN

DBNAME='kraken_dog_bear_human_Di.db'

# download taxonomy
kraken2-build --download-taxonomy --db $DBNAME


# add ref genome from dog to kraken database
kraken2-build --add-to-library /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/01_REF/GCF_014441545.1_ROS_Cfam_1.0_genomic.fna --no-masking --db $DBNAME


# download American black bear
cd REFS
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/020/975/775/GCF_020975775.1_gsc_jax_bbear_1.0/GCF_020975775.1_gsc_jax_bbear_1.0_genomic.fna.gz
if [ "$(md5sum GCF_020975775.1_gsc_jax_bbear_1.0_genomic.fna.gz | awk '{ print $1 }')" == "38e099b357366455579be84d88deec8f" ]; then echo "MD5 checksum matches."; else echo "MD5 checksum does not match."; fi
## MD5 checksum matches.
# unzip
gunzip GCF_020975775.1_gsc_jax_bbear_1.0_genomic.fna.gz
# add to kraken database
cd ..
kraken2-build --add-to-library REFS/GCF_020975775.1_gsc_jax_bbear_1.0_genomic.fna --no-masking --db $DBNAME


# download human
cd REFS
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
if [ "$(md5sum GCF_000001405.40_GRCh38.p14_genomic.fna.gz  | awk '{ print $1 }')" == "c30471567037b2b2389d43c908c653e1" ]; then echo "MD5 checksum matches."; else echo "MD5 checksum does not match."; fi
## MD5 checksum matches.
# unzip
gunzip GCF_000001405.40_GRCh38.p14_genomic.fna.gz
# add to kraken database
cd ..
kraken2-build --add-to-library REFS/GCF_000001405.40_GRCh38.p14_genomic.fna --no-masking --db $DBNAME


# add Di
#We have dowloanded it from NCBI databes, since it includes the taxonomy info from NCBI
cd REFS
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/077/395/GCA_001077395.1_ASM107739v1/GCA_001077395.1_ASM107739v1_genomic.fna.gz
if [ "$(md5sum GCA_001077395.1_ASM107739v1_genomic.fna.gz | awk '{ print $1 }')" == "32025403ce2abb28e69d754ed6110156" ]; then echo "MD5 checksum matches."; else echo "MD5 checksum does not match."; fi
## MD5 checksum matches.
# unzip
gunzip GCA_001077395.1_ASM107739v1_genomic.fna.gz
# add to kraken database
cd ..
kraken2-build --add-to-library REFS/GCA_001077395.1_ASM107739v1_genomic.fna --no-masking --db $DBNAME




# build database
bsub.py --threads 8 50 kraken_build "kraken2-build --threads 8 --build --db $DBNAME"
# this can take a while

# create loop
while read line; do
    kraken2 --db $DBNAME --report $line\.kraken --paired ../$line\/$line\_val_1.fq.gz ../$line\/$line\_val_2.fq.gz;
done < fastq.txt
```




## How well does D. repens map to D. immitis?

- Can make fake reads from a D. repens assembly, then map them to D. immitis. This allows us to see how well D. repens maps to D. immitis in a controlled way.

I can try using BBMap randomreads.sh:
This generates random synthetic reads from a reference genome: https://github.com/BioInfoTools/BBMap/blob/master/sh/randomreads.sh

D. repens reference genome: https://doi.org/10.1016/j.ijpara.2019.04.008 

```bash
# load modules
module load bbtools/39.01
module load samtools/1.14--hb421002_0

cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/05_ANALYSIS/OUTGROUPS/DREPENS/DATA

# download ref genome for D. repens
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/008/729/115/GCA_008729115.1_FGCZ_Drep_1.0/*.fna.gz
if [ "$(md5sum GCA_008729115.1_FGCZ_Drep_1.0_genomic.fna.gz | awk '{ print $1 }')" == "5187f349044c7ca68c7305e566c0c577" ]; then echo "MD5 checksum matches."; else echo "MD5 checksum does not match."; fi
## MD5 checksum matches.
# unzip
gunzip GCA_008729115.1_FGCZ_Drep_1.0_genomic.fna.gz
# index
bsub.py 2 repens_index "samtools faidx DATA/GCA_008729115.1_FGCZ_Drep_1.0_genomic.fna"

cd ..
bsub.py 10 repens_randomreads "randomreads.sh ref=DATA/GCA_008729115.1_FGCZ_Drep_1.0_genomic.fna out1=Drep_1.fq.gz out2=Drep_2.fq.gz simplenames=t length=150 coverage=35 paired=t midq=32"


# Not working bc the read names are too long - try another method
module load art/2016.06.05--h869255c_2

bsub.py 4 repens_art "art_illumina --paired --in DATA/GCA_008729115.1_FGCZ_Drep_1.0_genomic.fna --len 150 --fcov 30 --mflen 200 --sdev 10 --noALN --out repens_fake"
```
This worked better, the read names are much shorter.


### Map fake reads to combined D. immitis & Wol & dog genome

```bash
# load modules
module load mapping-helminth/v1.0.9

# Set variables
WORKING_DIR=/lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA
REF=${WORKING_DIR}/01_REF/reference_di_wol_dog.fa

# Run nextflow mapping pipeline
bsub.py 10 repens_fake_mapping "mapping-helminth --input repens_fake_wgs.mapping.manifest --reference ${REF}"
```

Needed to run a bit longer than the 12h limit. Results are in: /lustre/scratch125/pam/teams/team333/sd21/dirofilaria_immitis/ROSE_AUS_HW/FAKEREADS

![](images/Drepens_FAKEREADS_mapping.PNG)

~18% of reads mapped, which is still pretty low. So there is clearly low mapping due to divergence.


### Try variant calling with the mapped fake reads

```bash
cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/05_ANALYSIS/OUTGROUPS/DREPENS/VARIANTS

# create bam list
ls /lustre/scratch125/pam/teams/team333/sd21/dirofilaria_immitis/ROSE_AUS_HW/FAKEREADS/results/Drepens/*.bam > original.bamlist

bsub.py 1 run_gatk_hc "./run_gatk_hc.sh"
```

run_gatk_hc.sh:
```bash
#!/bin/bash

#-------------------------------------------------------------------------------
# run_gatk_hc.sh
#-------------------------------------------------------------------------------

# stephen doyle
# Jan 2023

# Export environment variables
export PREFIX=DIMMITIS_POPGEN  # prefix for output files
export REFERENCE=/lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/01_REF/dimmitis_WSI_2.2.fa  # path to reference genome
export BAM_LIST=/lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/05_ANALYSIS/OUTGROUPS/DREPENS/VARIANTS/original.bamlist  # path to list of BAM files

# Load GATK module
module load gatk/4.1.4.1

# Load samtools and fastaq modules for reference step
module load samtools/1.14--hb421002_0
module load fastaq/3.17.0-docker3

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
### 01. Prepare reference files
#-------------------------------------------------------------------------------

func_build_reference() {
    # Check if the reference genome file already exists
    if [ -f "${REFERENCE_FILES}/REF.fa" ]; then
        echo "Reference is already setup. Moving on."
        exit 0
    else
        # Copy the reference genome file to the REFERENCE_FILES directory
        cp "${REFERENCE}" "${REFERENCE_FILES}/REF.fa"
        # Create an index file for the reference genome
        samtools faidx "${REFERENCE_FILES}/REF.fa"
        # Create a dictionary file for the reference genome
        samtools dict "${REFERENCE_FILES}/REF.fa" > "${REFERENCE_FILES}/REF.dict"

        # Append the BAM_LIST file to the bam.list file in the REFERENCE_FILES directory
        cat "${BAM_LIST}" >> "${REFERENCE_FILES}/bam.list"

        # Split the reference genome file into chunks of approximately 10 Mb in size
        fastaq split_by_base_count "${REFERENCE_FILES}/REF.fa" "${REFERENCE_FILES}/REFsplit" 10000000

        # For each chunked genome section, create a list of contig/scaffold names
        for i in "${REFERENCE_FILES}/REFsplit"* ; do
            # Extract the name of the chunked genome section
            NAME=$( echo "${i}" | awk -F '/' '{print $NF}' )
            # Extract the contig/scaffold names from the chunked genome section
            grep ">" "${i}" | sed 's/>//g' > "${REFERENCE_FILES}/${NAME}.list"
        done
    fi
}

export -f func_build_reference


#-------------------------------------------------------------------------------
### 02. Make GVCF per sample
#-------------------------------------------------------------------------------
func_make_gvcf() { 

# make jobs
COUNT=0
while read BAM; do
	n=1

    SAMPLE=$( echo ${BAM} | awk -F '/' '{print $NF}' | sed -e 's/.bam//g' )

    # check if the sample directory exists already - if yes, stop and move on
    if [ -d "${GATK_HC_GVCFs}/${SAMPLE}_GATK_HC_GVCF_complete" ] ; then
        echo -e "\nThere is already a run started/completed with this sample name. Rename and start again, or move on to the enxt sample\n"
        continue
	else

    # make sample directories
    mkdir ${GATK_HC_GVCFs}/${SAMPLE}_GATK_HC_GVCF_started
	mkdir ${GATK_HC_GVCFs}/${SAMPLE}_GATK_HC_GVCF_started/LOGFILES

     echo "gatk GatherVcfsCloud \\" > ${GATK_HC_GVCFs}/${SAMPLE}_GATK_HC_GVCF_started/run_gather_${SAMPLE}_gvcf

    for SEQUENCE in ${REFERENCE_FILES}/REFsplit*list; do
        SEQUENCE=$( echo ${SEQUENCE} | awk -F '/' '{print $NF}' )
	    echo -e "gatk HaplotypeCaller --input ${BAM} --output ${GATK_HC_GVCFs}/${SAMPLE}_GATK_HC_GVCF_started/${n}.${SAMPLE}.tmp.g.vcf.gz --reference ${REFERENCE_FILES}/REF.fa --intervals ${REFERENCE_FILES}/${SEQUENCE} --emit-ref-confidence GVCF " > ${GATK_HC_GVCFs}/${SAMPLE}_GATK_HC_GVCF_started/run_hc_${SAMPLE}.tmp.job_${n};
	    echo -e "--input ${GATK_HC_GVCFs}/${SAMPLE}_GATK_HC_GVCF_started/${n}.${SAMPLE}.tmp.g.vcf.gz \\" >> ${GATK_HC_GVCFs}/${SAMPLE}_GATK_HC_GVCF_started/run_gather_${SAMPLE}_gvcf;
	    let "n+=1";
        done;

	echo -e "--output ${GATK_HC_GVCFs}/${SAMPLE}_GATK_HC_GVCF_started/${SAMPLE}.g.vcf.gz; tabix -p vcf ${GATK_HC_GVCFs}/${SAMPLE}_GATK_HC_GVCF_started/${SAMPLE}.g.vcf.gz" >> ${GATK_HC_GVCFs}/${SAMPLE}_GATK_HC_GVCF_started/run_gather_${SAMPLE}_gvcf;

	echo -e "rm ${GATK_HC_GVCFs}/${SAMPLE}_GATK_HC_GVCF_started/*.tmp.* && mv ${GATK_HC_GVCFs}/${SAMPLE}_GATK_HC_GVCF_started ${GATK_HC_GVCFs}/${SAMPLE}_GATK_HC_GVCF_complete" > ${GATK_HC_GVCFs}/${SAMPLE}_GATK_HC_GVCF_started/run_clean_${SAMPLE};

	chmod a+x ${GATK_HC_GVCFs}/${SAMPLE}_GATK_HC_GVCF_started/run_*

	# setup job conditions
	JOBS=$( ls -1 ${GATK_HC_GVCFs}/${SAMPLE}_GATK_HC_GVCF_started/run_hc_* | wc -l )
	ID="U$(date +%s)"

	#submit job array to call variants put scaffold / contig
	bsub -q long -R'span[hosts=1] select[mem>15000] rusage[mem=15000]' -n 6 -M15000 -J "gatk_make_gvcf_${ID}_[1-$JOBS]%100" -e "${GATK_HC_GVCFs}/${SAMPLE}_GATK_HC_GVCF_started/LOGFILES/gatk_make_gvcf_${ID}_[1-$JOBS].e" -o "${GATK_HC_GVCFs}/${SAMPLE}_GATK_HC_GVCF_started/LOGFILES/gatk_make_gvcf_${ID}_[1-$JOBS].o" "${GATK_HC_GVCFs}/${SAMPLE}_GATK_HC_GVCF_started/run_hc_${SAMPLE}.tmp.job_\$LSB_JOBINDEX"

	#submit job to gather gvcfs into a single, per sample gvcf
	bsub -q normal -w "done(gatk_make_gvcf_${ID}_[1-$JOBS])" -R'span[hosts=1] select[mem>500] rusage[mem=500]' -n 1 -M500 -J "gatk_gather_gvcf_${ID}" -e "${GATK_HC_GVCFs}/${SAMPLE}_GATK_HC_GVCF_started/LOGFILES/gatk_gather_gvcf_${ID}.e" -o "${GATK_HC_GVCFs}/${SAMPLE}_GATK_HC_GVCF_started/LOGFILES/gatk_gather_gvcf_${ID}s.o" "${GATK_HC_GVCFs}/${SAMPLE}_GATK_HC_GVCF_started/run_gather_${SAMPLE}_gvcf"

	# clean up
	bsub -q normal -w "done(gatk_gather_gvcf_${ID})" -R'span[hosts=1] select[mem>500] rusage[mem=500]' -n 1 -M500 -J "gatk_clean_gvcf_${ID}" -e "${GATK_HC_GVCFs}/${SAMPLE}_GATK_HC_GVCF_started/LOGFILES/gatk_clean_gvcf_${ID}.e" -o "${GATK_HC_GVCFs}/${SAMPLE}_GATK_HC_GVCF_started/LOGFILES/gatk_clean_gvcf_${ID}.o" "${GATK_HC_GVCFs}/${SAMPLE}_GATK_HC_GVCF_started/run_clean_${SAMPLE}"

	sleep 1
	fi
done < ${BAM_LIST}


# check that GVCF directories are complete before finishing this step
while true; do
  found=0
  count=$(cat $BAM_LIST | wc -l)
while read -r NAME; do 
    if [ -d "${GATK_HC_GVCFs}/${NAME}_GATK_HC_GVCF_complete" ]; then
      found=$((found + 1))
    fi
    done  < <(cat "$BAM_LIST" | awk -F '/' '{print $NF}' | sed 's/.bam//g' )  
  if [ ${found} -eq ${count} ]; then
    echo "Directories are complete. Moving on."
    break
  fi
  echo "Directories not found, waiting..."
  sleep 20
done

}

export -f func_make_gvcf



#-------------------------------------------------------------------------------
### 03. Merge GVCFs
#-------------------------------------------------------------------------------
func_merge_gvcf() { 

ls -1 ${GATK_HC_GVCFs}/*complete/*gz > ${GATK_HC_MERGED}/gvcf.list

[ -d ${GATK_HC_MERGED}/LOGFILES ] || mkdir -p ${GATK_HC_MERGED}/LOGFILES


n=1
for SEQUENCE in ${REFERENCE_FILES}/REFsplit*list; do
    SEQUENCE=$( echo ${SEQUENCE} | awk -F '/' '{print $NF}' )
    echo -e "gatk CombineGVCFs -R ${REFERENCE_FILES}/REF.fa --intervals ${REFERENCE_FILES}/${SEQUENCE} \\" > ${GATK_HC_MERGED}/run_merge_gvcfs.tmp.job_${n}
    while read SAMPLE; do
        echo -e "--variant ${SAMPLE} \\" >> ${GATK_HC_MERGED}/run_merge_gvcfs.tmp.job_${n};
   done < ${GATK_HC_MERGED}/gvcf.list
   echo -e "--output ${GATK_HC_MERGED}/${n}.cohort.tmp.g.vcf.gz" >> ${GATK_HC_MERGED}/run_merge_gvcfs.tmp.job_${n};
   let "n+=1"; 
done

chmod a+x ${GATK_HC_MERGED}/run_merge_gvcfs.tmp.job_*

# setup job conditions
JOBS=$( ls -1 ${GATK_HC_MERGED}/run_merge_gvcfs.tmp.job_* | wc -l )
ID="U$(date +%s)"

#submit job array to call variants put scaffold / contig
bsub -q long -R'span[hosts=1] select[mem>30000] rusage[mem=30000]' -n 10 -M30000 -J "gatk_merge_gvcf_[1-$JOBS]%100" -e "${GATK_HC_MERGED}/LOGFILES/gatk_merge_gvcf_[1-$JOBS].e" -o "${GATK_HC_MERGED}/LOGFILES/gatk_merge_gvcf_[1-$JOBS].o" "${GATK_HC_MERGED}/run_merge_gvcfs.tmp.job_\$LSB_JOBINDEX"

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
for SEQUENCE in ${REFERENCE_FILES}/REFsplit*list; do
    SEQUENCE=$( echo ${SEQUENCE} | awk -F '/' '{print $NF}' )
    echo -e "gatk GenotypeGVCFs \
    -R ${REFERENCE_FILES}/REF.fa \
    -V ${GATK_HC_MERGED}/${n}.cohort.tmp.g.vcf.gz \
    --intervals ${REFERENCE_FILES}/${SEQUENCE} \
    -O ${GATK_HC_MERGED}/${n}.cohort.tmp.vcf.gz -G StandardAnnotation -G AS_StandardAnnotation" > ${GATK_HC_MERGED}/run_hc_genotype.tmp.job_${n};
    let "n+=1"; 
done

chmod a+x ${GATK_HC_MERGED}/run_hc_genotype*

# setup job conditions
JOBS=$( ls -1 ${GATK_HC_MERGED}/run_hc_genotype* | wc -l )
ID="U$(date +%s)"

bsub -q long -R'span[hosts=1] select[mem>20000] rusage[mem=20000]' -n 6 -M20000 -J "gatk_genotype_cohort_gvcf_[1-$JOBS]" -e "${GATK_HC_MERGED}/LOGFILES/gatk_genotype_cohort_gvcf_[1-$JOBS].e" -o "${GATK_HC_MERGED}/LOGFILES/gatk_genotype_cohort_gvcf_[1-$JOBS].o" "${GATK_HC_MERGED}/run_hc_genotype.tmp.job_*\$LSB_JOBINDEX"

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
    rm ${GATK_HC_MERGED}/*tmp*

}

export -f func_finish_vcf





#-------------------------------------------------------------------------------
# running the pipeline
#-------------------------------------------------------------------------------

# func_build_reference
bsub -E 'test -e /nfs/users/nfs_r/rp24' -R "select[mem>1000] rusage[mem=1000]" -M1000 -o ${LOG_FILES}/gatk_01_build_reference.o -e ${LOG_FILES}/gatk_01_build_reference.e -J gatk_01_build_reference_${PREFIX} func_build_reference

# func_make_gvcf
bsub -w "done(gatk_01_build_reference_${PREFIX})" -E 'test -e /nfs/users/nfs_r/rp24' -R "select[mem>50000] rusage[mem=50000]" -q long -M50000 -n20 -o ${LOG_FILES}/gatk_02_make_gvcf.o -e ${LOG_FILES}/gatk_02_make_gvcf.e -J gatk_02_make_gvcf_${PREFIX} func_make_gvcf

# func_merge_gvcf
bsub -w "done(gatk_02_make_gvcf_${PREFIX})" -E 'test -e /nfs/users/nfs_r/rp24' -R "select[mem>50000] rusage[mem=50000]" -q long -M50000 -n20 -o ${LOG_FILES}/gatk_03_merge_gvcf.o -e ${LOG_FILES}/gatk_03_merge_gvcf.e -J gatk_03_merge_gvcf_${PREFIX} func_merge_gvcf

# func_genotype_gvcfs
bsub -w "done(gatk_03_merge_gvcf_${PREFIX})" -E 'test -e /nfs/users/nfs_r/rp24' -R "select[mem>50000] rusage[mem=50000]" -q long -M50000 -n20 -o ${LOG_FILES}/gatk_04_genotype_gvcfs.o -e ${LOG_FILES}/gatk_04_genotype_gvcfs.e -J gatk_04_genotype_gvcfs_${PREFIX} func_genotype_gvcfs

# func_finish_vcf
bsub -w "done(gatk_04_genotype_gvcfs_${PREFIX})" -E 'test -e /nfs/users/nfs_r/rp24' -R "select[mem>1000] rusage[mem=1000]" -q long -M1000 -n1 -o ${LOG_FILES}/gatk_05_finish_vcf.o -e ${LOG_FILES}/gatk_05_finish_vcf.e -J gatk_05_finish_vcf_${PREFIX} func_finish_vcf
```
Successfully completed.


Now run the filter1 parameters and see how many SNPs we have.

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
bsub.py 4 run_snps_qc "run_snps_qc.sh"
```

```bash
#!/bin/bash

# Load modules
module load gatk/4.1.4.1

WORKING_DIR=/lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/05_ANALYSIS/OUTGROUPS/DREPENS/VARIANTS/FILTER1

# set reference, vcf, and mitochondrial and Wb contig
REFERENCE=/lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/01_REF/dimmitis_WSI_2.2.fa
VCF=/lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/05_ANALYSIS/OUTGROUPS/DREPENS/VARIANTS/gatk_hc_DIMMITIS_POPGEN/GATK_HC_MERGED/DIMMITIS_POPGEN.cohort.2024-07-21.vcf.gz
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
--output DIMMITIS_POPGEN.cohort.2024-07-21.nuclearSNPs.vcf

# select nuclear INDELs
gatk SelectVariants \
--reference ${REFERENCE} \
--variant ${VCF} \
--select-type-to-include INDEL \
--exclude-intervals ${MIT_CONTIG} \
--exclude-intervals ${WB_CONTIG} \
--output DIMMITIS_POPGEN.cohort.2024-07-21.nuclearINDELs.vcf

# select mitochondrial SNPs
gatk SelectVariants \
--reference ${REFERENCE} \
--variant ${VCF} \
--select-type-to-include SNP \
--intervals ${MIT_CONTIG} \
--output DIMMITIS_POPGEN.cohort.2024-07-21.mitoSNPs.vcf

# select mitochondrial INDELs
gatk SelectVariants \
--reference ${REFERENCE} \
--variant ${VCF} \
--select-type-to-include INDEL \
--intervals ${MIT_CONTIG} \
--output DIMMITIS_POPGEN.cohort.2024-07-21.mitoINDELs.vcf

# select WB SNPs
gatk SelectVariants \
--reference ${REFERENCE} \
--variant ${VCF} \
--select-type-to-include SNP \
--intervals ${WB_CONTIG} \
--output DIMMITIS_POPGEN.cohort.2024-07-21.WbSNPs.vcf

# select WB INDELs
gatk SelectVariants \
--reference ${REFERENCE} \
--variant ${VCF} \
--select-type-to-include INDEL \
--intervals ${WB_CONTIG} \
--output DIMMITIS_POPGEN.cohort.2024-07-21.WbINDELs.vcf

# won't make tables for the various SNP/indel data because I'm just going to stick with the same filtering parameters used in 'filter1'.
```

Won't be making density plots this time, just use the ones obtained from filter1.


### Attending to the quantiles, thresholds for specific parameters are established

```bash
# load modules
module load bsub.py/0.42.1
module load gatk/4.1.4.1

WORKING_DIR=/lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/05_ANALYSIS/OUTGROUPS/DREPENS/VARIANTS/FILTER1
cd ${WORKING_DIR}

# set reference
REFERENCE=/lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/01_REF/dimmitis_WSI_2.2.fa

#Nuclear
bsub.py 1 filter_nuclearSNPs "gatk VariantFiltration \
--reference ${REFERENCE} \
--variant DIMMITIS_POPGEN.cohort.2024-07-21.nuclearSNPs.vcf \
--filter-expression 'QUAL < 40 || DP < 1837 || DP > 11972 || MQ < 21.52 || SOR > 7.112 || QD < 0.780 || FS > 8.275 || MQRankSum < -5.870 || ReadPosRankSum < -2.792 || ReadPosRankSum > 2.240' \
--filter-name "SNP_filtered" \
--output DIMMITIS_POPGEN.cohort.2024-07-21.nuclearSNPs.filtered.vcf"

bsub.py 1 filter_nuclearINDELS "gatk VariantFiltration \
--reference ${REFERENCE} \
--variant DIMMITIS_POPGEN.cohort.2024-07-21.nuclearINDELs.vcf \
--filter-expression 'QUAL < 46 || DP < 787 || DP > 12600 || MQ < 22.83 || SOR > 6.736 || QD < 0.990 || FS > 7.324 || MQRankSum < -4.797 || ReadPosRankSum < -4.178 || ReadPosRankSum > 1.960' \
--filter-name "INDEL_filtered" \
--output DIMMITIS_POPGEN.cohort.2024-07-21.nuclearINDELs.filtered.vcf"

#Mitochondrial
bsub.py 1 filter_mitoSNPs "gatk VariantFiltration \
--reference ${REFERENCE} \
--variant DIMMITIS_POPGEN.cohort.2024-07-21.mitoSNPs.vcf \
--filter-expression 'QUAL < 77 || DP < 199695 || DP > 659990 || MQ < 20.92 || SOR > 12.402 || QD < 3.320 || FS > 41.869 || MQRankSum < -5.622 || ReadPosRankSum < -3.003 || ReadPosRankSum > 5.452' \
--filter-name "SNP_filtered" \
--output DIMMITIS_POPGEN.cohort.2024-07-21.mitoSNPs.filtered.vcf"

bsub.py 1 filter_mitoINDELS "gatk VariantFiltration \
--reference ${REFERENCE} \
--variant DIMMITIS_POPGEN.cohort.2024-07-21.mitoINDELs.vcf \
--filter-expression 'QUAL < 61 || DP < 200280 || DP > 655589 || MQ < 21.36 || SOR > 14.646 || QD < 0.288 || FS > 82.622 || MQRankSum < -7.928 || ReadPosRankSum < -7.416 || ReadPosRankSum > 6.170' \
--filter-name "INDEL_filtered" \
--output DIMMITIS_POPGEN.cohort.2024-07-21.mitoINDELs.filtered.vcf"

#Wolbachia
bsub.py 1 filter_WbSNPs "gatk VariantFiltration \
--reference ${REFERENCE} \
--variant DIMMITIS_POPGEN.cohort.2024-07-21.WbSNPs.vcf \
--filter-expression 'QUAL < 51 || DP < 93526 || DP > 140940 || MQ < 22.00 || SOR > 11.661 || QD < 1.790 || FS > 46.426 || MQRankSum < -15.506 || ReadPosRankSum < -6.165 || ReadPosRankSum > 13.108' \
--filter-name "SNP_filtered" \
--output DIMMITIS_POPGEN.cohort.2024-07-21.WbSNPs.filtered.vcf"

bsub.py 1 filter_WbINDELS "gatk VariantFiltration \
--reference ${REFERENCE} \
--variant DIMMITIS_POPGEN.cohort.2024-07-21.WbINDELs.vcf \
--filter-expression 'QUAL < 42 || DP < 93456 || DP > 145612 || MQ < 22.017 || SOR > 11.268 || QD < 1.210 || FS > 59.220 || MQRankSum < -6.716 || ReadPosRankSum < -6.413 || ReadPosRankSum > 6.014' \
--filter-name "INDEL_filtered" \
--output DIMMITIS_POPGEN.cohort.2024-07-21.WbINDELs.filtered.vcf"

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














- Can use raw reads from public D. repens data & see how well that maps to D. immitis

There's raw data for the D. repens ref genome (SRR8742586, SRR9613458-65), but it's pacbio single reads and they'll probably contain extra errors. Didn't find any other potential D. repens raw data.
