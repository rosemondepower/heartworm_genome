# Data processing for 2 final samples of the collection

## Check md5

```bash
cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/02_FASTQ/EXTRA/X201SC24040556-Z01-F003/X201SC24040556-Z01-F003
md5sum -c MD5.txt # all ok
```

## Merge fastq files for same samples

module load bsub.py/0.42.1
cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/00_SCRIPTS/LOGS
bsub.py 4 merge_fastq_extra "../merge_fastq_extra.sh"

```bash
cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/02_FASTQ/EXTRA/merged

# JS6766
zcat /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/02_FASTQ/EXTRA/X201SC24040556-Z01-F003/X201SC24040556-Z01-F003/01.RawData/JS6766/JS6766_DKDN240018872-1A_HH22TDSXC_L2_1.fq.gz /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/02_FASTQ/EXTRA/X201SC24040556-Z01-F003/X201SC24040556-Z01-F003/01.RawData/JS6766/JS6766_DKDN240018872-1A_HH2FVDSXC_L2_1.fq.gz | gzip > JS6766_1.fq.gz

zcat /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/02_FASTQ/EXTRA/X201SC24040556-Z01-F003/X201SC24040556-Z01-F003/01.RawData/JS6766/JS6766_DKDN240018872-1A_HH22TDSXC_L2_2.fq.gz /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/02_FASTQ/EXTRA/X201SC24040556-Z01-F003/X201SC24040556-Z01-F003/01.RawData/JS6766/JS6766_DKDN240018872-1A_HH2FVDSXC_L2_2.fq.gz | gzip > JS6766_2.fq.gz
```

# Check that the files merged correctly

```bash
# F reads
cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/02_FASTQ/EXTRA/X201SC24040556-Z01-F003/X201SC24040556-Z01-F003/01.RawData/JS6766
# Raw data files
for f in *_1.fq.gz; do echo $f;zcat $f|wc -l ; done > /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/02_FASTQ/EXTRA/COUNT/count_raw_1.txt
cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/02_FASTQ/EXTRA/COUNT
# Every 2nd line moved to new column
sed 'N;s/\n/ /g' count_raw_1.txt | column -t > raw_1.txt
# Get list of unique sample names. Make sample list file.
awk '{print $1}' OFS="\t" raw_1.txt | cut -c1-6 | uniq > samples_1.txt
# Finds all files for each individual sample. Saves to new file for each sample.
for f in $(cat samples_1.txt); do
grep $f raw_1.txt > ${f}_raw_1.txt
done
# prints sample ID and total at the bottom. Extract last line.
for f in JS*_raw_1.txt; do awk '{sum+=$2;print $1" "$2} END {print "'$f'", sum}' $f | tail -1 > total_$f; done
# combine files
cat total_JS*_raw_1.txt > total_raw_1.txt
# remove files I don't need anymore
rm *JS*.txt


# Merged data files
cd ../merged
# Forward reads
for f in *_1.fq.gz; do echo $f;zcat $f|wc -l ; done > ../COUNT/count_merged_1.txt
cd ../COUNT
# Every 2nd line moved to new column
sed 'N;s/\n/ /g' count_merged_1.txt | column -t > merged_1.txt
# Join the raw & merged stats for FORWARD reads
paste total_raw_1.txt merged_1.txt | column -s $'\t' -t > total_both_1.txt
# Make txt file into csv file
mv total_both_1.txt total_both_1.csv


# R reads
cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/02_FASTQ/EXTRA/X201SC24040556-Z01-F003/X201SC24040556-Z01-F003/01.RawData/JS6766
# Raw data files
for f in *_2.fq.gz; do echo $f;zcat $f|wc -l ; done > /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/02_FASTQ/EXTRA/COUNT/count_raw_2.txt
cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/02_FASTQ/EXTRA/COUNT
# Every 2nd line moved to new column
sed 'N;s/\n/ /g' count_raw_2.txt | column -t > raw_2.txt
# Get list of unique sample names. Make sample list file.
awk '{print $1}' OFS="\t" raw_2.txt | cut -c1-6 | uniq > samples_2.txt
# Finds all files for each individual sample. Saves to new file for each sample.
for f in $(cat samples_2.txt); do
grep $f raw_2.txt > ${f}_raw_2.txt
done
# prints sample ID and total at the bottom. Extract last line.
for f in JS*_raw_2.txt; do awk '{sum+=$2;print $1" "$2} END {print "'$f'", sum}' $f | tail -1 > total_$f; done
# combine files
cat total_JS*_raw_2.txt > total_raw_2.txt
# remove files I don't need anymore
rm *JS*.txt


# Merged data files
cd ../merged
# Forward reads
for f in *_2.fq.gz; do echo $f;zcat $f|wc -l ; done > ../COUNT/count_merged_2.txt
cd ../COUNT
# Every 2nd line moved to new column
sed 'N;s/\n/ /g' count_merged_2.txt | column -t > merged_2.txt
# Join the raw & merged stats for FORWARD reads
paste total_raw_2.txt merged_2.txt | column -s $'\t' -t > total_both_2.txt
# Make txt file into csv file
mv total_both_2.txt total_both_2.csv
```
Yep the numbers match.

## FastQC


We want to get some stats on the raw data.

```bash
cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/02_FASTQ/EXTRA
cp /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/02_FASTQ/EXTRA/X201SC24040556-Z01-F003/X201SC24040556-Z01-F003/01.RawData/JS6765/*.gz .
cp /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/02_FASTQ/EXTRA/X201SC24040556-Z01-F003/X201SC24040556-Z01-F003/01.RawData/JS6766/*.gz .
md5sum *.gz
`3860497aaddd1dcbabc71c4cad7d0526  JS6765_DKDN240018871-1A_HH22TDSXC_L2_1.fq.gz
b5002ecd1171ffc310269e7a594c847b  JS6765_DKDN240018871-1A_HH22TDSXC_L2_2.fq.gz
a259bf2af121a4ae9fa15d0de6bc83d4  JS6766_DKDN240018872-1A_HH22TDSXC_L2_1.fq.gz
376e11df03d3bc170bc736f35948f5da  JS6766_DKDN240018872-1A_HH22TDSXC_L2_2.fq.gz
bc0959724242813fac7175f6d0b881b8  JS6766_DKDN240018872-1A_HH2FVDSXC_L2_1.fq.gz
cc4a12e1a25dae55f59fcf33d0c5f2b4  JS6766_DKDN240018872-1A_HH2FVDSXC_L2_2.fq.gz
`
# all ok
```


```bash
# set variables
WORKING_DIR=/lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA
OUT_DIR=${WORKING_DIR}/03_ANALYSIS/01_PREP/FASTQC/RAW
cd ${OUT_DIR}
# make a list of the files
for file in ${WORKING_DIR}/02_FASTQ/EXTRA/*.gz; do
  basename ${file}
done > ${OUT_DIR}/fastq_extra.list
FASTQ_LIST=${OUT_DIR}/fastq_extra.list

# load modules
module load fastqc/0.12.1--hdfd78af_0
module load bsub.py/0.42.1

# set up run files
n=1
while read SAMPLE; do
echo -e "fastqc -t 1 -o ${OUT_DIR} ${WORKING_DIR}/02_FASTQ/EXTRA/${SAMPLE}" > run_fastqc_raw_extra_${SAMPLE}.tmp.job_${n};
let "n+=1";
done < ${FASTQ_LIST}

chmod a+x run_fastqc_raw_extra*

#run
for i in run_fastqc_raw_extra*; do
    bsub.py --threads 1 4 ${i} "./${i}";
done

# clean up
mv run_fastqc_raw_extra*.e run_fastqc_raw_extra*.o LOGS
rm run_fastqc_raw_extra*

# check for any errors
cd LOGS
grep -i "Exited" *extra*.o
grep -i "Successfully completed" *extra*.o | wc -l
grep -i "Error" *extra*.e
# All ok

# multiqc
# Load module
module load multiqc/1.17--pyhdfd78af_1
module load bsub.py/0.42.1
bsub.py 4 multiqc_raw "multiqc ${OUT_DIR} -o ${OUT_DIR}"
```


## Trimmomatic

```bash
# set variables
WORKING_DIR=/lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA
IN_DIR=${WORKING_DIR}/02_FASTQ/EXTRA
cd ${WORKING_DIR}/03_ANALYSIS/01_PREP/TRIMMOMATIC/EXTRA
OUT_DIR=${WORKING_DIR}/03_ANALYSIS/01_PREP/TRIMMOMATIC/EXTRA

# load modules
module load trimmomatic/0.39--1

# make a list of sample names
for file in ${WORKING_DIR}/02_FASTQ/EXTRA/*_1.f*q.gz; do
  basename ${file}
done > ${OUT_DIR}/samples_extra.list

# set up run files
n=1
while read LONG_SAMPLE_1; do
LONG_SAMPLE_2=${LONG_SAMPLE_1/_1./_2.} 
SHORT_SAMPLE=$(basename "$LONG_SAMPLE_1" | cut -d'_' -f1) 

JOB_FILE="run_trimmomatic_extra_${SHORT_SAMPLE}.tmp.job_${n}"

echo -e "trimmomatic PE \
-threads 1 \
-phred33 \
${IN_DIR}/${LONG_SAMPLE_1} \
${IN_DIR}/${LONG_SAMPLE_2} \
${OUT_DIR}/${SHORT_SAMPLE}_1_trimpaired.fq.gz ${OUT_DIR}/${SHORT_SAMPLE}_1_trimunpaired.fq.gz \
${OUT_DIR}/${SHORT_SAMPLE}_2_trimpaired.fq.gz ${OUT_DIR}/${SHORT_SAMPLE}_2_trimunpaired.fq.gz \
SLIDINGWINDOW:10:20 MINLEN:50" > ${JOB_FILE};
let "n+=1";
done < samples_extra.list

chmod a+x run_trimmomatic_extra*

#run
for i in run_trimmomatic_extra*; do
    bsub.py --threads 1 20 ${i} "./${i}";
done


# SLIDINGWINDOW:10:20 means it will scan the read with a 10-base wide sliding window, cutting when the average quality per base drops below 20.

# Instead of SLIDINGWINDOW, in my previous practice code I used 'AVGQUAL:30 MINLEN:150'.

# clean up
mkdir LOGS
mv run_trimmomatic_*.e run_trimmomatic_*.o LOGS
rm run_trimmomatic_*

# check for any errors
cd LOGS
grep -i "Exited" *.o
grep -i "Successfully completed." *.o | wc -l
grep -i "error" *.e
# no errors
```

## Map trimmed reads to combined D. immitis/Wol/dog genome

### Nextflow pipeline 

Steve ran mapping-helminth module on the data.

```bash
mapping-helminth --input wgs.mapping.manifest_extra --reference /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/01_REF/reference_di_wol_dog.fa
```

## Extract reads that mapped to the *D. immitis* genome

If I mapped to the *D. immitis* and dog genomes separately, there could be reads that mapped to both genomes. To avoid this, I mapped to the combined D. immitis/dog genome. I can now extract the reads that mapped to only the *D. immitis* genome and use this for downstream analyses.


```bash
WORKING_DIR=/lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA
OUT_DIR=${WORKING_DIR}/03_ANALYSIS/02_MAP/EXTRACT/EXTRA
REF=${WORKING_DIR}/01_REF/dimmitis_WSI_2.2.fa
SAMPLE_LIST=${WORKING_DIR}/03_ANALYSIS/02_MAP/EXTRA/samples_extra.list

# Load modules
module load bsub.py/0.42.1
module load samtools/1.14--hb421002_0

# We've already indexed the reference previously.
BED=${WORKING_DIR}/01_REF/dimmitis_WSI_2.2.bed

# Extract reads that only mapped to D. immitis.
cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/02_MAP/EXTRA

n=1
while read SAMPLE; do
echo "samtools view --threads 4 --bam --with-header --target-file ${BED} /lustre/scratch125/pam/teams/team333/sd21/dirofilaria_immitis/POPGEN/NEWDATA_2024/results/${SAMPLE}/${SAMPLE}.bam > ${OUT_DIR}/${SAMPLE}_extract.bam" > run_extract_${SAMPLE}.tmp.job_${n};
let "n+=1";
done < ${SAMPLE_LIST}
# Should still be in sorted form
# -b flag makes sure the output is bam
# -h flag includes the header in SAM output

chmod a+x run_extract_*

#run
for i in run_extract_*; do
    bsub.py --threads 4 20 ${i} "./${i}";
done

mkdir LOGS
mv run_extract_*.e run_extract_*.o LOGS
rm run_extract_*

# check for any errors
cd LOGS
grep -i "Exited" *.o
grep -i "Successfully completed." *.o | wc -l
grep -i "error" *.e
# All ok


# I do not have to sort the bam file again, it should still be sorted.

## QC
# How many D. immitis reads were extracted?

WORKING_DIR=/lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA
OUT_DIR=${WORKING_DIR}/03_ANALYSIS/02_MAP/EXTRACT/EXTRA
SAMPLE_LIST=${WORKING_DIR}/03_ANALYSIS/02_MAP/EXTRA/samples_extra.list

module load samtools/1.14--hb421002_0
module load multiqc/1.17--pyhdfd78af_1

while read SAMPLE; do
samtools flagstat ${OUT_DIR}/${SAMPLE}_extract.bam > ${OUT_DIR}/${SAMPLE}_extract_flagstat.txt;
done < ${SAMPLE_LIST}

# multiqc to combine flagstat files for all samples
multiqc ./*flagstat.txt EXTRA/*flagstat.txt -o ./
```


## Index the extracted bam files

```bash
cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/02_MAP/EXTRA

n=1
while read SAMPLE; do
echo "samtools index ${OUT_DIR}/${SAMPLE}_extract.bam" > run_extract_index_${SAMPLE}.tmp.job_${n};
let "n+=1";
done < ${SAMPLE_LIST}

chmod a+x run_extract_index_*

#run
for i in run_extract_index_*; do
    bsub.py --threads 2 10 ${i} "./${i}";
done

# check for any errors
grep -i "Exited" *.o
grep -i "Successfully completed." *.o | wc -l
grep -i "error" *.e
# All ok

mv run_extract_index*.e run_extract_index*.o LOGS
rm run_extract_index*

# get extracted bam list for variant calling
cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/02_MAP/EXTRACT/EXTRA
for file in *.bam; do
echo "$(pwd)/$file" >> extra.bamlist;
done

ls *.bam > extra.bamlist
```