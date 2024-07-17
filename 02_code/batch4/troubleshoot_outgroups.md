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

Needed to run a bit longer than the 12h limit.

![](images/Drepens_FAKEREADS_mapping.PNG)

~18% of reads mapped, which is still pretty low. So there is clearly low mapping due to divergence.




- Can use raw reads from public D. repens data & see how well that maps to D. immitis

Get the raw data for the D. repsn ref genome. It's pacbio single reads.

```bash
# load modules
module load cellgen/sratoolkit/3.0.10

# download raw data that made the ref genome
chmod a+x repens_sra_download.sh
bsub.py 20 repens_sra "repens_sra_download.sh"

# repens_sra_download.sh:
cd DATA
fastq-dump --origfmt --gzip SRR8742586
fastq-dump --origfmt --gzip SRR9613458
fastq-dump --origfmt --gzip SRR9613459
fastq-dump --origfmt --gzip SRR9613460
fastq-dump --origfmt --gzip SRR9613461
fastq-dump --origfmt --gzip SRR9613462
fastq-dump --origfmt --gzip SRR9613463
fastq-dump --origfmt --gzip SRR9613464
fastq-dump --origfmt --gzip SRR9613465
```


## FastQC

We want to get some stats on the raw D. repens data.

```bash
# set variables
WORKING_DIR=/lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA
OUT_DIR=${WORKING_DIR}/03_ANALYSIS/05_ANALYSIS/OUTGROUPS/DREPENS/FASTQC
cd ${WORKING_DIR}/03_ANALYSIS/05_ANALYSIS/OUTGROUPS/DREPENS

# make a list of sample names
for file in ${WORKING_DIR}/DATA/*.fastq.gz; do
  basename ${file}
done > samples.list

# load modules
module load fastqc/0.12.1--hdfd78af_0
module load bsub.py/0.42.1

# set up run files
n=1
while read SAMPLE; do
SAMPLE_BASE=$(basename "${SAMPLE}" .fastq.gz)
echo -e "fastqc -t 4 -o ${OUT_DIR} ${WORKING_DIR}/03_ANALYSIS/05_ANALYSIS/OUTGROUPS/DREPENS/DATA/${SAMPLE}" > run_fastqc_raw_${SAMPLE_BASE}.tmp.job_${n};
let "n+=1";
done < samples.list

chmod a+x run_fastqc_raw*

#run
for i in run_fastqc_raw*; do
    bsub.py --threads 4 20 ${i} "./${i}";
done

# clean up
mkdir LOGS
mv run_fastqc_raw_*.e run_fastqc_raw_*.o LOGS
rm run_fastqc_raw_*

# check for any errors
cd LOGS
grep -i "Exited" *.o
grep -i "Successfully completed" *.o | wc -l
grep -i "Error" *.e
# All ok

# multiqc
# Load module
module load multiqc/1.17--pyhdfd78af_1
module load bsub.py/0.42.1
bsub.py 4 multiqc_raw "multiqc ${OUT_DIR} -o ${OUT_DIR}"
```

There's only 0.2 M seqs per sample, but the median read length is much longer. They're all from the same BioSample, so you could merge them together.

```bash
module load bsub.py/0.42.1

cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/05_ANALYSIS/OUTGROUPS/DREPENS/DATA

bsub.py 10 repens_merge "repens_merge.sh"

zcat DATA/SRR8742586.fastq.gz DATA/SRR9613458.fastq.gz DATA/SRR9613459.fastq.gz DATA/SRR9613460.fastq.gz DATA/SRR9613461.fastq.gz DATA/SRR9613462.fastq.gz DATA/SRR9613463.fastq.gz DATA/SRR9613464.fastq.gz DATA/SRR9613465.fastq.gz | gzip > DATA/Drepens_public.fastq.gz

```


## Trimmomatic

```bash
# set variables
WORKING_DIR=/lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/05_ANALYSIS/OUTGROUPS/DREPENS
IN_DIR=${WORKING_DIR}/DATA
cd ${WORKING_DIR}
OUT_DIR=${WORKING_DIR}/TRIMMOMATIC

# load modules
module load trimmomatic/0.39--1

# set up run files
n=1
while read SAMPLE; do

SAMPLE_BASE=$(basename "${SAMPLE}" .fastq.gz)

JOB_FILE="run_trimmomatic_${SAMPLE_BASE}.tmp.job_${n}"

echo -e "trimmomatic SE \
-threads 1 \
-phred33 \
${IN_DIR}/${SAMPLE} \
${OUT_DIR}/${SAMPLE_BASE}_trim.fq.gz \
SLIDINGWINDOW:10:20 MINLEN:50" > ${JOB_FILE};
let "n+=1";
done < samples.list

chmod a+x run_trimmomatic*

#run
for i in run_trimmomatic*; do
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