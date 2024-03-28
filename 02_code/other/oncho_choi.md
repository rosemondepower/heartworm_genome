# Nuclear, mtDNA and Wolbachia of Onchocerca volvulus

Is low number of SNPs specific to Dirofilaria immitis, or is it common amongst filarial worms?


O. volvulus data from:
- https://doi.org/10.1038%2Fnmicrobiol.2016.207
- 27 worms from West Africa, Ecuador and Uganda
- They got 1.3M SNPs in Ov, 167 SNPs in mtDNA, 772 SNPs in wOv

## Download data

```bash
cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/OTHER/ONCHO_CHOI/01_SEQ

module load bsub.py/0.42.1
module load cellgen/sratoolkit/3.0.7


while read ACCESSION; do
echo -e "fastq-dump --gzip ${ACCESSION}" > ../${n}.run_oncho_choi.tmp.${ACCESSION};
let "n+=1";
done < ../accession.list
chmod a+x ../*.run_oncho_choi.tmp.*

# run
for i in ../*.run_oncho_choi.tmp.*; do
bsub.py --queue normal --threads 4 10 oncho_choi_data "${i}";
done
```

## Get reference sequences ready

- Human: used genome assembly GRCh38.p14
- Oncho: used Oncho genome on WormBase ParaSite: Bioproject PRJEB513

```bash
#!/bin/bash

# PBS directives 
#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N mapping_prep
#PBS -l select=1:ncpus=1:mem=15GB
#PBS -l walltime=01:30:00
#PBS -m e
#PBS -q defaultQ
#PBS -o mapping_prep.txt

# qsub ../mapping_prep.pbs

# Set working directory
cd /project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_ALL/data/analysis/mapping

# Load modules
module load bwa/0.7.17

# Combine the 2 references
cat dimmitis_WSI_2.2.fa GCA_014441545.1_ROS_Cfam_1.0_genomic.fna > reference_di_wol_dog.fa 
# Already did this step in the test run. Can just copy over the joined reference file.
# cp /project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_test/data/analysis/mapping/reference_di_wol_dog.fa .

# Also copy over the D. immitis/Wol reference
# cp /project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_test/data/analysis/mapping/dimmitis_WSI_2.2.fa .
# Already did these steps prior.

# index reference sequence
bwa index reference_di_wol_dog.fa

# Perform mapping, sam-to-bam conversion, filtering, and indexing. Need to have separate table with the relevant sample names I want. Then need to use variables to replace the sample names. This will reduce time/effort and minimise typos. I could make a loop, but since I will be running it on the Artemis server as a job, it may be best to write out a script.
```

