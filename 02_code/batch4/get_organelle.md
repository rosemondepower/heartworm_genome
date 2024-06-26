# Get Organelle to get consensus mitochondrial genomes

```bash
# Install getorganelles
conda create --name getorganelle
conda activate getorganelle
conda install -c bioconda getorganelle
conda list
# installed getorganelle v1.7.7.1
get_organelle_config.py -a all # download animal_mt

cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/05_MITO
mkdir GetOrganelle
WORKING_DIR=/lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/05_MITO/GetOrganelle
cd ${WORKING_DIR}

# Run on 1 test sample
module help bsub.py/0.42.1
bsub.py 4 getorganelle_test "get_organelle_from_reads.py -1 /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/02_FASTQ/JS6277_1.fq.gz -2 /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/02_FASTQ/JS6277_2.fq.gz -F animal_mt -o ${WORKING_DIR}/JS6277-mitogenome1 -R 10 --max-reads 8000000"

# Choose ~ 1 sample per population and run
n=1
SAMPLE_LIST=/lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/05_MITO/GetOrganelle/samples.list
while IFS=$'\t' read -r SAMPLE_OLD SAMPLE_NEW; do
    echo "get_organelle_from_reads.py -1 /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/02_FASTQ/${SAMPLE_OLD}_1.fq.gz -2 /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/02_FASTQ/${SAMPLE_OLD}_2.fq.gz -F animal_mt -o ${WORKING_DIR}/${SAMPLE_NEW} -R 10 --max-reads 8000000" > run_getorganelles_${SAMPLE_NEW}.tmp.job_${n};
    let "n+=1";
done < ${SAMPLE_LIST}

chmod a+x run_getorganelles_*

#run
for i in run_getorganelles_*; do
    bsub.py --threads 4 10 ${i} "./${i}";
done

```