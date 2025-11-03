# Assembling mitochondrial genomes

## Get Organelle

Assemble the mitochondrial genomes of all samples (exclude the MF samples).
- Downloaded get_organelle_from_reads.py v1.7.7.1 (https://github.com/Kinggerm/GetOrganelle)

```bash
# Install getorganelles
conda create --name getorganelle
conda activate getorganelle
conda install -c bioconda getorganelle
conda install -c bioconda matplotlib
conda list
# installed getorganelle v1.7.7.1
get_organelle_config.py -a all # download animal_mt

cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/05_MITO
mkdir GetOrganelle
WORKING_DIR=/lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/05_MITO/GetOrganelle
cd ${WORKING_DIR}

# Load modules
module load bsub.py/0.42.1

# Make job scripts
n=1
SAMPLE_LIST=/lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/05_MITO/GetOrganelle/samples.txt
while IFS=$'\t' read -r SAMPLE_OLD SAMPLE_NEW; do
    echo "get_organelle_from_reads.py -1 /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/02_FASTQ/${SAMPLE_OLD}_*1.f*q.gz -2 /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/02_FASTQ/${SAMPLE_OLD}_*2.f*q.gz -t 4 -F animal_mt -o ${WORKING_DIR}/${SAMPLE_NEW} -R 20 --reduce-reads-for-coverage 80" > run_getorganelles_${SAMPLE_NEW}.tmp.job_${n};
    let "n+=1";
done < ${SAMPLE_LIST}
## want 80-100x coverage

# extra samples
SAMPLE_EXTRA_LIST=/lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/05_MITO/GetOrganelle/samples_extra.txt
while IFS=$'\t' read -r SAMPLE_OLD SAMPLE_NEW; do
    echo "get_organelle_from_reads.py -1 /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/02_FASTQ/EXTRA/${SAMPLE_OLD}_*1.f*q.gz -2 /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/02_FASTQ/EXTRA/${SAMPLE_OLD}_*2.f*q.gz -t 4 -F animal_mt -o ${WORKING_DIR}/${SAMPLE_NEW} -R 20 --reduce-reads-for-coverage 80" > run_getorganelles_${SAMPLE_NEW}.tmp.job_${n};
    let "n+=1";
done < ${SAMPLE_EXTRA_LIST}

# Make job scripts executable
chmod a+x run_getorganelles_*

# Run GetOrganelle
for i in run_getorganelles_*; do
    bsub.py --threads 4 50 ${i} "./${i}";
done
```


## Re-running some samples

```bash
# ITA_NEA_AD_004 & ITA_NEA_AD_005 (D. repens) outgroup samples have low coverage, re-run these
rm -r ITA_NEA_AD_004 ITA_NEA_AD_005
SAMPLE_OUTGROUP_LIST=/lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/05_MITO/GetOrganelle/samples_outgroup.txt
while IFS=$'\t' read -r SAMPLE_OLD SAMPLE_NEW; do
    echo "get_organelle_from_reads.py -1 /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/02_FASTQ/${SAMPLE_OLD}_*1.f*q.gz -2 /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/02_FASTQ/${SAMPLE_OLD}_*2.f*q.gz -t 4 -F animal_mt -o ${WORKING_DIR}/${SAMPLE_NEW} -R 20 --reduce-reads-for-coverage 120" > run_getorganelles_${SAMPLE_NEW}_rerun.tmp.job_${n};
    let "n+=1";
done < ${SAMPLE_OUTGROUP_LIST}

# Make job scripts executable
chmod a+x run_getorganelles_*rerun*

# Run GetOrganelle
for i in run_getorganelles_*rerun*; do
    bsub.py --threads 4 20 ${i} "./${i}";
done
# better now


## Some GRC samples having issues, try re-running with adjusted parameters
rm -r GRC_XAN_AD_010 GRC_XAN_AD_011 GRC_XAN_AD_012
SAMPLE_RERUN_LIST=/lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/05_MITO/GetOrganelle/samples_rerun.txt
while IFS=$'\t' read -r SAMPLE_OLD SAMPLE_NEW; do
    echo "get_organelle_from_reads.py -1 /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/02_FASTQ/${SAMPLE_OLD}_*1.f*q.gz -2 /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/02_FASTQ/${SAMPLE_OLD}_*2.f*q.gz -t 8 -F animal_mt -o ${WORKING_DIR}/${SAMPLE_NEW} -R 20 --reduce-reads-for-coverage 500" > run_getorganelles_${SAMPLE_NEW}_rerun.tmp.job_${n};
    let "n+=1";
done < ${SAMPLE_RERUN_LIST}

# Make job scripts executable
chmod a+x run_getorganelles_*rerun*

# Run GetOrganelle
for i in run_getorganelles_GRC*rerun*; do
    bsub.py --threads 8 50 ${i} "./${i}";
done
# still low coverage and some warnings


# Get average coverage for each sample
output_file="coverage_summary.txt"
for sample_dir in */; do
    sample_dir=${sample_dir%/}
    if [ -f "$sample_dir/get_org.log.txt" ]; then
        coverage_line=$(grep "Average animal_mt base-coverage =" "$sample_dir/get_org.log.txt")
        if [ ! -z "$coverage_line" ]; then
            echo "$sample_dir" >> $output_file
            echo "$coverage_line" >> $output_file
            echo "" >> $output_file 
        fi
    fi
done

# Clean up
mv run_getorganelles_*.e run_getorganelles_*.o LOGS
rm run_getorganelles_*

# Check
cd LOGS
grep -i "Exited" *.o
grep -i "Successfully completed" *.o | wc -l
grep -i "Error" *.e
# All ok

# Add sample IDs to filenames and edit header of fasta file
cd ASSEMBLIES
while IFS=$'\t' read -r sequence sequence_id; do
cp ../${sequence} ${sequence_id}
new_header=$(basename "${sequence_id}" .fasta)
sed -i "1s/.*/>${new_header}/" "${sequence_id}"
done < ../sequence_id.txt

"
cp: cannot stat '../GRC_XAN_AD_010/animal_mt*.graph1.1.path_sequence.fasta': No such file or directory
sed: can't read GRC_XAN_AD_010.fasta: No such file or directory
cp: cannot stat '../GRC_XAN_AD_011/animal_mt*.graph1.1.path_sequence.fasta': No such file or directory        
sed: can't read GRC_XAN_AD_011.fasta: No such file or directory
cp: cannot stat '../GRC_XAN_AD_012/animal_mt*.graph1.1.path_sequence.fasta': No such file or directory        
sed: can't read GRC_XAN_AD_012.fasta: No such file or directory
"

# do the GRC ones now that they're finished running
cd ASSEMBLIES
while IFS=$'\t' read -r sequence sequence_id; do
cp ../${sequence} ${sequence_id}
new_header=$(basename "${sequence_id}" .fasta)
sed -i "1s/.*/>${new_header}/" "${sequence_id}"
done < ../sequence_id_rerun.txt


## CLC

**Mitochondrial references obtained**

D. immitis
- (NC_005305): https://doi.org/10.1017/S0031182003003275
- (mDi_Athens_2.1): https://doi.org/10.1096%2Ffj.12-205096
- (mDi_Pavia_2.1): https://doi.org/10.1096%2Ffj.12-205096

D. sp. 'Thailand II'
- (MH823370-72): https://doi.org/10.1111/tbed.13033

D. sp. 'hongkongensis'
- (KX265050): https://doi.org/10.1371/journal.pntd.0005028

D. repens
- (KX265047-49): https://doi.org/10.1371/journal.pntd.0005028
- (KR071802): unpublished

O. volvulus
- (KT599912): https://doi.org/10.1590/0074-02760150350
- (AF015193): https://doi.org/10.1016/S0166-6851(98)00102-9
- (AP017695): unpublished

O. ochengi
- (KX181289): unpublished
- (PRJEB1809): https://doi.org/10.1101/gr.138420.112 ### didnt download this
- (KX181290): unpublished
- (AP017694): unpublished
- (AP017693): unpublished

O. flexuosa
- (AP017692): unpublished
- (HQ214004): https://doi.org/10.1186/1471-2164-13-145

All sequences (including ref seqs) were manually edited to start at cox1.
First, all my D. immitis sequences (including the D. immitis ref seqs) were aligned and cleaned. Then, the outgroup samples (including outgroup ref seqs) were added to the alignment, which was manually inspected and cleaned once more.

Genes were extracted and separate alignments created. Now, I will make a tree based on the entire mtDNA genome, and the individual genes.

## IQ-TREE

```bash
# load modules
module load bsub.py/0.42.1
module load iqtree/1.6.12--he513fc3_1

cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/05_MITO/IQTREE/LOGS

# run iqtree on whole mtDNA alignment
bsub.py 4 iqtree_mito "iqtree -s ../INPUT/Mito_align_cleaned_outgroups.fa -bb 1000 -redo"

## At the beginning of each run, iqtree performs a composition chi-square test for every sequence in the alignment. The purpose is to test for homogeneity of character composition (e.g., nucleotide for DNA, amino-acid for protein sequences). A sequence is denoted failed if its character composition significantly deviates from the average composition of the alignment. This test should be regarded as an explorative tool which might help to nail down problems in a dataset. One would typically not remove failing sequences by default. But if the tree shows unexpected topology the test might point in direction of the origin of the problem.

## quite a few samples failed the chi-square test, but all of these samples were the outgroups.

# run iqtree on ATP6
bsub.py 2 iqtree_ATP6 "iqtree -s ../INPUT/Mito_align_cleaned_outgroups_ATP6.fa -bb 1000 -redo"

# run iqtree on COX1
bsub.py 2 iqtree_COX1 "iqtree -s ../INPUT/Mito_align_cleaned_outgroups_COX1.fa -bb 1000 -redo"

# run iqtree on COX2
bsub.py 2 iqtree_COX2 "iqtree -s ../INPUT/Mito_align_cleaned_outgroups_COX2.fa -bb 1000 -redo"

# run iqtree on COX3
bsub.py 2 iqtree_COX3 "iqtree -s ../INPUT/Mito_align_cleaned_outgroups_COX3.fa -bb 1000 -redo"

# run iqtree on CYTB
bsub.py 2 iqtree_CYTB "iqtree -s ../INPUT/Mito_align_cleaned_outgroups_CYTB.fa -bb 1000 -redo"

# run iqtree on ND1
bsub.py 2 iqtree_ND1 "iqtree -s ../INPUT/Mito_align_cleaned_outgroups_ND1.fa -bb 1000 -redo"

# run iqtree on ND2
bsub.py 2 iqtree_ND2 "iqtree -s ../INPUT/Mito_align_cleaned_outgroups_ND2.fa -bb 1000 -redo"

# run iqtree on ND3
bsub.py 2 iqtree_ND3 "iqtree -s ../INPUT/Mito_align_cleaned_outgroups_ND3.fa -bb 1000 -redo"

# run iqtree on ND4
bsub.py 2 iqtree_ND4 "iqtree -s ../INPUT/Mito_align_cleaned_outgroups_ND4.fa -bb 1000 -redo"

# run iqtree on ND4L
bsub.py 2 iqtree_ND4L "iqtree -s ../INPUT/Mito_align_cleaned_outgroups_ND4L.fa -bb 1000 -redo"

# run iqtree on ND5
bsub.py 2 iqtree_ND5 "iqtree -s ../INPUT/Mito_align_cleaned_outgroups_ND5.fa -bb 1000 -redo"

# run iqtree on ND6
bsub.py 2 iqtree_ND6 "iqtree -s ../INPUT/Mito_align_cleaned_outgroups_ND6.fa -bb 1000 -redo"
```

Visualised trees in FigTree & MEGA.