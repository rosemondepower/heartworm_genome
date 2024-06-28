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

# Choose ~ 1 sample per population and run
n=1
SAMPLE_LIST=/lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/03_ANALYSIS/05_MITO/GetOrganelle/samples.list
while IFS=$'\t' read -r SAMPLE_OLD SAMPLE_NEW; do
    echo "get_organelle_from_reads.py -1 /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/02_FASTQ/${SAMPLE_OLD}_*1.f*q.gz -2 /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/02_FASTQ/${SAMPLE_OLD}_*2.f*q.gz -t 4 -F animal_mt -o ${WORKING_DIR}/${SAMPLE_NEW} -R 30 --reduce-reads-for-coverage 1000 --max-reads 20000000 -w 95" > run_getorganelles_${SAMPLE_NEW}.tmp.job_${n};
    let "n+=1";
done < ${SAMPLE_LIST}

chmod a+x run_getorganelles_*

#run
for i in run_getorganelles_*; do
    bsub.py --threads 4 10 ${i} "./${i}";
done

# clean up
mv run_getorganelles_*.e run_getorganelles_*.o LOGS
rm run_getorganelles_*

cd LOGS
grep -i "Exited" *.o
grep -i "Successfully completed" *.o | wc -l
grep -i "Error" *.e
# All ok

# Add sample ID to sequences
cd ..
while IFS=, read -r sequence sequence_id; do
cp ${sequence} ${sequence_id}
done < sequence_id.csv

# Combine all samples into 1 fasta file
for file in *.fasta; do
    [ "$file" = "heartworm.fasta" ] && continue
    header="${file%.fasta}"
    second_line=$(sed -n '2p' "$file")
    echo ">${header}" >> heartworm.fasta
    echo "${second_line}" >> heartworm.fasta;
done
mitos

```

Questions
- should I align by MUSCLE or clustalW and does it matter?
- I have 120+ samples from around the world (show map) - should i subsample to make a tree?
- Where to put mutation rate in BEAST?  - mutation rate (C.elegans) = 2.7e-9 (https://doi.org/10.1073/pnas.0904895106)
- What nucleotide substitution model would I do?
-How many times do i run BEAST? So i physically run at least 2?

Notes/advice from Simon:
- dont need to do a subset based on location - can easily assemble 120 samples. can use the entire mitochondrial genome, align it, then run it through IQtree pretty quickly. If nodes are short, there isnt necessarily anything wrong they could just be closely related?
- Can split up based on genes if I want (but entire mtDNA is the first step). Can get annotated mtDNA from genbank, use genius to get certain genes, and do trees based on those.
- Don't need to run astral to get a species tree if it's all just 1 species.
- You should run BEAST multiple time and the likelihood scores should be similar

Notes from demo:
- probably shouldnt run across whole mito genome because there's diff rates of evolution, there's the D loop that has a lot of mutations
- best to do the genius method to pull out certain genes of itnerest.
- but you need an annotated genome to see where the genes are - can use web server MITOS to annotate genome -> you just click "nematode" etc and that's it
-Whether you use MUSCLE or clustalW doesnt really matter in most cases
-you would want to partition it based on codon position 1,2,3, esp if looking at protein coding genes. Substitution rate at 2nd position is far lower than the 3rd codon position bc it doesn't change the amino acid and get passed onto the next generation. Especially important if I'm just looking at the 1 species and there isn't much variation.
- need mitochondrial worm mutation rate (not genome-wide bc that's totally different). If we can't find anything email, and we can still find a number it just might be less confident.
- species model is usually birth0death that makes the most sense and everyone uses
- JModel test package find the best model for you and you dont have to fiddle with anything - leave as default
- there's a model selection package which spits out the maximum likelihood or bayes factor? Choose the one that's the highest and go with that one. Can run a few.
- Looking at tracer is the most important thing - checking that it stabilises, generally just look at the first few parameters, they dont look at all

- need outgroups in the other trees but dont need outgroup for BEAST if it's just the 1 species
- using mutation rate is best -put it in priors --> clockRate --> Normal distribution, enter the mean and st deviation (sigma) that the paper reports
- if you use mito-wide mutation rate then you need to use it on the entire mito. If you're only looking at cox1 then you need to use cox1-specific mutation rate.
- you can have cox 1 and other genes, partition it between the cox 1 and others, then apply the cox1 mutation rate as a prior, but don't have a prior for the others (it will use the cox1 relatively so estimate the other genes).



# Mitos

Annotating mitochondrial genome.

```bash
module load samtools/1.14--hb421002_0

cd /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/01_REF

# Get mito reference
samtools faidx dimmitis_WSI_2.2.fa "dirofilaria_immitis_chrMtDNA" > dimmitis_WSI_2.2_chrMtDNA.fa
# index it
samtools faidx dimmitis_WSI_2.2_chrMtDNA.fa
```
