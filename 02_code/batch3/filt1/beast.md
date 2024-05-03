# BEAST

## Make tree using BEAST with nuclear data
I have played around with BEAST using the whole mitochondrial genomes. Now, I want to do the same with my nuclear data. Choose some representative nuclear genes.

Try "5.8S-ITS2-28S rRNA"
- Primer F: 5′ AGTGCGAATTGCAGACGCATTGAG 3′
- Primer R: 5′ AGCGGGTAATCACGACTGAGTTGA 3′
- ~542 bp
- https://doi.org/10.1186/s12917-023-03803-0

```bash
# First, I need to get the genome positions of this region

# Get primers in FASTA format
File: "primers1.fasta"
>Forward_Primer
AGTGCGAATTGCAGACGCATTGAG
>Reverse_Primer
AGCGGGTAATCACGACTGAGTTGA

# Run blastn
module load blast+/2.12.0
cd /scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping/BEAST/nuclear
blastn -query primers1.fasta -subject /project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_ALL/batch1/analysis/mapping/dimmitis_WSI_2.2.fa -out primer1_results.txt
```
Result: no blast matches


Try "ND1 rDNA gene"
- Primer F: 5′ ATGGCCTAAAGGGGCGTAAG 3′
- Primer R: 5′ TAGACTGAAGCCCCTGATGC 3′
- https://doi.org/10.4269/ajtmh.13-0579

```bash
# First, I need to get the genome positions of this region

# Get primers in FASTA format
File: "primers2.fasta"
>Forward_Primer
ATGGCCTAAAGGGGCGTAAG
>Reverse_Primer
TAGACTGAAGCCCCTGATGC

# Run blastn
blastn -query primers2.fasta -subject /project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_ALL/batch1/analysis/mapping/dimmitis_WSI_2.2.fa -out primer2_results.txt
```
Result: no blast matches


Try "5s rDNA intergenic spacer region"
- Primer F: 5′ GTTAAGCAACGTTGGGCCTGG 3′
- Primer R: 5′ TTGACAGATCGGACGAGATG 3′
- https://doi.org/10.1645/GE-2208.1

```bash
# First, I need to get the genome positions of this region

# Get primers in FASTA format
File: "primers3.fasta"
>Forward_Primer
GTTAAGCAACGTTGGGCCTGG
>Reverse_Primer
TTGACAGATCGGACGAGATG

# Run blastn
blastn -word_size 7 -query primers3.fasta -subject /project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_ALL/batch1/analysis/mapping/dimmitis_WSI_2.2.fa -out primer3_results.txt
```

This doesn't seem to be working, perhaps the primer sequences are too sho9rt for blast to find. Try bowtie instead.

```bash
module load bowtie2/2.4.4
bowtie2-build /project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_ALL/batch1/analysis/mapping/dimmitis_WSI_2.2.fa dimmitis_WSI_2.2_index
bowtie2 -x dimmitis_WSI_2.2_index -f cox1.fasta -S cox1_alignment.sam

```

