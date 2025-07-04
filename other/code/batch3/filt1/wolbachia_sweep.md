# Investigate a Wolbachia selective sweep in D. immitis

First, I will download the Wolbachia data from: (https://doi.org/10.1645/GE-2208.1) and re-create their phylogenetic tree based on the ftsZ gene of Wolbachia. 

Fetch these ftsz sequences from GenBank:
GQ217523
AY523519.1
AJ010273.1
AJ010268.1
AJ276501.1
AJ010266.1
AJ010267.1
AJ415416.1
DQ842341.1
AY583315.1
DQ093835.1
AJ010271.1
AJ628414.1
AJ292345.2
AJ344216.1
DQ842340.1
CP001391.1
DQ842337.1
DQ842309.1


Where is the ftsz region located in my Wolbachia data? Need to extract it for each sample.
- get wolbachia reference genome by itself
- find ftsz region




```bash
module load samtools/1.9
module load bcftools/1.11
module load tabix/0.2.6

DIR=/scratch/RDS-FSC-Heartworm_MLR-RW/batch3/fastq/Wol/my_data
REF_DIR=/project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_ALL/batch1/analysis/mapping
VCF_DIR=/scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping/filter/FINAL_SETS
cd $DIR

# get Wolbachia from reference sequence
samtools faidx $REF_DIR/dimmitis_WSI_2.2.fa "dirofilaria_immitis_chrWb" > $REF_DIR/dimmitis_WSI_2.2_chrWb.fa
# index it
samtools faidx $REF_DIR/dimmitis_WSI_2.2_chrWb.fa
```

Now transfer this Wolbachia reference sequence into CLC. 
The primers used for the ftsZ region: "(ftsZunif 5–3 GG(CT)AA(AG) GGTGC(AG)GCAGAAGA; ftsZunir 5–3 ATC(AG) AT(AG)CCAGTTGCAAG) were designed based on an alignment of all arthropod and filarial Wolbachia ftsZ sequences. These primers enabled amplification of a 775bp fragment (737 bp excluding primers) from F. candida and each termite specimen." from https://doi.org/10.1093/oxfordjournals.molbev.a004087 (Nathan Lo's paper).

Nucleotides in brackets mean that it can anneal to either one. Don't know which one works for D. immitis, so tested all 8 combinations (forward primer) and all 4 combinations (reverse primer) to see which one annealed without any mismatches. Result: ftsZunif_TAA and ftsZunir_AA binded to the D. immitis Wolbachia ref genome with 0 mismatches, producing a ~714 bp fragment which seems sort of similar to what they mentioned for termites.

Now I can extract that particular region from the reference and use this to extract the region in my data.


```bash
# Extracting the specific ftsZ region in the reference
samtools faidx $REF_DIR/dimmitis_WSI_2.2.fa "dirofilaria_immitis_chrWb:851873-852624" > $REF_DIR/dimmitis_WSI_2.2_chrWb_851873-852624.fa
# index it
samtools faidx $REF_DIR/dimmitis_WSI_2.2_chrWb_851873-852624.fa


# zip my Wol vcf file
bgzip -c $VCF_DIR/wb_samples3x_missing0.9.recode.RENAMED.vcf > $VCF_DIR/wb_samples3x_missing0.9.recode.RENAMED.vcf.gz
# index it

# try 1 sample
bcftools consensus -f $REF_DIR/dimmitis_WSI_2.2_chrWb_851873-852624.fa -s "JS6349" $VCF_DIR/wb_samples3x_missing0.9.recode.RENAMED.vcf.gz > JS6349_wb_851873-852624.fasta
# this worked, now loop through each sample
```

```bash
#!/bin/bash

# PBS directives 
#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N ftsz_seq
#PBS -l select=1:ncpus=2:mem=20GB
#PBS -l walltime=00:30:00
#PBS -m e
#PBS -q defaultQ
#PBS -o ftsz_seq.txt
#PBS -M rosemonde.power@sydney.edu.au

# Submit job
## qsub ../ftsz_seq.sh

module load bcftools/1.11

DIR=/scratch/RDS-FSC-Heartworm_MLR-RW/batch3/fastq/Wol/my_data
WOL_REF=/project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_ALL/batch1/analysis/mapping/dimmitis_WSI_2.2_chrWb_851873-852624.fa
WOL_VCF=/scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping/filter/FINAL_SETS/wb_samples3x_missing0.9.recode.RENAMED.vcf.gz

cd ${DIR}

bcftools query -l ${WOL_VCF} | while read -r SAMPLE; do
echo ">${SAMPLE}" > ${DIR}/${SAMPLE}_Wb_851873-852624.fasta 
bcftools consensus -f ${WOL_REF} -s ${SAMPLE} ${WOL_VCF} | tail -n +2 >> ${DIR}/${SAMPLE}_Wb_851873-852624.fasta 
echo "FASTA sequence generated for ${SAMPLE}"
done

# write sample name to first line
# append sequence starting from 2nd line
```

Transfer my data and all the Genbank sequences into CLC and align to the reference. Make a tree to see whether all my Wolbachia are all just 1 strain at that ftsz region (they should be, there were no variants in my data when I was making consensus sequences).



## Get whole Wolbachia genome from my data

```bash
#!/bin/bash

# PBS directives 
#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N wol_seq
#PBS -l select=1:ncpus=2:mem=20GB
#PBS -l walltime=00:30:00
#PBS -m e
#PBS -q defaultQ
#PBS -o wol_seq.txt
#PBS -M rosemonde.power@sydney.edu.au

# Submit job
## qsub ../wol_seq.sh

module load bcftools/1.11

DIR=/scratch/RDS-FSC-Heartworm_MLR-RW/batch3/fastq/Wol/my_data/whole
WOL_REF=/project/RDS-FSC-Heartworm_MLR-RW/HW_WGS_ALL/batch1/analysis/mapping/dimmitis_WSI_2.2_chrWb.fa
WOL_VCF=/scratch/RDS-FSC-Heartworm_MLR-RW/batch3/analysis/mapping/filter/FINAL_SETS/wb_samples3x_missing0.9.recode.RENAMED.vcf.gz

cd ${DIR}

bcftools query -l ${WOL_VCF} | while read -r SAMPLE; do
echo ">${SAMPLE}" > ${DIR}/${SAMPLE}_Wb.fasta 
bcftools consensus -f ${WOL_REF} -s ${SAMPLE} ${WOL_VCF} | tail -n +2 >> ${DIR}/${SAMPLE}_Wb.fasta 
echo "FASTA sequence generated for ${SAMPLE}"
done

# write sample name to first line
# append sequence starting from 2nd line


```