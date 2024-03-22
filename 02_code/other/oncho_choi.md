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

