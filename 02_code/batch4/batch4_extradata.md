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

# Check that the files merged correctly

# F reads
cd ..
mkdir COUNT
cd COUNT
# Raw data files
for f in /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/02_FASTQ/EXTRA/X201SC24040556-Z01-F003/X201SC24040556-Z01-F003/01.RawData/JS6766/*_1.fq.gz; do echo $f;zcat $f|wc -l ; done > count_raw_1.txt

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
cd ../COUNT

# Raw data files
for f in /lustre/scratch125/pam/teams/team333/rp24/DIRO/DATA/02_FASTQ/EXTRA/X201SC24040556-Z01-F003/X201SC24040556-Z01-F003/01.RawData/JS6766/*_2.fq.gz; do echo $f;zcat $f|wc -l ; done > count_raw_2.txt

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
