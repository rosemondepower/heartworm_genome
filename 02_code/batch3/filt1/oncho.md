# Add in Onchocerca volvulus data as an outgroup

I would ideally like to use my own D. repens data once I get it, but for now I will download random public data from Genbank (SRR16526445) - note that it is not published and may not be reliable. However, this is just for testing purposes for the fstats.

## Download fastq

```bash
#!/bin/bash

#PBS -P RDS-FSC-Heartworm_MLR-RW
#PBS -N fastq
#PBS -l select=1:ncpus=1:mem=30GB 
#PBS -l walltime=10:00:00
#PBS -q defaultQ
#PBS -o fastq.txt 

module load sratoolkit/3.0.3

cd /scratch/RDS-FSC-Heartworm_MLR-RW/Oncho/fastq

fastq-dump --split-files --origfmt --gzip SRR16526445
```