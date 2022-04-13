#!/bin/bash

#SBATCH -p batch
#SBATCH -J counts
#SBATCH -n 8

#Read count quantification
#Script to counts the number of reads aligned to the reference genome using HTSeq.

module load htseq

cd /home/emurungi/gitau/marion/raw

# process paired-end data

BAM_DIR="/home/emurungi/gitau/marion/ACE/bam"

bam_file_name=$(basename "$BAM_DIR" .sorted.bam)
GTF_FILE=GCF_000001405.39_GRCh38.p13_genomic.gtf

for bam_file in ${BAM_DIR}/*.sorted.bam; do

        htseq-count \
            -f bam \
            -r pos \
            -s no \
            -t exon \
            -i gene \
            $bam_file \
            $GTF_FILE \
            > ${bam_file}.counts.txt

done

