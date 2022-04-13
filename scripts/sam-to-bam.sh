#!/bin/bash

#SBATCH -p batch
#SBATCH -J sam2bam
#SBATCH -n 8

# load the blast module
module load samtools

mkdir -p /home/emurungi/gitau/marion/ACE/bam

cd /home/emurungi/gitau/marion/ACE/sam

SAM_DIR="/home/emurungi/gitau/marion/ACE/sam"

BAM_DIR="/home/emurungi/gitau/marion/ACE/bam"


# convert sam file to sorted bam files
for sam_file in ${SAM_DIR}/*.sam; do
        sam_file_name=$(basename "$sam_file" .sam)
        samtools view -S -b $sam_file | \
        samtools sort -n -o ${BAM_DIR}/${sam_file_name}.sorted.bam
done
