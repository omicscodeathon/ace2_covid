#!/bin/bash

#SBATCH -p batch
#SBATCH -J alignment
#SBATCH -n 8

# load the blast module
module load hisat2

mkdir -p /home/emurungi/gitau/marion/ACE/sam

cd /home/emurungi/gitau/marion/ACE/trimmed

FNA_DIR="/home/emurungi/gitau/marion/raw"

SAM_DIR="/home/emurungi/gitau/marion/ACE/sam"

#Align fastqs to indexed reference genome
SAMPLES="SRR15042704 SRR15042705 SRR15042706 SRR15042707"

for SAMPLE in $SAMPLES; do

        hisat2 \
                 -x ${FNA_DIR}/GCF_000001405.39_GRCh38.p13_genomic.fna_index_hisat2 \
                 -1 ${SAMPLE}_1.trimmed.fastq \
                 -2 ${SAMPLE}_2.trimmed.fastq \
                 -S ${SAM_DIR}/${SAMPLE}.sam \
                 -p 6 \
                --summary-file ${SAMPLE}.txt \
                --new-summary

done

