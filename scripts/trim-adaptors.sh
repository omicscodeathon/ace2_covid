#!/bin/bash

#SBATCH -p batch
#SBATCH -J trimming
#SBATCH -n 8

module load cutadapt

#Script to trim adaptors


mkdir -p /home/emurungi/gitau/marion/ACE/trimmed


FASTQ_DIR="/home/emurungi/gitau/marion/ACE"

OUT_DIR="/home/emurungi/gitau/marion/ACE/trimmed"


SAMPLES="SRR15042704 SRR15042705 SRR15042706 SRR15042707"

for SAMPLE in $SAMPLES; do

cutadapt -a AGATCGGAAGAG -A AGATCGGAAGAG -m15 -o ${OUT_DIR}/${SAMPLE}_1.trimmed.fastq -p ${OUT_DIR}/${SAMPLE}_2.trimmed.fastq ${FASTQ_DIR}/${SAMPLE}_1.fastq ${FASTQ_DIR}/${SAMPLE}_2.fastq

done


