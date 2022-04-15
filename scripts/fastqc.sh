#!/bin/bash

cd /home/nyamarim/qualitycontrol

#Script to generate FastQC reports using the FastQC tool

# fastq files directory
FASTQ_DIR="/home/nyamarim/data/healthy"

# fastqc reports directory
REPORT_DIR="/home/nyamarim/qualitycontrol"

for file in $FASTQ_DIR/*.fastq; do
   fastqc ${file} -o ${REPORT_DIR}
done

