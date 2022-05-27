#!/bin/bash

cd /srv/data/my_shared_data_folder/ace2covid/results/qualitycontrol

#Script to generate FastQC reports using the FastQC tool

# fastq files directory
FASTQ_DIR="/srv/data/my_shared_data_folder/ace2covid/raw"

# fastqc reports directory
REPORT_DIR="/srv/data/my_shared_data_folder/ace2covid/results/qualitycontrol"

for file in $FASTQ_DIR/*.fastq; do
   fastqc ${file} -o ${REPORT_DIR}
done

