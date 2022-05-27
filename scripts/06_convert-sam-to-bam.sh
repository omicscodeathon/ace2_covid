#!/bin/bash

mkdir -p /srv/data/my_shared_data_folder/ace2covid/results/bam

cd /srv/data/my_shared_data_folder/ace2covid/results/sam

SAM_DIR="/srv/data/my_shared_data_folder/ace2covid/results/sam"

BAM_DIR="/srv/data/my_shared_data_folder/ace2covid/results/bam"


# convert sam file to sorted bam files
for sam_file in ${SAM_DIR}/*.sam; do
        sam_file_name=$(basename "$sam_file" .sam)
        samtools view -S -b $sam_file | \
        samtools sort -n -o ${BAM_DIR}/${sam_file_name}.sorted.bam
done
