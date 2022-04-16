#!/bin/bash

#Read count quantification
#Script to counts the number of reads aligned to the reference genome using HTSeq.

#module load htseq

cd /srv/data/my_shared_data_folder/ace2covid/results/feature-counts

# process paired-end data

BAM_DIR="/srv/data/my_shared_data_folder/ace2covid/results/bam-files"

bam_file_name=$(basename "$BAM_DIR" .sorted.bam)
GTF_FILE=/srv/data/my_shared_data_folder/ace2covid/data/ref-index/GCF_000001405.39_GRCh38.p13_genomic.gtf

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

