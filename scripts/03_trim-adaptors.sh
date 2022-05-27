#!/usr/bin/bash

cd /srv/data/my_shared_data_folder/ace2covid/raw


mkdir -p /srv/data/my_shared_data_folder/ace2covid/results/trimmed


FASTQ_DIR="/srv/data/my_shared_data_folder/ace2covid/raw"

OUT_DIR="/srv/data/my_shared_data_folder/ace2covid/results/trimmed"

SAMPLES="SRR13081305 SRR13081307 SRR13081308 SRR13081315 SRR13081350 SRR13081368 SRR15042705 SRR15042706 SRR15042707 SRR15042708 SRR15042710 SRR15042714 SRR15042716 SRR15042724 SRR15042728 SRR15042729 SRR15042749 SRR15042750 SRR15042751 SRR15042752 SRR15042753 SRR15042754 SRR15042755 SRR15042756 SRR15042757 SRR15042758 SRR15042759 SRR15042760 SRR15042761 SRR15042762 SRR15042763 SRR15042764 SRR15042765 SRR15042766 SRR15042767 SRR15042768 SRR15042769 SRR15042770 SRR15042771"

for SAMPLE in $SAMPLES; do

cutadapt -a AGATCGGAAGAG -A AGATCGGAAGAG -m15 -o ${OUT_DIR}/${SAMPLE}_1.trimmed.fastq -p ${OUT_DIR}/${SAMPLE}_2.trimmed.fastq ${FASTQ_DIR}/${SAMPLE}_1.fastq ${FASTQ_DIR}/${SAMPLE}_2.fastq

done


