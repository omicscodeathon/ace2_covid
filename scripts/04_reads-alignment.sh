#!/usr/bin/bash

cd /home/nyamarim/raw/ref

# load the blast module
module load hisat2

mkdir -p /srv/data/my_shared_data_folder/ace2covid/results/sam/

cd /srv/data/my_shared_data_folder/ace2covid/raw/results/trimmed/

FNA_DIR="/srv/data/my_shared_data_folder/ace2covid/raw/ref/"

SAM_DIR="/srv/data/my_shared_data_folder/ace2covid/results/sam/"

#Align fastqs to indexed reference genome
SAMPLES="SRR13081305 SRR13081307 SRR13081308 SRR13081315 SRR13081350 SRR13081368 SRR15042705 SRR15042706 SRR15042707 SRR15042708 SRR15042710 SRR15042714 SRR15042716 SRR15042724 SRR15042728 SRR15042729 SRR15042749 SRR15042750 SRR15042751 SRR15042752 SRR15042753 SRR15042754 SRR15042755 SRR15042756 SRR15042757 SRR15042758 SRR15042759 SRR15042760 SRR15042761 SRR15042762 SRR15042763 SRR15042764 SRR15042765 SRR15042766 SRR15042767 SRR15042768 SRR15042769 SRR15042770 SRR15042771"

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

