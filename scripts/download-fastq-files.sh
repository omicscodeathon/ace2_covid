#!/usr/bin/bash

#SBATCH -p batch
#SBATCH -J ascp
#SBATCH -n 8

cd /home/emurungi/gitau/marion/ACE

fasterq-dump --split-files SRR15042704
fasterq-dump --split-files SRR15042705
fasterq-dump --split-files SRR15042706
fasterq-dump --split-files SRR15042707

