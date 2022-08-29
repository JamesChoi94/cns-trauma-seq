#!/bin/bash
#BSUB -J cns-trauma-seq_download-data-curl
#BSUB -P lemmon
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -W 4:00
#BSUB -q general
#BSUB -n 8
#BSUB -B
#BSUB -N
#BSUB -u jsc228@miami.edu

####################################################
# Set directory variables and conda env
####################################################

cd /scratch/projects/lemmon/jsc228/cns-trauma-seq
mkdir data
cd data
mkdir bams
cd /scratch/projects/lemmon/jsc228/cns-trauma-seq

####################################################
# Job
####################################################

cat milich2021-ftp-urls.sh | parallel -j 8 "{}"

####################################################
# Script End
####################################################