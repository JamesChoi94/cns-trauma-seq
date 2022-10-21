#!/bin/bash
#BSUB -J cns-trauma-seq_download-data
#BSUB -P lemmon
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -W 4:00
#BSUB -q general
#BSUB -n 4
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

cat milich2021-ftp-urls.txt | xargs -n 1 -P 8 wget

####################################################
# Script End
####################################################