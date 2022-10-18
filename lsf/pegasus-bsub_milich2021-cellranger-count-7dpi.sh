#!/bin/bash
#BSUB -J cns-trauma-seq_cellranger-count-7dpi
#BSUB -P lemmon
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -W 10:00
#BSUB -q bigmem
#BSUB -n 10
#BSUB -R "rusage[mem=4096]"
#BSUB -B
#BSUB -N
#BSUB -u jsc228@miami.edu

####################################################
# Set directory variables and conda env
####################################################

export PATH=/nethome/jsc228/cellranger-6.0.2:$PATH

# ! NOTE: in-run $HOME is not cns-trauma-seq/
cd /scratch/projects/lemmon/jsc228/cns-trauma-seq/data
mkdir cellranger
cd cellranger

####################################################
# Job
####################################################

# 7dpi sample 1
cellranger count \
  --id=Lindsay-Milich-78 \
  --transcriptome=/scratch/projects/lemmon/jsc228/cns-trauma-seq/refdata-gex-mm10-2020-A \
  --fastqs=/scratch/projects/lemmon/jsc228/cns-trauma-seq/data/bamtofastq/7dpi_sample1_possorted_genome_bam.bam.1/Lindsay-Milich-78_MissingLibrary_1_HV5F5BGX7
# 7dpi sample 2
cellranger count \
  --id=2019-05-23-Lindsay-Milich-7 \
  --transcriptome=/scratch/projects/lemmon/jsc228/cns-trauma-seq/refdata-gex-mm10-2020-A \
  --fastqs=/scratch/projects/lemmon/jsc228/cns-trauma-seq/data/bamtofastq/7dpi_sample2_possorted_genome_bam.bam.1/2019-05-23-Lindsay-Milich-7_0_1_HHTVNBGXB

####################################################
# Script End
####################################################