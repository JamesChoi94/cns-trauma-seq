#!/bin/bash
#BSUB -J cns-trauma-seq_bamtofastq
#BSUB -P lemmon
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -W 6:00
#BSUB -q general
#BSUB -n 8
#BSUB -B
#BSUB -N
#BSUB -u jsc228@miami.edu

####################################################
# Set directory variables and conda env
####################################################

export PATH=/nethome/jsc228/cellranger-6.0.2:$PATH
cd /scratch/projects/lemmon/jsc228/cns-trauma-seq
mkdir data/bamtofastq/

####################################################
# Job
####################################################

/nethome/jsc228/cellranger-6.0.2/lib/bin/bamtofastq \
  --traceback \
  --nthreads=8 \
  data/bams/3dpi_sample1_possorted_genome_bam.bam.1 \
  data/bamtofastq/3dpi_sample1_possorted_genome_bam.bam.1/

  /nethome/jsc228/cellranger-6.0.2/lib/bin/bamtofastq \
  --traceback \
  --nthreads=8 \
  data/bams/3dpi_sample2_possorted_genome_bam.bam.1 \
  data/bamtofastq/3dpi_sample2_possorted_genome_bam.bam.1/
  
####################################################
# Script End
####################################################