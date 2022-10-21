#!/bin/bash 
#BSUB -J cns-trauma-seq_velocyto
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

conda activate cns-trauma-seq

export PATH=/nethome/jsc228/cellranger-6.0.2:$PATH

# ! NOTE: in-run $HOME is not cns-trauma-seq/
cd /scratch/projects/lemmon/jsc228/cns-trauma-seq/


####################################################
# Job
####################################################

velocyto run10x \
  data/cellranger/Lindsay-Milich-14/outs/possorted_genome_bam.bam \
  refdata-gex-mm10-2020-A/genes/genes.gtf

velocyto run10x \
  data/cellranger/2019-05-22-Lindsay-Milich-3/outs/possorted_genome_bam.bam \
  refdata-gex-mm10-2020-A/genes/genes.gtf

# fin