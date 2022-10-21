#!/bin/bash
#BSUB -J cns-trauma-seq_cellranger-count-uninj
#BSUB -P lemmon
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -W 12:00
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


# uninj sample1
cellranger count \
  --id=Lindsay-Milich-95 \
  --transcriptome=/scratch/projects/lemmon/jsc228/cns-trauma-seq/refdata-gex-mm10-2020-A \
  --fastqs=/scratch/projects/lemmon/jsc228/cns-trauma-seq/data/bamtofastq/uninj_sample1_possorted_genome_bam.bam.1/Lindsay-Milich-95_MissingLibrary_1_HMVKHBGX7
# uninj sample2
cellranger count \
  --id=Lindsay-Milich-2-1 \
  --transcriptome=/scratch/projects/lemmon/jsc228/cns-trauma-seq/refdata-gex-mm10-2020-A \
  --fastqs=/scratch/projects/lemmon/jsc228/cns-trauma-seq/data/bamtofastq/uninj_sample2_possorted_genome_bam.bam.1/Lindsay-Milich-2-1_0_1_HM22WBGX9
# uninj sample3
cellranger count \
  --id=Lee-U-2-2020-09-23 \
  --transcriptome=/scratch/projects/lemmon/jsc228/cns-trauma-seq/refdata-gex-mm10-2020-A \
  --fastqs=/scratch/projects/lemmon/jsc228/cns-trauma-seq/data/bamtofastq/uninj_sample3_possorted_genome_bam.bam.1/Lee-U-2-2020-09-23_0_1_HTYVNBGXG


####################################################
# Script End
####################################################