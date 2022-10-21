#!/bin/bash
#BSUB -J cns-trauma-seq_cellranger-count
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

# 1dpi sample1
cellranger count \
  --id=Lindsay-Milich \
  --transcriptome=/scratch/projects/lemmon/jsc228/cns-trauma-seq/refdata-gex-mm10-2020-A \
  --fastqs=/scratch/projects/lemmon/jsc228/cns-trauma-seq/data/bamtofastq/1dpi_sample1_possorted_genome_bam.bam.1/Lindsay-Milich_MissingLibrary_1_HMVM5BGX7
# 1dpi sample 2
cellranger count \
  --id=2019-05-21-Lindsay-Milich-1 \
  --transcriptome=/scratch/projects/lemmon/jsc228/cns-trauma-seq/refdata-gex-mm10-2020-A \
  --fastqs=/scratch/projects/lemmon/jsc228/cns-trauma-seq/data/bamtofastq/1dpi_sample2_possorted_genome_bam.bam.1/2019-05-21-Lindsay-Milich-1_0_1_HGK7YBGXB
# 1dpi sample 3
cellranger count \
  --id=Lee-1-2-2020-09-22 \
  --transcriptome=/scratch/projects/lemmon/jsc228/cns-trauma-seq/refdata-gex-mm10-2020-A \
  --fastqs=/scratch/projects/lemmon/jsc228/cns-trauma-seq/data/bamtofastq/1dpi_sample3_possorted_genome_bam.bam.1/Lee-1-2-2020-09-22_0_1_HNTGLBGXF
# 3dpi sample 1
cellranger count \
  --id=Lindsay-Milich-14 \
  --transcriptome=/scratch/projects/lemmon/jsc228/cns-trauma-seq/refdata-gex-mm10-2020-A \
  --fastqs=/scratch/projects/lemmon/jsc228/cns-trauma-seq/data/bamtofastq/3dpi_sample1_possorted_genome_bam.bam.1/Lindsay-Milich-14_MissingLibrary_1_HV5CFBGX7
# 3dpi sample 2
cellranger count \
  --id=2019-05-22-Lindsay-Milich-3 \
  --transcriptome=/scratch/projects/lemmon/jsc228/cns-trauma-seq/refdata-gex-mm10-2020-A \
  --fastqs=/scratch/projects/lemmon/jsc228/cns-trauma-seq/data/bamtofastq/3dpi_sample2_possorted_genome_bam.bam.1/2019-05-22-Lindsay-Milich-3_0_1_HHTVCBGXB
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