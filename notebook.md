---
title: "cns-trauma-seq notebook"
author: "James Choi"
date: 'Last compiled: 2022-08-31'
output:
  html_document:
    keep_md: yes
  pdf_document: default
  word_document: default
editor_options:
  chunk_output_type: console
---




## About

Single cell sequencing resource of data from central nervous system traumatic injury.


## 2022-08-29

Workflow for this study and single cell model were separated into two tasks where the first task was to preprocess sequencing data (e.g. `.fastq`/`.fq` files) into gene-count numeric matrices and the second task was to perform cluster analysis and differential expression tests. 

One of the first issues is that I only have access to `.bam` and not `.fastq` files for certain samples. 10X Genomics has a software tool `bamtofastq` that converts 10X Genomics `.bam` files back into `.fastq` so that they can be used as input for analysis re-runs. I specifically need the raw `.bam` files, and not the files decompressed from `.sra`, due to [a known issue](https://support.10xgenomics.com/docs/bamtofastq#issues). To use these tools I'll need to use the Pegasus cluster. 

I requested the following data from our SRA submission (Milich2021):

  * BioProject: PRJNA682392
  * https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA682392&o=acc_s%3Aa
  * `data/Milich2021_SRR_Acc_List.txt`

I saved the *Accession List* from under *Download* column of the table under the *Select* tab.

By continuing to explore the data deposited at SRA, I saw the data available like so: https://trace.ncbi.nlm.nih.gov/Traces/index.html?view=run_browser&acc=SRR13192645&display=data-access

The `.sra` and `.bam` file format options were found here.

So I needed a method to take these accessions and get the full list of download links. Using WSL2 on a Windows machine, I set up a `conda` environment named `cns-trauma-seq`. I installed the tool [`ffq` developed by Pachter Lab at UC Berkeley](https://github.com/pachterlab/ffq) which allows me to access SRA metadata using only accession codes. 

I confirmed that `ffq` is working properly using SRR13192645. The info I needed was under `"{"SRR13192645":"files":{ftp:["url":...]}}`. 

```{}
>> ffq SRR13192645 | head
[2022-08-29 15:34:27,261]    INFO Parsing run SRR13192645
{
    "SRR13192645": {
        "accession": "SRR13192645",
        "experiment": "SRX9627035",
        "study": "SRP295673",
        "sample": "SRS7828004",
        "title": "NextSeq 500 paired end sequencing; GSM4955359: uninj_sample1; Mus musculus; RNA-Seq",
        "attributes": {
            "assembly": "GRCm38",
            "ENA-FIRST-PUBLIC": "2021-03-02",
            
>> ffq -o ffq-output.json SRR13192645 SRR13192646 SRR13192647 SRR13192648 SRR13192649 SRR13192650 SRR13192651 SRR13192652 SRR13192653 SRR13192654
```


I will eventually need to reorganize the `json` output into a batch process command. I could do this programmatically or by hand but first I need to download a test dataset. 

`ssh` to Pegasus and tested `wget https://sra-pub-src-2.s3.amazonaws.com/SRR13192645/uninj_sample1_possorted_genome_bam.bam.1`, which was ETA of 20 minutes. For all 10 samples, that's just over 3 hours. First I parse `ffq-output.json` to generate a `txt` file of just links.


```r
links <- jsonlite::fromJSON(txt = "ffq-output.json")
links <- lapply(links, function(x) x[["files"]][["ftp"]])
links <- sapply(links, function(x) x$url[x$filetype == "bam"][1])
tmp.names <- sapply(strsplit(links, "/"), function(x) rev(x)[])
writeLines(text = links, con = "milich2021-ftp-urls.txt")
# writeLines(text = paste0('curl -L ', links, ' -o ', tmp.names), con =
# 'milich2021-ftp-urls.txt')
```


## 2022-08-31

Download of bams was successful. Each bam was between 35-40Gb. 

Next was test `bamtofastq` tool from `cellranger` by 10X Genomics. I used an interactive job session on Pegasus. I followed the instructions here: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/6.0/installation

I'm downloading `cellranger-6.0.2` because that is the version Dr Griswold at HIHG used for the aged/proliferation SCI scRNAseq study. That way, I can re-align two Milich2021 samples for the aged/proliferation study and I can test-run `bamtofastq`. 

In pegasus interactive mode:

```
# Download cellranger binary
$ curl -o cellranger-6.0.2.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-6.0.2.tar.gz?Expires=1662008119&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC1leHAvY2VsbHJhbmdlci02LjAuMi50YXIuZ3oiLCJDb25kaXRpb24iOnsiRGF0ZUxlc3NUaGFuIjp7IkFXUzpFcG9jaFRpbWUiOjE2NjIwMDgxMTl9fX1dfQ__&Signature=bdlpEB~od9rl6LMWCIZqvP-ESlKuyg~Upc45K0nQlHgprSBPviOZWKj6uuRA4UZK4IjQfqR94TFmgdkztg-Tv8kV3eFK0Cojyrxmokhqrqwnp0Zr5Rj4XjU66Z8jlQEtQfhjdoBfkqu34J7ZBreknFLxs9XEzygPMxew-~CwIGdVzNaIqVw9E0Vm0rAdQD44gEFdAcEgJcjyccBHJZuHD7yddfi7R58osZ~0LDTq4eaGNlJQaCTbElLgEPKoRolKWqBsOMwU1hn28vKUZjcX~h3XJyYzovqZahavghQ7tXYf6-epzB~vGjY43XlsvMFFSCbGSxI5gertXfWhg4ToLg__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"

$ tar -xvzf cellranger-6.0.2.tar.gz

# Download reference data files (genome annotations, etc.)
$ curl -O https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-mm10-2020-A.tar.gz

$ tar -xzvf refdata-gex-GRCh38-2020-A.tar.gz

# Also added to PATH in .bashrc
$ export PATH=/opt/cellranger-6.0.2:$PATH
$ cellranger sitecheck > sitecheck.txt
$ cellranger testrun --id=tiny

$ cd cellranger-6.0.2
$ cd /scratch/projects/lemmon/jsc228/cns-trauma-seq
$ /nethome/jsc228/cellranger-6.0.2/lib/bin/bamtofastq --help

bamtofastq v1.3.1
10x Genomics BAM to FASTQ converter.

    Tool for converting 10x BAMs produced by Cell Ranger or Long Ranger back to
    FASTQ files that can be used as inputs to re-run analysis. The FASTQ files
    emitted by the tool should contain the same set of sequences that were
    input to the original pipeline run, although the order will not be
    preserved.  The FASTQs will be emitted into a directory structure that is
    compatible with the directories created by the 'mkfastq' tool.

    10x BAMs produced by Long Ranger v2.1+ and Cell Ranger v1.2+ contain header
    fields that permit automatic conversion to the correct FASTQ sequences.

    Older 10x pipelines require one of the arguments listed below to indicate
    which pipeline created the BAM.
```
