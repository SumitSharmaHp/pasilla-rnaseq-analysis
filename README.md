# Pasilla RNA-seq Analysis

## Overview

This project demonstrates an end-to-end RNA-seq differential expression analysis workflow using the Drosophila melanogaster Pasilla dataset.

The analysis was performed on Ubuntu Linux using standard bioinformatics tools and R/Bioconductor packages.

## Workflow

1. Quality Control (FastQC)
2. Read Trimming (fastp)
3. Genome Alignment (HISAT2)
4. BAM Processing (SAMtools)
5. Read Quantification (featureCounts)
6. Differential Expression Analysis (DESeq2)

## Dataset

Drosophila melanogaster Pasilla RNA-seq dataset.

Samples analyzed:

* SRR031714
* SRR031716
* SRR031724
* SRR031726

## Tools Used

* FastQC
* fastp
* HISAT2
* SAMtools
* featureCounts
* R
* DESeq2

## Results

After filtering low-expression genes:

* Total genes analyzed: 9868
* Upregulated genes: 364
* Downregulated genes: 404
* Total significant genes: 768

Generated outputs:

* Gene count matrix
* Read assignment summary
* Differential expression results

## Repository Structure

results/
├── DESeq2_results.csv
├── pasilla_counts.txt
└── pasilla_counts.txt.summary

## Author

Sumit Sharma

M.Sc. Bioinformatics
