---
layout: project
title: TitanCNA/scripts/
project: TitanCNA
repo: gavinha/TitanCNA
permalink: /:path/:basename:output_ext
---

# Scripts and Pipeline code for TitanCNA analysis
## Contents
* [Standard Whole Genome/Exome Sequencing Analysis](#wgs)
* [10X Genomics Whole Genome Sequencing Analysis](#tenx)

## 1. Standard Whole Genome/Exome Sequencing Analysis
a. We provide a Snakemake implementation of the pipeline to perform TITAN analysis starting from BAM files to TitanCNA results.
Please see [snakemake/](snakemake) for more details.

b. R scripts for running TitanCNA analysis  
The `R_scripts` directory contains the actual R script to run the (R component of the) TitanCNA analysis as well as an R script to select the optimal solution.  
While the Snakemake pipeline will make use of R script within `R_scripts`, you can learn about more details of the R script here [R_scripts/](R_scripts/).  

## 2. 10X Genomics Whole Genome Sequencing Analysis
These are scripts to perform TITAN analysis from 10X Genomics whole genome sequencing data.  Please see [TenX_scripts/](TenX_scripts/) for more details.  
A Snakemake pipeline will be provided soon...
