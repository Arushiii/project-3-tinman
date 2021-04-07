# Project Description

This repository contains the scripts used to replicate parts of the results in the paper "A comprehensive study design reveals treatment- and transcript abundance–dependent concordance between RNA-seq and microarray data". This project is divided into three parts: data curator, programmer, and analyst. It contains the bash code and R code that was used to pre-process and analyze the down-stream data.

Wang, Charles, Binsheng Gong, Pierre R. Bushel, Jean Thierry-Mieg, Danielle Thierry-Mieg, Joshua Xu, Hong Fang, et al. 2014. “A comprehensive study design reveals treatment- and transcript abundance–dependent concordance between RNA-seq and microarray data” Nature Biotechnology 32 (9): 926–32. PMID: 4243706


# Contributors

Data Curator: Teresa Rice - tpillars@bu.edu

Programmer: Maha Naim - mnaim21@bu.edu

Analyst: Arushi Shrivastava - arushi08@bu.edu

# Repository Contents
Prepare for severly uncreative file names

## commands.txt
Copy toxgroup 4 samples into working directory and run fastqc on the fastq.gz files. You need to manually create 2 directories(out_fastqc and fastq), manually move files to respective folders.

## commands2.txt
Runs STAR analysis using rat genome reference from BU, and puts output into directory star_results2 (need to mkdir star_results2 manually).

## commands3.txt
Runs MultiQC on out_fastqc files produced from commands.txt

## datacuratorRcommands.R
Utilizes libraries: tidyverse, dplyr, and gridExtra. Individually loads in each star log.final.out file into a table (modify the object to include the sample ID) full join all STAR output tables. Creates STAR summary table of select QC metrics, creates png from this table.
