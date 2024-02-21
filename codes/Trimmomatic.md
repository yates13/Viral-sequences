---
title: Bash Script Execution Report
---

# Bash Script Execution Report

## Introduction

This report documents the execution of the Bash script for processing Illumina DNA reads using Trimmomatic.

## Script Description

The Bash script processes the specified directory of Illumina DNA reads using Trimmomatic. It trims adapter sequences and performs quality filtering on the reads.

## Directory to Process

The script processes the following directory:

```bash
directory="102021_178"
echo "Processing directory: $directory"
```

## Execution
The script executes the Trimmomatic command with the following parameters:

# Input files:
 - 1_S178_R1_001.fastq.gz
 - 1_S178_R2_001.fastq.gz

# Output files:
 - ${directory}_trimmed_R1.fastq.gz
 - ${directory}_trimmed_R1_unpaired.fastq.gz
 - ${directory}_trimmed_R2.fastq.gz
 - ${directory}_trimmed_R2_unpaired.fastq.gz

```bash
trimmomatic PE -threads 8 -phred33 "/home/aubsdy002/Illumina DNA Reads/$directory/1_S178_R1_001.fastq.gz" "/home/aubsdy002/Illumina DNA Reads/$directory/1_S178_R2_001.fastq.gz" \
               "/home/aubsdy002/trimmed/$directory/${directory}_trimmed_R1.fastq.gz" "/home/aubsdy002/trimmed/$directory/${directory}_trimmed_R1_unpaired.fastq.gz" \
               "/home/aubsdy002/trimmed/$directory/${directory}_trimmed_R2.fastq.gz" "/home/aubsdy002/trimmed/$directory/${directory}_trimmed_R2_unpaired.fastq.gz" \
               LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:50
```
## Conclusion
The Bash script successfully processed the specified directory of Illumina DNA reads using Trimmomatic, trimming adapter sequences and performing quality filtering.
