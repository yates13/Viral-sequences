---
title: Viral Sequences Project
---

## Introduction

This project aims to analyze and study viral sequences using various bioinformatics tools and techniques.

## Table of Contents

- [Introduction](#introduction)
- [Installation](#installation)
- [Usage](#usage)
- [Results](#results)
- [Contributing](#contributing)
  

## Installation

Before running the viral sequence analysis, ensure that the following software and packages are installed:

- BWA (Burrows-Wheeler Aligner)
- Samtools
- BCFtools
- GATK (Genome Analysis Toolkit)
- FastQC

Make sure these tools are installed and accessible in your system's PATH.

## Usage

1. Navigate to the project directory:

   ```bash
   cd viral-sequencUpdate the paths and filenames in the provided scripts according to your dataset and reference genomes.
   ```

 ## Perform fastqc to determine the quality of your data (refer to fastqc in codes file)
 
  ```bash
# Find 'your_file.fastq' for your sample
fastqc "your_file.fastq"
```

## Run the main script for read mapping:

   ```bash
./mapping_k181.sh
```

## Variant call files
  ```bash
bcftools mpileup -Ou -f k181.fasta sorted_wt.bam | bcftools call -mv > wt_variants.vcf
```

The output files, including alignment files, and directories.

## Results
The results of this project include:

- Aligned BAM files (sorted and indexed)
- Variant files (VCF format)
- Consensus sequences
- These outputs can be further analyzed and visualized using various tools such as IGV or custom scripts.

## Contributing
Contributions to this project are welcome. If you encounter any issues, have suggestions for improvements, or want to contribute new features, please submit an issue or a pull request.es
