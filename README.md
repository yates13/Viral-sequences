---
title: "Viral Sequences Project"
---

## Table of Contents

- [Introduction](#introduction)
- [Installation](#installation)
- [Usage](#usage)
- [Results](#results)
- [Contributing](#contributing)
  
## Introduction

This project aims to analyze and study viral sequences using various bioinformatics tools and techniques.

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
   cd viral-sequence "Update the paths and filenames in the provided scripts according to your dataset and reference genomes".
   ```

 ## Perform fastqc to determine the quality of your data (refer to fastqc in codes file)

 2. Run fatqc on your sample
    
  ```bash
# Find 'your_file.fastq' for your sample
fastqc "your_file.fastq"
```

## Run the main script for read mapping:

3. Trimming with Trimmomatic
   
   Change fq1 and fq2 to take in your sample in "mapping_k181.sh" in the codes folder. Also change the refernce genome to the genome of yopur choosing.
   
  ```bash
#!/bin/bash

# STEP 1: Read mapping to reference genome using BWA
module load bwa/0.7.12

# Indexing the reference genome
bwa index k181.fasta

# STEP 2: Mapping with BWA

# Input fastq files
fq1="102021_178_trimmed_R1.fastq.gz"
fq2="102021_178_trimmed_R2.fastq.gz"

echo "Working with file $fq1"

# Perform read mapping using BWA and save the output to SAM file
bwa mem k181.fasta $fq1 $fq2 > output_k181_wt.sam

# Convert the SAM file to BAM format
samtools view -bS output_k181_wt.sam > output_k181_wt.bam

# Sort the BAM file and save the sorted output to a new BAM file
samtools sort -o sorted_wt.bam -O BAM -T temp output_k181_wt.bam

# Index the sorted BAM file
samtools index sorted_wt.bam

echo "Alignment completed."
```


## Variant call files

4. Create variant files to look for mutations

  ```bash

#!/bin/bash

# STEP 4: Variant calling using samtools and bcftools

# Define input files
reference_genome="GCF_000859805.1_ViralProj15181_genomic.fna"
mapped_reads="wt_sorted.bam"
variants_vcf="variants.vcf"

# Load samtools and bcftools modules
module load samtools
module load bcftools

# Generate VCF file of variants
samtools mpileup -uf $reference_genome $mapped_reads | bcftools call -mv > $variants_vcf

echo "Variant calling completed."

# seperate thing 
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
