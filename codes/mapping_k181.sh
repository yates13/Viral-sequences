#!/bin/bash

# STEP 1: Read mapping to reference genome using BWA
module load bwa/0.7.12

# Indexing the reference genome
bwa index k181.fasta

# STEP 2: Mapping with BWA
module load bwa/0.7.12

# Input fastq files
fq1="102021_178_trimmed_R1.fastq.gz"
fq2="102021_178_trimmed_R2.fastq.gz"

echo "Working with file $fq1"

bwa mem k181.fasta $fq1 $fq2 > output_k181_wt.sam

samtools view -bS output_k181_wt.sam > output_k181_wt.bam

samtools sort -o sorted_wt.bam -O BAM -T temp output_k181_wt.bam
samtools index sorted_wt.bam

echo "Alignment completed."

