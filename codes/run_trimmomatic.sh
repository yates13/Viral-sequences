#!/bin/bash

# Directory to process
directory="102021_178"

echo "Processing directory: $directory"

# Run Trimmmomatic on the files within the directory
trimmomatic PE -threads 8 -phred33 "/home/aubsdy002/Illumina DNA Reads/$directory/1_S178_R1_001.fastq.gz" "/home/aubsdy002/Illumina DNA Reads/$directory/1_S178_R2_001.fastq.gz" \
               "/home/aubsdy002/trimmed/$directory/${directory}_trimmed_R1.fastq.gz" "/home/aubsdy002/trimmed/$directory/${directory}_trimmed_R1_unpaired.fastq.gz" \
               "/home/aubsdy002/trimmed/$directory/${directory}_trimmed_R2.fastq.gz" "/home/aubsdy002/trimmed/$directory/${directory}_trimmed_R2_unpaired.fastq.gz" \
               LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:50

echo "Finished processing directory: $directory"

