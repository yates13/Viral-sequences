This Bash script is designed to perform read mapping to a reference genome using BWA (Burrows-Wheeler Aligner) and post-processing steps with Samtools. Let's break down the script step by step:

## STEP 1: Indexing the Reference Genome

```bash
module load bwa/0.7.12
bwa index k181.fasta
```
The module load command is used to load the BWA module.
bwa index is used to create an index of the reference genome file k181.fasta.

## STEP 2: Mapping Reads with BWA

### Input fastq files
fq1="102021_178_trimmed_R1.fastq.gz"

fq2="102021_178_trimmed_R2.fastq.gz"

### Perform read mapping using BWA and save the output to a SAM file
```bash
bwa mem k181.fasta $fq1 $fq2 > output_k181_wt.sam
```
fq1 and fq2 store the paths to the input FASTQ files.
bwa mem is used to perform read mapping to the reference genome, generating a SAM file (output_k181_wt.sam).
## Post-processing with Samtools
### Convert the SAM file to BAM format
```bash
samtools view -bS output_k181_wt.sam > output_k181_wt.bam
```
### Sort the BAM file and save the sorted output to a new BAM file
```bash
samtools sort -o sorted_wt.bam -O BAM -T temp output_k181_wt.bam
```
### Index the sorted BAM file
samtools index sorted_wt.bam
samtools view converts the SAM file to BAM format (output_k181_wt.bam).
samtools sort sorts the BAM file (output_k181_wt.bam) and saves the sorted output to a new BAM file (sorted_wt.bam).
samtools index creates an index for the sorted BAM file (sorted_wt.bam).
Finally, the script concludes with a message indicating the completion of the alignment process.

```bash
echo "Alignment completed."
```
