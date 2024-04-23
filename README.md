# Acropora_cryptic_species

This repository contains scripts to process and analyse genomic data of Acropora hyacinthus on the GBR. Data and code support the results show in:

Naugle, M. S., Denis, H., Mocellin, V. J., Laffy, P. W., Popovic, I., Bay, L. K., & Howells, E. J. (2024). Environmental, host, and symbiont drivers of heat tolerance in a species complex of reef-building corals. bioRxiv, 2024-01.

A version of this manuscript is avalible at: https://www.biorxiv.org/content/10.1101/2024.01.31.575130v1

# fastq_toBam

The following scripts are used to process raw fastq files to BAM files for whole genome sequence data

1) **fastqc-array.sh**
  - All demultiplexed fastq files are visually summarised to explore number of reads, read quality and adapter contamination. 

2) **multiqc.sh**
  - Visual comparison and collective summary of FASTQC reports. Multiqc has various table output options for fastq statistics, which are useful for plotting.

3) **trimmomatic-array.sh**
  - Basic quality filtering of FASTQC files to remove poor-quality reads and adapter contamination
  - Requies list of fastq files (R1 only)
  - Keep bases with phred-score quality > 20 in sliding window of 4 bp (average)
  - Remove adapter sequences in user specified file (**adapters.txt**) with recommended ILLUMINACLIP parameters
  - Remove reads with length < 50bp following trimming
  - keepBothReads option = True (retains paired reads with adapter read-through)

4) **bwa2samtools-array.sh**
  - Map paired reads to the Acropora hyacinthus reference genome `GCA_020536085.1_Ahyacinthus.chrsV1_genomic` available on NCBI:
    https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_020536085.1/
  - Requires a bwa index from the referenc genome
  - Requies a list of R1 clean (paired reads) fastq files

```
# create bwa index
bwa index GCA_020536085.1_Ahyacinthus.chrsV1_genomic.fna

# create list of fastq files
ls *_R1_paired.fastq.gz > fastqClean

```

5) **picard2bam-array.sh**
  - Assign read group information (required for variant calling)
  - Mark and remove duplicates using PICARD

```
# check Read Groups

module load samtools/1.13-gcc-10.3.0
for FILENAME in *markedRG.bam; do
samtools view -H ${FILENAME} | grep '^@RG'
done
```

6) **extract-bam-unmapped-array.sh** [OPTIONAL]
  - Extract unmapped reads from BAM files
  - Unmapped reads will contain all reads not mapped to the reference genome, including symbionts and microbes

7) **samtools_bamStats.sh**
  - Examine bam statistics (e.g., read counts, primary and secondary mapping, coverage) using samtools
  - BAMs can also be individually examined using the Integrative Genomics Viewer 

Additionals Scripts
   - **count-reads.sh** : count reads in a fastq file
   - **samtools-index.sh** : create index file (.bai) for each BAM file
 

