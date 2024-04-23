#!/bin/bash --login
#SBATCH --job-name="bwa_map"
#SBATCH --nodes=1             
#SBATCH --ntasks-per-node=1  
#SBATCH --cpus-per-task=2		
#SBATCH --mem=20G			
#SBATCH --time=24:00:00	
#SBATCH --account=a_riginos	
#SBATCH --partition=general	
#SBATCH -o bwa_%A_%a.o        
#SBATCH -e bwa_%A_%a.e      
#SBATCH --array=1-684    

LIST=/scratch/project/rrap_ahya/analysis/lists/fastqClean
INDIR=/scratch/project/rrap_ahya/analysis/trimmomatic/fastq_cleanUSE
OUTDIR=/scratch/project/rrap_ahya/analysis/bwa/markedRG_bamUSE
REF=/scratch/project/rrap_ahya/genome/GCA_020536085.1_Ahyacinthus.chrsV1/GCA_020536085.1_Ahyacinthus.chrsV1_genomic.fna

# Remember to first make bwa index from reference genome!
# bwa index GCA_020536085.1_Ahyacinthus.chrsV1_genomic.fna

# file name variable is associated Array index
FILENAME=`cat ${LIST} | tail -n +${SLURM_ARRAY_TASK_ID} | head -1`
BASE=`basename ${FILENAME} _R1_paired.fastq.gz`

# execute bwa using paired reads
module load bwa/0.7.17-gcc-10.3.0   
module load samtools/1.13-gcc-10.3.0

# map paired reads to reference genome
bwa mem -t 2 ${REF} ${INDIR}/${BASE}_R1_paired.fastq.gz ${INDIR}/${BASE}_R2_paired.fastq.gz > ${OUTDIR}/${BASE}.sam

# create sorted bam
samtools view -b ${OUTDIR}/${BASE}.sam > ${OUTDIR}/${BASE}_UNSORTED.bam
samtools sort ${OUTDIR}/${BASE}_UNSORTED.bam -O BAM -o ${OUTDIR}/${BASE}_UNDEDUP.bam

# removed temp files
rm ${OUTDIR}/${BASE}_UNSORTED.bam
rm ${OUTDIR}/${BASE}.sam

# Notes: 
# -t : number of threads
# -O : BAM output format  
# -o : output file or can use '>'
