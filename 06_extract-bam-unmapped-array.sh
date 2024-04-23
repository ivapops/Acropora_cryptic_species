#!/bin/bash --login
#SBATCH --job-name="unmap"     
#SBATCH --nodes=1              
#SBATCH --ntasks-per-node=1    
#SBATCH --cpus-per-task=2		
#SBATCH --mem=80G				
#SBATCH --time=20:00:00		
#SBATCH --account=a_riginos	
#SBATCH --partition=general		
#SBATCH -o bwa_%A_%a.o        
#SBATCH -e bwa_%A_%a.e  
#SBATCH --array=1-684        

REF=/scratch/project/rrap_ahya/genome/GCA_020536085.1_Ahyacinthus.chrsV1/GCA_020536085.1_Ahyacinthus.chrsV1_genomic.fna
INDIR=/scratch/project/rrap_ahya/analysis/bwa/markedRG_bamUSE
OUTDIR=/scratch/project/rrap_ahya/analysis/bwa/unmapped_bamUSE

# Pull in all the files within directory. This will include path to file!
FILES=($(ls -1 ${INDIR}/RRAP-ECT01*.bam))

# Allows slurm to enter individual files into the job array
FILENAME=${FILES[$SLURM_ARRAY_TASK_ID]} 
BASE=`basename ${FILENAME} _markedRG.bam`

echo "My input file is ${FILENAME}"

module load samtools/1.13-gcc-10.3.0

echo ${FILENAME} >> ${OUTDIR}/sample_name.txt		# generate sample names
samtools view -c ${FILENAME} >> ${OUTDIR}/allCounts.txt		# generate read counts
samtools view -b -f 12 -F 256 ${FILENAME} > ${OUTDIR}/${BASE}.unmapped.bam		# generate bam
samtools view -c ${OUTDIR}/${BASE}.unmapped.bam >> ${OUTDIR}/unmappedCounts.txt		# generate unmapped read counts

# Notes:
# -f 4: extract all unmapped reads
# -f12: extract only the reads where read 1 is unmapped AND read 2 is unmapped (= both mates are unmapped). Applies only to paired reads.


