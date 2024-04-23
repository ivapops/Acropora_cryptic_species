#!/bin/bash --login
#SBATCH --job-name="bwa_stat"
#SBATCH --nodes=1              
#SBATCH --ntasks-per-node=1  
#SBATCH --cpus-per-task=1	
#SBATCH --mem=20G		
#SBATCH --time=4:00:00	
#SBATCH --account=a_riginos		
#SBATCH --partition=general	
#SBATCH -o stat_%A_%a.o     
#SBATCH -e stat_%A_%a.e 
#SBATCH --array=1-684    

module load samtools/1.13-gcc-10.3.0

LIST=RRAP-list
INDIR=/scratch/project/rrap_ahya/analysis/bwa/markedRG_bamUSE
OUTDIR=/scratch/project/rrap_ahya/analysis/bwa/bam_stats

# file name variable is associated Array index
FILENAME=`cat ${LIST} | tail -n +${SLURM_ARRAY_TASK_ID} | head -1`
BASE=`basename ${FILENAME} _markedRG.bam`

cd ${INDIR}

# Generate index bai and statistics from a BAM file: aligned reads
# refer to samtools-index.sh
# samtools index ${BASE}_markedRG.bam

samtools idxstats ${BASE}_markedRG.bam > ${OUTDIR}/${BASE}-index_stats.txt
samtools coverage ${BASE}_markedRG.bam > ${OUTDIR}/${BASE}-coverage.txt
samtools flagstat -O tsv ${BASE}_markedRG.bam > ${OUTDIR}/${BASE}-flagstats.txt
