#!/bin/bash --login
#SBATCH --job-name="index"
#SBATCH --nodes=1             
#SBATCH --ntasks-per-node=1    
#SBATCH --cpus-per-task=1		
#SBATCH --mem=8G			
#SBATCH --time=10:00:00			
#SBATCH --account=a_riginos		
#SBATCH --partition=general	
#SBATCH -o indx_%A_%a.o      
#SBATCH -e indx_%A_%a.e        
#SBATCH --array=1-684      

module load samtools/1.13-gcc-10.3.0

LIST=myBAMS
# e.g., RRAP-ECO3-2021-Ahya-CBHE-218_L3

INDIR=/analysis/bwa/markedRG_bamUSE
REF=/genome/GCA_020536085.1_Ahyacinthus.chrsV1/GCA_020536085.1_Ahyacinthus.chrsV1_genomic.fna

# file name variable is associated Array index
FILENAME=`cat ${LIST} | tail -n +${SLURM_ARRAY_TASK_ID} | head -1`

samtools index ${INDIR}/${FILENAME}
