#!/bin/bash --login
#SBATCH --job-name="fastqc"  
#SBATCH --nodes=1               
#SBATCH --ntasks-per-node=1		
#SBATCH --cpus-per-task=1		
#SBATCH --mem=10G				
#SBATCH --time=10:00:00			
#SBATCH --account=a_riginos		
#SBATCH --partition=general		
#SBATCH -o fqc_%A_%a.o          
#SBATCH -e fqc_%A_%a.e	        
#SBATCH --array=1-684        	  

module load fastqc/0.11.9-java-11 

LIST=/analysis/fastqc/fastqClean_ECT
INDIR=/analysis/fastqc/fastq_cleanUSE

# file name variable is associated ARRAY index
FILENAME=`cat ${LIST} | tail -n +${SLURM_ARRAY_TASK_ID} | head -1`

echo "My interval is ${FILENAME}"

fastqc --extract ${INDIR}/${FILENAME} -f fastq -o output

