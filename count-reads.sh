#!/bin/bash --login
#SBATCH --job-name="count"   
#SBATCH --nodes=1              
#SBATCH --ntasks-per-node=1    
#SBATCH --cpus-per-task=1		
#SBATCH --mem=10G			
#SBATCH --time=4:00:00		
#SBATCH --account=a_riginos		
#SBATCH --partition=general		
#SBATCH -o reads_%A_%a.o      
#SBATCH -e reads_%A_%a.e	   
#SBATCH --array=1-684

LIST=/analysis/fastqc/ECT/fastq_cleanList
INDIR=/genomic_datasets/fastq_cleanUSE

FILENAME=`cat ${LIST} | tail -n +${SLURM_ARRAY_TASK_ID} | head -1`

echo "My interval is ${FILENAME}"

NAME=`echo ${FILENAME} `
COUNT=`echo $(zcat ${INDIR}/${FILENAME} | wc -l)/4|bc`

echo $NAME $COUNT  >> cleanCounts.txt

