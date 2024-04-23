#!/bin/bash --login
#SBATCH --job-name="fastqc"     # job name
#SBATCH --nodes=1               # use 1 node
#SBATCH --ntasks-per-node=1		# use 1 for single and multi core jobs
#SBATCH --cpus-per-task=1		# number of cores per job
#SBATCH --mem=10G				# RAM per job given in megabytes (M), gigabytes (G), or terabytes (T)
#SBATCH --time=10:00:00			# walltime
#SBATCH --account=a_riginos		# group account name
#SBATCH --partition=general		# queue name
#SBATCH -o fqc_%A_%a.o          # standard output
#SBATCH -e fqc_%A_%a.e	        # standard error
#SBATCH --array=1-684        	# job array     

module load fastqc/0.11.9-java-11 

LIST=/scratch/project/rrap_ahya/analysis/fastqc/ECT/fastqClean_ECT
INDIR=/QRISdata/Q4020/Iva_Popovic/genomic_datasets/fastq_cleanUSE

# file name variable is associated ARRAY index
FILENAME=`cat ${LIST} | tail -n +${SLURM_ARRAY_TASK_ID} | head -1`

echo "My interval is ${FILENAME}"

fastqc --extract ${INDIR}/${FILENAME} -f fastq -o output

