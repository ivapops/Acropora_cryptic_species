#!/bin/bash --login
#SBATCH --job-name="index"      # job name
#SBATCH --nodes=1               # use 1 node
#SBATCH --ntasks-per-node=1     # use 1 for single and multi core jobs
#SBATCH --cpus-per-task=1		# number of cores per job
#SBATCH --mem=8G				# RAM per job given in megabytes (M), gigabytes (G), or terabytes (T)
#SBATCH --time=10:00:00			# walltime
#SBATCH --account=a_riginos		# group account name
#SBATCH --partition=general		# queue name
#SBATCH -o indx_%A_%a.o         # standard output
#SBATCH -e indx_%A_%a.e         # standard error
#SBATCH --array=1-684        	# job array

module load samtools/1.13-gcc-10.3.0

LIST=myBAMS
# e.g., RRAP-ECO3-2021-Ahya-CBHE-218_L3

INDIR=/scratch/project/rrap_ahya/analysis/bwa/markedRG_bamUSE
REF=/scratch/project/rrap_ahya/genome/GCA_020536085.1_Ahyacinthus.chrsV1/GCA_020536085.1_Ahyacinthus.chrsV1_genomic.fna

# file name variable is associated Array index
FILENAME=`cat ${LIST} | tail -n +${SLURM_ARRAY_TASK_ID} | head -1`

samtools index ${INDIR}/${FILENAME}
