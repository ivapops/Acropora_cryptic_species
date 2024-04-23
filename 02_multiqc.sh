#!/bin/bash --login
#SBATCH --job-name="multiqc"    # job name
#SBATCH --nodes=1				# use 1 node
#SBATCH --ntasks-per-node=1     # use 1 for single and multi core jobs
#SBATCH --cpus-per-task=1		# number of cores per job
#SBATCH --mem=10G				# RAM per job given in megabytes (M), gigabytes (G), or terabytes (T)
#SBATCH --time=10:00:00			# walltime
#SBATCH --account=a_riginos		# group account name
#SBATCH --partition=general		# queue name
#SBATCH -o mfqc_%A.o            # standard output
#SBATCH -e mfqc_%A.e	        # standard error

module load anaconda3/2022.05
source activate multiqc

# use FASTQC output directory as input
multiqc output/*_fastqc.zip --interactive


