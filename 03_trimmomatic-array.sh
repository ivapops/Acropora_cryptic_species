#!/bin/bash --login
#SBATCH --job-name="trim"       # job name
#SBATCH --nodes=1               # use 1 node
#SBATCH --ntasks-per-node=1     # use 1 for single and multi core jobs
#SBATCH --cpus-per-task=1		# number of cores per job
#SBATCH --mem=80G				# RAM per job given in megabytes (M), gigabytes (G), or terabytes (T)
#SBATCH --time=20:00:00			# walltime
#SBATCH --account=a_riginos		# group account name
#SBATCH --partition=general		# queue name
#SBATCH --array=1-94        	# job array
#SBATCH -e trim_%A_%a.e      	# standard error
#SBATCH -o trim_%A_%a.o

# text fille including fastq file names (R1 only)
LIST=fastqList
INDIR=/QRISdata/Q4020/Iva_Popovic/genomic_datasets/raw_data/all_fastq
OUTDIR=/scratch/project/rrap_ahya/analysis/trimmomatic/fastq_clean

# file name variable is associated ARRAY index
FILENAME=`cat ${LIST} | tail -n +${SLURM_ARRAY_TASK_ID} | head -1`
BASE=`basename ${FILENAME} _1.fq.gz`

module load anaconda3/2022.05
source activate trimmomatic

trimmomatic PE -phred33 ${INDIR}/${BASE}_1.fq.gz ${INDIR}/${BASE}_2.fq.gz \
${OUTDIR}/${BASE}_R1_paired.fastq.gz ${OUTDIR}/${BASE}_R1_unpaired.fastq.gz \
${OUTDIR}/${BASE}_R2_paired.fastq.gz ${OUTDIR}/${BASE}_R2_unpaired.fastq.gz \
ILLUMINACLIP:adapters.txt:2:30:10:4:true SLIDINGWINDOW:4:20 MINLEN:50

# Notes:
# trimmomatic-0.39
# Keep bases with phred-score quality > 20 in sliding window of 4 bp (average)
# Remove adapter sequences in user specified file
# Remove reads with length < 50bp following trimming
# **keepBothReads** option = True
