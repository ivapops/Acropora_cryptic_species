#!/bin/bash --login
#SBATCH --job-name="pc2bam"    # job name
#SBATCH --nodes=1               # use 1 node
#SBATCH --ntasks-per-node=1     # use 1 for single and multi core jobs
#SBATCH --cpus-per-task=10		# number of cores per job
#SBATCH --mem=20G				# RAM per job given in megabytes (M), gigabytes (G), or terabytes (T)
#SBATCH --time=10:00:00			# walltime
#SBATCH --account=a_riginos		# group account name
#SBATCH --partition=general		# queue name
#SBATCH -o pcd_%A_%a.o          # standard output
#SBATCH -e pcd_%A_%a.e          # standard error
#SBATCH --array=1-684        	# job array

module load anaconda3/2022.05
source activate picard

LIST=/scratch/project/rrap_ahya/analysis/lists/bamList

# Name of sample e.g., RRAP-ECO3-2021-Ahya-CBHE-218_L3

INDIR=/scratch/project/rrap_ahya/analysis/bwa/markedRG_bamUSE_2023
OUTDIR=/scratch/project/rrap_ahya/analysis/bwa/markedRG_bamUSE_2023
REF=/scratch/project/rrap_ahya/genome/GCA_020536085.1_Ahyacinthus.chrsV1/GCA_020536085.1_Ahyacinthus.chrsV1_genomic.fna

# file name variable is associated Array index
FILENAME=`cat ${LIST} | tail -n +${SLURM_ARRAY_TASK_ID} | head -1`
BASE=`basename ${FILENAME} _UNDEDUP.bam`

# assign read group information
# Note: modify 'sed' string as needed to match sample name
picard AddOrReplaceReadGroups \
        INPUT=${INDIR}/${BASE}_UNDEDUP.bam OUTPUT=${OUTDIR}/${BASE}_RG.UNDEDUP.bam \
        RGID=$(echo $BASE | cut -d"_" -f1 | sed "s/RRAP-.*-202.*-A/A/") \
        RGLB=$(echo $BASE | cut -d"-" -f1,2,3) \
        RGPL=Illumina \
        RGPU=$(echo $BASE | cut -d"_" -f2) \
        RGSM=$(echo $BASE | cut -d"_" -f1 | sed "s/RRAP-.*-202.*-A/A/") \


# mark and remove PCR duplicates
picard MarkDuplicates \
        VALIDATION_STRINGENCY=LENIENT \
        TMP_DIR=${TMPDIR}/MarkDUP \
        INPUT=${OUTDIR}/${BASE}_RG.UNDEDUP.bam \
        OUTPUT=${OUTDIR}/${BASE}_markedRG.bam \
        REMOVE_DUPLICATES=true \
        METRICS_FILE=${OUTDIR}/${BASE}-markDup_metrics.txt 

# remove temp files
# rm *_RG.UNDEDUP.bam
