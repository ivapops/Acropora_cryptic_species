#!/bin/bash --login
#SBATCH --job-name="pc2bam" 
#SBATCH --nodes=1           
#SBATCH --ntasks-per-node=1  
#SBATCH --cpus-per-task=10	
#SBATCH --mem=20G			
#SBATCH --time=10:00:00		
#SBATCH --account=a_riginos	
#SBATCH --partition=general		
#SBATCH -o pcd_%A_%a.o      
#SBATCH -e pcd_%A_%a.e       
#SBATCH --array=1-684   

module load anaconda3/2022.05
source activate picard

LIST=/analysis/lists/bamList

# Name of sample e.g., RRAP-ECO3-2021-Ahya-CBHE-218_L3

INDIR=/analysis/bwa/markedRG_bamUSE_2023
OUTDIR=/analysis/bwa/markedRG_bamUSE_2023
REF=/genome/GCA_020536085.1_Ahyacinthus.chrsV1/GCA_020536085.1_Ahyacinthus.chrsV1_genomic.fna

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
