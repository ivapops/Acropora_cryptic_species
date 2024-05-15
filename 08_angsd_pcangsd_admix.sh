#!/bin/bash --login
#SBATCH --job-name="angsd"  
#SBATCH --nodes=1        
#SBATCH --ntasks-per-node=1    
#SBATCH --cpus-per-task=4	
#SBATCH --mem=400G				
#SBATCH --time=300:00:00	
#SBATCH --account=a_riginos
#SBATCH --partition=general
#SBATCH -o beagle_%A.o          
#SBATCH -e beagle_%A.e	     

# 1) ANGSD v0.934: create beagle file

module load htslib/1.12-gcc-10.3.0

export PATH="/software/angsd:$PATH"
export PATH="/software/angsd/misc/$PATH"

LIST=/analysis/angsd/ECT_only/bamList_625_Ahya
NAME=ECT_625

REF=/genome/GCA_020536085.1_Ahyacinthus.chrsV1/GCA_020536085.1_Ahyacinthus.chrsV1_genomic.fna
FAI=/genome/GCA_020536085.1_Ahyacinthus.chrsV1/GCA_020536085.1_Ahyacinthus.chrsV1_genomic.fna.fai

# < 5% missing individuals
Q5="-minMapQ 30 -minQ 30 -doMajorMinor 1 -doMaf 1 -SNP_pval 1e-6 -minInd 595 -doCounts 1 -dumpCounts 2 -setMinDepthInd 3"

# beagle output (maf=0)
angsd -gl 1 -doGlf 2 -nThreads 4 ${Q5} -rf chrs_only -bam ${LIST} -ref ${REF} -fai ${FAI} -out outQ5_${NAME}_chrs


# 2) PCAngsd: individual covariance matrix

module load anaconda3/2022.05

for BEAGLE in *chrs.beagle.gz; do
BASE=`basename ${BEAGLE} .beagle.gz`

# PCA with admix option

for K in 2 3 4 5; do
pcangsd -b ${BEAGLE} --minMaf 0.05 --admix --admix_auto 10000 --admix_K ${K} --sites_save --pi_save -o pca_${BASE}_${K}

done

# Notes:
#  --pi_save             Save individual allele frequencies
#  --dosage_save         Save genotype dosages
#  --sites_save          Save boolean vector of used sites
