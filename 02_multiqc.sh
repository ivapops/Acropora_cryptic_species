#!/bin/bash --login
#SBATCH --job-name="multiqc"    
#SBATCH --nodes=1				
#SBATCH --ntasks-per-node=1     
#SBATCH --cpus-per-task=1		
#SBATCH --mem=10G				
#SBATCH --time=10:00:00			
#SBATCH --account=a_riginos		
#SBATCH --partition=general		
#SBATCH -o mfqc_%A.o            
#SBATCH -e mfqc_%A.e	        

module load anaconda3/2022.05
source activate multiqc

# use FASTQC output directory as input
multiqc output/*_fastqc.zip --interactive


