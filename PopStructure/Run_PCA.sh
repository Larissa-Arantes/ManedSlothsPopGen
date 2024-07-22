#!/bin/bash
#SBATCH -J PCA
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=larissasa@zedat.fu-berlin.de
#SBATCH --partition=begendiv,main
#SBATCH --qos=standard
#SBATCH --cpus-per-task=10
#SBATCH --mem=1G
#SBATCH --time=24:00:00
#SBATCH -o /scratch/larissasa/bradypus/09_Resequencing/x_scripts/logs/PCA.%j.out
#SBATCH -e /scratch/larissasa/bradypus/09_Resequencing/x_scripts/logs/PCA.%j.err

export PATH="/scratch/larissasa/miniconda/bin:${PATH}"
source activate psmc

VCF="/scratch/larissasa/bradypus/09_Resequencing/04_VCFcombined/AllInd.Genot.PASS.biallelic.maxmiss0.8.LDprunning.vcf.gz"
INPATH="/scratch/larissasa/bradypus/09_Resequencing/04_VCFcombined"
OUT="/scratch/larissasa/bradypus/09_Resequencing/06_PopStructure"

# prune and create pca
plink --vcf ${VCF} --double-id --allow-extra-chr --set-missing-var-ids @:# --memory 300000 --threads 10 --make-bed --pca --out ${OUT}/Bradypus_pca
