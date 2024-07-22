#!/bin/bash
#SBATCH -J snpEff
#SBATCH --mail-type=END,FAIL
...

export LMOD_DISABLE_SAME_NAME_AUTOSWAP="no"

export PATH="/scratch/larissasa/miniconda/bin:${PATH}"
source activate psmc
module load BCFtools/1.11-GCC-10.2.0

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#         Prepare VCF             #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

VCF="/scratch/larissasa/bradypus/09_Resequencing/11_snpEff/01_snpeff/01_IdentifyAncestralAllele/AllInd_Genot.PASS.SNPID.biallellic.polarized.vcf"
VCFmiss="/scratch/larissasa/bradypus/09_Resequencing/11_snpEff/01_snpeff/01_IdentifyAncestralAllele/AllInd.Genot.PASS.SNPID.biallellic.polarized.maxmiss0.8"
prefix="BraTor_pop_maxmiss_0.8"
OUT="/scratch/larissasa/bradypus/09_Resequencing/11_snpEff/02_snpeff_maxmiss_0.8"

cd ${OUT}
vcftools --vcf ${VCF} --max-missing 0.8 --recode --out ${VCFmiss}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#          Run snpEff             #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

cd /scratch/larissasa/programs/snpEff
java -Xmx16g -jar snpEff.jar BtorquatusToga ${VCFmiss}.recode.vcf -stats ${OUT}/${prefix}_summary.html -csvStats ${OUT}/${prefix}_summary.csv > ${OUT}/${prefix}.ann.vcf

#SnpEff added functional annotations in the ANN info field (eigth column in the VCF output file)
