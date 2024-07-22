#!/bin/bash
#SBATCH -J snpEff
#SBATCH --mail-type=END,FAIL
...

cd /scratch/larissasa/bradypus/09_Resequencing/11_snpEff/02_snpeff_maxmiss_0.8
VCF="/scratch/larissasa/bradypus/09_Resequencing/11_snpEff/02_snpeff_maxmiss_0.8/BraTor_pop_maxmiss_0.8.ann.vcf"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#    Remove sites with warnings   #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

grep "WARNING_SEQUENCE_NOT_AVAILABLE\|WARNING_TRANSCRIPT_INCOMPLETE\|WARNING_TRANSCRIPT_MULTIPLE_STOP_CODONS\|WARNING_TRANSCRIPT_NO_START_CODON" ${VCF} | sed -e '1d' | cut -f1,2 > warning_sites.txt

wc -l warning_sites.txt
# X warning_sites.txt

# Filter from annotated vcf file

vcftools --vcf ${VCF} --recode --recode-INFO-all --exclude-positions warning_sites.txt --out BraTor_pop.ann_filtered

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#        Get site categories        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

grep "LOF" BraTor_pop.ann_filtered.recode.vcf | sed -e '1d' | cut -f1,2 > lof_sites.txt

grep "missense_variant" BraTor_pop.ann_filtered.recode.vcf | sed -e '1d' | cut -f1,2 > missense_sites.txt

grep "synonymous_variant" BraTor_pop.ann_filtered.recode.vcf | sed -e '1d' | cut -f1,2 > synonymous_sites.txt

grep "intergenic_region" BraTor_pop.ann_filtered.recode.vcf | sed -e '1d' | cut -f1,2 > intergenic_sites.txt
  
# Subsample intergenic sites:

get_seeded_random()
{
  seed="$1"
  openssl enc -aes-256-ctr -pass pass:"$seed" -nosalt \
    </dev/zero 2>/dev/null
}

shuf -n 100000 --random-source=<(get_seeded_random 42) intergenic_sites.txt > intergenic_subsampled_sites.txt 

wc -l *_sites.txt


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#        Extract genotypes          #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Get vcf file for each site category

for i in lof missense synonymous intergenic intergenic_subsampled
do
	vcftools --vcf BraTor_pop.ann_filtered.recode.vcf --recode --recode-INFO-all --positions ${i}_sites.txt --out BraTor_pop.ann_filtered_${i}_snps
done

# Convert vcf to plink
# Get plink traw file of 012 genotypes for all inds and each SNP category

for i in missense synonymous lof intergenic intergenic_subsampled
do
	plink2 --vcf BraTor_pop.ann_filtered_${i}_snps.recode.vcf --export A-transpose --allow-extra-chr --out BraTor_pop.ann_filtered_${i}_snps
done


#~~ Get allele frequencies within populations

for i in missense synonymous lof intergenic intergenic_subsampled
do
	plink2 --vcf BraTor_pop.ann_filtered_${i}_snps.recode.vcf --keep Bahia.txt --freq --out BraTor_pop.ann_filtered_${i}_snps_Bahia --allow-extra-chr --debug
done

for i in missense synonymous lof intergenic intergenic_subsampled
do
	plink2 --vcf BraTor_pop.ann_filtered_${i}_snps.recode.vcf --keep RiodeJaneiro.txt --freq --out BraTor_pop.ann_filtered_${i}_snps_RiodeJaneiro --allow-extra-chr --debug
done
