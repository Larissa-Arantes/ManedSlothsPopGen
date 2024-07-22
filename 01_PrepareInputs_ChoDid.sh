#!/bin/bash
#SBATCH -J snapp
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=larissasa@zedat.fu-berlin.de
#SBATCH --partition=begendiv,main
#SBATCH --qos=standard
#SBATCH --constraint=no_gpu
#SBATCH --cpus-per-task=20
#SBATCH --mem=2G
#SBATCH --time=2-00:00:00
#SBATCH -o /scratch/larissasa/bradypus/09_Resequencing/12_Snapp/ChoDid/logs/snapp.%j.out
#SBATCH -e /scratch/larissasa/bradypus/09_Resequencing/12_Snapp/ChoDid/logs/snapp.%j.err

export LMOD_DISABLE_SAME_NAME_AUTOSWAP="no"

######################################################################################
#  1) Choose genomes and run mapping and SNP calling using output                    
######################################################################################

#Run jATG SNP calling for all species against ChoDid
#jATG: /scratch/larissasa/programs/jATG_3/4.snp_calling
#BraVar GCA_004027775
#Cho Did GCF_015220235
#Cho Hof GCA_000164785.2
#Bra Tor RJ BT04
#Bra Tor BA BTPF27
#Bra Tor BA M184BA

######################################################################################
#  2) Prepare VCF
######################################################################################

# 2.1) Merge VCFs
# VCF should NOT have masked regions, sex chromosome and mitogenome -> already done by jATG

cd /scratch/larissasa/bradypus/09_Resequencing/12_Snapp/ChoDidOnly

#module load BCFtools/1.11-GCC-10.2.0
#bcftools merge --file-list vcf_list_ChoDidRef -O z -o Combined_ChoDidRef.vcf.gz

# 2.2) Filter VCF (following https://github.com/ForBioPhylogenomics/tutorials/blob/main/species_tree_inference_with_snp_data/README.md)
# Exclude all sites at which no alternative alleles are called for any of the individuals ("AC==0"), 
# all sites at which only alternative alleles are called ("AC==AN"), 
# and missing sites ("F_MISSING > 0").
# ensure that only bi-allelic SNPs are included in the reduced dataset by using the BCFtools options -m2 and -M2 which set both the minimum (-m) and maximum (-M) number of alleles to 2.

#TRY NO MISSING ALLOWED
#bcftools view -e 'AC==0 || AC==AN || F_MISSING > 0' -m2 -M2 -O z -o Combined_ChoDidRef_filtered.vcf.gz Combined_ChoDidRef.vcf.gz

# 2.3) Thin SNPs
#The autosomal SNPs were thinned with a 10-kb window size to minimize the influence of linkage imbalance between loci.

#plink2 --vcf Combined_ChoDidRef_filtered.vcf.gz --double-id --allow-extra-chr --set-all-var-ids @:# --threads 10 --parallel 1 10 --bad-ld --indep-pairwise 50 5 0.01 --out LDpruning_ChoDid

# edit plink output file for use in vcftools
#awk '{gsub(":", "\t");print}' LDpruning_ChoDid.prune.in > LDpruning_ChoDid_edited.prune.in

# remove SNPs under linkage disequilibrium from all SNPs VCF file
#vcftools --gzvcf Combined_ChoDidRef_filtered.vcf.gz --positions LDpruning_ChoDid_edited.prune.in --recode --out Combined_ChoDidRef_filtered_LDpruning

#Test pruning with vcftools thin
#vcftools --gzvcf Combined_ChoDidRef_filtered.vcf.gz --thin 50000 --recode --out Combined_ChoDidRef_pruned

#Generate vcf without ChoHof
#vcftools --vcf Combined_ChoDidRef_pruned100000.recode.vcf --remove-indv ChoHof --recode --out Combined_ChoDidRef_pruned100000_NoCHoHof

######################################################################################
#  3) Prepare xml file for Beast
######################################################################################

#Following https://github.com/mmatschiner/tutorials/tree/master/divergence_time_estimation_with_snp_data#q1

#mkdir Test1
#module add Ruby/3.3.0-GCCcore-12.3.0
#ruby snapp_prep.rb -v Combined_ChoDidRef_pruned100000_NoCHoHof.recode.vcf -t samples.txt -c constraints.txt -s phylogeny.ChoDid.newick -l 10000000 -o Test1/Test1

#    -v, --vcf FILENAME               File with SNP data in vcf format (default: none).
#    -t, --table FILENAME             File with table linking species and individuals (default: example.spc.txt).
#    -c, --constraints FILENAME       File with age constraint information (default: example.con.txt).
#    -s, --starting-tree FILENAME     File with starting tree in Nexus or Newick format (default: none).
#    -l, --length LENGTH              Number of MCMC generations (default: 500000).
#    -m, --max-snps NUMBER            Maximum number of SNPs to be used (default: no maximum).

# do snapp_prep.rb filter by missing data 0???

#To limit the dataset to 1,000 randomly selected SNPs and to set a chain length of 500,000 MCMC iterations,
#You may notice that the chain length of 100,000 MCMC iterations is extremely short compared to those used in the BEAST2 analyses of other tutorials. 
#Using much shorter chain lengths with SNAPP than BEAST2 is quite common, given that the SNAPP model has far fewer model parameters than most models 
#used in other BEAST2 analyses, and that each individual iteration is much slower with SNAPP due to the high computational demand of the integration over all possible gene trees at each SNP.

#Starting tree
#Note that besides the requirement that the tree should be compatible with all constraints, the topology and branch lengths chosen for the starting tree 
#should not have any effect on the outcome of the SNAPP analysis and can therefore be chosen arbitrarily.

######################################################################################
#  4) Run Beast                                                
######################################################################################

module load Beast/2.7.3-GCC-11.3.0
#packagemanager -add SNAPP
#packagemanager -add snapper
#examples files here: /trinity/shared/easybuild/software/Beast/2.7.3-GCC-11.3.0/examples/
beast -threads 20 snapp.xml 

