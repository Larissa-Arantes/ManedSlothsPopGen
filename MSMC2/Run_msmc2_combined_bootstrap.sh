#!/bin/bash
#SBATCH -J msmc2
#SBATCH --mail-type=END,FAIL
#...

export PATH="/scratch/larissasa/miniconda/bin:${PATH}"
source activate psmc
module load BCFtools/1.11-GCC-10.2.0

#################################################################################################################################
####################################        Installation instructions        ####################################################
#################################################################################################################################

#Installation
#wget https://github.com/stschiff/msmc2/releases/download/v2.1.4/msmc2_Linux
#chmod 775 msmc2_Linux
#Get following scripts:
#https://github.com/stschiff/msmc-tools/blob/master/vcfAllSiteParser.py
#https://github.com/stschiff/msmc-tools/blob/master/generate_multihetsep.py

#################################################################################################################################
######################## 1) Generate input using vcf files and mask files for each chromosome ###################################
#################################################################################################################################

ind="BA"
inputs="/scratch/larissasa/bradypus/09_Resequencing/08_msmc/combinedBA_4ind_onlyindpair/input"
MASK_GENOME="/scratch/larissasa/bradypus/09_Resequencing/08_msmc/mappabilityMask"
path="/scratch/larissasa/bradypus/09_Resequencing/08_msmc"
outPath="/scratch/larissasa/bradypus/09_Resequencing/08_msmc/combinedBA_4ind_onlyindpair/output"
outBoot="/scratch/larissasa/bradypus/09_Resequencing/08_msmc/combinedBA_4ind_onlyindpair/bootstrap"

mkdir /scratch/larissasa/bradypus/09_Resequencing/08_msmc/combinedBA_4ind_onlyindpair
mkdir ${inputs}

for C in {1..25};do
	python3 /scratch/larissasa/programs/msmc2/generate_multihetsep.py --chr ${C} --mask=${MASK_GENOME}/chromosome_${C}.bed.gz --mask=${path}/BTPF17/inputs/BTPF17_mask_${C}.bed.gz --mask=${path}/BTPF22/inputs/BTPF22_mask_${C}.bed.gz --mask=${path}/BTPF26/inputs/BTPF26_mask_${C}.bed.gz  --mask=${path}/BTPF27/inputs/BTPF27_mask_${C}.bed.gz  ${path}/BTPF17/inputs/BTPF17_VCFout_${C}.vcf.gz ${path}/BTPF22/inputs/BTPF22_VCFout_${C}.vcf.gz ${path}/BTPF26/inputs/BTPF26_VCFout_${C}.vcf.gz ${path}/BTPF27/inputs/BTPF27_VCFout_${C}.vcf.gz > ${inputs}/MSMC2_input_${C}.txt
done

#################################################################################################################################
############################################ 2) Run bootstrap ###################################################################
#################################################################################################################################

mkdir ${outBoot}
python3 /scratch/larissasa/programs/msmc2/multihetsep_bootstrap.py -n 50 -s 20000000 --chunks_per_chromosome 10 --nr_chromosomes 24 ${outBoot}/${ind}.bootstrap ${inputs}/MSMC2_input_*.txt

#################################################################################################################################
################################################### 3) Run MSMC2 ################################################################
#################################################################################################################################

cd ${outBoot}
for C in {1..50};do
	echo "running msmc2 on bootstraps for ${ind} rep${C}"
	/scratch/larissasa/programs/msmc2/msmc2_Linux -t 16 -p 1*2+25*1+1*2+1*3 -i 20 -o ${outBoot}/${ind}_MSMC_output_boot_rep${C} -I 0-1,2-3,4-5,6-7 ${outBoot}/${ind}.bootstrap_${C}/*txt
done

#-p, --timeSegmentPattern=<string> : pattern of fixed time segments [default=1*2+25*1+1*2+1*3]
#-i, --maxIterations=<size_t> : number of EM-iterations [default=20]
#-t, --nrThreads=<size_t> : nr of threads to use (defaults to nr of CPUs)
#-I, --pairIndices: a single comma-separated list like this "-I 0,1,4,5", the program will run over all pairs of haplotypes within #this set of indices.  A list of pairs, like this: "-I 0-1,2-3,4-5", program will run only those specified pairs.

#wget https://github.com/jessicarick/msmc2_scripts/blob/master/msmc_2_generateInput_singleInd.sh
#chmod 775 /scratch/larissasa/programs/msmc_2_generateInput_singleInd.sh

