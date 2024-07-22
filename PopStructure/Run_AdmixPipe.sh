

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#          Installation           #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#https://github.com/stevemussmann/admixturePipeline
#module load singularity
#singularity pull docker://mussmann/admixpipe:3.2
#git clone https://github.com/stevemussmann/admixturePipeline.git
#cd admixturePipeline
#singularity run -B /scratch/larissasa/programs -B /scratch/larissasa/bradypus /scratch/larissasa/programs/admixpipe_3.2.sif
#-B to bind correct folders is essential for recognizing softwares from the container (it should read python from the container and not my python)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#     Iniciate screen             #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

screen
srun -c 10 --mem 10G -p begendiv,main -q prio -t 24:00:00 --pty /bin/bash
module load singularity
singularity run -B /scratch/larissasa/programs/admixturePipeline/ -B /localscratch/tmp/larissasa -B /scratch/larissasa/bradypus /scratch/larissasa/programs/admixpipe_3.2.sif
export PATH=$PATH:/app/bin/

cd /scratch/larissasa/bradypus/09_Resequencing/06_PopStructure

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#         Run admixture           #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#Load AdmixPipe v3
/scratch/larissasa/programs/admixturePipeline/admixturePipeline.py -m /scratch/larissasa/bradypus/09_Resequencing/x_scripts/popmap.txt -v /scratch/larissasa/bradypus/09_Resequencing/04_VCFcombined/AllInd.Genot.PASS.biallelic.maxmiss0.8.LDprunning.vcf -k 1 -K 6 -n 10 -c 10 -R 10

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#    Run CLUMPAK                  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#CLUMPAK identifies multimodal population structure results from replicates performed in Admixture (i.e. major and minor clusters),

/scratch/larissasa/programs/admixturePipeline/submitClumpak.py -p AllInd.Genot.PASS.biallelic.maxmiss0.8.LDprunning -M

#-p / --prefix: Specify the prefix from your admixpipe run. This will either be the prefix of your plink file (i.e., without .bed or .map extension), or same name as your VCF file without the .vcf file extension.
#-M / --mainpipeline: Run the CLUMPAK main pipeline locally in the Docker container (or on your computer).

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#         Run distruct            #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

/scratch/larissasa/programs/admixturePipeline/distructRerun.py -a /scratch/larissasa/bradypus/09_Resequencing/06_PopStructure -d /scratch/larissasa/bradypus/09_Resequencing/06_PopStructure/clumpakOutput/ -k 1 -K 6

#-a is used to provide the path to the directory of results produced by admixturePipeline.py, 
#-d is used to give the name of the directory that the program will use as input, 
#-k (lower case) specifies the lowest clustering value that you tested in Admixture, and 
#-K (upper case) specifies the highest clustering value you tested

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#      Cross-validation           #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#Summarize the variability of cross-validation values and log likelihood scores across multiple runs of admixture.

/scratch/larissasa/programs/admixturePipeline/cvSum.py

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#      Choose best K              #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#This module aids in evaluating appropriateness of K values by running the evalAdmix program on each Admixture replicate, as well as the summarized .Q outputs from CLUMPAK for all major and minor clusters
#EvalAdmix algorithm calculates the correlation of residuals between actual and predicted genotypes for a given Admixture model, with an optimal K-value represented by correlations approximating 0 

/scratch/larissasa/programs/admixturePipeline/runEvalAdmix.py -p AllInd.Genot.PASS.biallelic.maxmiss0.8.LDprunning -k 1 -K 6 -m /scratch/larissasa/bradypus/09_Resequencing/x_scripts/popmap.txt -n 5

#-p / --prefix: Specify your .ped or .bed file prefix from your initial admixturePipeline.py run. If you input a VCF file, this will be the name of that VCF file, except without the .vcf extension (required).
#-k / --minK: Specify the minimum K value (required).
#-K / --maxK: Specify the maximum K value (required).
#-n / --np: Provide the number of processors to use for evalAdmix (optional; default=1).
