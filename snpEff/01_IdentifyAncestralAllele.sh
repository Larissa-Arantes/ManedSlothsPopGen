#!/bin/bash
#SBATCH -J IndAncAllele
#SBATCH --mail-type=END,FAIL
...

export LMOD_DISABLE_SAME_NAME_AUTOSWAP="no"

cd /scratch/larissasa/bradypus/09_Resequencing/11_snpEff/01_snpeff/01_IdentifyAncestralAllele

#######################################################################################
# We used the bedtools getFasta command to generate a FASTA file made of short (70 bp) 
# subsequences covering the entirety of each of these three reference genomes. For each 
# chromosome, we extracted 70 bp fragments with fragment starting points separated by 10 bp 
# (adjacent fragments were tiled such that they overlapped by 60 bp) to increase the   
# proportion of sites successfully mapped to the reference genome
#######################################################################################

# 1) Create a fasta only with chromosome scaffolds
# Genomes should be hard masked

seqkit grep -n -f ChoDidChrom.txt /scratch/larissasa/Choloepus/jATG/Choloepus_didactylus/mChoDid/2.masking/3_masker/GCF_015220235.1_mChoDid1.pri_genomic.HM.fa > mChoDid1_mainChrom_masked.fasta
seqkit grep -n -f DasNovChrom.txt /scratch/larissasa/bradypus/03_OtherGenomes/Dasnov/jATG/Dasypus_novemcinctus/mDasNov/2.masking/3_masker/GCF_030445035.1_mDasNov1.hap2_genomic.HM.fa > mDasNov1_mainChrom_masked.fasta
seqkit grep -n -f TamTetChrom.txt /scratch/larissasa/bradypus/03_OtherGenomes/TamTet/jATG/Tamandua_tetradactyla/mTamTet/2.masking/3_masker/GCA_023851605.1_mTamTet1.pri_genomic.HM.fa > mTamTet1_mainChrom_masked.fasta

# 2) Linearize fasta file

seqkit seq -w 0 mChoDid1_mainChrom_masked.fasta > mChoDid1_mainChrom_masked_linear.fasta
seqkit seq -w 0 mDasNov1_mainChrom_masked.fasta > mDasNov1_mainChrom_masked_linear.fasta
seqkit seq -w 0 mTamTet1_mainChrom_masked.fasta > mTamTet1_mainChrom_masked_linear.fasta
rm mChoDid1_mainChrom_masked.fasta 
rm mDasNov1_mainChrom_masked.fasta
rm mTamTet1_mainChrom_masked.fasta

# 3) Remove masked sequence (Ns) by breaking non-masked regions in new sequences

python RemoveMaskedRegions.py mChoDid1_mainChrom_masked_linear.fasta mChoDid1_mainChrom_masked_pruned.fa
python RemoveMaskedRegions.py mDasNov1_mainChrom_masked_linear.fasta mDasNov1_mainChrom_masked_pruned.fa
python RemoveMaskedRegions.py mTamTet1_mainChrom_masked_linear.fasta mTamTet1_mainChrom_masked_pruned.fa

# 4) find and replace lower case (soft masked) sites in each fasta file

sed '/^>/!s/g/G/g; /^>/!s/a/A/g; /^>/!s/t/T/g; /^>/!s/c/C/g' mChoDid1_mainChrom_masked_pruned.fa > mChoDid1_mainChrom_masked_pruned_cap.fa
sed '/^>/!s/g/G/g; /^>/!s/a/A/g; /^>/!s/t/T/g; /^>/!s/c/C/g' mDasNov1_mainChrom_masked_pruned.fa > mDasNov1_mainChrom_masked_pruned_cap.fa
sed '/^>/!s/g/G/g; /^>/!s/a/A/g; /^>/!s/t/T/g; /^>/!s/c/C/g' mTamTet1_mainChrom_masked_pruned.fa > mTamTet1_mainChrom_masked_pruned_cap.fa
rm mChoDid1_mainChrom_masked_pruned.fa
rm mDasNov1_mainChrom_masked_pruned.fa
rm mTamTet1_mainChrom_masked_pruned.fa

# 5) index the fasta files

module load SAMtools/1.17-GCC-12.2.0
samtools faidx mChoDid1_mainChrom_masked_pruned_cap.fa
samtools faidx mDasNov1_mainChrom_masked_pruned_cap.fa
samtools faidx mTamTet1_mainChrom_masked_pruned_cap.fa

# 6) get chromosome sizes

cut -f1,2  mChoDid1_mainChrom_masked_pruned_cap.fa.fai > chromSizes.mChoDid1_mainChrom_masked_pruned_cap.fasta
cut -f1,2  mDasNov1_mainChrom_masked_pruned_cap.fa.fai > chromSizes.mDasNov1_mainChrom_masked_pruned_cap.fasta
cut -f1,2  mTamTet1_mainChrom_masked_pruned_cap.fa.fai > chromSizes.mTamTet1_mainChrom_masked_pruned_cap.fasta

# 7) use Rscript makeReadBed_X.R to create a bed file with coordinates for short fragments (70 bp / 10 bp overlapping)
# create one Rscript for each genome

module load R
Rscript makeReadBed_ChoDid.R
Rscript makeReadBed_TamTet.R
Rscript makeReadBed_DasNov.R

# 8) Run bedtools to fasta fasta short reads.

module load BEDTools/2.30.0-GCC-11.3.0
bedtools getfasta -fi mChoDid1_mainChrom_masked_pruned_cap.fa -bed ChoDid_read.bed > ChoDid_reads.fasta
bedtools getfasta -fi mDasNov1_mainChrom_masked_pruned_cap.fa -bed DasNov_read.bed > DasNov_reads.fasta
bedtools getfasta -fi mTamTet1_mainChrom_masked_pruned_cap.fa -bed TamTet_read.bed > TamTet_reads.fasta


#######################################################################################
# Convert these FASTA files to FASTQ format 
#######################################################################################

python FastaToFastq.py ChoDid_reads.fasta ChoDid_reads.fastq
python FastaToFastq.py DasNov_reads.fasta DasNov_reads.fastq
python FastaToFastq.py TamTet_reads.fasta TamTet_reads.fastq

#######################################################################################
# Align the sequence data to our new maned sloth reference genome using BWA
#######################################################################################

module load BWA/0.7.17-GCCcore-11.2.0
bwa index /scratch/larissasa/bradypus/10_RefereceGenome/mBraTor.curated.fasta mBraTor.curated                     
bwa mem -B 3 /scratch/larissasa/bradypus/10_RefereceGenome/mBraTor.curated.fasta ChoDid_reads.fastq > ChoDid.sam
bwa mem -B 3 /scratch/larissasa/bradypus/10_RefereceGenome/mBraTor.curated.fasta DasNov_reads.fastq > DasNov.sam
bwa mem -B 3 /scratch/larissasa/bradypus/10_RefereceGenome/mBraTor.curated.fasta TamTet_reads.fastq > TamTet.sam

#######################################################################################
# Convert the aligned reads from SAM to BAM format and then sorted the BAM files using SAMtools
#######################################################################################

module load SAMtools/1.17-GCC-12.2.0
samtools view -S -b ChoDid.sam > ChoDid.bam
samtools sort ChoDid.bam -o ChoDid_sorted.bam
rm ChoDid.sam 

samtools view -S -b DasNov.sam > DasNov.bam
samtools sort DasNov.bam -o DasNov_sorted.bam
rm DasNov.sam

samtools view -S -b TamTet.sam > TamTet.bam
samtools sort TamTet.bam -o TamTet_sorted.bam
rm TamTet.sam

#######################################################################################
# Use SAMtools to merge bam files
#######################################################################################

samtools merge xenarthra_merged.bam ChoDid_sorted.bam DasNov_sorted.bam TamTet_sorted.bam

#######################################################################################
# Get consensus fasta from merged bam 
#######################################################################################

angsd -i xenarthra_merged.bam -P 4 -doFasta 2 -doCounts 1 -out outgroup_consensus.fasta
gunzip outgroup_consensus.fasta.fa.gz

#######################################################################################
# Back to the VCF of population data
# Extract SNP positions from filtered VCF file and convert to bed file
#######################################################################################

# Add SNP ID to VCF
python add_SNP_ID_to_VCF.py

# Convert vcf to ped/map
plink2 --vcf /scratch/larissasa/bradypus/09_Resequencing/11_snpEff/VCF_PASS/AllInd_Genot.PASS.vcf.SNPID.vcf --allow-extra-chr --recode ped --out /scratch/larissasa/bradypus/09_Resequencing/11_snpEff/VCF_PASS/AllInd_Genot.PASS.vcf.SNPID

# Keep only biallellic snsp
plink2 --vcf /scratch/larissasa/bradypus/09_Resequencing/11_snpEff/VCF_PASS/AllInd_Genot.PASS.vcf.SNPID.vcf --max-alleles 2 --allow-extra-chr --export vcf --threads 10 --out /scratch/larissasa/bradypus/09_Resequencing/11_snpEff/VCF_PASS/AllInd_Genot.PASS.vcf.SNPID.biallellic

# Convert VCF to bed
module load R
Rscript vcf_to_bed.R

#######################################################################################
# Extract ancestral alleles from outgroup consensus
#######################################################################################

module load BEDTools/2.30.0-GCC-11.3.0
bedtools getfasta -fi outgroup_consensus.fasta.fa -bed /scratch/larissasa/bradypus/09_Resequencing/11_snpEff/VCF_PASS/AllInd_Genot.PASS.bed -fo ancestral_alleles.out

#######################################################################################
# Format ancestral allele file (from fasta format to a two-column file with SNP ID and reference allele)
#######################################################################################

python format_ancestral_alleles_v4.py

#######################################################################################
# Use Plink to polarize alleles - remove SNPs with Ns in alignment (and therefore can't infer) 
# Switch ref alleles when they are wrong
# Extract SNPs that could be identified in the alignment and rotate
# allele codes where necessary
#######################################################################################

awk '{print $1}' ancestral_alleles.snpid.txt > ancestral_positions.snpid.txt

plink2 --vcf /scratch/larissasa/bradypus/09_Resequencing/11_snpEff/VCF_PASS/AllInd_Genot.PASS.vcf.SNPID.biallellic.vcf.gz --threads 10 --extract ancestral_positions.snpid.txt --ref-allele force ancestral_alleles.snpid.txt 2 1 --recode vcf --allow-extra-chr --out AllInd_Genot.PASS.SNPID.biallellic.polarized.temp

#######################################################################################
# Remove these SNPs with ref allele mismatches i.e. ancestral allele
# from alignment matched neither of the SNPs
#######################################################################################

sed -i '/variant$/{N;s/\n/ /}' AllInd_Genot.PASS.SNPID.biallellic.polarized.temp.log
grep 'Warning' AllInd_Genot.PASS.SNPID.biallellic.polarized.temp.log | awk {'print $7'} | sed 's/.//;s/.$//' | sed 's/.$//' > mismatches.txt

plink2 --vcf AllInd_Genot.PASS.SNPID.biallellic.polarized.temp.vcf --exclude mismatches.txt --allow-extra-chr --export vcf-4.2 --threads 10 --out AllInd_Genot.PASS.SNPID.biallellic.polarized --allow-misleading-out-arg

rm AllInd_Genot.PASS.polarized.temp.vcf
