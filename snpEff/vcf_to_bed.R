library(dplyr)
library(purrr)
library(tidyr)
library(data.table)

# Set options to avoid scientific notation
options(scipen = 999)

#~~ Create bed file of SNPs to extract from consensus fasta

map <- fread("/scratch/larissasa/bradypus/09_Resequencing/11_snpEff/VCF_PASS/AllInd_Genot.PASS.vcf.SNPID.biallellic.maxmiss0.8.map")

bed <- map %>%
	select(V1, V4) %>%
   	mutate(start = V4 - 1,
	              name = "SNP") %>%
        select(V1, start, V4, name)

write.table(bed, "/scratch/larissasa/bradypus/09_Resequencing/11_snpEff/VCF_PASS/AllInd_Genot.PASS.vcf.SNPID.biallellic.maxmiss0.8.bed",
	    quote = F, row.names = F, col.names = F, sep = "\t")

