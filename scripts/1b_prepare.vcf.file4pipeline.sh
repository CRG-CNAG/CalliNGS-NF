#!/bin/bash
#$ -N prepare_vcf
#$ -cwd
#$ -q rg-el7
#$ -o prepare_vcf.out
#$ -e prepare_vcf.err
#$ -M anna.vlasova@crg.eu
#$ -m e
#$ -l virtual_free=40G

module load Java/1.8.0_74

# This script will prepare vcf file for allele read counter 
# It will filter out regions with low mapability score, those fall into blacklisted regions (Anshul's one) and those with too hight read coverage.

variantFile=/no_backup_isis/rg/avlasova/NGS_2017/snv/NA12878.chr22.vcf.gz
blacklisted=/no_backup_isis/rg/avlasova/NGS_2017/refs/ENCFF001TDO.sorted.bed

genome=/no_backup_isis/rg/avlasova/NGS_2017/refs/Homo_sapiens.GRCh37.chromosomes.chr22.fa

#1. filter out those SNPs from blacklisted regions

#I can use either vcftools here
/software/rg/el6.3/vcftools-0.1.12a/bin/vcftools --gzvcf $variantFile --out $variantFile.filtered --exclude-bed $blacklisted --recode

#After filtering, kept XX out of a possible XX Sites

#or by using GATK tools -- for this I need genome reference file
#java -jar -XX:ParallelGCThreads=4 /software/rg/el6.3/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T VariantFiltration -R $genome --out $variantFile.filtered1.vcf --mask $blacklisted --variant $variantFile --maskName 'blacklisted'
#XX snp was masked

#2. Filter homozygous SNPs -  I will not do it so far, lets see how many of them I'll get with RNAseq data

#java -jar -XX:ParallelGCThreads=4 /software/rg/el6.3/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T VariantFiltration -R $genome --out $variantFile.filtered2.vcf --variant $variantFile.filtered1.vcf --genotypeFilterExpression "isHomRef==1||isHomVar==1"  --genotypeFilterName  'homozygous'

# within non-overlaping exonic regions 
#java -jar -XX:ParallelGCThreads=4 /software/rg/el6.3/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T VariantFiltration -R $genome --out $variantFile.filtered.vcf --mask $exonic_regions --variant $variantFile.filtered2.vcf --maskName 'notExonic' --filterNotInMask 

#remove intermediate files

#rm $variantFile.filtered1.vcf* $variantFile.filtered2.vcf*

#select only variants that pass all filters 
#this step requires a lot of memory, increased from 20G to 40G
#java -jar /software/rg/el6.3/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T  SelectVariants -R $genome -V $variantFile.filtered.vcf -o $variantFile.final.vcf --excludeFiltered

#(grep '#' $variantFile.final.vcf ; grep -v 'homozygous' $variantFile.final.vcf) > $variantFile.final_nohomiz.vcf


