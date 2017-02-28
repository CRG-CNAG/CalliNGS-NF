#!/bin/bash
#$ -N AE_GATK
#$ -cwd
#$ -q rg-el7,long-sl7
#$ -o ae_gatk.out
#$ -e ae_gatk.err
#$ -M anna.vlasova@crg.eu
#$ -m e
#$ -l virtual_free=62G
#$ -l h_rt=24:00:00
#$ -pe smp 4

module load STAR/2.5.1b-goolf-1.4.10-no-OFED
module load Java/1.8.0_74

GATK=/software/rg/el6.3/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar 
genome=/no_backup/rg/avlasova/NGS_2017/refs/Homo_sapiens.GRCh37.chromosomes.chr22.fa

bam1=/no_backup/rg/avlasova/NGS_2017/rep1/2pass/final.sorted.bam
vcfFile=/no_backup/rg/avlasova/NGS_2017/out.recode.vcf
resultFile=/no_backup/rg/avlasova/NGS_2017/results/ASER.out

echo "java -jar -XX:ParallelGCThreads=8 $GATK -R $genome  -T ASEReadCounter -o $resultFile -I $bam1 -sites $vcfFile"

time java -jar -XX:ParallelGCThreads=8 $GATK -R $genome  -T ASEReadCounter -o $resultFile -I $bam1 -sites $vcfFile
#real    1m11.729s
#user    1m14.085s
#sys     0m1.860s

#resulting file contains following fields
#refAllele: reference allele
#altAllele: alternate allele
#refCount: number of reads that support the reference allele
#altCount: number of reads that support the alternate allele
#totalCount: total number of reads at the site that support both reference and alternate allele and any other alleles present at the site
#lowMAPQDepth: number of reads that have low mapping quality
#lowBaseQDepth: number of reads that have low base quality
#rawDepth: total number of reads at the site that support both reference and alternate allele and any other alleles present at the site
#otherBases: number of reads that support bases other than reference and alternate bases
#improperPairs: number of reads that have malformed pairs

