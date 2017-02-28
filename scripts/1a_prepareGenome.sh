#!/bin/bash
#$ -N prepareGenome
#$ -cwd
#$ -q rg-el7
#$ -o rnaSNPs.prepareGenome.out
#$ -e rnaSNPs.prepareGenome.out
#$ -M anna.vlasova@crg.eu
#$ -m e
#$ -l virtual_free=62G
#$ -l h_rt=24:00:00
#$ -pe smp 4

# 4/11/2016  - pipeline to call variants from RNAseq data. this is the first part of allele specific expression pipeline. The second part will be in the ae_pipeline.sh file.
module load STAR/2.5.1b-goolf-1.4.10-no-OFED
module load Java/1.8.0_74

picard=/software/rg/el6.3/picard-tools-2.7.1/picard.jar

#path to the newest GATK toolkit /software/rg/el6.3/GenomeAnalysisTK-3.6
#path to the newest pickard tools /software/rg/el6.3/picard-tools-2.7.1/

#this pipeline reproduces steps from the GATK best practics of SNP calling with RNAseq data
#https://software.broadinstitute.org/gatk/guide/article?id=3891

#1. set up variables - genome, transcriptome and reads
genome=/no_backup/rg/avlasova/NGS_2017/refs/Homo_sapiens.GRCh37.chromosomes.chr22.fa
genomeDir=/no_backup/rg/avlasova/NGS_2017/genome/

mkdir -p $genomeDir

#2. create STAR genome index file.

echo "STAR --runMode genomeGenerate --genomeDir $genomeDir --genomeFastaFiles $genome  --runThreadN 4"
time STAR --runMode genomeGenerate --genomeDir $genomeDir --genomeFastaFiles $genome  --runThreadN 4

#real    0m33.646s
#user    0m53.266s
#sys     0m4.483s

#3. For the Split and Trim GATK tool I need index file *.fai for genome and *.dict file for genome
# to create dictionary and index files for genome
echo "samtools faidx $genome"
time samtools faidx $genome

#real    0m0.563s
#user    0m0.332s
#sys     0m0.011s

echo "java -jar  -XX:ParallelGCThreads=8 $picard CreateSequenceDictionary R= $genome O= $genomeBase.dict"

time java -jar  -XX:ParallelGCThreads=8 $picard CreateSequenceDictionary R= $genome O= $genomeBase.dict

#real    0m3.162s
#user    0m0.602s
#sys     0m0.079s


