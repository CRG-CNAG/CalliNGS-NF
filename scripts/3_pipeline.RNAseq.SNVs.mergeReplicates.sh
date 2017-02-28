#!/bin/bash
#$ -N RNAseq_SNVs
#$ -cwd
#$ -q rg-el7,long-sl7
#$ -o rnaSNPs.mergeRep.out
#$ -e rnaSNPs.mergeRep.out
#$ -M anna.vlasova@crg.eu
#$ -m e
#$ -l virtual_free=62G
#$ -l h_rt=24:00:00
#$ -pe smp 8

module load STAR/2.5.1b-goolf-1.4.10-no-OFED
module load Java/1.8.0_74

#path to the newest GATK toolkit /software/rg/el6.3/GenomeAnalysisTK-3.6
#path to the newest pickard tools /software/rg/el6.3/picard-tools-2.7.1/

GATK=/software/rg/el6.3/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar 
genome=/no_backup/rg/avlasova/NGS_2017/refs/Homo_sapiens.GRCh37.chromosomes.chr22.fa

#For each replicate do individual mappings first, till variant calling step.

folderAnalysis=/no_backup/rg/avlasova/NGS_2017/results/
listBams=$folderAnalysis/bam.list


initVCF=$folderAnalysis/output.samtools.vcf.gz
#echo "java -jar -XX:ParallelGCThreads=8 $GATK -T HaplotypeCaller -R $genome -I $listBams -dontUseSoftClippedBases -stand_call_conf 20.0 -stand_emit_conf 20.0 -o $initVCF"

#time java -jar -XX:ParallelGCThreads=8 $GATK -T HaplotypeCaller -R $genome -I $listBams -dontUseSoftClippedBases -stand_call_conf 20.0 -stand_emit_conf 20.0 -o $initVCF
#real    28m43.085s
#user    30m5.581s
#sys     0m17.918s

# with '-L chr22:17549845-17700284' option 
#real    0m29.495s
#user    0m28.436s
#sys 0m0.512s



echo "samtools mpileup"
#time samtools mpileup -uf $genome /no_backup/rg/avlasova/NGS_2017/rep1/2pass/final.sorted.bam  /no_backup/rg/avlasova/NGS_2017/rep2/2pass/final.sorted.bam | bcftools call -vmO z -o $initVCF

time samtools mpileup -uf  $genome /no_backup/rg/avlasova/NGS_2017/rep1/2pass/final.sorted.bam  | bcftools  view -  -vmO z   > results/output.samtools.vcf
real    3m39.305s
user    3m38.451s
sys 0m0.557s

#11061 SNPs

#7.  Variant filtering
finalVCF=$folderAnalysis/final.samtools.vcf
echo "java -jar -XX:ParallelGCThreads=4 /software/rg/el6.3/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T VariantFiltration -R $genome -V $initVCF -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o $finalVCF"

#As in DNA-seq, we recommend filtering based on Fisher Strand values (FS > 30.0) and Qual By Depth values (QD < 2.0).
#We recommend that you filter clusters of at least 3 SNPs that are within a window of 35 bases between them by adding -window 35 -cluster 3 to your command. This filter recommendation is specific for RNA-seq data. 

time java -jar -XX:ParallelGCThreads=8 /software/rg/el6.3/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T VariantFiltration -R $genome -V $initVCF -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o $finalVCF

#real    0m6.147s
#user    0m12.662s
#sys     0m0.544s


