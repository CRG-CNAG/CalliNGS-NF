#!/bin/bash
#$ -N RNAseq_SNVs
#$ -cwd
#$ -q rg-el7,long-sl7
#$ -t 1-2
#$ -o rnaSNPs.rep$TASK_ID.newOptions.out
#$ -e rnaSNPs.rep$TASK_ID.newOptions.out
#$ -M anna.vlasova@crg.eu
#$ -m e
#$ -l virtual_free=62G
#$ -l h_rt=24:00:00
#$ -pe smp 8

# pipeline to call variants from RNAseq data. This is the first part of allele specific expression pipeline. 

module load STAR/2.5.1b-goolf-1.4.10-no-OFED
module load Java/1.8.0_74

GATK=/software/rg/el6.3/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar 
picard=/software/rg/el6.3/picard-tools-2.7.1/picard.jar

SEEDFILE=/no_backup/rg/avlasova/NGS_2017/input_data.txt
fileId=$(cat $SEEDFILE | head -n $SGE_TASK_ID | tail -n 1)

#fileId='rep1'
#this pipeline reproduces steps from the GATK best practics of SNP calling with RNAseq data
#https://software.broadinstitute.org/gatk/guide/article?id=3891

#1. setup variables - genome, transcriptome and reads
genome=/no_backup/rg/avlasova/NGS_2017/refs/Homo_sapiens.GRCh37.chromosomes.chr22.fa
genomeBase=/no_backup/rg/avlasova/NGS_2017/refs/Homo_sapiens.GRCh37.chromosomes.chr22
genomeDir=/no_backup/rg/avlasova/NGS_2017/genome/

read1=/no_backup/rg/avlasova/NGS_2017/rnaseq/$fileId\_1.fq.gz
read2=/no_backup/rg/avlasova/NGS_2017/rnaseq/$fileId\_2.fq.gz

# setup folders
mkdir -p $fileId
runDir_pass1=/no_backup/rg/avlasova/NGS_2017/$fileId/1pass/
mkdir -p $runDir_pass1
genomeDir2=/no_backup/rg/avlasova/NGS_2017/$fileId/genome2/
mkdir -p $genomeDir2
runDir_pass2=/no_backup/rg/avlasova/NGS_2017/$fileId/2pass/
mkdir -p $runDir_pass2

#2. Align reads to genome
cd $runDir_pass1
echo "STAR --genomeDir $genomeDir --readFilesIn $read1 $read2 --runThreadN 8 --readFilesCommand zcat"

time STAR --genomeDir $genomeDir --readFilesIn $read1 $read2 --runThreadN 8 --readFilesCommand zcat --outFilterType BySJout --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999

#real    1m28.599s
#user    8m33.088s
#sys     0m3.822s

# 2d pass STAR - improve alignmnets using table of splice junctions. create a new index

echo "STAR --runMode genomeGenerate --genomeDir $genomeDir2 --genomeFastaFiles $genome --sjdbFileChrStartEnd $runDir_pass1/SJ.out.tab --sjdbOverhang 75 --runThreadN 8"

time STAR --runMode genomeGenerate --genomeDir $genomeDir2 --genomeFastaFiles $genome --sjdbFileChrStartEnd $runDir_pass1/SJ.out.tab --sjdbOverhang 75 --runThreadN 8
#real    0m50.026s
#user    1m5.731s
#sys     0m3.748s

# Final alignments:
cd $runDir_pass2
echo "STAR --genomeDir $genomeDir2 --readFilesIn $read1 $read2 --runThreadN 8  --readFilesCommand zcat"

time STAR --genomeDir $genomeDir2 --readFilesIn $read1 $read2 --runThreadN 8  --readFilesCommand zcat --outFilterType BySJout --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999
#real    1m2.729s
#user    6m32.568s
#sys     0m2.217s

#3. Add groups. Note, no mark duplicates! 

inputFile=$runDir_pass2/Aligned.out.sam
output_groupFile=$runDir_pass2/groups_added_sorted.bam
echo "java -jar -XX:ParallelGCThreads=8 $picard AddOrReplaceReadGroups I=$inputFile O=$output_groupFile SO=coordinate RGID=groupid RGLB=library RGPL=illumina RGPU=machine RGSM=GM12878"

time java -jar -XX:ParallelGCThreads=8 $picard AddOrReplaceReadGroups I=$inputFile O=$output_groupFile SO=coordinate RGID=$fileId RGLB=library RGPL=illumina RGPU=machine RGSM=GM12878

#real    1m1.741s
#user    1m9.555s
#sys     0m4.562s

echo "samtools index $output_groupFile"
time samtools index $output_groupFile
#real    0m3.468s
#user    0m2.986s
#sys     0m0.056s

#4.  Split'N'Trim and reassign mapping qualities
splitFile=$runDir_pass2/split.bam
# this script requires an index file *.fai for genome and *.dict file for genome

#to split and reassign mapping qualities
echo "java -jar -XX:ParallelGCThreads=8  $GATK -T SplitNCigarReads -R $genome -I $output_groupFile -o $splitFile -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS"

time java -jar -XX:ParallelGCThreads=8  $GATK -T SplitNCigarReads -R $genome -I $output_groupFile -o $splitFile -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS --fix_misencoded_quality_scores
#real    2m16.599s
#user    2m59.034s
#sys     0m10.206s

echo "samtools index $splitFile"
time samtools index $splitFile
#real    0m3.372s
#user    0m3.284s
#sys     0m0.058s

#5. Other GATK steps
#  Indel Realignment (optional, effect is marginal + requires list of known SNPs and indels - i.e. from 1000 genomes or dbSNPs)
#  Base Recalibration (optional, effect is marginal on a good quality data, requires information about known sites similar to the tool above)

baseRecFile1=$runDir_pass2/final.rnaseq.grp
baseRecFile2=$runDir_pass2/final.bam
DBSNP=/no_backup/rg/avlasova/NGS_2017/snv/NA12878.chr22.vcf.gz.filtered.recode.vcf

time java -jar $GATK -T BaseRecalibrator -nct 8 --default_platform illumina -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -knownSites $DBSNP -cov ContextCovariate -R $genome -I $splitFile --downsampling_type NONE -o $baseRecFile1
#real    1m20.128s
#user    7m35.777s
#sys     0m11.232s

time java -jar $GATK -T PrintReads -R $genome -I $splitFile -BQSR $baseRecFile1 -o $baseRecFile2
#real    5m33.240s
#user    6m8.974s
#sys     0m15.185s

#6.select only unique alignments, no multimappings.
finalBam=$runDir_pass2/final.sorted.bam
time (samtools view -H $baseRecFile2; (samtools view $baseRecFile2|awk '$14~/NH\:i\:1/'))|samtools view -Sb -  > $finalBam
#real    1m23.224s
#user    1m46.059s
#sys     0m5.591s

time samtools index $finalBam
#real    0m6.097s
#user    0m5.933s
#sys     0m0.125s

#7. Variant calling -- this part will be in separate script, I will call variants on the merged bams
#In total this part will take ~37,5 min. 

#initVCF=output.vcf
#echo "java -jar -XX:ParallelGCThreads=8 $GATK -T HaplotypeCaller -R $genome -I $finalBam -dontUseSoftClippedBases -stand_call_conf 20.0 -stand_emit_conf 20.0 -o $initVCF"

#time java -jar -XX:ParallelGCThreads=8 $GATK -T HaplotypeCaller -R $genome -I $finalBam -dontUseSoftClippedBases -stand_call_conf 20.0 -stand_emit_conf 20.0 -o $initVCF
#7.  Variant filtering
#finalVCF=final.vcf
#echo "java -jar -XX:ParallelGCThreads=4 /software/rg/el6.3/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T VariantFiltration -R $genome -V $initVCF -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o $finalVCF"

#As in DNA-seq, we recommend filtering based on Fisher Strand values (FS > 30.0) and Qual By Depth values (QD < 2.0).
#We recommend that you filter clusters of at least 3 SNPs that are within a window of 35 bases between them by adding -window 35 -cluster 3 to your command. This filter recommendation is specific for RNA-seq data. 

#time java -jar -XX:ParallelGCThreads=4 /software/rg/el6.3/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T VariantFiltration -R $genome -V $initVCF -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o $finalVCF

