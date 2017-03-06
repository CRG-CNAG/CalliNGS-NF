/*
 * Copyright (c) 2017, Centre for Genomic Regulation (CRG).
 *
 *   This file is part of 'NGS17WS-NF'.
 *
 *   NGS17WS-NF is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   NGS17WS-NF is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with NGS17WS-NF.  If not, see <http://www.gnu.org/licenses/>.
 */
 
  
/* 
 * Calling variants in RNAseq
 * 
 * This pipeline reproduces steps from the GATK best practics of SNP calling with RNAseq data
 * https://software.broadinstitute.org/gatk/guide/article?id=3891
 * 
 * Anna Vlasova 
 * Emilio Palumbo 
 * Paolo Di Tommaso
 * Evan Floden 
 */
 
params.genome = "$baseDir/data/genome.fa"
params.variant = "$baseDir/data/NA12878.chr22.vcf.gz"
params.blacklist = "$baseDir/data/ENCFF001TDO.sorted.bed" 
params.reads = "$baseDir/data/*_{1,2}.fastq.gz"
params.gatk = 'GenomeAnalysisTK'

GATK = params.gatk
genome_file = file(params.genome)
variant_file = file(params.variant)
blacklist_file = file(params.blacklist)

reads_ch = Channel.fromFilePairs(params.reads)

/* 
 * create STAR genome index file.
 */
process '1a_prepape_genome_star' {
  input: file(genome) from genome_file 
  output: file(genome_dir) into genome_index_ch
  """
  mkdir genome_dir
  STAR --runMode genomeGenerate --genomeDir genome_dir --genomeFastaFiles $genome  --runThreadN ${task.cpus}
  """
}

/* 
 * For the Split and Trim GATK tool I need index file *.fai for genome and *.dict file for genome
 * to create dictionary and index files for genome
 */
process '1a_prepape_genome_samtools' {
  input: file genome from genome_file 
  output: file "${genome}.fai" into genome_index, genome_index1  
  
  """
  samtools faidx $genome
  """
}

process '1a_prepare_genome_picard' {
  input: file genome from genome_file 
  output: file "${genome.baseName}.dict" into genome_dict, genome_dict1
  
  """
  PICARD=`which picard.jar`
  java -jar -XX:ParallelGCThreads=8 \$PICARD CreateSequenceDictionary R= $genome O= ${genome.baseName}.dict
  """
}

process '1b_prepare_vcf_file' {
  input: 
  file(variantFile) from variant_file
  file(blacklisted) from blacklist_file

  output:
  set file("${variantFile.baseName}.filtered.recode.vcf.gz"), file("${variantFile.baseName}.filtered.recode.vcf.gz.tbi") into prepared_vcf
    
  """
  vcftools --gzvcf $variantFile -c --exclude-bed $blacklisted --recode | bgzip -c > ${variantFile.baseName}.filtered.recode.vcf.gz
  tabix ${variantFile.baseName}.filtered.recode.vcf.gz
  """
}

/* 
 * Align reads to genome
 */
process '2_rnaseq_star' {
  input: 
  file genome from genome_file 
  file genomeDir from genome_index_ch.first()
  set pairId, file(reads) from reads_ch 

  output: 
  set file('Aligned.sortedByCoord.out.bam'), file('Aligned.sortedByCoord.out.bam.bai') into output_groupFile

  """
  # Align reads to genome
  STAR --genomeDir $genomeDir --readFilesIn $reads --runThreadN ${task.cpus} --readFilesCommand zcat --outFilterType BySJout --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999
    
  # 2d pass STAR - improve alignmnets using table of splice junctions. create a new index  
  mkdir genomeDir  
  STAR --runMode genomeGenerate --genomeDir genomeDir --genomeFastaFiles $genome --sjdbFileChrStartEnd SJ.out.tab --sjdbOverhang 75 --runThreadN ${task.cpus}  
    
  # Final alignments  
  STAR --genomeDir genomeDir --readFilesIn $reads --runThreadN ${task.cpus} --readFilesCommand zcat --outFilterType BySJout --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outSAMtype BAM SortedByCoordinate --outSAMattrRGline ID:$pairId LB:library PL:illumina PU:machine SM:GM12878

  # Index the BAM file
  samtools index Aligned.sortedByCoord.out.bam
  """
}

process '2_rnaseq_gatk_split_n_cigar' {
  container 'biodckrdev/gatk'
  
  input: 
  file genome from genome_file 
  file index from genome_index
  set file(bam), file(index) from output_groupFile
  file genome_dict

  output:
  set file('split.bam'), file('split.bam.bai') into output_split
  
  """
  # Split'N'Trim and reassign mapping qualities
  
  $GATK -T SplitNCigarReads -R $genome -I $bam -o split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS --fix_misencoded_quality_scores
  """
}

process '2_rnaseq_gatk_recalibrate' {
  container 'biodckrdev/gatk'
  
  input: 
  file genome from genome_file 
  file index from genome_index1
  set file(bam), file(index) from output_split
  file genome_dict1
  set file(variant_file), file(variant_file_index) from prepared_vcf

  output:
  set file('final.bam'), file('final.bam.bai') into output_final
  
  """
  #  Indel Realignment and Base Recalibration

  $GATK -T BaseRecalibrator -nct 8 --default_platform illumina -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -knownSites $variant_file -cov ContextCovariate -R $genome -I $bam --downsampling_type NONE -o final.rnaseq.grp

  $GATK -T PrintReads -R $genome -I $bam -BQSR final.rnaseq.grp -o final.bam
  """
}

