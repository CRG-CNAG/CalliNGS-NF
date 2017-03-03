/* 
 * this pipeline reproduces steps from the GATK best practics of SNP calling with RNAseq data
 * https://software.broadinstitute.org/gatk/guide/article?id=3891
 */
 
params.genome = "$baseDir/data/genome.fa"
params.variant = "$baseDir/data/NA12878.chr22.vcf.gz"
params.blacklist = "$baseDir/data/ENCFF001TDO.sorted.bed" 
params.reads = "$baseDir/data/*_{1,2}.fastq.gz"

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
  input: file(genome) from genome_file 
  
  """
  samtools faidx $genome
  """
}

process '1a_prepare_genome_picard' {
  input: file(genome) from genome_file 
  output: file("${genome.baseName}.dict") into genome
  
  """
  PICARD=`which picard.jar`
  java -jar -XX:ParallelGCThreads=8 \$PICARD CreateSequenceDictionary R= $genome O= ${genome.baseName}.dict
  """
}

process '1b_prepare_vcf_file' {
  input: 
    file(variantFile) from variant_file
    file(blacklisted) from blacklist_file
    
  """
  vcftools --gzvcf $variantFile --out ${variantFile}.filtered --exclude-bed $blacklisted --recode
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
  file 'Aligned.sortedByCoord.out.bam*'

  """
  # Align reads to genome
  STAR \
    --genomeDir $genomeDir \
    --readFilesIn $reads \
    --runThreadN ${task.cpus} \
    --readFilesCommand zcat \
    --outFilterType BySJout \
    --alignSJoverhangMin 8 \
    --alignSJDBoverhangMin 1 \
    --outFilterMismatchNmax 999
    
  # 2d pass STAR - improve alignmnets using table of splice junctions. create a new index  
  mkdir genomeDir  
  STAR \
    --runMode genomeGenerate \
    --genomeDir genomeDir \
    --genomeFastaFiles $genome \
    --sjdbFileChrStartEnd SJ.out.tab \
    --sjdbOverhang 75 \
    --runThreadN ${task.cpus}  
    
  # Final alignments  
  STAR \
    --genomeDir genomeDir \
    --readFilesIn $reads \
    --runThreadN ${task.cpus} \
    --readFilesCommand zcat \
    --outFilterType BySJout \
    --alignSJoverhangMin 8 \
    --alignSJDBoverhangMin 1 \
    --outFilterMismatchNmax 999 \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMattrRGline ID:$pairId LB:library PL:illumina PU:machine SM:GM12878 

  samtools index Aligned.sortedByCoord.out.bam

  """
}


