params.gatk       = '/usr/local/bin/GenomeAnalysisTK.jar'
params.gatk_launch = "java -jar $params.gatk" 
params.results    = "results"

GATK = params.gatk_launch

/**********
 * PART 1: Data preparation
 *
 * Process 1A: Create a FASTA genome index (.fai) with samtools for GATK
 */


process PREPARE_GENOME_SAMTOOLS { 
  tag "$genome.baseName"
 
  input: 
    path genome  
 
  output: 
    path "${genome}.fai"   
  
  script:
  """
  samtools faidx ${genome}
  """
}

/*
 * Process 1B: Create a FASTA genome sequence dictionary with Picard for GATK
 */

process PREPARE_GENOME_PICARD {
  tag "$genome.baseName"
  label 'mem_xlarge'

  input:
    path genome
  output:
    path "${genome.baseName}.dict"

  script:
  """
  PICARD=`which picard.jar`
  java -jar \$PICARD CreateSequenceDictionary R= $genome O= ${genome.baseName}.dict
  """
}

/*
 * Process 1C: Create STAR genome index file.
 */

process PREPARE_STAR_GENOME_INDEX {
  tag "$genome.baseName"

  input:
    path genome 
  
  output:
    path "genome_dir"

  script:
  """
  mkdir genome_dir

  STAR --runMode genomeGenerate \
       --genomeDir genome_dir \
       --genomeFastaFiles ${genome} \
       --runThreadN ${task.cpus}
  """
}

/*
 * Process 1D: Create a file containing the filtered and recoded set of variants
 */

process PREPARE_VCF_FILE {
  tag "$variantsFile.baseName"

  input: 
    path variantsFile
    path blacklisted

  output:
    tuple path("${variantsFile.baseName}.filtered.recode.vcf.gz"), path("${variantsFile.baseName}.filtered.recode.vcf.gz.tbi")
  
  script:  
  """
  vcftools --gzvcf $variantsFile -c \
           --exclude-bed ${blacklisted} \
           --recode | bgzip -c \
           > ${variantsFile.baseName}.filtered.recode.vcf.gz

  tabix ${variantsFile.baseName}.filtered.recode.vcf.gz
  """
}

/**********
 *
 * Process 2: Align RNA-Seq reads to the genome with STAR
 */

process RNASEQ_MAPPING_STAR {
  tag "$replicateId"

  input: 
    path genome  
    path genomeDir
    tuple replicateId, path(reads) 

  output: 
    tuple val(replicateId), path('Aligned.sortedByCoord.out.bam'), path('Aligned.sortedByCoord.out.bam.bai')

  script:    
  """
  # ngs-nf-dev Align reads to genome
  STAR --genomeDir $genomeDir \
       --readFilesIn $reads \
       --runThreadN ${task.cpus} \
       --readFilesCommand zcat \
       --outFilterType BySJout \
       --alignSJoverhangMin 8 \
       --alignSJDBoverhangMin 1 \
       --outFilterMismatchNmax 999
    
  # 2nd pass (improve alignmets using table of splice junctions and create a new index)  
  mkdir genomeDir  
  STAR --runMode genomeGenerate \
       --genomeDir genomeDir \
       --genomeFastaFiles $genome \
       --sjdbFileChrStartEnd SJ.out.tab \
       --sjdbOverhang 75 \
       --runThreadN ${task.cpus}  
    
  # Final read alignments  
  STAR --genomeDir genomeDir \
       --readFilesIn $reads \
       --runThreadN ${task.cpus} \
       --readFilesCommand zcat \
       --outFilterType BySJout \
       --alignSJoverhangMin 8 \
       --alignSJDBoverhangMin 1 \
       --outFilterMismatchNmax 999 \
       --outSAMtype BAM SortedByCoordinate \
       --outSAMattrRGline ID:$replicateId LB:library PL:illumina PU:machine SM:GM12878

  # Index the BAM file
  samtools index Aligned.sortedByCoord.out.bam
  """
}

/**********
 *
 * Process 3: Split reads that contain Ns in their CIGAR string.
 *            Creates k+1 new reads (where k is the number of N cigar elements) 
 *            that correspond to the segments of the original read beside/between 
 *            the splicing events represented by the Ns in the original CIGAR.
 */

process RNASEQ_GATK_SPLITNCIGAR {
  tag "$replicateId"
  label 'mem_large'
  
  input: 
    path genome
    path index
    path genome_dict
    tuple val(replicateId), path(bam), path(index)

  output:
    tuple val(replicateId), path('split.bam'), path('split.bai')
  
  script:
  """
  # SplitNCigarReads and reassign mapping qualities
  $GATK -T SplitNCigarReads \
          -R $genome -I $bam \
          -o split.bam \
          -rf ReassignOneMappingQuality \
          -RMQF 255 -RMQT 60 \
          -U ALLOW_N_CIGAR_READS \
          --fix_misencoded_quality_scores
  """
}

/***********
 *
 * Process 4: Base recalibrate to detect systematic errors in base quality scores, 
 *            select unique alignments and index
 *             
 */

process RNASEQ_GATK_RECALIBRATE {
  tag "$replicateId"
  label 'mem_large'    

  input: 
    path genome  
    path index
    path dict
    tuple val(replicateId), path(bam), path(index)
    tuple path(variants_file), path(variants_file_index)

  output:
    tuple val(sampleId), path("${replicateId}.final.uniq.bam"), path("${replicateId}.final.uniq.bam.bai")
  
  script: 
  sampleId = replicateId.replaceAll(/[12]$/,'')
  """
  # Indel Realignment and Base Recalibration
  $GATK -T BaseRecalibrator \
          --default_platform illumina \
          -cov ReadGroupCovariate \
          -cov QualityScoreCovariate \
          -cov CycleCovariate \
          -knownSites ${variants_file} \
          -cov ContextCovariate \
          -R ${genome} -I ${bam} \
          --downsampling_type NONE \
          -nct ${task.cpus} \
          -o final.rnaseq.grp 

  $GATK -T PrintReads \
          -R ${genome} -I ${bam} \
          -BQSR final.rnaseq.grp \
          -nct ${task.cpus} \
          -o final.bam

  # Select only unique alignments, no multimaps
  (samtools view -H final.bam; samtools view final.bam| grep -w 'NH:i:1') \
  |samtools view -Sb -  > ${replicateId}.final.uniq.bam

  # Index BAM files
  samtools index ${replicateId}.final.uniq.bam
  """
}

/***********
 * PART 5: GATK Variant Calling
 *
 * Process 5: Call variants with GATK HaplotypeCaller.
 *            Calls SNPs and indels simultaneously via local de-novo assembly of 
 *            haplotypes in an active region.
 *            Filter called variants with GATK VariantFiltration.    
 */

process RNASEQ_CALL_VARIANTS {
  tag "$sampleId"
  label 'mem_large'

  input:
    path genome
    path index 
    path dict 
    tuple val(sampleId), path(bam), path(bai) 
  
  output: 
    tuple val(sampleId), path('final.vcf')

  script:
  """
  # fix absolute path in dict file
  sed -i 's@UR:file:.*${genome}@UR:file:${genome}@g' $dict
  echo "${bam.join('\n')}" > bam.list
  
  # Variant calling
  $GATK -T HaplotypeCaller \
          -R $genome -I bam.list \
          -dontUseSoftClippedBases \
          -stand_call_conf 20.0 \
          -o output.gatk.vcf.gz

  # Variant filtering
  $GATK -T VariantFiltration \
          -R $genome -V output.gatk.vcf.gz \
          -window 35 -cluster 3 \
          -filterName FS -filter "FS > 30.0" \
          -filterName QD -filter "QD < 2.0" \
          -o final.vcf
  """
}

/***********
 * PART 6: Post-process variants file and prepare for Allele-Specific Expression and RNA Editing Analysis
 *
 * Process 6A: Post-process the VCF result  
 */

process POST_PROCESS_VCF {
  tag "$sampleId"
  publishDir "$params.results/$sampleId" 

  input:
    tuple val(sampleId), path('final.vcf')
    tuple path('filtered.recode.vcf.gz'), path('filtered.recode.vcf.gz.tbi')
  output: 
    tuple val(sampleId), path('final.vcf'), path('commonSNPs.diff.sites_in_files') 
  
  script:
  '''
  grep -v '#' final.vcf | awk '$7~/PASS/' |perl -ne 'chomp($_); ($dp)=$_=~/DP\\=(\\d+)\\;/; if($dp>=8){print $_."\\n"};' > result.DP8.vcf
  
  vcftools --vcf result.DP8.vcf --gzdiff filtered.recode.vcf.gz  --diff-site --out commonSNPs
  '''
}

/* 
 * Process 6B: Prepare variants file for allele specific expression (ASE) analysis
 */

process PREPARE_VCF_FOR_ASE {
  tag "$sampleId"
  publishDir "$params.results/$sampleId" 

  input: 
    tuple val(sampleId), path('final.vcf'), path('commonSNPs.diff.sites_in_files')

  output: 
    tuple val(sampleId), path('known_snps.vcf') 
    path 'AF.histogram.pdf'

  script:
  '''
  awk 'BEGIN{OFS="\t"} $4~/B/{print $1,$2,$3}' commonSNPs.diff.sites_in_files  > test.bed
    
  vcftools --vcf final.vcf --bed test.bed --recode --keep-INFO-all --stdout > known_snps.vcf

  grep -v '#'  known_snps.vcf | awk -F '\\t' '{print $10}' \
               |awk -F ':' '{print $2}'|perl -ne 'chomp($_); \
               @v=split(/\\,/,$_); if($v[0]!=0 ||$v[1] !=0)\
               {print  $v[1]/($v[1]+$v[0])."\\n"; }' |awk '$1!=1' \
               >AF.4R

  gghist.R -i AF.4R -o AF.histogram.pdf
  '''
}

/* 
 * Process 6C: Allele-Specific Expression analysis with GATK ASEReadCounter.
 *             Calculates allele counts at a set of positions after applying 
 *             filters that are tuned for enabling allele-specific expression 
 *             (ASE) analysis
 */
 
process ASE_KNOWNSNPS {
  tag "$sampleId"
  publishDir "$params.results/$sampleId" 
  label 'mem_large'  

  input:
    path genome 
    path index
    path dict 
    tuple val(sampleId), path(vcf),  path(bam), path(bai)
  
  output:
    path "ASE.tsv"
  
  script:
  """
  echo "${bam.join('\n')}" > bam.list
    
  $GATK -R ${genome} \
          -T ASEReadCounter \
          -o ASE.tsv \
          -I bam.list \
          -sites ${vcf}
  """
}


/* 
 * Group data for allele-specific expression.
 * 
 * The `bam_ch` emites tuples having the following structure, holding the final BAM/BAI files:
 *  
 *   ( sample_id, file_bam, file_bai )
 * 
 * The `vcf_ch` channel emits tuples having the following structure, holding the VCF file:
 *  
 *   ( sample_id, output.vcf ) 
 * 
 * The BAMs are grouped together and merged with VCFs having the same sample id. Finally 
 * it creates a channel named `grouped_vcf_bam_bai_ch` emitting the following tuples: 
 *  
 *   ( sample_id, file_vcf, List[file_bam], List[file_bai] )
 */

def group_per_sample(bam_ch, vcf_ch) {
  bam_ch
    .groupTuple()
    .phase(vcf_ch)
    .map{ left, right -> 
      def sampleId = left[0]
      def bam = left[1]
      def bai = left[2]
      def vcf = right[1]
      tuple(sampleId, vcf, bam, bai)  
    }
}