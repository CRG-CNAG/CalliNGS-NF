
/*
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
  gatk CreateSequenceDictionary -R $genome -O ${genome.baseName}.dict
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
    path denylisted

  output:
    tuple \
      path("${variantsFile.baseName}.filtered.recode.vcf.gz"), \
      path("${variantsFile.baseName}.filtered.recode.vcf.gz.tbi")
  
  script:  
  """
  vcftools --gzvcf $variantsFile -c \
           --exclude-bed ${denylisted} \
           --recode | bgzip -c \
           > ${variantsFile.baseName}.filtered.recode.vcf.gz

  tabix ${variantsFile.baseName}.filtered.recode.vcf.gz
  """
}

/*
 * Process 2: Align RNA-Seq reads to the genome with STAR
 */

process RNASEQ_MAPPING_STAR {
  tag "$replicateId"

  input: 
    path genome
    path genomeDir
    tuple val(replicateId), path(reads) 

  output: 
    tuple \
      val(replicateId), \
      path('Aligned.sortedByCoord.uniq.bam'), \
      path('Aligned.sortedByCoord.uniq.bam.bai')

  script:
  """
  # ngs-nf-dev Align reads to genome
  STAR --genomeDir $genomeDir \
       --readFilesIn $reads \
       --runThreadN $task.cpus \
       --readFilesCommand zcat \
       --outFilterType BySJout \
       --alignSJoverhangMin 8 \
       --alignSJDBoverhangMin 1 \
       --outFilterMismatchNmax 999

  # Run 2-pass mapping (improve alignmets using table of splice junctions and create a new index)  
  STAR --genomeDir $genomeDir \
       --readFilesIn $reads \
       --runThreadN $task.cpus \
       --readFilesCommand zcat \
       --outFilterType BySJout \
       --alignSJoverhangMin 8 \
       --alignSJDBoverhangMin 1 \
       --outFilterMismatchNmax 999 \
       --sjdbFileChrStartEnd SJ.out.tab \
       --outSAMtype BAM SortedByCoordinate \
       --outSAMattrRGline ID:$replicateId LB:library PL:illumina PU:machine SM:GM12878

  # Select only unique alignments, no multimaps
  (samtools view -H Aligned.sortedByCoord.out.bam; samtools view Aligned.sortedByCoord.out.bam| grep -w 'NH:i:1') \
  |samtools view -Sb - > Aligned.sortedByCoord.uniq.bam
  
  # Index the BAM file
  samtools index Aligned.sortedByCoord.uniq.bam
  """
}


/*
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
  gatk SplitNCigarReads \
            -R $genome \
            -I $bam \
            --refactor-cigar-string \
            -O split.bam
  """
}

/*
 * Process 4: Base recalibrate to detect systematic errors in base quality scores, 
 *            select unique alignments and index
 *             
 */

process RNASEQ_GATK_RECALIBRATE {
  tag "$replicateId"
  label "mem_large"

  input: 
    path genome
    path index
    path dict
    tuple val(replicateId), path(bam), path(index)
    tuple path(variants_file), path(variants_file_index)

  output:
    tuple \
      val(sampleId), \
      path("${replicateId}.final.uniq.bam"), \
      path("${replicateId}.final.uniq.bam.bai")
  
  script: 
  sampleId = replicateId.replaceAll(/[12]$/,'')
  """
  # Indel Realignment and Base Recalibration
  gatk BaseRecalibrator \
          -R $genome \
          -I $bam \
          --known-sites $variants_file \
          -O final.rnaseq.grp 

  gatk ApplyBQSR \
          -R $genome -I $bam \
          --bqsr-recal-file final.rnaseq.grp \
          -O ${replicateId}.final.uniq.bam

  # Index BAM files
  samtools index ${replicateId}.final.uniq.bam
  """
}

/*
 * Process 5: Call variants with GATK HaplotypeCaller.
 *            Calls SNPs and indels simultaneously via local de-novo assembly of 
 *            haplotypes in an active region.
 *            Filter called variants with GATK VariantFiltration.    
 */

process RNASEQ_CALL_VARIANTS {
  tag "$sampleId"
  label "mem_xlarge"

  input:
    path genome
    path index
    path dict
    tuple val(sampleId), path(bam), path(bai)
 
  output: 
    tuple val(sampleId), path('final.vcf')

  script:
  def bam_params = bam.collect{ "-I $it" }.join(' ')
  """
  # fix absolute path in dict file
  sed -i 's@UR:file:.*${genome}@UR:file:${genome}@g' $dict
  
  # Variant calling
  gatk HaplotypeCaller \
          --native-pair-hmm-threads ${task.cpus} \
          --reference ${genome} \
          --output output.gatk.vcf.gz \
          ${bam_params} \
          --standard-min-confidence-threshold-for-calling 20.0 \
          --dont-use-soft-clipped-bases 

  # Variant filtering
  gatk VariantFiltration \
          -R ${genome} -V output.gatk.vcf.gz \
          --cluster-window-size 35 --cluster-size 3 \
          --filter-name FS --filter-expression \"FS > 30.0\" \
          --filter-name QD --filter-expression \"QD < 2.0\" \
          -O final.vcf
  """
}

/*
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
    tuple val(sampleId), path('known_snps.vcf.gz'), path('known_snps.vcf.gz.tbi')
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

  # gghist.R -i AF.4R -o AF.histogram.pdf
  # Known SNPs have to be zipped and indexed for being used
  bgzip -c known_snps.vcf  > known_snps.vcf.gz
  tabix -p vcf known_snps.vcf.gz
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
  label "mem_large"

  input:
    path genome
    path index
    path dict
    tuple val(sampleId), path(vcf), path(tbi), path(bam), path(bai)
  
  output:
    path "ASE.tsv"
  
  script:
  def bam_params = bam.collect{ "-I $it" }.join(' ')
  """
  gatk ASEReadCounter \
          -R ${genome} \
          -O ASE.tsv \
          ${bam_params} \
          -V ${vcf}
  """
}

/* 
 * Group data for allele-specific expression.
 * 
 * The `bam_for_ASE_ch` emites tuples having the following structure, holding the final BAM/BAI files:
 *  
 *   ( sample_id, file_bam, file_bai )
 * 
 * The `vcf_for_ASE` channel emits tuples having the following structure, holding the VCF file:
 *  
 *   ( sample_id, output.vcf ) 
 * 
 * The BAMs are grouped together and merged with VCFs having the same sample id. Finally 
 * it returns a channel emitting the following tuples: 
 *  
 *   ( sample_id, vcf_file, tbi_file, List[file_bam], List[file_bai] )
 */

def group_per_sample(bam_ch, vcf_ch) {
  bam_ch
    .groupTuple()
    .join(vcf_ch)
    .map{ left, right -> 
      def sampleId = left[0]
      def bam = left[1]
      def bai = left[2]
      def vcf = right[1]
      def tbi = right[2]
      tuple(sampleId, vcf, tbi, bam, bai)
    }
}
