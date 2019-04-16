/*
 * Copyright (c) 2017-2019, Centre for Genomic Regulation (CRG).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. 
 * 
 * This Source Code Form is "Incompatible With Secondary Licenses", as
 * defined by the Mozilla Public License, v. 2.0.
 */
 
  
/* 
 * 'CalliNGS-NF' - A Nextflow pipeline for variant calling with NGS data
 * 
 * This pipeline that reproduces steps from the GATK best practics of SNP 
 * calling with RNAseq data procedure:
 * https://software.broadinstitute.org/gatk/guide/article?id=3891
 * 
 * Anna Vlasova 
 * Emilio Palumbo 
 * Paolo Di Tommaso
 * Evan Floden 
 * Luca Cozzuto
 */


/*
 * Define the default parameters
 */ 

params.genome     = "$baseDir/data/genome.fa"
params.variants   = "$baseDir/data/known_variants.vcf.gz"
params.blacklist  = "$baseDir/data/blacklist.bed" 
params.reads      = "$baseDir/data/reads/*_{1,2}.fq.gz"
params.results    = "results"

log.info """\
C A L L I N G S  -  N F    v 1.1 
================================
genome   : $params.genome
reads    : $params.reads
variants : $params.variants
blacklist: $params.blacklist
results  : $params.results
"""

/*
 *  Parse the input parameters
 */

genome_file     = file(params.genome)
variants_file   = file(params.variants)
blacklist_file  = file(params.blacklist)
reads_ch        = Channel.fromFilePairs(params.reads)


/**********
 * PART 1: Data preparation
 *
 * Process 1A: Create a FASTA genome index (.fai) with samtools for GATK
 */

process '1A_prepare_genome_samtools' { 
  tag "$genome.baseName"
 
  input: 
      file genome from genome_file 
 
  output: 
      file "${genome}.fai" into genome_index_ch  
  
  script:
  """
  samtools faidx ${genome}
  """
}


/*
 * Process 1B: Create a FASTA genome sequence dictionary with Picard for GATK
 */

process '1B_prepare_genome_picard' {
  tag "$genome.baseName"
  label 'mem_xlarge'

  input:
      file genome from genome_file
  output:
      file "${genome.baseName}.dict" into genome_dict_ch

  script:
  """
  PICARD=`which picard.jar`
  java -jar \$PICARD CreateSequenceDictionary R= $genome O= ${genome.baseName}.dict
  """
}


/*
 * Process 1C: Create STAR genome index file.
 */

process '1C_prepare_star_genome_index' {
  tag "$genome.baseName"

  input:
      file genome from genome_file
  output:
      file "genome_dir" into genome_dir_ch

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

process '1D_prepare_vcf_file' {
  tag "$variantsFile.baseName"

  input: 
      file variantsFile from variants_file
      file blacklisted from blacklist_file

  output:
      set file("${variantsFile.baseName}.filtered.recode.vcf.gz"), file("${variantsFile.baseName}.filtered.recode.vcf.gz.tbi") into prepared_vcf_ch
  
  script:  
  """
  vcftools --gzvcf $variantsFile -c \
           --exclude-bed ${blacklisted} \
           --recode | bgzip -c \
           > ${variantsFile.baseName}.filtered.recode.vcf.gz

  tabix ${variantsFile.baseName}.filtered.recode.vcf.gz
  """
}

/*
 *  END OF PART 1
 *********/



/**********
 * PART 2: STAR RNA-Seq Mapping
 *
 * Process 2: Align RNA-Seq reads to the genome with STAR
 */

process '2_rnaseq_mapping_star' {
  tag "$replicateId"

  input: 
      file genome from genome_file 
      file genomeDir from genome_dir_ch
      set replicateId, file(reads) from reads_ch 

  output: 
      set replicateId, file('Aligned.sortedByCoord.uniq.bam'), file('Aligned.sortedByCoord.uniq.bam.bai') into aligned_bam_ch

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

  # Select only unique alignments, no multimaps
  (samtools view -H Aligned.sortedByCoord.out.bam; samtools view Aligned.sortedByCoord.out.bam| grep -w 'NH:i:1') \
  |samtools view -Sb -  > Aligned.sortedByCoord.uniq.bam
  
  # Index the BAM file
  samtools index Aligned.sortedByCoord.uniq.bam
  """
}

/*
 *  END OF PART 2
 ******/


/**********
 * PART 3: GATK Prepare Mapped Reads
 *
 * Process 3: Split reads that contain Ns in their CIGAR string.
 *            Creates k+1 new reads (where k is the number of N cigar elements) 
 *            that correspond to the segments of the original read beside/between 
 *            the splicing events represented by the Ns in the original CIGAR.
 */

process '3_rnaseq_gatk_splitNcigar' {
  tag "$replicateId"
  label 'mem_large'
  
  input: 
      file genome from genome_file 
      file index from genome_index_ch
      file genome_dict from genome_dict_ch
      set replicateId, file(bam), file(index) from aligned_bam_ch

  output:
      set replicateId, file('split.bam'), file('split.bai') into splitted_bam_ch
  
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
 *  END OF PART 3
 ******/


/***********
 * PART 4: GATK Base Quality Score Recalibration Workflow
 *
 * Process 4: Base recalibrate to detect systematic errors in base quality scores, 
 *            select unique alignments and index
 *             
 */

process '4_rnaseq_gatk_recalibrate' {
  tag "$replicateId"
  label 'mem_large'    

  input: 
      file genome from genome_file 
      file index from genome_index_ch
      file dict from genome_dict_ch
      set replicateId, file(bam), file(index) from splitted_bam_ch
      set file(variants_file), file(variants_file_index) from prepared_vcf_ch

  output:
      set sampleId, file("${replicateId}.final.uniq.bam"), file("${replicateId}.final.uniq.bam.bai") into (final_output_ch, bam_for_ASE_ch)
  
  script: 
  sampleId = replicateId.replaceAll(/[12]$/,'')
  """
  # Indel Realignment and Base Recalibration
  gatk BaseRecalibrator \
          -R ${genome} \
          -I ${bam} \
          --known-sites ${variants_file} \
          -O final.rnaseq.grp 

  gatk ApplyBQSR \
          -R ${genome} -I ${bam} \
          --bqsr-recal-file final.rnaseq.grp \
          -O ${replicateId}.final.uniq.bam

  # Index BAM files
  samtools index ${replicateId}.final.uniq.bam
  """
}

/*
 *  END OF PART 4
 ******/



/***********
 * PART 5: GATK Variant Calling
 *
 * Process 5: Call variants with GATK HaplotypeCaller.
 *            Calls SNPs and indels simultaneously via local de-novo assembly of 
 *            haplotypes in an active region.
 *            Filter called variants with GATK VariantFiltration.    
 */


process '5_rnaseq_call_variants' {
  tag "$sampleId"
  label 'mem_xlarge'

  input:
      file genome from genome_file
      file index from genome_index_ch
      file dict from genome_dict_ch
      set sampleId, file(bam), file(bai) from final_output_ch.groupTuple()
 
  output: 
      set sampleId, file('final.vcf') into vcf_files

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
          -R $genome -V output.gatk.vcf.gz \
          --cluster-window-size 35 --cluster-size 3 \
          --filter-name FS --filter-expression \"FS > 30.0\" \
          --filter-name QD --filter-expression \"QD < 2.0\" \
          -O final.vcf
  """
}

/*
 *  END OF PART 5
 ******/


/***********
 * PART 6: Post-process variants file and prepare for Allele-Specific Expression and RNA Editing Analysis
 *
 * Process 6A: Post-process the VCF result  
 */

process '6A_post_process_vcf' {
  tag "$sampleId"
  publishDir "$params.results/$sampleId" 

  input:
      set sampleId, file('final.vcf') from vcf_files
      set file('filtered.recode.vcf.gz'), file('filtered.recode.vcf.gz.tbi') from prepared_vcf_ch 
  output: 
      set sampleId, file('final.vcf'), file('commonSNPs.diff.sites_in_files') into vcf_and_snps_ch
  
  script:
  '''
  grep -v '#' final.vcf | awk '$7~/PASS/' |perl -ne 'chomp($_); ($dp)=$_=~/DP\\=(\\d+)\\;/; if($dp>=8){print $_."\\n"};' > result.DP8.vcf
  
  vcftools --vcf result.DP8.vcf --gzdiff filtered.recode.vcf.gz  --diff-site --out commonSNPs
  '''
}

/* 
 * Process 6B: Prepare variants file for allele specific expression (ASE) analysis
 */

process '6B_prepare_vcf_for_ase' {
  tag "$sampleId"
  publishDir "$params.results/$sampleId" 

  input: 
      set sampleId, file('final.vcf'), file('commonSNPs.diff.sites_in_files') from vcf_and_snps_ch
  output: 
      set sampleId, file('known_snps.vcf.gz'), file('known_snps.vcf.gz.tbi') into vcf_for_ASE
      file('AF.histogram.pdf') into gghist_pdfs

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
  # Known SNPs have to be zipped and indexed for being used
  bgzip -c known_snps.vcf  > known_snps.vcf.gz
  tabix -p vcf known_snps.vcf.gz
  '''
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
 * it creates a channel named `grouped_vcf_bam_bai_ch` emitting the following tuples: 
 *  
 *   ( sample_id, file_vcf, List[file_bam], List[file_bai] )
 */

bam_for_ASE_ch
    .groupTuple()
    .phase(vcf_for_ASE)
    .map{ left, right -> 
      def sampleId = left[0]
      def bam = left[1]
      def bai = left[2]
      def vcf = right[1]
      def tbi = right[2]
      tuple(sampleId, vcf, tbi, bam, bai)  
    }
    .set { grouped_vcf_bam_bai_ch }


/* 
 * Process 6C: Allele-Specific Expression analysis with GATK ASEReadCounter.
 *             Calculates allele counts at a set of positions after applying 
 *             filters that are tuned for enabling allele-specific expression 
 *             (ASE) analysis
 */

process '6C_ASE_knownSNPs' {
  tag "$sampleId"
  publishDir "$params.results/$sampleId" 
  label 'mem_large'  

  input:
      file genome from genome_file 
      file index from genome_index_ch
      file dict from genome_dict_ch
      set sampleId, file(vcf), file(tbi), file(bam), file(bai) from grouped_vcf_bam_bai_ch
  
  output:
      file "ASE.tsv"
  
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
 *  END OF PART 6
 ******/

