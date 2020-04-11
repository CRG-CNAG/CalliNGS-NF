/*
 * Copyright (c) 2017-2018, Centre for Genomic Regulation (CRG).
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
 */


/*
 * Define the default parameters
 */ 

params.genome     = "$baseDir/data/genome.fa"
params.variants   = "$baseDir/data/known_variants.vcf.gz"
params.blacklist  = "$baseDir/data/blacklist.bed" 
params.reads      = "$baseDir/data/reads/rep1_{1,2}.fq.gz"
params.results    = "results"
params.gatk       = '/usr/local/bin/GenomeAnalysisTK.jar'
params.gatk_launch = "java -jar $params.gatk"
params.tmp        = "/home/fyan0011/mx82_scratch/tmp/" 

log.info """\
C A L L I N G S  -  N F    v 1.0 
================================
genome   : $params.genome
reads    : $params.reads
variants : $params.variants
blacklist: $params.blacklist
results  : $params.results
gatk     : $params.gatk
tmp      : $params.tmp
"""

/*
 *  Parse the input parameters
 */

GATK            = params.gatk_launch
genome_file     = file(params.genome)
variants_file   = file(params.variants)
blacklist_file  = file(params.blacklist)
reads_ch        = Channel.fromFilePairs(params.reads)
tmp_path        = file(params.tmp)

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
      file "chrfilter.bed" into chr_filter  

  script:
  """
  samtools faidx ${genome}
  ## need to escape for dollar sign
  awk 'BEGIN{OFS="\t"} /^([0-9]*|MT|X|Y)\t/ {print \$1,0,\$2}' "${genome}.fai" > chrfilter.bed
  """
}

/*
 * Process 1B: Create a FASTA genome sequence dictionary with Picard for GATK
 */

process '1B_prepare_genome_picard' {
  tag "$genome.baseName"
  label 'star'

  input:
      file genome from genome_file
  output:
      file "${genome.baseName}.dict" into genome_dict_ch

  script:
  """
  picard CreateSequenceDictionary R= $genome O= ${genome.baseName}.dict
  """
}


/*
 * Process 1C: Create STAR genome index file.
 */

process '1C_prepare_star_genome_index' {
  tag "$genome.baseName"
  label 'star'

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
  label 'star'

  input:
      path tmp from tmp_path 
      file chr_filter from chr_filter
      file genome from genome_file 
      file genomeDir from genome_dir_ch
      set replicateId, file(reads) from reads_ch 

  output: 
      set replicateId, file('Aligned.sorted.mdups.bam'), file('Aligned.sorted.mdups.bai') into aligned_bam_ch

  script:    
  """
  # use new 2pass mode of STAR
 
  STAR --twopassMode Basic \
       --genomeDir $genomeDir \
       --readFilesIn $reads \
       --runThreadN ${task.cpus} \
       --readFilesCommand zcat \
       --outFilterType BySJout \
       --alignSJoverhangMin 8 \
       --alignSJDBoverhangMin 1 \
       --outFilterMismatchNmax 999 \
       --outSAMtype BAM SortedByCoordinate \
       --outSAMattrRGline ID:$replicateId LB:library PL:illumina PU:machine SM:$replicateId
  
  # Index the BAM file
  # samtools index Aligned.sortedByCoord.out.bam
  
  # Filter chromosome
  samtools view -bh -L $chr_filter -o Aligned.sort.filter.bam Aligned.sortedByCoord.out.bam 

  # Markduplicates  
  picard MarkDuplicates TMP_DIR=${tmp} \
                        VALIDATION_STRINGENCY=LENIENT \
                        CREATE_INDEX=true \
                        INPUT=Aligned.sort.filter.bam \
                        OUTPUT=Aligned.sorted.mdups.bam \
                        METRICS_FILE=dup_metrics.txt 
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
  label 'star'
  
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
  $GATK -T SplitNCigarReads \
          -R $genome -I $bam \
          -o split.bam \
          -rf ReassignOneMappingQuality \
          -RMQF 255 -RMQT 60 \
          -U ALLOW_N_CIGAR_READS
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
  label 'star'    

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
  label 'star'
  publishDir "$params.results/$sampleId", mode: 'copy'

  input:
      file variants from variants_file
      file genome from genome_file
      file index from genome_index_ch
      file dict from genome_dict_ch
      set sampleId, file(bam), file(bai) from final_output_ch.groupTuple()
  
  output:
      file 'final.vcf'

  script:
  """
  # fix absolute path in dict file
  sed -i 's@UR:file:.*${genome}@UR:file:${genome}@g' $dict
  echo "${bam.join('\n')}" > bam.list
  
  # VCF indexing
  tabix -f -p vcf $variants
 
  # Variant calling
  $GATK -T HaplotypeCaller \
          -R $genome -I bam.list \
          -dontUseSoftClippedBases \
          -stand_call_conf 20.0 \
          --dbsnp ${variants} \
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

/*
 *  END OF PART 5
 ******/

