/*
 * Copyright (c) 2017, Centre for Genomic Regulation (CRG).
 *
 *   This file is part of 'CalliNGS-NF': 
 *   A Nextflow pipeline for Variant Calling with NGS data
 *
 *   CalliNGS-NF is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   CalliNGS-NF is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with CalliNGS-NF.  If not, see <http://www.gnu.org/licenses/>.
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
params.bams       = ''

log.info "C A L L I N G S  -  N F    v 1.0" 
log.info "================================"
log.info "genome   : $params.genome"
log.info "reads    : $params.reads"
log.info "bams     : $params.bams"
log.info "variants : $params.variants"
log.info "blacklist: $params.blacklist"
log.info "results  : $params.results" 
log.info "gatk     : $params.gatk"
log.info ""

/*
 *  Parse the input parameters
 */

GATK            = params.gatk
genome_file     = file(params.genome)
variants_file   = file(params.variants)
blacklist_file  = file(params.blacklist)

/*
 * Create reads channel if --map=true, else create bam channel.
 */ 
if (!params.bams) 
  Channel
        .fromFilePairs(params.reads)
        .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
        .set { reads_ch }

else 
  Channel
        .fromPath(params.bams)
        .ifEmpty { error "Cannot find any BAM files matching: ${params.bams}" }
        .map { it -> [ it.baseName, it, it+'.bai'] }
        .set { aligned_bam_ch } 

     
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
if (!params.bams) {
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
}


/*
 * Process 1D: Create a file containing the filtered and recoded set of variants
 */

process '1D_prepare_vcf_file' {
  tag "$variantsFile.baseName"

  input: 
      file variantsFile from variants_file
      file blacklisted from blacklist_file
      file genome from genome_file
      file genone_index from genome_index_ch
      file genome_dict from genome_dict_ch

  output:
      set file("${variantsFile.baseName}.finalPrep.vcf.gz"), \
          file("${variantsFile.baseName}.finalPrep.vcf.gz.tbi") \
          into prepared_vcf_ch
  
   shell:
   '''
   # Step 0
   # Change R/N/Y/W/S/B characters

   zcat !{variantsFile} | \
       awk -F '\t' '$4!~/^[NRYWSB]/&&$5!~/^[NRYWSB]/||$1~/^\\#/' \
       > !{variantsFile.baseName}.noRNY.vcf

   # Step 4
   # Remove those variants that have exactly the same position 
   grep -v '#'  !{variantsFile.baseName}.noRNY.vcf \
                |awk '{print $1" "$2}' \
                | sort|uniq -c |awk '$1>1' \
                |awk 'BEGIN{OFS="\t"}{print $2,$3,$3}' \
                > problematic.sites.bed

   bedtools intersect -header -v -a  !{variantsFile.baseName}.noRNY.vcf \
            -b problematic.sites.bed > \
            !{variantsFile.baseName}.dupsFiltered.vcf

   # Step 1 
   # Select only those variants that passed all filters
    java -jar !{GATK} \
         -T SelectVariants \
         -R !{genome} \
         -V !{variantsFile.baseName}.dupsFiltered.vcf \
         -o !{variantsFile.baseName}.allFiltered.vcf.gz --excludeFiltered

   # Step 2
   # Select only Heterozygous variants
    zcat !{variantsFile.baseName}.allFiltered.vcf.gz \
         | grep -v '0|0\\|1|1' > \
         !{variantsFile.baseName}.heterozygous.vcf
    gzip !{variantsFile.baseName}.heterozygous.vcf

   # Step 3 
   # Remove the blacklisted regions from the variants file
   bedtools intersect -header \
            -v -a !{variantsFile.baseName}.heterozygous.vcf.gz \
            -b !{blacklisted} > \
            !{variantsFile.baseName}.finalPrep.vcf

   # Step 5
   # Compress + index 
   # NB: tabix is working only with bgzipped version of the files
   bgzip !{variantsFile.baseName}.finalPrep.vcf
   tabix !{variantsFile.baseName}.finalPrep.vcf.gz
   '''
}

/*
 *  END OF PART 1
 *********/



/**********
 * PART 2: STAR RNA-Seq Mapping
 *
 * Process 2: Align RNA-Seq reads to the genome with STAR
 */


if (!params.bams) {
  process '2_rnaseq_mapping_star' {
    tag "$replicateId"

    input: 
      file genome from genome_file 
      file genomeDir from genome_dir_ch
      set replicateId, file(reads) from reads_ch 

    output: 
      set replicateId, \
          file('Aligned.sortedByCoord.out.bam'), \
          file('Aligned.sortedByCoord.out.bam.bai') \
          into aligned_bam_ch

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
  java -jar $GATK -T SplitNCigarReads \
                  -R $genome -I $bam \
                  -o split.bam \
                  -rf ReassignOneMappingQuality \
                  -RMQF 255 -RMQT 60 \
                  -U ALLOW_N_CIGAR_READS \
                  --fix_misencoded_quality_scores
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
    
  input: 
      file genome from genome_file 
      file index from genome_index_ch
      file dict from genome_dict_ch
      set replicateId, file(bam), file(index) from splitted_bam_ch
      set file(variants_file), file(variants_file_index) from prepared_vcf_ch

  output:
      set sampleId, file("${replicateId}.final.uniq.bam")  into final_output_ch
      set sampleId, file("${replicateId}.final.uniq.bam"), file("${replicateId}.final.uniq.bam.bai") into bam_for_ASE_ch
  
  script: 
  sampleId = replicateId.replaceAll(/[12]$/,'')
  """
  # Indel Realignment and Base Recalibration
  java -jar $GATK -T BaseRecalibrator \
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

  java -jar $GATK -T PrintReads \
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

  input:
      file genome from genome_file
      file index from genome_index_ch
      file dict from genome_dict_ch
      set sampleId, file(bam) from final_output_ch.groupTuple()
  
  output: 
      set sampleId, file('final.vcf') into vcf_files

  script:
  """
  echo "${bam.join('\n')}" > bam.list
  
  # Variant calling
  java -jar $GATK -T HaplotypeCaller \
                  -R $genome -I bam.list \
                  -dontUseSoftClippedBases \
                  -stand_call_conf 20.0 \
                  -o output.gatk.vcf.gz

  # Variant filtering
  java -jar $GATK -T VariantFiltration \
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


/***********
 * PART 6: Post-process variants file and prepare for Allele-Specific Expression and RNA Editing Analysis
 *
 * Process 6A: Post-process the VCF result  
 */

process '6A_post_process_vcf' {
  tag "$sampleId"
  publishDir "$params.results/$sampleId" 
  
  input:
      set sampleId, \
          file('final.vcf') from vcf_files
      set file('filtered.recode.vcf.gz'), \
          file('filtered.recode.vcf.gz.tbi') from prepared_vcf_ch 

  output: 
      set sampleId, file('final.vcf'), \
          file('commonSNPs.diff.sites_in_files') into vcf_and_snps_ch
  
  script:
  '''
  grep -v '#' final.vcf | awk '$7~/PASS/' |\
      perl -ne 'chomp($_); ($dp)=$_=~/DP\\=(\\d+)\\;/; \
      if($dp>=8){print $_."\\n"};' \
      > result.DP8.vcf
  
  vcftools --vcf result.DP8.vcf \
           --gzdiff filtered.recode.vcf.gz  \
           --diff-site --out commonSNPs
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
      set sampleId, file('known_snps.vcf') into vcf_for_ASE
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
      tuple(sampleId, vcf, bam, bai)  
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
  
  input:
      file genome from genome_file 
      file index from genome_index_ch
      file dict from genome_dict_ch
      set sampleId, file(vcf),  file(bam), file(bai) from grouped_vcf_bam_bai_ch
  
  output:
      file "ASE.tsv"
  
  script:
  """
  echo "${bam.join('\n')}" > bam.list
    
  java -jar $GATK -R ${genome} \
                  -T ASEReadCounter \
                  -o ASE.tsv \
                  -I bam.list \
                  -sites ${vcf}
  """
}

/*
 *  END OF PART 6
 ******/

