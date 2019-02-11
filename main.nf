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


/* 
 * Import modules 
 */

require 'modules.nf', params:[gatk: params.gatk, results: params.results]

log.info """\
C A L L I N G S  -  N F    v 1.0 
================================
genome   : $params.genome
reads    : $params.reads
variants : $params.variants
blacklist: $params.blacklist
results  : $params.results
gatk     : $params.gatk
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

genome_index_ch = '1A_prepare_genome_samtools'(genome_file)

/*
 * Process 1B: Create a FASTA genome sequence dictionary with Picard for GATK
 */

genome_dict_ch = '1B_prepare_genome_picard'(genome_file)


/*
 * Process 1C: Create STAR genome index file.
 */

genome_dir_ch = '1C_prepare_star_genome_index'(genome_file)


/*
 * Process 1D: Create a file containing the filtered and recoded set of variants
 */

prepared_vcf_ch = '1D_prepare_vcf_file'(variants_file, blacklist_file)




/**********
 *
 * Process 2: Align RNA-Seq reads to the genome with STAR
 */

aligned_bam_ch = '2_rnaseq_mapping_star'( genome_file, genome_dir_ch, reads_ch)



/**********
 *
 * Process 3: Split reads that contain Ns in their CIGAR string.
 *            Creates k+1 new reads (where k is the number of N cigar elements) 
 *            that correspond to the segments of the original read beside/between 
 *            the splicing events represented by the Ns in the original CIGAR.
 */

splitted_bam_ch = '3_rnaseq_gatk_splitNcigar'(genome_file, genome_index_ch, genome_dict_ch, aligned_bam_ch)


/***********
 *
 * Process 4: Base recalibrate to detect systematic errors in base quality scores, 
 *            select unique alignments and index
 *             
 */

'4_rnaseq_gatk_recalibrate'(
            genome_file, genome_index_ch, 
            genome_dict_ch, 
            splitted_bam_ch, 
            prepared_vcf_ch)
    . into { final_output_ch; bam_for_ASE_ch }


/***********
 * PART 5: GATK Variant Calling
 *
 * Process 5: Call variants with GATK HaplotypeCaller.
 *            Calls SNPs and indels simultaneously via local de-novo assembly of 
 *            haplotypes in an active region.
 *            Filter called variants with GATK VariantFiltration.    
 */


vcf_files = '5_rnaseq_call_variants'( 
              genome_file, 
              genome_index_ch, 
              genome_dict_ch, 
              final_output_ch.groupTuple())



/***********
 * PART 6: Post-process variants file and prepare for Allele-Specific Expression and RNA Editing Analysis
 *
 * Process 6A: Post-process the VCF result  
 */

vcf_and_snps_ch = '6A_post_process_vcf'( vcf_files, prepared_vcf_ch )

/* 
 * Process 6B: Prepare variants file for allele specific expression (ASE) analysis
 */

(vcf_for_ASE, gghist_pdfs) = '6B_prepare_vcf_for_ase'( vcf_and_snps_ch )


/* 
 * Process 6C: Allele-Specific Expression analysis with GATK ASEReadCounter.
 *             Calculates allele counts at a set of positions after applying 
 *             filters that are tuned for enabling allele-specific expression 
 *             (ASE) analysis
 */

'6C_ASE_knownSNPs'(
      genome_file, 
      genome_index_ch, 
      genome_dict_ch, 
      group_per_sample(bam_for_ASE_ch, vcf_for_ASE) )

/*
 *  END OF PART 6
 ******/

