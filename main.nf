/*
 * Copyright (c) 2020, Seqera Labs.
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
 */

/* 
 * Enable DSL 2 syntax
 */
nextflow.enable.dsl = 2

/*
 * Define the default parameters
 */ 
params.genome     = "$baseDir/data/genome.fa"
params.variants   = "$baseDir/data/known_variants.vcf.gz"
params.blacklist  = "$baseDir/data/blacklist.bed" 
params.reads      = "$baseDir/data/reads/rep1_{1,2}.fq.gz"
params.results    = "results"
params.gatk       = '/usr/local/bin/GenomeAnalysisTK.jar'

log.info """\
C A L L I N G S  -  N F    v 2.1 
================================
genome   : $params.genome
reads    : $params.reads
variants : $params.variants
blacklist: $params.blacklist
results  : $params.results
gatk     : $params.gatk
"""

/* 
 * Import modules 
 */
include { 
  PREPARE_GENOME_SAMTOOLS;
  PREPARE_GENOME_PICARD; 
  PREPARE_STAR_GENOME_INDEX;
  PREPARE_VCF_FILE;
  RNASEQ_MAPPING_STAR;
  RNASEQ_GATK_SPLITNCIGAR; 
  RNASEQ_GATK_RECALIBRATE;
  RNASEQ_CALL_VARIANTS;
  POST_PROCESS_VCF;
  PREPARE_VCF_FOR_ASE;
  ASE_KNOWNSNPS;
  group_per_sample } from './modules.nf' 

/* 
 * main pipeline logic
 */
workflow {
      reads_ch = Channel.fromFilePairs(params.reads)

      // PART 1: Data preparation
      PREPARE_GENOME_SAMTOOLS(params.genome)
      PREPARE_GENOME_PICARD(params.genome)
      PREPARE_STAR_GENOME_INDEX(params.genome)
      PREPARE_VCF_FILE(params.variants, params.blacklist)

      // PART 2: STAR RNA-Seq Mapping
      RNASEQ_MAPPING_STAR( 
            params.genome, 
            PREPARE_STAR_GENOME_INDEX.out, 
            reads_ch )

      // PART 3: GATK Prepare Mapped Reads
      RNASEQ_GATK_SPLITNCIGAR(
            params.genome, 
            PREPARE_GENOME_SAMTOOLS.out, 
            PREPARE_GENOME_PICARD.out, 
            RNASEQ_MAPPING_STAR.out )

      // PART 4: GATK Base Quality Score Recalibration Workflow
      RNASEQ_GATK_RECALIBRATE(
                  params.genome, 
                  PREPARE_GENOME_SAMTOOLS.out, 
                  PREPARE_GENOME_PICARD.out, 
                  RNASEQ_GATK_SPLITNCIGAR.out, 
                  PREPARE_VCF_FILE.out)

      // PART 5: GATK Variant Calling
      RNASEQ_CALL_VARIANTS( 
            params.genome, 
            PREPARE_GENOME_SAMTOOLS.out, 
            PREPARE_GENOME_PICARD.out, 
            RNASEQ_GATK_RECALIBRATE.out.groupTuple() )

      // PART 6: Post-process variants file and prepare for 
      // Allele-Specific Expression and RNA Editing Analysis
      POST_PROCESS_VCF( 
            RNASEQ_CALL_VARIANTS.out, 
            PREPARE_VCF_FILE.out )

      PREPARE_VCF_FOR_ASE( POST_PROCESS_VCF.out )

      ASE_KNOWNSNPS(
            params.genome, 
            PREPARE_GENOME_SAMTOOLS.out, 
            PREPARE_GENOME_PICARD.out, 
            group_per_sample(
                  RNASEQ_GATK_RECALIBRATE.out, 
                  PREPARE_VCF_FOR_ASE.out[0]) )
}
