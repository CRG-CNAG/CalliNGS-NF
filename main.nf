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
nextflow.preview.dsl = 2
include 'modules.nf' params(gatk: params.gatk, results: params.results)

log.info """\
C A L L I N G S  -  N F    v 2.0 
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

workflow {

      PREPARE_GENOME_SAMTOOLS(genome_file)

      PREPARE_GENOME_PICARD(genome_file)

      PREPARE_STAR_GENOME_INDEX(genome_file)

      PREPARE_VCF_FILE(variants_file, blacklist_file)

      RNASEQ_MAPPING_STAR( 
            genome_file, 
            PREPARE_STAR_GENOME_INDEX.out, 
            reads_ch)

      RNASEQ_GATK_SPLITNCIGAR(
            genome_file, 
            PREPARE_GENOME_SAMTOOLS.out, 
            PREPARE_GENOME_PICARD.out, 
            RNASEQ_MAPPING_STAR.out)

      RNASEQ_GATK_RECALIBRATE(
                  genome_file, PREPARE_GENOME_SAMTOOLS.out, 
                  PREPARE_GENOME_PICARD.out, 
                  RNASEQ_GATK_SPLITNCIGAR.out, 
                  PREPARE_VCF_FILE.out)

      RNASEQ_CALL_VARIANTS( 
            genome_file, 
            PREPARE_GENOME_SAMTOOLS.out, 
            PREPARE_GENOME_PICARD.out, 
            RNASEQ_GATK_RECALIBRATE.out.groupTuple())

      POST_PROCESS_VCF( 
            RNASEQ_CALL_VARIANTS.out, 
            PREPARE_VCF_FILE.out )

      PREPARE_VCF_FOR_ASE( POST_PROCESS_VCF.out )

      ASE_KNOWNSNPS(
            genome_file, 
            PREPARE_GENOME_SAMTOOLS.out, 
            PREPARE_GENOME_PICARD.out, 
            group_per_sample(   
                  RNASEQ_GATK_RECALIBRATE.out, 
                  PREPARE_VCF_FOR_ASE.out[0]) )

}