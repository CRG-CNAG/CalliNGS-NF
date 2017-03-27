# NGS 2017 Post-Conference Workshop

## CalliNGS-NF
A Nextflow pipeline for basic Variant Calling Analysis with NGS RNA-Seq data

[![nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.24.0-brightgreen.svg)](http://nextflow.io)
![CircleCI status](https://circleci.com/gh/CRG-CNAG/CalliNGS-NF.png?style=shield)

## Quickstart 

Install Nextflow by using the following commands: 

    curl get.nextflow.io | bash 
    
Download the Docker image with this command (optional) : 

    docker pull cbcrg/callings-nf@sha256:b65a7d721b9dd2da07d6bdd7f868b04039860f14fa514add975c59e68614c310
    
Note: the Docker image contains all the required dependencies except GATK which 
cannot be included due to license restriction. 

Download the `GenomeAnalysisTK.jar` from [this link](https://software.broadinstitute.org/gatk/download/).    

Launch the pipeline execution with the following command: 

    nextflow run CRG-CNAG/CalliNGS-NF --gatk <path/to/GenomeAnalysisTK.jar>


## Pipeline Description

The RNA sequencing (RNA-seq) data, in additional to the expression information, can be used to obtain somatic variants present in the genes of the analysed organism. The CalliNGS-NF pipeline process RNAseq data to obtain small variants(SNVs), single polymorphisms (SNPs) and small INDELs (insertions, deletions). The pipeline is an implementation of the GATK Best Practices for variant calling on RNAseq and include all major steps of the analysis, [link](http://gatkforums.broadinstitute.org/gatk/discussion/3892/the-gatk-best-practices-for-variant-calling-on-rnaseq-in-full-detail). 

In additional to the GATK best practics, the pipeline includes steps to compare obtained SNVs with known variants and to calculate allele specific counts for the overlapped SNVs.

## Input files

The CalliNGS-NF pipeline needs as the input following files:
* RNAseq reads, <i>*.fastq</i>
* Genome assembly, <i>*.fa</i>
* Known variants, <i>*.vcf</i>

## Pipeline results

For each sample with `sampleID` the pipeline creates a number of output files inside a current working folder.
Here is a brief description of output files:
* <i>`sampleID`.final.vcf</i>,  somatic SNVs called from the RNAseq data
* <i>`sampleID`.diff.sites_in_filesM</i>, comparison of the SNVs from RNAseq data with the set of known variants
* <i>`sampleID`.known.vcf</i>, SNVs that are common between RNAseq calls and known variants
* <i>`sampleID`.ASE.tsv</i>, allele counts at a positions of SNVs (only for common SNVs)
* <i>`sampleID`.FA.hisotgram.pdf</i>, a histogram plot for allele frequency (only for common SNVs)


## Schematic Outline
![Image](../callings-nf-dev/figures/workflow.png?raw=true)

## Requirements 

* Java 7/8
* [Docker](https://www.docker.com/) 1.10 (or higher) or [Sigularity](http://singularity.lbl.gov) engine
* [GATK](https://software.broadinstitute.org/gatk/) 3.7 
