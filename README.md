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

The RNA sequencing (RNA-seq) data, in additional to the expression information, can be used to obtain information about somatic variants present in the genes of the organism. The CalliNGS-NF pipeline process RNAseq data to obtain small variants(SNVs), single polymorphisms (SNPs) and small INDELs (insertions, deletions). This pipeline is an implementation of the GATK Best Practices for variant calling on RNAseq and include all major steps of the analysis, [link](http://gatkforums.broadinstitute.org/gatk/discussion/3892/the-gatk-best-practices-for-variant-calling-on-rnaseq-in-full-detail). 

In additional to the GATK best practics, the pipeline includes steps to compare obtained SNVs with known variants and to calculate allele specific counts for the overlapped SNVs.

## Input files

The CalliNGS-NF pipeline needs as the imput following files:
* RNAseq reads, *.fastq format
* Genome assembly, *.fa format
* Known variants, *.vcf format

## Pipeline results

The pipeline creates a number of output files inside a current working folder. 

## Schematic Outline
![Image](../callings-nf-dev/figures/workflow.png?raw=true)

## Requirements 

* Java 7/8
* [Docker](https://www.docker.com/) 1.10 (or higher) or [Sigularity](http://singularity.lbl.gov) engine
* [GATK](https://software.broadinstitute.org/gatk/) 3.7 
