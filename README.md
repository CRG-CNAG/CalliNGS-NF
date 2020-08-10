# CalliNGS-NF
A Nextflow pipeline for Variant Calling Analysis with NGS RNA-Seq data based on GATK best practices.

[![nextflow](https://img.shields.io/badge/nextflow-%E2%89%A520.01.0-brightgreen.svg)](http://nextflow.io)
[![Build Status](https://travis-ci.org/CRG-CNAG/CalliNGS-NF.svg?branch=master)](https://travis-ci.org/CRG-CNAG/CalliNGS-NF)

## Quickstart 

Install Nextflow by using the following command: 

    curl -s https://get.nextflow.io | bash 
    
Download the Docker image with this command (optional) : 

    docker pull cbcrg/callings-nf:gatk4


Launch the pipeline execution with the following command: 

    nextflow run CRG-CNAG/CalliNGS-NF -profile docker

Note: the Docker image contains all the required dependencies. Add the `-profile docker` 
  to enable the containerised execution to the example command line shown below. 

## Pipeline Description

The RNA sequencing (RNA-seq) data, in additional to the expression information, can be used to obtain somatic variants present in the genes of the analysed organism. The CalliNGS-NF pipeline processes RNAseq data to obtain small variants(SNVs), single polymorphisms (SNPs) and small INDELs (insertions, deletions). The pipeline is an implementation of the GATK best practices for variant calling on RNAseq and includes all major steps of the analysis, [link](http://gatkforums.broadinstitute.org/gatk/discussion/3892/the-gatk-best-practices-for-variant-calling-on-rnaseq-in-full-detail). 

In addition to the GATK best practics, the pipeline includes steps to compare obtained SNVs with known variants and to calculate allele specific counts for the overlapped SNVs.

## Input files

The CalliNGS-NF pipeline needs as the input following files:
* RNAseq reads, `*.fastq`
* Genome assembly, `*.fa`
* Known variants, `*.vcf`
* Denylisted regions of the genome, `*.bed`

The RNAseq read file names should match the following naming convention:  *sampleID{1,2}_{1,2}.extension* 

where: 
* *sampleID* is the identifier of the sample;
* the first number **1** or **2** is the replicate ID;
* the second number **1** or **2** is the read pair in the paired-end samples;
* *extension* is the read file name extension eg. `fq`, `fq.gz`, `fastq.gz`, etc. 

example: `ENCSR000COQ1_2.fastq.gz`.

## Pipeline parameters

#### `--reads` 
   
* Specifies the location of the reads FASTQ file(s).
* Multiple files can be specified using the usual wildcards (*, ?), in this case make sure to surround the parameter string
  value by single quote characters (see the example below)
* By default it is set to the CalliNGS-NF's location: `$baseDir/data/reads/rep1_{1,2}.fq.gz`
* See above for naming convention of samples, replicates and pairs read files.

Example: 

    $ nextflow run CRG-CNAG/CalliNGS-NF --reads '/home/dataset/*_{1,2}.fq.gz'


#### `--genome`

* The location of the genome fasta file.
* It should end in `.fa`.
* By default it is set to the CalliNGS-NF's location: `$baseDir/data/genome.fa`.

Example:

    $ nextflow run CRG-CNAG/CalliNGS-NF --genome /home/user/my_genome/human.fa
    

#### `--variants`

* The location of the known variants VCF file.
* It should end in `.vcf` or `vcf.gz`.
* By default it is set to the CalliNGS-NF's location: `$baseDir/data/known_variants.vcf.gz`.

Example:

    $ nextflow run CRG-CNAG/CalliNGS-NF --variants /home/user/data/variants.vcf


#### `--denylist` (formely `--blacklist`)

* The location of the denylisted genome regions in bed format.
* It should end in `.bed`.
* By default it is set to the CalliNGS-NF's location: `$baseDir/data/denylist.bed`.

Example:

    $ nextflow run CRG-CNAG/CalliNGS-NF --denylist /home/user/data/denylisted_regions.bed


#### `--results` 
   
* Specifies the folder where the results will be stored for the user.  
* It does not matter if the folder does not exist.
* By default is set to CalliNGS-NF's folder: `results` 

Example: 

    $ nextflow run CRG-CNAG/CalliNGS-NF --results /home/user/my_results
    

    
    
## Pipeline results

For each sample the pipeline creates a folder named `sampleID` inside the directory specified by using the `--results` command line option (default: `results`).
Here is a brief description of output files created for each sample:

file | description 
---- | ----
`final.vcf` | somatic SNVs called from the RNAseq data
`diff.sites_in_files` | comparison of the SNVs from RNAseq data with the set of known variants
`known_snps.vcf` | SNVs that are common between RNAseq calls and known variants
`ASE.tsv` | allele counts at a positions of SNVs (only for common SNVs)
`AF.histogram.pdf` | a histogram plot for allele frequency (only for common SNVs)


## Schematic Outline
![Image](../master/figures/workflow.png?raw=true)

## Requirements 

* [Nextflow](https://www.nextflow.io) 20.07.1 (or later)
* Java 8 or later
* [Docker](https://www.docker.com/) 1.10 (or later) or [Singularity](http://singularity.lbl.gov) engine
* [GATK](https://gatk.broadinstitute.org/) 4.1.x 

Note: CalliNGS-NF can be used without a container engine by installing in your system all the 
required software components reported in the following section. See the included 
[Dockerfile](docker/Dockerfile) for the configuration details.
 

## Components 

CalliNGS-NF uses the following software components and tools: 

* Java 8 
* Samtools 1.3.1
* Vcftools 0.1.14
* STAR 2.5.2b
* GATK 4.1
* R 3.1.1 
* Awk
* Perl
* Grep
