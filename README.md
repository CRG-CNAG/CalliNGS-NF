# NGS 2017 Post-Conference Workshop

## CalliNGS-NF
A Nextflow pipeline for basic Variant Calling Analysis with NGS RNA-Seq data


## Quickstart 

Clone the repository in your computer and change to the project root directory: 

    https://github.com/CRG-CNAG/ngs2017ws-nf.git
    cd CRG-CNAG/ngs2017ws-nf.git

Download the Docker image with this command: 

    docker pull cbcrg/ngs2017ws-nf
    
Note: the Docker image contains all the required dependencies except GATK with 
cannot be included due to its license restriction. 

Download the `GenomeAnalysisTK.jar` from [this link](https://software.broadinstitute.org/gatk/download/)
 and copy it in the project root directory.      


Launch the pipeline execution with the following command: 

    NXF_VER=0.24.0-SNAPSHOT nextflow run main.nf 

## Pipeline Description

## Schematic Outline
![alt tag](https://raw.githubusercontent.com/CRG-CNAG/ngs2017ws-nf/callings-nf/figures/workflow.png)


