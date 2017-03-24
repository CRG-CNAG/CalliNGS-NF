# NGS 2017 Post-Conference Workshop

## CalliNGS-NF
A Nextflow pipeline for basic Variant Calling Analysis with NGS RNA-Seq data

## Quickstart 

Clone the repository in your computer and change to the project root directory: 

    git clone https://github.com/CRG-CNAG/CalliNGS-NF.git
    cd CalliNGS-NF

Download the Docker image with this command: 

    docker pull cbcrg/callings-nf@sha256:b65a7d721b9dd2da07d6bdd7f868b04039860f14fa514add975c59e68614c310
    
Note: the Docker image contains all the required dependencies except GATK which 
cannot be included due to its license restriction. 

Download the `GenomeAnalysisTK.jar` from [this link](https://software.broadinstitute.org/gatk/download/)
 and copy it in the project root directory.      


Launch the pipeline execution with the following command: 

    nextflow run main.nf 


## Pipeline Description


## Schematic Outline
![Image](../callings-nf-dev/figures/workflow.png?raw=true)
