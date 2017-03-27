#!/bin/bash 
#
# basic output files check 
# improve it for a stronger validation check
# 
[[ -s results/rep.final.vcf ]] || false
[[ -s results/rep.diff.sites_in_files ]] || false
[[ -s results/rep.known_snps.vcf ]] || false
[[ -s results/rep.AF.histogram.pdf ]] || false
[[ -s results/rep.ASE.tsv ]] || false


