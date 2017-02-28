#!/bin/bash
#$ -N processVcf
#$ -cwd

time grep -v '#'  results/final.vcf | awk '$7~/PASS/' |perl -ne 'chomp($_); ($dp)=$_=~/DP\=(\d+)\;/; if($dp>=8){print $_."\n"};' > results/final.DP8.vcf
#real    0m0.173s
#user    0m0.038s
#sys     0m0.018s

#5763 SNPs with GATK and 3826 with samtools

#compare consensus SNPs vs one called by Illumina platinum genomes - select those overlap (AE expression) and those unique (RNAediting)
time vcftools --vcf results/final.DP8.vcf --diff snv/NA12878.chr22.vcf.gz.filtered.recode.vcf  --diff-site --out results/commonSNPs
#Found 3936 sites common to both files. 68%  (1979, 52% with samtools)
#Found 1819 sites only in main file.
#Found 48538 sites only in second file.
#Found 7 non-matching overlapping sites. (97 with samtools)
#After filtering, kept 5762 out of a possible 5762 Sites
#Run Time = 0.00 seconds

#real    0m0.432s
#user    0m0.226s
#sys     0m0.088s

#run AE expression quantification for both replicates independently, then I will compare results

#select only common SNPs
time awk 'BEGIN{OFS="\t"} $4~/B/{print $1,$2,$3}' results/commonSNPs.diff.sites_in_files  > test.bed
#real    0m0.453s
#user    0m0.036s
#sys     0m0.008s

time vcftools --vcf  results/final.vcf --bed test.bed --recode --keep-INFO-all
#real    0m1.041s
#user    0m0.912s
#sys     0m0.011s

#then AE_pipeline... per each replicate? or to merge bam?

#allele frequency plot
time grep -v '#'  out.recode.vcf|awk -F '\t' '{print $10}' |awk -F ':' '{print $2}'|perl -ne 'chomp($_); @v=split(/\,/,$_); if($v[0]!=0 ||$v[1] !=0){print  $v[1]/($v[1]+$v[0])."\n"; }' |awk '$1!=1' >AF.4R
#real    0m0.103s
#user    0m0.071s
#sys     0m0.011s

 time ~abreschi/R/scripts/gghist.R -i AF.4R
#real    0m1.902s
#user    0m0.873s
#sys     0m0.100s

# Those SNPs that are different from Illumina SNPs -> RNAediting. Calculate frequency of events 
time awk '$4~/1/{print $5"_"$7}' results/commonSNPs.diff.sites_in_files |sort|uniq -c|sort -k1,1nr|head
#real    0m0.057s
#user    0m0.041s
#sys     0m0.010s

#    688 A_G  A->I,  I behaves as if it is G both in translation and when forming secondary structures. (wikipedia)
#    578 T_C  A->I , just reversed
#     74 G_T 	
#     65 G_A
#     64 C_T  C->T 
#     47 C_A
#     33 C_CA
#     23 T_A
#     22 G_C
#...



