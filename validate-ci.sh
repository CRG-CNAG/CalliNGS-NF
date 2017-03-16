#
# basic output files check 
# improve it for a stronger validation check
# 
[[ -s results/ASER.out ]] && false
[[ -s results/final.vcf ]] && false
[[ -s results/gghist.out.pdf ]] && false
[[ -s results/out.recode.vcf ]] && false
[[ -s results/rep.final.vcf ]] && false
[[ -s results/result.commonSNPs.diff.sites_in_files ]] && false


