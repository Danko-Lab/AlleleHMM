#Fig 5C
# how many genes in EACH hmm blocks?
# GM
cat gencode.v28lift37.annotation.bed |grep chr| cut -f 1-8 |awk 'BEGIN {OFS="\t"} ($8=="gene"){print $0}' |LC_ALL=C sort --temporary-directory=/workdir/sc2457/tmp/ --parallel=10 -V  > gencode.v28lift37.annotation_gene.bed
bedtools merge -i gencode.v28lift37.annotation_gene.bed -s -o collapse,distinct,distinct  -c 4,5,6 |awk 'BEGIN {OFS="\t"; a="111"} {print $1,$2,$3,$4,a,$6}' |LC_ALL=C sort --temporary-directory=/workdir/sc2457/tmp/ --parallel=10 -V > gencode.v28lift37.annotation_geneMerged.bed
Gene_bed=gencode.v28lift37.annotation_geneMerged.bed
intersectBed -sorted -a ${Gene_bed} -b counts_plus_hmm_regions_t1e-0${T}_interestingHets_IGV.bed -s -wb| cut -f 7-10 |uniq -c > counts_plus_hmm_regions_t1e-0${T}_interestingHets_IGV_genePerBlock.txt
intersectBed -sorted -a ${Gene_bed} -b counts_minus_hmm_regions_t1e-0${T}_interestingHets_IGV.bed -s -wb| cut -f 7-10 |uniq -c > counts_minus_hmm_regions_t1e-0${T}_interestingHets_IGV_genePerBlock.txt

wc -l counts_*.bed
# 2015 counts_minus_hmm_regions_t1e-05_interestingHets_IGV.bed
#  2013 counts_plus_hmm_regions_t1e-05_interestingHets_IGV.bed
#  4028 total


# how many genes in EACH hmm blocks?
# F1
Gene_bed=gencode.vM17.annotation_geneMerged.bed
intersectBed -sorted -a ${Gene_bed} -b counts_plus_hmm_regions_t1e-0${T}_interestingHets_IGV.bed -s -wb| cut -f 7-10 |uniq -c > counts_plus_hmm_regions_t1e-0${T}_interestingHets_IGV_genePerBlock.txt
intersectBed -sorted -a ${Gene_bed} -b counts_minus_hmm_regions_t1e-0${T}_interestingHets_IGV.bed -s -wb| cut -f 7-10 |uniq -c > counts_minus_hmm_regions_t1e-0${T}_interestingHets_IGV_genePerBlock.txt

wc -l counts_*.bed
#  1769 counts_minus_hmm_regions_t1e-05_interestingHets_IGV.bed
#  1716 counts_plus_hmm_regions_t1e-05_interestingHets_IGV.bed
#  3485 total


R CMD BATCH Stats_of_HMM_blocks.R