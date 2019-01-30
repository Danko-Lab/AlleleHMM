#Fig4C,D,E and Sup_fig7

#GM
cd /workdir/sc2457/GM_GroSeq_AlleleDB_20170606/test_in_total
Max_Pvalue=1
#locations of pipeline
PL=/workdir/sc2457/tools/After_AlleleDB_pipeline

# filter input files based on Min reads count and Max P-value
R --vanilla --slave --args $(pwd) interestingHets_plus.txt ${Min_count} ${Max_Pvalue} < ${PL}/filter_counts_file.R 
R --vanilla --slave --args $(pwd) interestingHets_minus.txt ${Min_count} ${Max_Pvalue} < ${PL}/filter_counts_file.R 
# make bed files
Input_counts_plus=interestingHets_plus_MinCount1_MaxPvalue1
Input_counts_minus=interestingHets_minus_MinCount1_MaxPvalue1

cat ${Input_counts_plus}.txt |\
awk 'BEGIN{OFS="\t"} {print $1,$2,$19, $20, $21,$14,$15,$16}' |\
awk 'BEGIN{OFS="\t"} (NR>1) {print $1,$2-1, $2 ,$3, $4, $5, $6} ; (NR==1){print "#"$1,$2-1, $2 ,$3, $4, $5, $6} ' > ${Input_counts_plus}.bed

cat ${Input_counts_minus}.txt |\
awk 'BEGIN{OFS="\t"} {print $1,$2,$19, $20, $21,$14,$15,$16}' |\
awk 'BEGIN{OFS="\t"} (NR>1) {print $1,$2-1, $2 ,$3, $4, $5, $6} ; (NR==1){print "#"$1,$2-1, $2 ,$3, $4, $5, $6} ' > ${Input_counts_minus}.bed

#Figure4C
# how many AlleleDB significant SNPs overlap with AlleleHMM blocks?
T=5
#GM
#   44577 ../interestingHets_minus.txt
#   48073 ../interestingHets_plus.txt
intersectBed -sorted -a <(cat ../interestingHets_plus.txt | awk '{OFS="\t"} (NR>1){print "chr"$1, $2-1, $2}') -b counts_plus_hmm_regions_t1e-0${T}_interestingHets_IGV.bed |wc -l #6654
intersectBed -sorted -a <(cat ../interestingHets_minus.txt | awk '{OFS="\t"} (NR>1){print "chr"$1, $2-1, $2}') -b counts_minus_hmm_regions_t1e-0${T}_interestingHets_IGV.bed |wc -l #6432
#F1
intersectBed -sorted -a ../interestingHets_plus_MinCount1_MaxPvalue1.bed -b counts_plus_hmm_regions_t1e-0${T}_interestingHets_IGV.bed |wc -l
intersectBed -sorted -a ../interestingHets_minus_MinCount1_MaxPvalue1.bed -b counts_minus_hmm_regions_t1e-0${T}_interestingHets_IGV.bed |wc -l



#Figure4D
## seperate AlleleDB results to (overlap with AlleleHMM) & (Not overlap with AlleleHMM)
#GM
T=5
Input_counts_plus=interestingHets_plus_MinCount1_MaxPvalue1
Input_counts_minus=interestingHets_minus_MinCount1_MaxPvalue1
intersectBed -sorted -s -a <(cat ../${Input_counts_plus}.bed|awk 'BEGIN{OFS="\t"} (NR>1) {print "chr"$1, $2, $3, $4, $5, "+", $7}') -b counts_plus_hmm_regions_t1e-0${T}_interestingHets_IGV.bed > ${Input_counts_plus}_inAlleleHMM.bed
intersectBed -sorted -s -a <(cat ../${Input_counts_minus}.bed|awk 'BEGIN{OFS="\t"} (NR>1) {print "chr"$1, $2, $3, $4, $5, "-", $7}') -b counts_minus_hmm_regions_t1e-0${T}_interestingHets_IGV.bed > ${Input_counts_minus}_inAlleleHMM.bed
intersectBed -sorted -v -s -a <(cat ../${Input_counts_plus}.bed|awk 'BEGIN{OFS="\t"} (NR>1) {print "chr"$1, $2, $3, $4, $5, "+", $7}') -b counts_plus_hmm_regions_t1e-0${T}_interestingHets_IGV.bed > ${Input_counts_plus}_NotInAlleleHMM.bed
intersectBed -sorted -v -s -a <(cat ../${Input_counts_minus}.bed|awk 'BEGIN{OFS="\t"} (NR>1) {print "chr"$1, $2, $3, $4, $5, "-", $7}') -b counts_minus_hmm_regions_t1e-0${T}_interestingHets_IGV.bed > ${Input_counts_minus}_NotInAlleleHMM.bed


# how many allele-specific SNPs consistent with the gene allele-specificity?
j=gencode.v28lift37.annotation_geneMerged
for h in ${Input_counts_plus} ${Input_counts_minus}
do
# overlap between AlleleDB, AlleleHMM, and allele-specific gene annotation
intersectBed -sorted -s -wb -a ${h}_inAlleleHMM.bed -b ${j}_interestingHets.bed > ${h}_inAlleleHMM_geneMP.bed

# overlap between AlleleDB, AlleleHMM, and all gene annotation
intersectBed -sorted -s -u -a ${h}_inAlleleHMM.bed -b ${j}.bed > ${h}_inAlleleHMM_inGene.bed

# overlap between AlleleDB, AlleleHMM, and symmetric gene annotation
intersectBed -sorted -s -u -a <(intersectBed -sorted -s -wb -v -a ${h}_inAlleleHMM.bed -b ${j}_interestingHets.bed) -b gencode.v28lift37.annotation_geneMerged.bed > ${h}_inAlleleHMM_geneS.bed


intersectBed -sorted -s -wb -a ${h}_NotInAlleleHMM.bed -b ${j}_interestingHets.bed > ${h}_NotInAlleleHMM_geneMP.bed
intersectBed -sorted -s -u -a ${h}_NotInAlleleHMM.bed -b ${j}.bed > ${h}_NotInAlleleHMM_inGene.bed
intersectBed -sorted -s -u -a <(intersectBed -sorted -s -wb -v -a ${h}_NotInAlleleHMM.bed -b ${j}_interestingHets.bed) -b gencode.v28lift37.annotation_geneMerged.bed > ${h}_NotInAlleleHMM_geneS.bed

done

for f in *_geneMP.bed
do 
echo $f
#printf "# of SNPs:"
#cat $f |wc -l

printf "Concordant with Gene:"
cat $f | awk 'BEGIN{OFS="\t"} ($7==$11) {print $0}' |wc -l
printf "Discordant with Gene:"
cat $f | awk 'BEGIN{OFS="\t"} ($7!=$11) {print $0}' |wc -l

fh=`echo $f|rev|cut -d _ -f 2-|rev`
printf "In gene:"
cat ${fh}_inGene.bed |wc -l
done



# SNPs (with at least 1 read mapped) covered by counts_plus_hmm_regions_t1e-0${T}_interestingHets_IGV.bed and counts_minus_hmm_regions_t1e-0${T}_interestingHets_IGV.bed
cd /workdir/sc2457/GM_GroSeq_AlleleDB_20170606/test_in_total/AlleleHMM_result
T=5
#GM
##AlleleHMM but not in AlleleDB
#j=gencode.v28lift37.annotation_geneMerged
for s in plus minus
do echo -e '#chr\tchromStart\tchromEnd\tBlockReadCount\tSNPreadCount\tStrand' > counts_${s}_noX_MinCount1_inAlleleHMM_t1e-0${T}_interestingHets.bed

# counts.txt overlap with AlleleHMM blocks
intersectBed -sorted -wb -a <(cat counts_${s}_noX_MinCount1_MaxPvalue1_IGV.bed| awk '{print "chr"$0}') -b counts_${s}_hmm_regions_t1e-0${T}_interestingHets_IGV.bed | awk '
BEGIN{OFS="\t"} {print $1, $2, $3, $10, $4, $6}' >> counts_${s}_noX_MinCount1_inAlleleHMM_t1e-0${T}_interestingHets.bed

# remove the SNPs overlap with AlleleDB (AlleleHMM and AlleleDB intersect)
head counts_${s}_noX_MinCount1_inAlleleHMM_t1e-0${T}_interestingHets.bed -n 1 > counts_${s}_noX_MinCount1_inAlleleHMM_t1e-0${T}_interestingHets_outAlleleDB.bed
intersectBed -sorted -s -v -a counts_${s}_noX_MinCount1_inAlleleHMM_t1e-0${T}_interestingHets.bed -b interestingHets_${s}_MinCount1_MaxPvalue1_inAlleleHMM.bed >> counts_${s}_noX_MinCount1_inAlleleHMM_t1e-0${T}_interestingHets_outAlleleDB.bed

# overlap AlleleHMM (exclude AlleleDB), and allele-specific gene annotation
intersectBed -sorted -s -wb -a counts_${s}_noX_MinCount1_inAlleleHMM_t1e-0${T}_interestingHets_outAlleleDB.bed -b ${j}_interestingHets.bed > counts_${s}_noX_MinCount1_inAlleleHMM_t1e-0${T}_interestingHets_outAlleleDB_geneMP.bed

# overlap AlleleHMM (exclude AlleleDB), and all gene annotation
intersectBed -sorted -s -u -a counts_${s}_noX_MinCount1_inAlleleHMM_t1e-0${T}_interestingHets_outAlleleDB.bed -b ${j}.bed > counts_${s}_noX_MinCount1_inAlleleHMM_t1e-0${T}_interestingHets_outAlleleDB_inGene.bed

# overlap AlleleHMM (exclude AlleleDB), and symmetric gene annotation
intersectBed -sorted -s -u -a <(intersectBed -sorted -s -v -a counts_${s}_noX_MinCount1_inAlleleHMM_t1e-0${T}_interestingHets_outAlleleDB.bed -b ${j}_interestingHets.bed) -b ${j}.bed > counts_${s}_noX_MinCount1_inAlleleHMM_t1e-0${T}_interestingHets_outAlleleDB_geneS.bed

done


for f in *_outAlleleDB_geneMP.bed
do 
echo $f
#printf "# of SNPs:"
#cat $f |wc -l

printf "Concordant with Gene:"
cat $f | awk 'BEGIN{OFS="\t"} (substr($4,1,1) == $10) {print $0}' |wc -l
printf "Discordant with Gene:"
cat $f | awk 'BEGIN{OFS="\t"} (substr($4,1,1) != $10) {print $0}' |wc -l

fh=`echo $f|rev|cut -d _ -f 2-|rev`
printf "In gene:"
cat ${fh}_inGene.bed |wc -l
done

#Sup_fig7
# the read counts at Concordant. Discordant, and Symmetric regions
for f in *AlleleHMM_geneMP.bed
do echo $f
cat $f | awk 'BEGIN{OFS="\t"} ($7==$11) {print $4+$5}' > ${f}_Concordant_counts.txt
cat $f | awk 'BEGIN{OFS="\t"} ($7!=$11) {print $4+$5}' > ${f}_Discordant_counts.txt
fh=`echo $f|rev|cut -d _ -f 2-|rev`
cat ${fh}_geneS.bed | awk 'BEGIN{OFS="\t"} {print $4+$5}' > ${f}_Symmetric_counts.txt
done

for f in *_outAlleleDB_geneMP.bed
do echo $f
cat $f | awk 'BEGIN{OFS="\t"} (substr($4,1,1) == $10) {print substr($5,3,1)+substr($5,5,1)}' > ${f}_Concordant_counts.txt
cat $f | awk 'BEGIN{OFS="\t"} (substr($4,1,1) != $10) {print substr($5,3,1)+substr($5,5,1)}' > ${f}_Discordant_counts.txt
fh=`echo $f|rev|cut -d _ -f 2-|rev`
cat ${fh}_geneS.bed | awk 'BEGIN{OFS="\t"} {print substr($5,3,1)+substr($5,5,1)}' > ${f}_Symmetric_counts.txt
done

cat counts_*_noX_MinCount1_inAlleleHMM_t1e-05_interestingHets_outAlleleDB_geneMP.bed_Concordant_counts.txt > H_Concordant_counts.txt
cat counts_*_noX_MinCount1_inAlleleHMM_t1e-05_interestingHets_outAlleleDB_geneMP.bed_Discordant_counts.txt > H_Discordant_counts.txt
cat counts_*_noX_MinCount1_inAlleleHMM_t1e-05_interestingHets_outAlleleDB_geneMP.bed_Symmetric_counts.txt  > H_Symmetric_counts.txt
cat interestingHets_*_MinCount1_MaxPvalue1_inAlleleHMM_geneMP.bed_Concordant_counts.txt > A_Concordant_counts.txt
cat interestingHets_*_MinCount1_MaxPvalue1_inAlleleHMM_geneMP.bed_Discordant_counts.txt > A_Discordant_counts.txt
cat interestingHets_*_MinCount1_MaxPvalue1_inAlleleHMM_geneMP.bed_Symmetric_counts.txt > A_Symmetric_counts.txt
cat interestingHets_*_MinCount1_MaxPvalue1_NotInAlleleHMM_geneMP.bed_Concordant_counts.txt > D_Concordant_counts.txt
cat interestingHets_*_MinCount1_MaxPvalue1_NotInAlleleHMM_geneMP.bed_Discordant_counts.txt > D_Discordant_counts.txt
cat interestingHets_*_MinCount1_MaxPvalue1_NotInAlleleHMM_geneMP.bed_Symmetric_counts.txt > D_Symmetric_counts.txt


# how many MP switches in gene annotaiton of AlleleDB SNPs that overlape or not overlap with AlleleHMM?
j=gencode.v28lift37.annotation_geneMerged
cat ${j}.merged_cov_binomtest.bed | 
awk 'BEGIN {OFS="\t"} NR==1 {print $1, $2, $3, "state", $9, "strand", $6, $7, $8}; 
(NR>1 && $6>$7 && $9 <= thresh) {print $1, $2, $3, "M", $9, $4, $6, $7, $8}; 
(NR>1 && $7>$6 && $9 <= thresh) {print $1, $2, $3, "P", $9, $4, $6, $7, $8};
(NR>1 && $9 > thresh)           {print $1, $2, $3, "S", $9, $4, $6, $7, $8};
' thresh=$(awk 'END {print $6}' ${j}.merged_cov_binomtest_FDR.txt)   < ${j}.merged_cov_binomtest.bed > ${j}.merged_cov_binomtest_IGV.bed


bed_file=${j}.merged_cov_binomtest_IGV.bed 

intersectBed -sorted -s -wb -a ${Input_counts_plus}_inAlleleHMM.bed -b ${bed_file} | awk 'BEGIN {OFS="\t"} {print $7,$8,$9,$10,$11}' |uniq | awk 'BEGIN {OFS="\t"} {print $2,$3,$4,$5}'|uniq -c > interestingHets_AlleleDB_in_AlleleHMM_MP_switch_counts.txt
intersectBed -sorted -s -wb -a ${Input_counts_plus}_NotInAlleleHMM.bed -b ${bed_file}| awk 'BEGIN {OFS="\t"} {print $7,$8,$9,$10,$11}' |uniq | awk 'BEGIN {OFS="\t"} {print $2,$3,$4,$5}'|uniq -c > interestingHets_AlleleDB_out_AlleleHMM_MP_switch_counts.txt
intersectBed -sorted -s -wb -a ${Input_counts_minus}_inAlleleHMM.bed -b ${bed_file} | awk 'BEGIN {OFS="\t"} {print $7,$8,$9,$10,$11}' |uniq | awk 'BEGIN {OFS="\t"} {print $2,$3,$4,$5}'|uniq -c >> interestingHets_AlleleDB_in_AlleleHMM_MP_switch_counts.txt
intersectBed -sorted -s -wb -a ${Input_counts_minus}_NotInAlleleHMM.bed -b ${bed_file}| awk 'BEGIN {OFS="\t"} {print $7,$8,$9,$10,$11}' |uniq | awk 'BEGIN {OFS="\t"} {print $2,$3,$4,$5}'|uniq -c >> interestingHets_AlleleDB_out_AlleleHMM_MP_switch_counts.txt

intersectBed -sorted -s -wb -a counts_plus_noX_MinCount1_inAlleleHMM_t1e-0${T}_interestingHets_outAlleleDB.bed -b ${bed_file}| awk 'BEGIN {OFS="\t"} {print substr($4,1,1),$7,$8,$9,$10}' |uniq | awk 'BEGIN {OFS="\t"} {print $2,$3,$4,$5}'|uniq -c > counts_noX_MinCount1_inAlleleHMM_t1e-0${T}_interestingHets_outAlleleDB_switch_counts.txt
intersectBed -sorted -s -wb -a counts_minus_noX_MinCount1_inAlleleHMM_t1e-0${T}_interestingHets_outAlleleDB.bed -b ${bed_file}| awk 'BEGIN {OFS="\t"} {print substr($4,1,1),$7,$8,$9,$10}' |uniq | awk 'BEGIN {OFS="\t"} {print $2,$3,$4,$5}'|uniq -c >> counts_noX_MinCount1_inAlleleHMM_t1e-0${T}_interestingHets_outAlleleDB_switch_counts.txt

R CMD BATCH Compare_allele-specificity_between_gene_AlleleHMM_AlleleDB.R

