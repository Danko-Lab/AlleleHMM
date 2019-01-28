#### Redo H3K27me3 reads in GRO-seq regions
cd /workdir/sc2457/GM_GroSeq_AlleleDB_20170606/test_in_total/AlleleHMM_result
for s in plus minus
do cat counts_${s}_hmm_regions_t1e-0${T}.merged_cov_binomtest.bed |LC_ALL=C sort -k1,1V -k2,2n | awk 'BEGIN {OFS="\t"; t=","; s=" ";c="chr"} NR==1 { print $1,$2,$3,$4t$6t$7, 111, s} 
NR>1 {print $1,$2,$3,$4t$6t$7, 111, s}' > counts_${s}_hmm_regions_t1e-0${T}.merged_cov_binomtest_IGV.bed &
done

T=5
FDR_SIMS=10; FDR_CUTOFF=0.1
PL=/workdir/sc2457/alleleDB/alleledb_pipeline
for s in plus minus
do input_f=/workdir/sc2457/GM_GroSeq_AlleleDB_20170606/test_in_total/AlleleHMM_result/counts_${s}_hmm_regions_t1e-0${T}.merged_cov_binomtest_IGV.bed

j=H3K27me3_${s}
PATBOWTIE=H3K27me3.pat.bowtie ; MATBOWTIE=H3K27me3.mat.bowtie
  bedtools coverage -a ${input_f} -b ${MATBOWTIE}_AMBremoved_sorted_specific.map2ref.sorted.bed | awk 'BEGIN {OFS="\t"; t="_"} {print $1t$2t$3, $0}'  > ${j}_counts_hmm_regions_t1e-0${T}.mat_cov.bed &
  bedtools coverage -a ${input_f} -b ${PATBOWTIE}_AMBremoved_sorted_specific.map2ref.sorted.bed | awk 'BEGIN {OFS="\t"; t="_"} {print $1t$2t$3, $0}'  > ${j}_counts_hmm_regions_t1e-0${T}.pat_cov.bed &
wait
join -t $'\t' -j 1 -o 1.1,1.2,1.3,1.4,1.5,1.8,2.8 ${j}_counts_hmm_regions_t1e-0${T}.mat_cov.bed ${j}_counts_hmm_regions_t1e-0${T}.pat_cov.bed > ${j}_counts_hmm_regions_t1e-0${T}.temp_cov.bed &
  bedtools coverage -a ${input_f} -b ${PATBOWTIE}_AMBremoved_sorted_identical.map2ref.sorted.bed | awk 'BEGIN {OFS="\t"; t="_"} {print $1t$2t$3, $0}'  > ${j}_counts_hmm_regions_t1e-0${T}.iden_cov.bed &
wait
  join -t $'\t' -j 1 -o 1.2,1.3,1.4,1.5,1.6,1.7,2.8 ${j}_counts_hmm_regions_t1e-0${T}.temp_cov.bed ${j}_counts_hmm_regions_t1e-0${T}.iden_cov.bed > ${j}_counts_hmm_regions_t1e-0${T}.merged_cov.bed &
wait
# here 
FDR_SIMS=10; FDR_CUTOFF=0.1
  python ${PL}/BinomialTestFor_merged_cov.bed.py ${j}_counts_hmm_regions_t1e-0${T}.merged_cov.bed ${j}_counts_hmm_regions_t1e-0${T}.merged_cov_binomtest.bed 
  mv ${j}_counts_hmm_regions_t1e-0${T}.merged_cov_binomtest.bed tmp.bed
  echo -e '#chrm\tchrmStart\tchrmEnd\tH3K27me3_BinomialTest\tGROseq_hmm_state\tmat_allele_count\tpat_allele_count\tidentical_reads_count\tBinom_p_value' > ${j}_counts_hmm_regions_t1e-0${T}.merged_cov_binomtest.bed 
  # DID NOT filter out those regions where mat reads counts is equal to pat reads counts
  cat tmp.bed |awk 'BEGIN{OFS="\t"; l=","} NR >1 && $6>$7 {print $1,$2,$3,"M,"$6l$7,$4,$6,$7,$8,$9} NR >1 && $6==$7 {print $1,$2,$3,"S,"$6l$7,$4,$6,$7,$8,$9}
  NR >1 && $6<$7 {print $1,$2,$3,"P,"$6l$7,$4,$6,$7,$8,$9}' >> ${j}_counts_hmm_regions_t1e-0${T}.merged_cov_binomtest.bed 
  rm tmp.bed

# output of BinomialTestFor_merged_cov.bed.py:(hmm+BinomialTest) if p-value <= 0.05, remain what it got from hmm (can ne M,P, or S), otherwise S.
  python ${PL}/FalsePosFor_merged_cov.bed.py ${j}_counts_hmm_regions_t1e-0${T}.merged_cov_binomtest.bed ${FDR_SIMS} ${FDR_CUTOFF} > ${j}_counts_hmm_regions_t1e-0${T}.merged_cov_binomtest_FDR.txt 
  #rm ${j}_counts_hmm_regions_t1e-0${T}.temp_cov.bed ${j}_counts_hmm_regions_t1e-0${T}.mat_cov.bed ${j}_counts_hmm_regions_t1e-0${T}.pat_cov.bed ${j}_counts_hmm_regions_t1e-0${T}.iden_cov.bed
  #rm ${j}_counts_hmm_regions_t1e-0${T}.merged_cov.bed
  awk 'NR==1 { print $0 } NR>1 && $9 <= thresh { print $0 }'  thresh=$(awk 'END {print $6}' ${j}_counts_hmm_regions_t1e-0${T}.merged_cov_binomtest_FDR.txt) < ${j}_counts_hmm_regions_t1e-0${T}.merged_cov_binomtest.bed > ${j}_counts_hmm_regions_t1e-0${T}_interestingHets.bed
done

R CMD BATCH H3K27me3_GROseq.R