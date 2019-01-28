#Sup_Fig2B
## How many SNPs per gene using GENECODE annotation (merged gene)
#F1
unfiltered_snp=/workdir/sc2457/mouse_AlleleSpecific/mouse_genome.sanger.ac.uk/REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.129S1_indelsNsnps_CAST.subsample.bam/P.CAST_M.129S1_indelsNsnps_CAST.subsample.bam.snp.unfiltered
cat ${unfiltered_snp} |awk '{OFS="\t"; c="chr"}{print c$1, $2-1, $2, $6 }' > temp.snp.bed
intersectBed -sorted -a temp.snp.bed -b gencode.vM17.annotation_geneMerged.bed -wb | cut -f 5-8 |LC_ALL=C  sort --parallel=30|uniq -c  > gencode.vM17.annotation_geneMerged_SNPCount.txt
intersectBed -sorted -a temp.snp.bed -b gencode.vM17.annotation_gene.bed -wb | cut -f 5-8 |LC_ALL=C  sort --parallel=30|uniq -c  > gencode.vM17.annotation_gene_SNPCount.txt

#GM
snp_bed=/workdir/sc2457/SNP/1000genome_vol1.ftp.release.20130502/snp.calls.bed
intersectBed -sorted -a ${snp_bed} -b gencode.v28lift37.annotation_geneMerged.bed -wb | cut -f 5-8 |LC_ALL=C  sort --parallel=30|uniq -c  > gencode.v28lift37.annotation_geneMerged_SNPCount.txt
intersectBed -sorted -a ${snp_bed} -b gencode.v28lift37.annotation_gene.bed -wb | cut -f 5-8 |LC_ALL=C  sort --parallel=30|uniq -c  > gencode.v28lift37.annotation_gene_SNPCount.txt

# run SNPs_per_gene.R to get Sup_Fig2B