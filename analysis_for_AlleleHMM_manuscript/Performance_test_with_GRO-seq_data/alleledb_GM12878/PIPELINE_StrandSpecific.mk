BASE=/workdir/sc2457/
PL:= $(BASE)/alleleDB/alleledb_pipeline
SNPS:=$(BASE)/SNP/1000genome_vol1.ftp.release.20130502/snp.call.all
CNVS:=$(BASE)/alleleseq.gersteinlab.org/NA12878_diploid_dec16.2012.alleleseq.input/rd_4327183snps_na12878_hg19.txt
BNDS:=hits.bed
MAPS:=$(BASE)/alleleseq.gersteinlab.org/NA12878_diploid/NA12878_diploid_2015_feb5_3versions/1kgp3-svs-pass_NA12878_hg19_150109_w_transcriptome/%s_NA12878.map
FDR_SIMS:=5
FDR_CUTOFF:=0.1
MIN_READ_COUNT:=1

PREFIX:=NULL

sourcefiles:=$(PREFIX)     
countfiles_plus:=$(PREFIX).plus.cnt  
countfiles_minus:=$(PREFIX).minus.cnt  

## bowtie files that contain aligned reads
MATBOWTIE:=$(PREFIX).mat.bowtie
PATBOWTIE:=$(PREFIX).pat.bowtie
MergedBOWTIE_PLUS:=$(PREFIX).merged.bowtie_plus
MergedBOWTIE_MINUS:=$(PREFIX).merged.bowtie_minus

## file that contains read IDs to be removed from bowtie files above
READS2FILTER:=originalmatpatreads.toremove.ids

target:interestingHets_plus.txt

#Generate Strand Specific countfiles output
$(MergedBOWTIE_PLUS): $(PATBOWTIE) $(MATBOWTIE) $(READS2FILTER)
	bash -c "python $(PL)/MergeBowtie.py \
           <(python $(PL)/filter_reads_out.py $(PATBOWTIE) - $(READS2FILTER)) \
           <(python $(PL)/filter_reads_out.py $(MATBOWTIE) - $(READS2FILTER)) \
           $(MAPS) | python ${PL}/seperate_strand_of_bowtie_output_alleleDB.py $(PREFIX).merged.bowtie"


$(countfiles_plus): $(MergedBOWTIE_PLUS)
	bash -c "python $(PL)/SnpCounts.py $(SNPS) $(MergedBOWTIE_PLUS) $(MAPS) $@"
 
$(MergedBOWTIE_MINUS): $(MergedBOWTIE_PLUS)
$(countfiles_minus): $(MergedBOWTIE_MINUS)
	bash -c "python $(PL)/SnpCounts.py $(SNPS) $(MergedBOWTIE_MINUS) $(MAPS) $@"


check:
	@echo $(sourcefiles)


counts_plus.txt: $(countfiles_plus) $(countfiles_minus)
	python $(PL)/CombineSnpCounts.py $(MIN_READ_COUNT) $(SNPS) $(BNDS) $(CNVS) counts_plus.txt counts_plus.log $(countfiles_plus)
	python $(PL)/CombineSnpCounts.py $(MIN_READ_COUNT) $(SNPS) $(BNDS) $(CNVS) counts_minus.txt counts_minus.log $(countfiles_minus)

# calculate false discovery rates
FDR_plus.txt: counts_plus.txt
	python $(PL)/FalsePos.py counts_plus.txt $(FDR_SIMS) $(FDR_CUTOFF) > FDR_plus.txt
	python $(PL)/FalsePos.py counts_minus.txt $(FDR_SIMS) $(FDR_CUTOFF) > FDR_minus.txt

interestingHets_plus.txt: counts_plus.txt FDR_plus.txt
	awk -f $(PL)/finalFilter.awk thresh=$(shell awk 'END {print $$6}' FDR_plus.txt) < counts_plus.txt > interestingHets_plus.txt
	awk -f $(PL)/finalFilter.awk thresh=$(shell awk 'END {print $$6}' FDR_minus.txt) < counts_minus.txt > interestingHets_minus.txt

clean:
	@rm -f FDR*.txt interestingHets*.txt counts*.txt

cleanall: clean
	@rm -f *.cnt

.DELETE_ON_ERROR:

