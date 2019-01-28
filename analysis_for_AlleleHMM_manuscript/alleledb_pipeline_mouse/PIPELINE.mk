BASE:=/gpfs/scratch/fas/gerstein/jc2296/personal_genomes/pgenomes/NA12891_pgenome_hg19
PL:=/gpfs/scratch/fas/gerstein/jc2296/software/AlleleSeq_pipeline_v1.2.rob.new
SNPS:=/gpfs/scratch/fas/gerstein/jc2296/personal_genomes/pgenomes/NA12891_pgenome_hg19/snp.calls
CNVS:=/gpfs/scratch/fas/gerstein/jc2296/personal_genomes/pgenomes/NA12891_pgenome_hg19/cnv_rd/rd.hiseq.cnvnator.NA12891.snp.calls
MAPS:=$(BASE)/%s_NA12891.map
BNDS:=hits.bed
FDR_SIMS:=6
FDR_CUTOFF:=0.05

PREFIX:=NULL

sourcefiles:=$(PREFIX)     
countfiles:=$(PREFIX).cnt  

## bowtie files that contain aligned reads
MATBOWTIE:=$(PREFIX).mat.bowtie
PATBOWTIE:=$(PREFIX).pat.bowtie

## file that contains read IDs to be removed from bowtie files above
READS2FILTER:=originalmatpatreads.toremove.ids



all: interestingHets.txt


$(countfiles): $(PATBOWTIE) $(MATBOWTIE)
	bash -c "python $(PL)/MergeBowtie.py \
           <(python $(PL)/filter_reads_out.py $(PATBOWTIE) - $(READS2FILTER)) \
           <(python $(PL)/filter_reads_out.py $(MATBOWTIE) - $(READS2FILTER)) \
           $(MAPS) | python $(PL)/SnpCounts.py $(SNPS) - $(MAPS) $@"





check:
	@echo $(sourcefiles)


counts.txt: $(countfiles)
	python $(PL)/CombineSnpCounts.py 6 $(SNPS) $(BNDS) $(CNVS) counts.txt counts.log $(countfiles)

# calculate false discovery rates
FDR.txt: counts.txt
	python $(PL)/FalsePos.py counts.txt $(FDR_SIMS) $(FDR_CUTOFF) > FDR.txt

interestingHets.txt: counts.txt FDR.txt
	awk -f $(PL)/finalFilter.awk thresh=$(shell awk 'END {print $$6}' FDR.txt) < counts.txt > interestingHets.txt

clean:
	@rm -f FDR.txt interestingHets.txt counts.txt

cleanall: clean
	@rm -f *.cnt

.DELETE_ON_ERROR:
