
# This makefile runs the Allele Specific Pipeline.  See the README for more information
#

BASE=/home/rdb9/working/Joel/041311

PL:=$(BASE)/Pipeline
SNPS:=$(BASE)/GM12878/snp.calls
CNVS:=$(BASE)/GM12878/CEU.SRP000032.2010_03.genotypes.vcf.cnv 
MAPS:=$(BASE)/GM12878/NA12878_diploid_genome_july30_2010/chr%s_NA12878.map

BINDINGSITES:=hits.bed
FDR_SIMS:=100
FDR_CUTOFF:=0.1
QVAL_CUTOFF:=0.05

bamfiles := $(wildcard *_bam.txt)
queryfiles := $(subst _bam.txt,_eland_query.txt,$(bamfiles))
filterfiles := $(subst _eland_query.txt,_filtered_query.txt, $(queryiles))
bowtiehitfiles = $(subst eland_query,AltRefMother,$(queryfiles)) $(subst eland_query,AltRefFather,$(queryfiles))

TAG:=$(shell date +%m%d%y_%H%M)

all: interestingHets.txt FDR.txt

vars:
	echo $(bamfiles)
	echo $(queryfiles)
	echo $(bowtiehitfiles)

%_eland_query.txt:%_bam.txt
	/usr/local/cluster/software/installation/samtools/samtools-0.1.14/samtools view $< | awk '{OFS="\t"; print ">"$$1"\n"$$10}' > $@

#convert fastq files to eland_query, if needed
#%_eland_query.txt:%_fastq.txt   # figure out how to do this conditionally
#	python $(PL)/fastq2result.py $< $@

# remove any reads containing Ns
%_filtered_query.txt:%_eland_query.txt
	python $(PL)/filter_query.py $< $@

# create joblist for sqPBS.py
job.list: $(filteredfiles)
	python $(PL)/MakeJobList.py $(PL) $(filteredfiles) > job.list

# setup parallel bowtie run file
bowtie.pbs: job.list
	python $(PL)/sqPBS.py d 4 Allele Pipeline job.list > bowtie.pbs # fix the 4 to use length of job.list, fix title

# This is the simple sequential version
#BOWTIEDONE: job.list
#	bash -v job.list
#	touch BOWTIEDONE

# bowtie all reads in parallel against both haplotype genomes
BOWTIEDONE: bowtie.pbs
	rm -f job.list.REMAINING
	qsub bowtie.pbs
	python $(PL)/waitfor.py job.list.REMAINING
	touch BOWTIEDONE

# Merge hits from both haplotypes, keeping best hit.  Also separate bowtie results by chromosome.  
MERGEDONE: BOWTIEDONE $(bowtiehitfiles)
	rm -f combined_AltRef.chr*
	#python $(PL)/MergeDriver.py $(MAPS) | awk '{ print $$0 >> "combined_AltRef."$$3}' # $$ escapes $
	python $(PL)/MergeDriver.py $(MAPS) > combined_AltRef
	awk '{ split($$3,a,"_") ; print $$0 >> "combined_AltRef."a[1]}' < combined_AltRef # $$ escapes $                                                 
	touch MERGEDONE

# process the mapped reads to generate counts at known snp locations
counts.txt: MERGEDONE
	python $(PL)/GetSnpCounts.py 5 $(SNPS) combined_AltRef.chr%s $(MAPS) $(BINDINGSITES) $(CNVS) counts.txt counts.log counts.ks

#filter the counts
counts_passed.txt: counts.txt
	awk -f $(PL)/filter.awk < counts.txt > counts_passed.txt

# augment the counts with a Benjamini Hochberg qvalue, and flag as asymmetric using cutoff 
# This was superceded by FDR calc below, but included in output file anyway.
counts_qval.txt: counts_passed.txt
	python $(PL)/NewAddQValue.py $(QVAL_CUTOFF) counts_passed.txt counts_qval.txt

# calculate false discovery rates
FDR.txt: counts_qval.txt
	python $(PL)/FalsePos.py counts_qval.txt $(FDR_SIMS) $(FDR_CUTOFF) > FDR.txt

# filter for counts that are significant
interestingHets.txt: counts_qval.txt FDR.txt
	awk -f $(PL)/filterPval.awk thresh=$(shell bash $(PL)/getFDR.sh FDR.txt) < counts_qval.txt > interestingHets.txt

# archive a tagged collection of files
archive: 
	mv counts.txt counts.txt.$(TAG)
	mv counts_qval.txt counts_qval.txt.$(TAG)
	mv counts_passed.txt counts_passed.txt.$(TAG)
	mv counts.log counts.log.$(TAG)
	mv interestingHets.txt interestingHets.txt.$(TAG)
	mv FDR.txt FDR.txt.$(TAG)

cleanpbs:
	rm -fr job.list* BOWTIEDONE bowtie.pbs DONE.* SQ_Files* PBS_*

cleanall: cleanpbs
	rm -f combined_AltRef.chr* counts*.txt FDR.txt interestingHets.txt *_AltRef* Merge.log MERGEDONE

.DELETE_ON_ERROR:
