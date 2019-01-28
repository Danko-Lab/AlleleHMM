#!/bin/bash
### alleledb.sh <indiv_category> <pgenome> <pgenome path> <hetsnp file; snp.calls.bed> <FASTQ_PATH> <fastq.gz/basename> <pipeline path> <PIPELINE.mk filename> <FDR threshold> <asb/ase> <flag>
### e.g. alleledb.sh NA12878_kilpinen-rnaseq alleledb_pg /gpfs/scratch/fas/gerstein/jc2296/personal_genomes/NA12878_pgenome_hg19/hiseq_NA12878_hg19 /gpfs/scratch/fas/gerstein/jc2296/personal_genomes/NA12878_pgenome_hg19/hiseq_NA12878_hg19/snp.calls.bed /scratch/fas/gerstein/jc2296/alleledb/alleleseq-runs/trios/NA12878/NA12878_hiseq/rnaseq/combined_ASE_POLYA_NA12878 kilpinen_NA12878_ERR356372_1.fastq.gz ~/software/AlleleSeq_pipeline_v1.2.rob.new PIPELINE_mod.mk.fdr0.05 0.05 ase 0
## note for paths no end '/' pls
## note also that the script grabs anything in the FASTQ_PATH with fastq.gz extention; hence the FASTQ_NAME field is only used as a basename
## requires the following scripts: intersectBed from bedtools, Bowtie1, map.back.ref.wrapper.sh.ori,  multis-and-unaligneds-wrapper.sh 
## note that the original fastqs are all post-modified to remove tabs and spaces; the original fastqs are not themselves changed but their output will have spaces replaced with '_', essentially: sed 's/[ \t]\+/_/g'



## this causes the script to exit if there's any error in any part of the pipeline 
## The shell does NOT exit if the command that fails is part of the command list 
## immediately following a while or until keyword, part of the test in an if statement, 
## part of a && or || list, or if the command's return value is being inverted via !
#set -e


NAME=$1  ##NA12878_kilpinen-rnaseq
PGENOME=$2  ##name of personal genome newSV
PGENOME_PATH=$3 ##/path/to/pgenome; no '/' at the end
SNPCALLS_PATH=$4 ##/path/to/snp.calls.bed; note this is the explicit file location

FASTQ_PATH=$5 ##/path/to/fastq
FASTQ=$6 ##basename/fastq.gz

PL=$7 ##/path/to/pipeline
PIPELINEFILE=$8 ## PIPELINE.mk filename
FDR_THRESH=$9

AS=${10} ##asb/ase
FLAG=${11} ##see below 

#########################
### FLAG CODE   #########
### 0:run everything
### 1:alignment1
### 2:map back to reference
### 3:intersectBed between reads and hetSNVs
### 4:flip alleles on reads that overlap >=1 hetSNV 
### 5:alignment2 of flipped reads and check if locations match
### 6:unaligned reads in a separate folder for analysis
### 7:multi reads in a separate folder for analysis
### 8:alleleseq run on bias-filtered fastq
###-8: if $1 -eq -8, clean up everything from module 8; creates a folder 'trash'
### 	 e.g. alleledb.sh -8

#############################################################
## -8 clean up everything from module 8
#############################################################
if [[ $1 -eq -8 ]] ; then

mkdir trash
mv *.R *.Rout .RData trash
mv betabinomial counts.* FDR* interestingHets.* original* trash
mv *.cnt *.bowtie trash
mv *.log trash

exit 0

fi

#############################
## print parameters to log
#############################
echo "NAME          =$1"  
echo "PGENOME       =$2"  
echo "PGENOME_PATH  =$3" 
echo "SNPCALLS_PATH =$4" 
echo "FASTQ_PATH    =$5"
echo "BASENAME/FASTQ=$6"
echo "PIPELINE/PL   =$7"
echo "PIPELINE_FNAME=$8"
echo "FDR_THRESHOLD =$9"
echo "ASB/ASE       =${10}"
echo "FLAG          =${11}"


#############################################################
## create allelic bias folder
#############################################################

##_ mkdir allelicbias-${PGENOME}-${NAME}

if [ ! -d allelicbias-${PGENOME}-${NAME} ]; then
	mkdir allelicbias-${PGENOME}-${NAME}
fi

cd allelicbias-${PGENOME}-${NAME}

#############################################################
## 1 alignment
#############################################################
if [[ ${FLAG} -eq 1 || ${FLAG} -eq 0 ]] ; then
	## setting up
	## make directories, go in directory, set up links
	##_ mkdir 1-alignment-${NAME}
	if [ -d 1-alignment-${NAME} ]; then
		rm -r 1-alignment-${NAME}
	fi
	mkdir 1-alignment-${NAME}
	cd 1-alignment-${NAME}
	ln -s ${FASTQ_PATH}
	
	## start log print
	echo "######################" > 1-alignment-${NAME}.log
	echo "## 1-alignment #######" >> 1-alignment-${NAME}.log
	echo "######################" >> 1-alignment-${NAME}.log
	
	date >> 1-alignment-${NAME}.log
	
	## align mat and pat
	echo ${FASTQ_PATH}/${FASTQ}.fastq.gz
	
	zcat ${FASTQ_PATH}/${FASTQ}.fastq.gz | gzip -c | ${PL}/alleledb_filter_input.sh ${PL} - | bowtie --best --strata -p 20 -v 2 -m 1 -f ${PGENOME_PATH}/AltRefMother/AltRefMother - > ${FASTQ}.mat.bowtie 2> ${FASTQ}.mat.log & 
	zcat ${FASTQ_PATH}/${FASTQ}.fastq.gz | gzip -c | ${PL}/alleledb_filter_input.sh ${PL} - | bowtie --best --strata -p 20 -v 2 -m 1 -f ${PGENOME_PATH}/AltRefFather/AltRefFather - > ${FASTQ}.pat.bowtie 2> ${FASTQ}.pat.log &
	
	#zcat ${FASTQ_PATH}/*.fastq.gz | gzip -c | ${PL}/alleledb_filter_input.sh ${PL} - | bowtie --best --strata -v 2 -m 1 -f ${PGENOME_PATH}/${NAME}_maternal - > ${FASTQ}.mat.bowtie 2> ${FASTQ}.mat.log & 
	#zcat ${FASTQ_PATH}/*.fastq.gz | gzip -c | ${PL}/alleledb_filter_input.sh ${PL} - | bowtie --best --strata -v 2 -m 1 -f ${PGENOME_PATH}/${NAME}_paternal - > ${FASTQ}.pat.bowtie 2> ${FASTQ}.pat.log &
	wait

	## postprocess
	echo -e "${FASTQ_PATH}/${FASTQ}.fastq.gz\n$(zcat ${FASTQ_PATH}/${FASTQ}.fastq.gz | wc -l)" >> 1-alignment-${NAME}.log
	for i in ${FASTQ_PATH}/${FASTQ}.fastq.gz; do echo -e "$i $(zcat $i | wc -l)" >> 1-alignment-${NAME}.log; done
	wc -l ${FASTQ}.mat.bowtie >> 1-alignment-${NAME}.log
	wc -l ${FASTQ}.pat.bowtie >> 1-alignment-${NAME}.log
	
	date >> 1-alignment-${NAME}.log

	cd ..
fi





#############################################################
## 2 map to ref
#############################################################
if [[ ${FLAG} -eq 2 || ${FLAG} -eq 0 ]] ; then
	## setting up
	## make directories, go in directory, set up links
	##_mkdir 2-map.back.ref-${NAME}
	if [ -d 2-map.back.ref-${NAME} ]; then
		rm -r 2-map.back.ref-${NAME}
	fi
	mkdir 2-map.back.ref-${NAME}
	cd 2-map.back.ref-${NAME}
	
	## start log print
	echo "#######################" > 2-map.back.ref-${NAME}.log
  echo "## 2-map to ref #######" >> 2-map.back.ref-${NAME}.log
  echo "#######################" >> 2-map.back.ref-${NAME}.log
  
  date >> 2-map.back.ref-${NAME}.log
  
	
	## chain files
	ln -s ${PGENOME_PATH}/mat2ref.chain
	ln -s ${PGENOME_PATH}/pat2ref.chain

	for i in ../1-alignment-${NAME}/*.bowtie
	do
	ln -s $i
	done

	## map
	## there is a lot of compute here
	${PL}/alleledb_map.back.ref.wrapper_forMouse.sh ${FASTQ}.mat.bowtie maternal MAT mat2ref.chain ${PL} &
	${PL}/alleledb_map.back.ref.wrapper_forMouse.sh ${FASTQ}.pat.bowtie paternal PAT pat2ref.chain ${PL} &
	
	wait

	## postprocess
	wc -l *.maternal.*.bed >> 2-map.back.ref-${NAME}.log
	wc -l *.paternal.*.bed >> 2-map.back.ref-${NAME}.log

	#### clean up/debugging code ####
	#cat \
	#${FASTQ}.[mp]at.bowtie.[1-9]_[mp]aternal.bed \
	#${FASTQ}.[mp]at.bowtie.1[0-9]_[mp]aternal.bed \
	#${FASTQ}.[mp]at.bowtie.2[0-2]_[mp]aternal.bed \
	#${FASTQ}.[mp]at.bowtie.[XY]_[mp]aternal.bed > ${FASTQ}.[mp]at.bowtie.${MATPATERNAL}.bed
	#mkdir toremove
    rm ${FASTQ}.mat.bowtie.*_maternal.bed ${FASTQ}.mat.bowtie.*_maternal.bowtie
	rm ${FASTQ}.pat.bowtie.*_paternal.bed ${FASTQ}.pat.bowtie.*_paternal.bowtie 
	
	date >> 2-map.back.ref-${NAME}.log
	
	cd ..
fi




#############################################################
## 3  intersectBed
#############################################################
if [[ ${FLAG} -eq 3 || ${FLAG} -eq 0 ]] ; then
	
	## set up
	##_mkdir 3-intersectBed-${NAME}
	if [ -d 3-intersectBed-${NAME} ] ; then
		rm -r 3-intersectBed-${NAME}
	fi
	mkdir 3-intersectBed-${NAME}
	cd 3-intersectBed-${NAME}
	ln -s ${SNPCALLS_PATH} snp.calls.bed
	
	
	## start log printing	
	echo "##########################" > 3-intersectBed-${NAME}.log
  echo "## 3-intersectBed  #######" >> 3-intersectBed-${NAME}.log
  echo "##########################" >> 3-intersectBed-${NAME}.log
  
  date >> 3-intersectBed-${NAME}.log

	for i in ../2-map.back.ref-${NAME}/*.map2ref.bed
	do
	ln -s $i
	done
	
	## intersect
	## note that my chain files do not have prefix "chr"; could do away if present
	##_intersectBed -a <(awk '{OFS="\t"}{FS="\t"}{print "chr"$0}' ${FASTQ}.mat.bowtie.maternal.map2ref.bed) -b snp.calls.bed -wa -wb > intersect.${FASTQ}.mat.snp.calls.txt &
	##_intersectBed -a <(awk '{OFS="\t"}{FS="\t"}{print "chr"$0}' ${FASTQ}.pat.bowtie.paternal.map2ref.bed) -b snp.calls.bed -wa -wb > intersect.${FASTQ}.pat.snp.calls.txt &
        intersectBed -a <(awk '{OFS="\t"}{FS="\t"}{print "chr"$0}' ${FASTQ}.mat.bowtie.maternal.map2ref.bed) -b snp.calls.bed -wa -wb > intersect.${FASTQ}.mat.snp.calls.txt &
	intersectBed -a <(awk '{OFS="\t"}{FS="\t"}{print "chr"$0}' ${FASTQ}.pat.bowtie.paternal.map2ref.bed) -b snp.calls.bed -wa -wb > intersect.${FASTQ}.pat.snp.calls.txt &
	wait
	
	## preprocess
	wc -l intersect.${FASTQ}.mat.snp.calls.txt >> 3-intersectBed-${NAME}.log
	wc -l intersect.${FASTQ}.pat.snp.calls.txt >> 3-intersectBed-${NAME}.log
	
	date >> 3-intersectBed-${NAME}.log
	
	cd ..
fi




#############################################################
## 4 flip the reads
#############################################################
if [[ ${FLAG} -eq 4 || ${FLAG} -eq 0 ]] ; then
	## set up
	##_mkdir 4-flip-${NAME}
	if [ -d 4-flip-${NAME} ] ; then
		rm -r 4-flip-${NAME}
	fi
	mkdir 4-flip-${NAME}
	cd 4-flip-${NAME} 
	
	## start log printing
	echo "############################" > 4-flip-${NAME}.log
  echo "## 4-flip the reads  #######" >> 4-flip-${NAME}.log
  echo "############################" >> 4-flip-${NAME}.log
  
  date >> 4-flip-${NAME}.log
  
	
	for i in ../3-intersectBed-${NAME}/intersect.*.snp.calls.txt
	do
	ln -s $i
	done
	
	## flip and convert to fastq
	sort -nk4,4 intersect.${FASTQ}.mat.snp.calls.txt | ${PL}/flipread flipread2fastq 0 5 intersect.${FASTQ}.mat.snp.calls stdin | gzip -c > intersect.${FASTQ}.mat.flipread.fastq.gz &
	sort -nk4,4 intersect.${FASTQ}.pat.snp.calls.txt | ${PL}/flipread flipread2fastq 0 5 intersect.${FASTQ}.pat.snp.calls stdin | gzip -c > intersect.${FASTQ}.pat.flipread.fastq.gz &
	
	wait
	
	## postprocess
	echo -e "$(zcat intersect.${FASTQ}.mat.flipread.fastq.gz | wc -l) intersect.${FASTQ}.mat.flipread.fastq.gz" >> 4-flip-${NAME}.log
	wc -l *.mat.*.ids >> 4-flip-${NAME}.log
	
	echo -e "$(zcat intersect.${FASTQ}.pat.flipread.fastq.gz | wc -l) intersect.${FASTQ}.pat.flipread.fastq.gz" >> 4-flip-${NAME}.log
	wc -l *.pat.*.ids >> 4-flip-${NAME}.log
	
	date >> 4-flip-${NAME}.log
	
	cd ..
fi






#############################################################
## 5 alignment2
#############################################################
if [[ ${FLAG} -eq 5 || ${FLAG} -eq 0 ]] ; then
	
	## set up
	##_mkdir 5-alignment2-${NAME}
	if [ -d 5-alignment2-${NAME} ] ; then
		rm -r 5-alignment2-${NAME}
	fi
	mkdir 5-alignment2-${NAME}
	cd 5-alignment2-${NAME}
	
	## start log printing
	echo "#####################################" > 5-alignment2-${NAME}.log
  echo "## 5-alignment2 flipped reads #######" >> 5-alignment2-${NAME}.log
  echo "#####################################" >> 5-alignment2-${NAME}.log
  
  date >> 5-alignment2-${NAME}.log


	for i in ../4-flip-${NAME}/*.fastq.gz
	do
	ln -s $i
	done

	## align flipped reads
	# mat flip to mat and pat genomes
	zcat intersect.${FASTQ}.mat.flipread.fastq.gz | bowtie -p 20 --un intersect.${FASTQ}.matflip2pat.flipread.unaligned --max intersect.${FASTQ}.matflip2pat.flipread.multi --best --strata -v 2 -m 1 -q ${PGENOME_PATH}/AltRefFather/AltRefFather - > intersect.${FASTQ}.matflip2pat.flipread.bowtie 2> intersect.${FASTQ}.matflip2pat.flipread.log &  

	# pat flip to mat and pat genomes
	zcat intersect.${FASTQ}.pat.flipread.fastq.gz | bowtie -p 20 --un intersect.${FASTQ}.patflip2mat.flipread.unaligned --max intersect.${FASTQ}.patflip2mat.flipread.multi --best --strata -v 2 -m 1 -q ${PGENOME_PATH}/AltRefMother/AltRefMother - > intersect.${FASTQ}.patflip2mat.flipread.bowtie 2> intersect.${FASTQ}.patflip2mat.flipread.log &


	# mat flip to mat and pat genomes
	#zcat intersect.${FASTQ}.mat.flipread.fastq.gz | bowtie --un intersect.${FASTQ}.matflip2pat.flipread.unaligned --max intersect.${FASTQ}.matflip2pat.flipread.multi --best --strata -v 2 -m 1 -q ${PGENOME_PATH}/${NAME}_paternal - > intersect.${FASTQ}.matflip2pat.flipread.bowtie 2> intersect.${FASTQ}.matflip2pat.flipread.log &  

	# pat flip to mat and pat genomes
	#zcat intersect.${FASTQ}.pat.flipread.fastq.gz | bowtie --un intersect.${FASTQ}.patflip2mat.flipread.unaligned --max intersect.${FASTQ}.patflip2mat.flipread.multi --best --strata -v 2 -m 1 -q ${PGENOME_PATH}/${NAME}_maternal - > intersect.${FASTQ}.patflip2mat.flipread.bowtie 2> intersect.${FASTQ}.patflip2mat.flipread.log &


	wait
	
	## make chain files
	ln -s ${PGENOME_PATH}/mat2ref.chain
	ln -s ${PGENOME_PATH}/pat2ref.chain
	
	## check if mat and pat locations of the same read match
	## map back to reference; mat and pat locations should match
	${PL}/alleledb_map.back.ref.wrapper_forMouse.sh intersect.${FASTQ}.matflip2pat.flipread.bowtie paternal PAT pat2ref.chain ${PL} &
	${PL}/alleledb_map.back.ref.wrapper_forMouse.sh intersect.${FASTQ}.patflip2mat.flipread.bowtie maternal MAT mat2ref.chain ${PL} &
	
	wait 
	
	## compare mat and pat reference coordinates
	join -t $'\t' <(sed 's/\#\*o\*\#/\t/g' intersect.${FASTQ}.matflip2pat.flipread.bowtie.paternal.map2ref.bed | awk '{OFS="\t"}{FS="\t"}{print $4,$1,$2,$3}' | sort ) \
	<(sed 's/\#\*o\*\#/\t/g' intersect.${FASTQ}.patflip2mat.flipread.bowtie.maternal.map2ref.bed | awk '{OFS="\t"}{FS="\t"}{print $4,$1,$2,$3}' | sort ) \
	| sort | uniq | awk '{OFS="\t"}{FS="\t"}{if($2!=$5 && (sqrt(($3-$6)^2)>10 || sqrt(($4-$7)^2)>10)){print $0}}' > matpat.remove.reads.location.not.matched.txt
	
	
	## remove the mapping intermediate files
	rm intersect.${FASTQ}.matflip2pat.flipread.bowtie.*_paternal.bed intersect.${FASTQ}.matflip2pat.flipread.bowtie.*_paternal.bowtie intersect.${FASTQ}.matflip2pat.flipread.bowtie.paternal.map2ref.bed intersect.${FASTQ}.matflip2pat.flipread.bowtie.paternal.unmap2ref.log
	rm intersect.${FASTQ}.patflip2mat.flipread.bowtie.*_maternal.bed intersect.${FASTQ}.patflip2mat.flipread.bowtie.*_maternal.bowtie intersect.${FASTQ}.patflip2mat.flipread.bowtie.maternal.map2ref.bed intersect.${FASTQ}.patflip2mat.flipread.bowtie.maternal.unmap2ref.log
	
	## postprocessing
	echo "Note that there are redundancies in read IDs of the all files from here on due to multiple flipped reads simulated from reads that overlap >1 hetSNVs; unless there are only reads with 1 SNV" >> 5-alignment2-${NAME}.log
	echo "Here, the matpat.remove.reads.~txt file is only sorted and unique by row, not by read IDs" >> 5-alignment2-${NAME}.log
	wc -l intersect.${FASTQ}.matflip2pat.flipread.* >> 5-alignment2-${NAME}.log
	wc -l intersect.${FASTQ}.patflip2mat.flipread.* >> 5-alignment2-${NAME}.log
	wc -l matpat.remove.reads.location.not.matched.txt >> 5-alignment2-${NAME}.log
	
	date >> 5-alignment2-${NAME}.log
	
	cd ..
fi





#############################################################
## 6 unaligned
#############################################################
if [[ ${FLAG} -eq 6 || ${FLAG} -eq 0 ]] ; then
	## set up
	##_mkdir 6-unaligned-${NAME}
	if [ -d 6-unaligned-${NAME} ] ; then
		rm -r 6-unaligned-${NAME}
	fi
	mkdir 6-unaligned-${NAME}
	cd 6-unaligned-${NAME}
	
	
	## start log printing
	echo "###########################################" > 6-unaligned-${NAME}.log
  echo "## 6-check unaligned flipped reads  #######" >> 6-unaligned-${NAME}.log
  echo "###########################################" >> 6-unaligned-${NAME}.log
  
  date >> 6-unaligned-${NAME}.log

	ln -s ../5-alignment2-${NAME}/intersect.${FASTQ}.matflip2pat.flipread.unaligned 
	ln -s ../5-alignment2-${NAME}/intersect.${FASTQ}.patflip2mat.flipread.unaligned 
	
	## original mat and pat
	ln -s ../3-intersectBed-${NAME}/intersect.${FASTQ}.mat.snp.calls.txt
	ln -s ../3-intersectBed-${NAME}/intersect.${FASTQ}.pat.snp.calls.txt

	## preprocess for unaligned analyses	
	echo "Note that there are redundancies in read IDs of the all files from here on due to multiple flipped reads simulated from reads that overlap >1 hetSNVs; unless there are only reads with 1 SNV"
	${PL}/alleledb_multis-and-unaligneds-wrapper.sh intersect.${FASTQ}.matflip2pat.flipread.unaligned intersect.${FASTQ}.mat.snp.calls.txt mat ${PL}
	${PL}/alleledb_multis-and-unaligneds-wrapper.sh intersect.${FASTQ}.patflip2mat.flipread.unaligned intersect.${FASTQ}.pat.snp.calls.txt pat ${PL}
	
	date >> 6-unaligned-${NAME}.log
	
	cd ..
fi





#############################################################
## 7 multi
#############################################################
if [[ ${FLAG} -eq 7 || ${FLAG} -eq 0 ]] ; then
	## set up
	##_mkdir 7-multi-${NAME}
	if [ -d 7-multi-${NAME} ] ; then
		rm -r 7-multi-${NAME}
	fi
	mkdir 7-multi-${NAME}
	cd 7-multi-${NAME}
	
	## start log printing
	echo "#######################################" > 7-multi-${NAME}.log
  echo "## 7-check multi flipped reads  #######" >> 7-multi-${NAME}.log
  echo "#######################################" >> 7-multi-${NAME}.log
  
  date >> 7-multi-${NAME}.log

	ln -s ../5-alignment2-${NAME}/intersect.${FASTQ}.matflip2pat.flipread.multi
	ln -s ../5-alignment2-${NAME}/intersect.${FASTQ}.patflip2mat.flipread.multi

	## original mat and pat
	ln -s ../3-intersectBed-${NAME}/intersect.${FASTQ}.mat.snp.calls.txt
	ln -s ../3-intersectBed-${NAME}/intersect.${FASTQ}.pat.snp.calls.txt
	
	## preprocess for multimappers
	echo "Note that there are redundancies in read IDs of the all files from here on due to multiple flipped reads simulated from reads that overlap >1 hetSNVs; unless there are only reads with 1 SNV"
	${PL}/alleledb_multis-and-unaligneds-wrapper.sh intersect.${FASTQ}.matflip2pat.flipread.multi intersect.${FASTQ}.mat.snp.calls.txt mat ${PL}
  ${PL}/alleledb_multis-and-unaligneds-wrapper.sh intersect.${FASTQ}.patflip2mat.flipread.multi intersect.${FASTQ}.pat.snp.calls.txt pat ${PL}
	
	date >> 7-multi-${NAME}.log
	
	cd ..
fi

## get out of allelicbias folder
#cd ..





#############################################################
## 8 fsieve original multi reads from original fastq
## then run alleleseq again on this filtered fastqs
#############################################################
if [[ ${FLAG} -eq 8  ]] ; then
#if [[ ${FLAG} -eq 8 ]] ; then
	echo "################################################" > 8-final-run.log
	echo "## 8 rerun alleleseq on bias filtered reads ####" >> 8-final-run.log
	echo "################################################" >> 8-final-run.log
	date >> /workdir/sc2457/GM_GroSeq_AlleleDB_20170606/ENCODE_AlleleDB_log_20170615.log
	pwd >> /workdir/sc2457/GM_GroSeq_AlleleDB_20170606/ENCODE_AlleleDB_log_20170615.log
	${FASTQ} >> /workdir/sc2457/GM_GroSeq_AlleleDB_20170606/ENCODE_AlleleDB_log_20170615.log


	date >> 8-final-run.log

	## 1) cut for BED and then sort and uniq to obtain ids
	## 2a) removes reads from original aligned bowtie files in folder 1
	## 2b) Run AlleleSeq-betabinomial on filtered bowtie
	
	#ln -s ./allelicbias-${PGENOME}-${NAME}/1-alignment-${NAME}/${FASTQ}.mat.bowtie
	#ln -s ./allelicbias-${PGENOME}-${NAME}/1-alignment-${NAME}/${FASTQ}.pat.bowtie
	mv ./1-alignment-${NAME}/${FASTQ}.mat.bowtie .
	mv ./1-alignment-${NAME}/${FASTQ}.pat.bowtie .

	## 1) create ID list of all reads that multimap	
    if [ ! -f originalmatpatreads.toremove.ids ]; then
    echo "originalmatpatreads.toremove.ids not found!"
	echo "Note that all redundancies in read IDs that are to be removed are made unique here"
	cat <(cat ./7-multi-${NAME}/originalpatreads.intersect.${FASTQ}.patflip2mat.flipread.multi.bed ./7-multi-${NAME}/originalmatreads.intersect.${FASTQ}.matflip2pat.flipread.multi.bed | sed 's/\#\*o\*\#/\t/g' | cut -f4) \
	./4-flip-${NAME}/intersect.${FASTQ}.mat.snp.calls.removedreads.ids \
	./4-flip-${NAME}/intersect.${FASTQ}.mat.snp.calls.taggedreads.ids \
	./4-flip-${NAME}/intersect.${FASTQ}.mat.snp.calls.tooManySNVsReads.ids \
	./4-flip-${NAME}/intersect.${FASTQ}.pat.snp.calls.removedreads.ids \
	./4-flip-${NAME}/intersect.${FASTQ}.pat.snp.calls.taggedreads.ids \
	./4-flip-${NAME}/intersect.${FASTQ}.pat.snp.calls.tooManySNVsReads.ids \
	<(cut -f1 ./5-alignment2-${NAME}/matpat.remove.reads.location.not.matched.txt) \
	 | sort | uniq > originalmatpatreads.toremove.ids
	fi
##_	./allelicbias-${PGENOME}-${NAME}/4-flip-${NAME}/intersect.${FASTQ}.mat.snp.calls.toomanySNVsReads.ids \
##_	./allelicbias-${PGENOME}-${NAME}/4-flip-${NAME}/intersect.${FASTQ}.pat.snp.calls.toomanySNVsReads.ids \



	## 2) run alleleseq with betabinomial calculations
	echo "##############################" >> 8-final-run.log
	echo "### ALLELESEQ-BINOMIAL-RUN ###" >> 8-final-run.log
	echo "##############################" >> 8-final-run.log
	date >> 8-final-run.log
	gzip -d *.bowtie.gz
	make -f ${PIPELINEFILE} PREFIX=${FASTQ} >> 8-final-run.log 2>&1
        #make -f ${PIPELINEFILE} BASE=${PGENOME_PATH} PL=${PL} SNPS=${PGENOME_PATH}/snp.calls CNVS=${PGENOME_PATH}/cnv_rd_${NAME}/rd.${NAME}.txt MAPS=${PGENOME_PATH}/%s_${NAME}.map FDR_CUTOFF=${FDR_THRESH} PREFIX=${FASTQ} >> 8-final-run.log
	
#	${PL}/alleledb_alleleseqOutput2betabinomFormat.sh ${NAME} ${AS} ${PL} counts
#	${PL}/alleledb_alleleseqOutput2betabinomFormat.sh ${NAME} ${AS} ${PL} interestingHets 
	
	# betabinomial
#	R --vanilla --slave --args $(pwd) < ${PL}/alleledb_calcOverdispersion_1.R > ./calcOverdispersion.Rout 
#    R --vanilla --slave --args $(pwd) $FDR_THRESH < ${PL}/alleledb_alleleseqBetabinomial_1.R > alleleseqBetabinomial.Rout

	
#	${PL}/alleledb_alleleseqOutput2betabinomFormat.sh ${NAME} ${AS} ${PL} interestingHets.betabinom
#	${PL}/alleledb_alleleseqOutput2betabinomFormat.sh ${NAME} ${AS} ${PL} counts.betabinom
	
	date >> 8-final-run.log
	mkdir ../../${FASTQ}_toremove
	mv [1234567]* ../../${FASTQ}_toremove/.
	#rm [1234567]* -r
	gzip *.bowtie&
	cp counts.txt ../${FASTQ}_counts.txt
	cp interestingHets.txt ../${FASTQ}_interestingHets.txt
fi


#############################################################
## 9 fsieve original multi reads from original fastq
## seperate the bowtie to plus and minus strand
## then run alleleseq again on this filtered fastqs
#############################################################
if [[ ${FLAG} -eq 9 || ${FLAG} -eq 0 ]] ; then
	echo "###################################################################" > 9-Strand-Specific-final-run.log
	echo "## 9 rerun alleleseq on bias filtered reads (Strand-Specific) ####" >> 9-Strand-Specific-final-run.log
	echo "##################################################################" >> 9-Strand-Specific-final-run.log
	
	date >> 9-Strand-Specific-final-run.log

	## 1) cut for BED and then sort and uniq to obtain ids
	## 2a) removes reads from original aligned bowtie files in folder 1
	## 2b) Run AlleleSeq-betabinomial on filtered bowtie with strand specificity
	
	#ln -s ./allelicbias-${PGENOME}-${NAME}/1-alignment-${NAME}/${FASTQ}.mat.bowtie
	#ln -s ./allelicbias-${PGENOME}-${NAME}/1-alignment-${NAME}/${FASTQ}.pat.bowtie
	ln -s ./1-alignment-${NAME}/${FASTQ}.mat.bowtie .
	ln -s ./1-alignment-${NAME}/${FASTQ}.pat.bowtie .

	## 1) create ID list of all reads that multimap	
    if [ ! -f originalmatpatreads.toremove.ids ]; then
    echo "originalmatpatreads.toremove.ids not found!"
	echo "Note that all redundancies in read IDs that are to be removed are made unique here"
	cat <(cat ./7-multi-${NAME}/originalpatreads.intersect.${FASTQ}.patflip2mat.flipread.multi.bed ./7-multi-${NAME}/originalmatreads.intersect.${FASTQ}.matflip2pat.flipread.multi.bed | sed 's/\#\*o\*\#/\t/g' | cut -f4) \
	./4-flip-${NAME}/intersect.${FASTQ}.mat.snp.calls.removedreads.ids \
	./4-flip-${NAME}/intersect.${FASTQ}.mat.snp.calls.taggedreads.ids \
	./4-flip-${NAME}/intersect.${FASTQ}.mat.snp.calls.tooManySNVsReads.ids \
	./4-flip-${NAME}/intersect.${FASTQ}.pat.snp.calls.removedreads.ids \
	./4-flip-${NAME}/intersect.${FASTQ}.pat.snp.calls.taggedreads.ids \
	./4-flip-${NAME}/intersect.${FASTQ}.pat.snp.calls.tooManySNVsReads.ids \
	<(cut -f1 ./5-alignment2-${NAME}/matpat.remove.reads.location.not.matched.txt) \
	 | sort | uniq > originalmatpatreads.toremove.ids
	fi
##_	./allelicbias-${PGENOME}-${NAME}/4-flip-${NAME}/intersect.${FASTQ}.mat.snp.calls.toomanySNVsReads.ids \
##_	./allelicbias-${PGENOME}-${NAME}/4-flip-${NAME}/intersect.${FASTQ}.pat.snp.calls.toomanySNVsReads.ids \



	## 2) run alleleseq with betabinomial calculations
	echo "###############################################" >> 9-Strand-Specific-final-run.log
	echo "### ALLELESEQ-BINOMIAL-RUN (Strand-Specific)###" >> 9-Strand-Specific-final-run.log
	echo "###############################################" >> 9-Strand-Specific-final-run.log
	date >> 9-Strand-Specific-final-run.log
	make -f ${PIPELINEFILE} PREFIX=${FASTQ} >> 9-Strand-Specific-final-run.log
        #make -f ${PIPELINEFILE} BASE=${PGENOME_PATH} PL=${PL} SNPS=${PGENOME_PATH}/snp.calls CNVS=${PGENOME_PATH}/cnv_rd_${NAME}/rd.${NAME}.txt MAPS=${PGENOME_PATH}/%s_${NAME}.map FDR_CUTOFF=${FDR_THRESH} PREFIX=${FASTQ} >> 9-Strand-Specific-final-run.log
	
	#seperate into plus amd minus folder
	mkdir plus
	cp counts_plus.txt plus/counts.txt
	cp interestingHets_plus.txt plus/interestingHets.txt
	mkdir minus
	cp counts_minus.txt minus/counts.txt
	cp interestingHets_minus.txt minus/interestingHets.txt

    cd plus
	${PL}/alleledb_alleleseqOutput2betabinomFormat.sh ${NAME} ${AS} ${PL} counts
	${PL}/alleledb_alleleseqOutput2betabinomFormat.sh ${NAME} ${AS} ${PL} interestingHets 
	
	# betabinomial
	#echo "folder=\"$(pwd)/\"; setwd(folder)" | cat - ${PL}/alleledb_calcOverdispersion.R > calcOverdispersion.R 
	#R CMD BATCH ./calcOverdispersion.R  
	R --vanilla --slave --args $(pwd) < ${PL}/alleledb_calcOverdispersion_1.R > ./calcOverdispersion.Rout  
	#echo "folder=\"$(pwd)/\"; setwd(folder)" | cat - ${PL}/alleledb_alleleseqBetabinomial.R > alleleseqBetabinomial.R 
	#R CMD BATCH "--args FDR.thresh=$FDR_THRESH" ./alleleseqBetabinomial.R  
	R --vanilla --slave --args $(pwd) $FDR_THRESH < ${PL}/alleledb_alleleseqBetabinomial_1.R > alleleseqBetabinomial.Rout

	#${PL}/alleledb_alleleseqOutput2betabinomFormat.sh ${NAME} ${AS} ${PL} interestingHets.betabinom
	#${PL}/alleledb_alleleseqOutput2betabinomFormat.sh ${NAME} ${AS} ${PL} counts.betabinom
	
	cd ../minus
	${PL}/alleledb_alleleseqOutput2betabinomFormat.sh ${NAME} ${AS} ${PL} counts
	${PL}/alleledb_alleleseqOutput2betabinomFormat.sh ${NAME} ${AS} ${PL} interestingHets 
	
	# betabinomial
	R --vanilla --slave --args $(pwd) < ${PL}/alleledb_calcOverdispersion_1.R > ./calcOverdispersion.Rout  
	R --vanilla --slave --args $(pwd) $FDR_THRESH < ${PL}/alleledb_alleleseqBetabinomial_1.R > alleleseqBetabinomial.Rout

	
	#${PL}/alleledb_alleleseqOutput2betabinomFormat.sh ${NAME} ${AS} ${PL} interestingHets.betabinom
	#${PL}/alleledb_alleleseqOutput2betabinomFormat.sh ${NAME} ${AS} ${PL} counts.betabinom



	date >> 9-Strand-Specific-final-run.log
#	cd ..
#	mkdir toremove
#	mv [1234567]* *.bowtie toremove/.
fi


##10 only runs the 2nd half of 9
if [[ ${FLAG} -eq 10 ]] ; then
	
	date >> 9-Strand-Specific-final-run.log

	## 2) run alleleseq with betabinomial calculations
	echo "###############################################" >> 9-Strand-Specific-final-run.log
	echo "### ALLELESEQ-BINOMIAL-RUN (Strand-Specific)###" >> 9-Strand-Specific-final-run.log
	echo "###############################################" >> 9-Strand-Specific-final-run.log
	date >> 9-Strand-Specific-final-run.log
	make -f ${PIPELINEFILE} PREFIX=${FASTQ} >> 9-Strand-Specific-final-run.log
        #make -f ${PIPELINEFILE} BASE=${PGENOME_PATH} PL=${PL} SNPS=${PGENOME_PATH}/snp.calls CNVS=${PGENOME_PATH}/cnv_rd_${NAME}/rd.${NAME}.txt MAPS=${PGENOME_PATH}/%s_${NAME}.map FDR_CUTOFF=${FDR_THRESH} PREFIX=${FASTQ} >> 9-Strand-Specific-final-run.log
	
	#seperate into plus amd minus folder
	mkdir plus
	cp counts_plus.txt plus/counts.txt
	cp interestingHets_plus.txt plus/interestingHets.txt
	mkdir minus
	cp counts_minus.txt minus/counts.txt
	cp interestingHets_minus.txt minus/interestingHets.txt

    cd plus
	${PL}/alleledb_alleleseqOutput2betabinomFormat.sh ${NAME} ${AS} ${PL} counts
	${PL}/alleledb_alleleseqOutput2betabinomFormat.sh ${NAME} ${AS} ${PL} interestingHets 
	
	# betabinomial
	#echo "folder=\"$(pwd)/\"; setwd(folder)" | cat - ${PL}/alleledb_calcOverdispersion.R > calcOverdispersion.R 
	#R CMD BATCH ./calcOverdispersion.R  
	R --vanilla --slave --args $(pwd) < ${PL}/alleledb_calcOverdispersion_1.R > ./calcOverdispersion.Rout  
	#echo "folder=\"$(pwd)/\"; setwd(folder)" | cat - ${PL}/alleledb_alleleseqBetabinomial.R > alleleseqBetabinomial.R 
	#R CMD BATCH "--args FDR.thresh=$FDR_THRESH" ./alleleseqBetabinomial.R  
	R --vanilla --slave --args $(pwd) $FDR_THRESH < ${PL}/alleledb_alleleseqBetabinomial_1.R > alleleseqBetabinomial.Rout

	#${PL}/alleledb_alleleseqOutput2betabinomFormat.sh ${NAME} ${AS} ${PL} interestingHets.betabinom
	#${PL}/alleledb_alleleseqOutput2betabinomFormat.sh ${NAME} ${AS} ${PL} counts.betabinom
	
	cd ../minus
	${PL}/alleledb_alleleseqOutput2betabinomFormat.sh ${NAME} ${AS} ${PL} counts
	${PL}/alleledb_alleleseqOutput2betabinomFormat.sh ${NAME} ${AS} ${PL} interestingHets 
	
	# betabinomial
	R --vanilla --slave --args $(pwd) < ${PL}/alleledb_calcOverdispersion_1.R > ./calcOverdispersion.Rout  
	R --vanilla --slave --args $(pwd) $FDR_THRESH < ${PL}/alleledb_alleleseqBetabinomial_1.R > alleleseqBetabinomial.Rout

	
	#${PL}/alleledb_alleleseqOutput2betabinomFormat.sh ${NAME} ${AS} ${PL} interestingHets.betabinom
	#${PL}/alleledb_alleleseqOutput2betabinomFormat.sh ${NAME} ${AS} ${PL} counts.betabinom



	date >> 9-Strand-Specific-final-run.log
fi





