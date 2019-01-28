#!/bin/bash
## USAGE 2-map.back.ref.sh <bowtiefile> <maternal/paternal> <MAT/PAT> <chain filename> <pipeline path>
## e.g. 2-map.back.ref.sh wgEncodeBroadHistoneGm12878CtcfStdAlnRep1.fastq.gz.mat.bowtie maternal MAT ~/alleleseq-v1.2
## path has no '/'
## requires liftOver from UCSC
## removes intermediate files

BOWTIEFILE=$1
MATPATERNAL=$2
MATPAT=$3
CHAINFILE=$4
PL=$5


##### split chromosomes by 3rd col
echo "2-map.back.ref: splitting chromosomes by 3rd col for ${BOWTIEFILE}..."

for i in ${BOWTIEFILE} 
do 
awk '{OFS="\t"}{FS="\t"}{ print >> "'$i'." $3 ".bowtie" }' $i &
done

wait

#### sort by start position
echo "2-map.back.ref: sorting files..."
for i in *${MATPATERNAL}.bowtie
do 
${PL}/alleledb_bowtie2bed $i | sort -nk1,1 -nk2,2 -nk3,3 > ${i%bowtie}bed &
done

wait

#### combine and liftOver
echo "2-map.back.ref: combining and liftOver ..."
liftOver <(cat ${BOWTIEFILE}.[1-9]_${MATPATERNAL}.bed ${BOWTIEFILE}.1[0-9]_${MATPATERNAL}.bed ${BOWTIEFILE}.2[0-2]_${MATPATERNAL}.bed ${BOWTIEFILE}.[XY]_${MATPATERNAL}.bed) ${CHAINFILE} ${BOWTIEFILE}.${MATPATERNAL}.map2ref.bed ${BOWTIEFILE}.${MATPATERNAL}.unmap2ref.log
