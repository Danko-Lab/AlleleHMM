## multis-and-unaligneds-wrapper.sh <1> <2> <original mat/pat> <pipeline path>
## find the original locations of multis and unaligned in their fastqs
## convert them to BED
## e.g.
## multis-and-unaligneds-wrapper.sh intersect.combined_pooled_NA12878_CTCF.fastq.gz.matflip2pat.flipread.unaligned intersect.combined_pooled_NA12878_CTCF.fastq.gz.mat.snps.calls.txt mat

FILE1=$1
FILE2=$2
MATPAT=$3
PL=$4

## add header to intersected.snp.txt
sed 's/\#\*o\*\#/\t/g' ${FILE2} | awk '{OFS="\t"}{FS="\t"}{print $4,$1,$2,$3,$5,$6,$7,$8,$9}' | sed '1i\id\tchr\tstart\tend\tsequence\tchrsnp\tstartsnp\tendsnp\tallele' > ${FILE2}.reads


## grep only read ids in unaligned file
#grep "^@" $FILE1 | sed 's/^@//g' | sed '1i\id' > ${FILE1}.reads
python ${PL}/fastq2result.py ${FILE1} - | grep "^>" | sed 's/^>//g' | sort | uniq | sed '1i\id' > ${FILE1}.reads 

## fsieve -- obtain where the original (unflipped) maternal reads are in the reference genome coordinates
${PL}/fsieve -s ${FILE1}.reads ${FILE2}.reads -o original${MATPAT}reads.${FILE1}

## rearrange and sort
sed 1d original${MATPAT}reads.${FILE1} | awk '{OFS="\t"}{FS="\t"}{print $2,$3,$4,$1"#*o*#"$5"#*o*#"$6"#*o*#"$7"#*o*#"$8"#*o*#"$9}' | ${PL}/sortByChr.sh - > original${MATPAT}reads.${FILE1}.bed

mkdir trash
mv original${MATPAT}reads.${FILE1} trash/


