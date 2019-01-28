## USAGE alleleCounts2betabinomFormat.sh <indiv_category> <asb/ase> <pipeline PATH> <counts/interestingHets/interestingHets.betabinom> 
## e.g. alleleCounts2betabinomFormat.sh NA12878_Pol2 asb interestingHets.betabinom

NAME=$1
ASBE=$2
ASOUTPUT=$4
PL=$3

sed 1d ${ASOUTPUT}.txt | awk '{OFS="\t"}{FS="\t"}{if(($10+$11+$12+$13) >= 6){print "chr"$0}}' | awk '{OFS="\t"}{FS="\t"}{print $1,$2-1,$2,$0}' | cut -f1-3,6- > ${ASOUTPUT}.min6.bed

awk '{OFS="\t"}{FS="\t"}{if($1 != "chrX" && $1 != "chrY"){print $0}}' ${ASOUTPUT}.min6.bed > ${ASOUTPUT}.min6.auto.bed



if [[ ${ASOUTPUT} == "counts" ]] ; then

	awk '{OFS=","}{FS="\t"}{print "chr"$1"\t"$2"\t"$3"\t"$4,$5,$6,$11,$12,$13,$14,$11+$12+$13+$14,$17,name}' name=${NAME} ${ASOUTPUT}.min6.bed | sed 's/^chr//g' > ${ASOUTPUT}.min6.mod.bed

	${PL}/alleledb_alleleBed2ratio_${ASBE} -a ${ASOUTPUT}.min6.mod.bed -b ${ASOUTPUT}.min6.mod.bed > ${ASOUTPUT}.min6.allelicRatio.mod.bed

	paste ${ASOUTPUT}.min6.bed ${ASOUTPUT}.min6.allelicRatio.mod.bed | awk '{OFS="\t"}{FS="\t"}{print $1,$2,$3,$4,$25,$5,$6,$11,$12,$13,$14,$11+$12+$13+$14,$17,$24}' | awk '{OFS="\t"}{FS="\t"}{if($1 != "chrX" && $1 != "chrY"){print $0}}' > ${ASOUTPUT}.min6.allelicRatio.mod.auto.bed


	sed '1i\chr\tstart\tend\tref\tallelicRatio\tmat_gtyp\tpat_gtyp\tcA\tcC\tcG\tcT\ttotal\tpval\tindiv_cat' ${ASOUTPUT}.min6.allelicRatio.mod.auto.bed > ${ASOUTPUT}.min6.allelicRatio.mod.auto.txt

fi

