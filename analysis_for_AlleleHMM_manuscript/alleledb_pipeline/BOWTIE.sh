
BASE=$1
REF=$2

bowtie --best --strata -v 2 -m 1 -f ${REF}/${REF} --un ${BASE}_${REF}_un.txt --max ${BASE}_${REF}_max.txt ${BASE}_filtered_query.txt ${BASE}_${REF}.txt

exit $?

