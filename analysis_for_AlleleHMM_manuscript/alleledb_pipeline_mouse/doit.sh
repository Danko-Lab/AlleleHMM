
PL=/home/fas/hpcprog/rdb9/working/Jieming/100813/AlleleSpecificPipeline
$PL/filter_input.sh $PL $1 | bowtie --best --strata -v 2 -m 1 -f AltRefFather/AltRefFather > bt.out