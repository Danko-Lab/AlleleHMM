#!/bin/env sh
# This script should take whatever input is provided and stream a fasta format suitable for bowtie

PL=$1
INFILE=$2

zcat ${INFILE} | python ${PL}/fastq2result.py - - |  python ${PL}/filter_query.py - - | python ${PL}/ConvertTags.py 