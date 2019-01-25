
'''
This script creates a SimpleQueue job list from the command line args

Each line will look something like:
cd /home2/rdb9/working/Joel/080609/DataSets/Fos; ./BOWTIE.sh s_1 AltRefFather; touch DONE.1

This file is an imput to sqPBS.py

'''

import os, sys, re

pat=re.compile('(\w+)_filtered_query.txt$')

tmplt = "cd %s; %s/BOWTIE.sh %s %s; touch DONE.%d"

cnt=0
dir=os.getcwd()

PL=sys.argv[1]

for f in sys.argv[2:]:
    # strip off _filtered_query.txt
    base=pat.match(f).group(1)
    for ref in ['AltRefFather', 'AltRefMother']:
        print tmplt % (dir, PL, base, ref, cnt)
        cnt+=1


