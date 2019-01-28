

'''

This script generates reads for testing the AllelePipeline.  Basically it will create
reads with the desired alleles.

We will start with a file that specifies some desired positions and allele counts:
chr1 8888 1:10:0:0
chr11 9999 10:5:0:0
...

we want to create fake reads that look like this:
@JOHNLENNON_0006:3:1:1605:1093#0/1
GGGGGTGTGACCCTTGAGGGCAGAAGNNNNNTCGAGAGGCCAGGCCTAACAGGGTTGGTAGCCGGAGTAAAGCCTG
+JOHNLENNON_0006:3:1:1605:1093#0/1
Yffdfccd\^dcdd^addaPZ_BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB

We will specify the desired read length.  The code will generate one read for each of the
allele counts, with the desired allele substituted into the correct position, and with the read starting at some
random location near the position.

'''

import Mapping2, fasta
import sys, os, random

def makeMappers(maptmplt):
    mappers={}
    cs=['chr%s' % str(c) for c in range(1,23)] + ['chrX', 'chrY', 'chrM']
 
    for c in cs:
        f=maptmplt % c
        if os.path.exists(f):
            mappers[c] = Mapping2.Mapping(f)
    return mappers

def parseAlleles(fp):
    ret=[]
    for l in fp:
        c, pos, alleles=l.rstrip().split()
        alleles=[int(e) for e in alleles.split(':')]
        ret.append([c, int(pos), alleles])
    return ret

if __name__=='__main__':
    allelefp=open(sys.argv[1])
    readlen=int(sys.argv[2])
    maptmplt=sys.argv[3]
    mappers=makeMappers(maptmplt)

    todo=parseAlleles(allelefp)
    fastacache={}

    idx=0
    for c, pos, alleles in todo:
        for allele, cnt in zip(['A', 'C', 'G', 'T'], alleles):
            
            for i in range(cnt):
                use=random.sample([1,2], 1)[0]
                startpos=pos-random.randint(0,readlen-1)
                hapstartpos=mappers[c].trans(0,use,startpos)
                hapendpos=hapstartpos+readlen
                happos=mappers[c].trans(0,use,pos)
                if happos<hapstartpos or happos>=hapendpos:
                    print "eh?"
                    continue
                
                fastafile="../GM12878/NA12878_diploid_genome_May3_2011/%s_NA12878_%s.fa" % (c, ['', 'paternal', 'maternal'][use])

                if fastafile not in fastacache:
                    fastacache[fastafile] = fasta.fastaReader(fastafile).next()[1].lower()

                tmp=[]
                for i in range(hapstartpos-1, hapendpos-1): # my seq is 0-based
                    if i == happos-1:
                        tmp.append(allele)
                    else:
                        tmp.append(fastacache[fastafile][i])
                
                tmp=''.join(tmp)
                tag="%d_%d_%d" % (idx, happos, hapstartpos)
                print "@%s\n%s\n+%s\n%s" % (tag, tmp, tag, "B"*len(tmp))
                idx+=1
    
