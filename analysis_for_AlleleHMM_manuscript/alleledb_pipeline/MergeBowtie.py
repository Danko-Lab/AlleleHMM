

'''
This script takes two bowtie outputs, which were produces by searching exactly the same reads against
two slightly different genomes.  The bowtie parameters used should return a maximum of one location per
read, with a maximum of 2 errors (--best --strata -m 1 -v 2) Here are the rules:

1) If neither matched, skip

2) If one matched, but the other didn't, use the match

3) If both matched:

  3a) If one matched with fewer mismatches than the other, use it

  3b) If both had the same number of mismatches
  
     3ba) If both matched to precisely the same location (chr, strand, position) use that location

     3bb) otherwise, skip

The reads should be in the same order, but there will be gaps if one reference produced a hit and the other didn't
 
'''

import sys, random, pdb, os
import utils

def dump(fp, l):
    print >>fp, l,

def score(e):
    if len(e) != 8:
        return 0
    else:
        return -len(e[7].split(','))

def choose(l1, l2, mappers):
    global ident, file1better, file2better, skipped
    e1=l1.rstrip().split('\t')
    e2=l2.rstrip().split('\t')

    assert e1[0]==e2[0]
    s1=score(e1)
    s2=score(e2)

    c1=utils.formatChrom(e1[2].split('_')[0])
    c2=utils.formatChrom(e2[2].split('_')[0])

    if s1>s2: 
        file1better+=1
        return l1
    if s2>s1: 
        file2better+=1
        return l2
    if e1[1] != e2[1] or c1 != c2:
        skipped +=1 
        return None

    if mappers:
        
        m=mappers[c1]
        if m.trans(1,0,int(e1[3])+1) == m.trans(2,0,int(e2[3])+1): # bowtie is 0-based, but the maps are 1-based
            ident+=1
            return random.sample([l1, l2], 1)[0]
        else:
            skipped +=1
            return None
    else:
        if e1[3] == e2[3]:
            ident+=1
            return random.sample([l1, l2], 1)[0]
        else:
            skipped += 1
            return None

def getLine(fp):
    l=fp.readline()
    if l:
        e=l.rstrip().split('\t')
        idx=int(e[0].split(':')[0])
    else:
        e, idx = None, None
    return l, e, idx

def process(f1, f2, mappers, of, logf):

    global ident, file1, file2, file1better, file2better, skipped

    # counts of different pairings:
    ident=0 # case 3ba
    file1=0 # case 2
    file2=0 # case 2
    file1better=0 # case 3a
    file2better=0 # case 3b
    skipped=0 # case 3bb

    l1, e1, idx1 = getLine(f1)
    l2, e2, idx2 = getLine(f2)
    
    while l1 or l2:

        if idx1==idx2:
            use=choose(l1, l2, mappers)
            if use: dump(of, use)
            l1, e1, idx1 = getLine(f1)
            l2, e2, idx2 = getLine(f2)
            continue
        elif idx1!=None and (idx2==None or idx1<idx2):
            file1+=1
            dump(of, l1)
            l1, e1, idx1 = getLine(f1)
        elif idx2!=None and (idx1==None or idx1>idx2):
            file2+=1
            dump(of, l2)
            l2, e2, idx2 = getLine(f2)
        else:
            raise Exception("shouldn't get here")
          
    print >> logf, '''ident %(ident)d
file1 %(file1)d
file2 %(file2)d
file1better %(file1better)d
file2better %(file2better)d
skipped %(skipped)d''' % globals()

    # simple sanity check
    total=ident+file1+file2+file1better+file2better+skipped
    if ident < total * 0.8:
      print >> logf, "WARNING, ident seems low.  Please check!!"

if __name__=='__main__':
    f1=open(sys.argv[1])
    f2=open(sys.argv[2])
    maptmplt=sys.argv[3]
    mappers=utils.makeMappers(maptmplt)
    process(f1, f2, mappers, sys.stdout, sys.stderr)
