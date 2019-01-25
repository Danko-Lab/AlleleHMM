
import sys, string

TABLE=string.maketrans('ACGTacgt', 'TGCAtgca')

def reverseComplement(seq):
    tmp=seq[::-1]
    return tmp.translate(TABLE)

for l in sys.stdin:
    id, flag, chrm, pos, mapq, cigar, rnext, pnext, tlen, seq, qual = l.split('\t')[:11]
    flag=int(flag)
    if flag & 16:
        seq = reverseComplement(seq)
    print ">%s\n%s" % (id, seq)
