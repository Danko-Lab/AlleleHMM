
import sys

def allindices(s, c):
    return [i for i,v in enumerate(s) if v==c]

if __name__=='__main__':
    infn=sys.argv[1]; outfn=sys.argv[2]
    if infn=='-':
        ifile=sys.stdin
    else:
        ifile=open(sys.argv[1])
    if outfn=='-':
        ofile=sys.stdout
    else:
        ofile=open(sys.argv[2], 'w')

    total=0
    ns=0
    fp = sys.stdin
    ofp = sys.stdout
    while True:
        header=fp.readline()
        if not header:
            break
        seq=fp.readline()
        total+=1
        nindex=seq.find('N')
        if nindex==-1:
            ofile.write(header)
            ofile.write(seq)
            #print >>ofp, header,
            #print >>ofp, seq,
        else:
            ns+=1

    print >> sys.stderr, "total %d ns %d" % (total, ns)


    
    
