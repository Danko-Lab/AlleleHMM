
'''
This script converts a fastq file to an eland_query file
'''

import sys

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

    while True:
        header=ifile.readline().rstrip()
        if header=="": break
        seq=ifile.readline().rstrip()
        header2=ifile.readline().rstrip()
        qual=ifile.readline().rstrip()

        assert(header[0]=='@')
        assert(header2[0]=='+')
        
        print >> ofile, '>' + header[1:]
        print >> ofile, seq

