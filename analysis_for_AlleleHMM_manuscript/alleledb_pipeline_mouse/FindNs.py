
import sys

if __name__=='__main__':
    
    cnts={}
    for l in open(sys.argv[1]):
        seq=l.rstrip().split('\t')[8]
        try:
            firstN=seq.index('N')
            val=cnts.setdefault(firstN, 0)
            cnts[firstN]=val+1
        except ValueError:
            pass

    ks=sorted(cnts.keys())
    for k in ks:
        print k, cnts[k]


    
    
