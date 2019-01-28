
import sys

def getSymCls(oldSymCls, qval, thresh):
    if not oldSymCls in ['Sym', 'Asym']:
        return oldSymCls
    return 'Asym' if qval<=thresh else 'Sym'

if __name__=='__main__':
    thresh=float(sys.argv[1])
    ifp = open(sys.argv[2])
    ofp = open(sys.argv[3], 'w')


    header=ifp.readline().rstrip()
    lines=[]
    vs=[]

    for i, l in enumerate(ifp):
        l=l.rstrip()
        lines.append(l)
        e=l.split('\t')
        pval=float(e[15])
        vs.append([pval, i]) # vs -> (pval, i)
    '''
    vs=[[0.002, 0],
        [0.012, 1],
        [0.060, 2],
        [0.020, 3],
        [0.030, 4],
        [0.040, 5]]
    '''

    vs.sort(reverse=True)
    n=len(vs)
    for pos, v in enumerate(vs):
        pval,i = v
        qval=pval*n/(n-pos)
        v.append(qval)

    vs2=[[v[1], v[0], v[2]] for v in vs] # vs2 -> (i, pval, qval)
    vs2.sort()

    print >>ofp, header
    for l, v in zip(lines, vs2):
        e=l.rstrip().split('\t')
        e[14]=getSymCls(e[14], v[2], thresh) # we reset the SymClas field here using the new thresh
        e[18] = "%f" % v[2]
        #e.append("%f"%v[2])
        print >>ofp, '\t'.join(e)

