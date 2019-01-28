
import gc, os, sys, string, re, pdb, scipy.stats, cPickle
import Mapping2, getNewer1000GSNPAnnotations, Bowtie, binom, GetCNVAnnotations, dictmerge, utils, InBindingSite

TABLE=string.maketrans('ACGTacgt', 'TGCAtgca')
USAGE="%s mindepth snpfile readfiletmplt maptmplt bindingsites cnvfile outfile logfile ksfile"

def reverseComplement(seq):
    tmp=seq[::-1]
    return tmp.translate(TABLE)

def makeMappers(maptmplt):
    mappers={}
    cs=['chr%s' % str(c) for c in range(1,23)] + ['chrX', 'chrY', 'chrM']
 
    for c in cs:
        f=maptmplt % c
        if os.path.exists(f):
            mappers[c] = Mapping2.Mapping(f)
    return mappers

THRESH1=0.90
THRESH2=0.05

SYMMETRIC="Sym"
ASYMMETRIC="Asym"
HOMOZYGOUS="Homo"
WEIRD="Weird"

tbl={
'a':('a','a'),
'c':('c','c'),
'g':('g','g'),
't':('t','t'),
'r':('a','g'),
'y':('c','t'),
's':('c','g'),
'w':('a','t'),
'k':('g','t'),
'm':('a','c')
}

def convert(a):
    return tbl[a.lower()]

def testCounts(counts, chrm, snprec):
    winningParent='?'
    ref_pos, mat_genotype, pat_genotype, child_genotype, mat_allele, pat_allele, typ, ref, hetSNP = snprec

    # first, make sure that the expected alleles are the bulk of the counts
    total = counts['a']+counts['c']+counts['g']+counts['t']
    a1,a2=convert(child_genotype)
    if a1==a2:
        allelecnts = counts[a1]
    else:
        allelecnts = counts[a1]+counts[a2]
        
    both=counts[a1]+counts[a2]

    sortedCounts=sorted([(counts['a'], 'a'), (counts['c'],'c'), (counts['g'], 'g'), (counts['t'], 't')], reverse=True)
    majorAllele=sortedCounts[0][1]

    smaller=min(counts[a1], counts[a2])
    #pval=binomialDist.cdf(smaller, both, 0.5)*2 # This had problems for large sample sizes.  Switched to using scipy
    pval = binom.binomtest(smaller, both, 0.5) # scipy.binom_test was unstable for large counts
    
    if float(allelecnts)/total < THRESH1:
        print >>LOGFP,  "WARNING %s:%d failed thresh 1 %d %d" % (chrm, ref_pos, allelecnts, total)
        return (WEIRD, pval, a1, a2, counts, winningParent)

    # if the snp was phased
    if mat_allele and pat_allele:
        if mat_allele.lower()==majorAllele.lower():
            winningParent='M'
        elif pat_allele.lower()==majorAllele.lower():
            winningParent='P'
        else:
            winningParent='?'

    if a1!=a2:
        # we expect roughly 50/50.  
        if pval < THRESH2:
            print >>LOGFP,  "NOTE %s:%d Looks interesting: failed thresh 2 %d %d %f" % (chrm, ref_pos, both, smaller, pval)
            print >>LOGFP,  "SNPS %s/%s, COUNTS a:%d c:%d g:%d t:%d" % (a1, a2, counts['a'], counts['c'], counts['g'], counts['t'])
            print >>LOGFP,  "Phasing P:%s M:%s D:%s" % (pat_allele, mat_allele, snprec)
            print >>LOGFP,  "\n"
            return (ASYMMETRIC, pval, a1, a2, counts, winningParent)
        else:
            return (SYMMETRIC, pval, a1, a2, counts, winningParent)
    else:
        return (HOMOZYGOUS, pval, a1, a2, counts, winningParent)

def process(chrm, snppos, counts, snprec, CNVhandler):
    ref_pos, mat_genotype, pat_genotype, child_genotype, mat_allele, pat_allele, typ, ref, hetSNP = snprec
    t, pval, a1, a2, counts, winningParent = testCounts(counts, chrm, snprec)
    
    #if t==ASYMMETRIC or t==SYMMETRIC:
    #    hetSnps+=1
    #if t==ASYMMETRIC:
    #    interestingSnps+=1
    if BShandler:
        inBS=1 if BShandler.check("chr%s"%chrm, snppos) else 0
    else:
        inBS=-1

    cnv=CNVhandler.getAnnotation("chr%s"%chrm, snppos)
    if cnv:
        cnv=cnv[2]
    else:
        cnv='1.0'

    #nd, np = scipy.stats.kstest(ksvals, 'uniform', (0.0, 1.0))  

    print >>OUTFP, utils.myFormat('\t', (chrm, snppos, ref, mat_genotype, pat_genotype, child_genotype, typ, mat_allele, pat_allele, counts['a'], counts['c'], counts['g'], counts['t'], winningParent, t, pval, inBS, cnv))
    OUTFP.flush()

# This is used to order the chromosomes 1,2,3,...,22,X,Y.  Tricky, eh?
def chrcmp(a, b):
    try:
        a=int(a)
    except:
        pass
    try:
        b=int(b)
    except:
        pass
    return cmp(a,b)

if __name__=='__main__':

    if len(sys.argv) < 7:
        print USAGE % sys.argv[0]
        sys.exit(-1)

    mindepth=int(sys.argv[1])
    snpfile=sys.argv[2]
    BindingSitefile=sys.argv[3]
    CNVFile=sys.argv[4]
    OUTFP = open(sys.argv[5], 'w')
    LOGFP = open(sys.argv[6], 'w')

    countfiles=sys.argv[7:]
    
    if os.access(BindingSitefile, os.R_OK):
        BShandler=InBindingSite.BSHandler(BindingSitefile)
    else:
        BShandler=None

    CNVhandler=GetCNVAnnotations.Handler(CNVFile)

    hetSnps=0
    interestingSnps=0

    gc.disable()

    pat=re.compile('chr(\w+)_([mp]aternal)')
    print >>OUTFP, utils.myFormat('\t', ['chrm', 'snppos   ', 'ref', 'mat_gtyp', 'pat_gtyp', 'c_gtyp', 'phase', 'mat_all', 'pat_all', 'cA', 'cC', 'cG', 'cT', 'winning', 'SymCls', 'SymPval', 'BindingSite', 'cnv'])

    ref_1000G=getNewer1000GSNPAnnotations.Handler(snpfile, None, 'PAT', hasHeader=True, onlyHets=True)

    counts={}
    for countfile in countfiles:
        temp=cPickle.load(open(countfile))
        dictmerge.accum(counts, temp, lambda : 0, lambda a, b: a+b)
    
    for chrm in sorted(counts.keys(), chrcmp):
        for pos in sorted(counts[chrm].keys()):
            total = sum(counts[chrm][pos].values())
            if total >= mindepth:
                process(chrm, pos, counts[chrm][pos], ref_1000G.getAnnotation(chrm, pos), CNVhandler)

