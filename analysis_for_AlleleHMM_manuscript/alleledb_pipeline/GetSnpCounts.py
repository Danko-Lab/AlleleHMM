
import pdb
import sys, bisect, scipy.stats, gc, os
from string import maketrans
import getNew1000GSNPAnnotations, InBindingSite, GetCNVAnnotations, Mapping2
import binom

MAXREADLEN=75

#tmp_trans={"paternal":0, "maternal":1}
# This is unfortunate; Alex's map files used REF, PAT, MAT, but everywhere else we use paternal, maternal.  
hap_trans={"paternal":"PAT", "maternal":"MAT"}

class chromRec(object):
    def __init__(self):
        self.reads=[]
        self.ends=[]

    def testOverlap(self, a1, a2, b):
        assert a1<a2
        return a1<=b and b<a2

    def getOverlappingReads(self, pos):
        readset = self.reads[bisect.bisect_left(self.ends, pos-MAXREADLEN):bisect.bisect_right(self.ends, pos+MAXREADLEN)]
        readset=[r for r in readset if self.testOverlap(int(r[1]), int(r[0]), pos)] # first part only approximate
        return readset

def parse_bowtie(f, mapper):
    d={}
    for l in open(f):
        mms=""
        vals=l.rstrip().split('\t')
        ident, strand, chrom, start, seq, qual, dummy = vals[0:7]
        chrom, hap = chrom.split('_')
        if len(vals)==8:
            mms=vals[7]

        nmm=len(mms.split(','))
        # end and start are hap coords!
        start=int(start)+1 # bowtie is 0-based.  
        end=start+len(seq)
        if not d.has_key(chrom):
            d[chrom]=chromRec()
        rec=d[chrom]
        # now we sort by end in hap coords.  Both hap and ref coords are included in the record.
        ref_start=mapper.trans(mapper.mt[hap_trans[hap]], mapper.mt["REF"], start)
        ref_end  =mapper.trans(mapper.mt[hap_trans[hap]], mapper.mt["REF"], end)
        if ref_start!=0 and ref_end!=0:
            rec.reads.append((ref_end, ref_start, end, start, nmm, seq, strand, ident, hap))

    for rec in d.itervalues():
        rec.reads.sort()
        rec.ends=[e[0] for e in rec.reads]

    return d

# FIX This is currently broken
def parse_eland(f):
    d={}
    for l in open(f):
        elems=l.rstrip().split()
        ident, seq, mapcode = elems[0:3]
        if mapcode in ['U0', 'U1', 'U2']:
            chrom, pos, strand = elems[6:9]
            pos=int(pos)
            end=pos+len(seq)
            if not d.has_key(chrom):
                d[chrom]=chromRec()
            rec=d[chrom]
            rec.reads.append((end, pos, mapcode, seq, strand, ident))

    for rec in d.itervalues():
        rec.reads.sort()
        rec.ends=[e[0] for e in rec.reads]

    return d

#NOTE: in the output P->Father, M->Mother, C->child

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

def testCounts(counts, snp):
    winningParent='?'
    snpchr, snppos, snprec = snp
    mat_genotype, pat_genotype, child_genotype, mat_allele, pat_allele, typ, ref, hetSNP = snprec

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
        print >>LOGFP,  "WARNING %s:%d failed thresh 1 %d %d" % (snpchr, snppos, allelecnts, total)
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
            print >>LOGFP,  "NOTE %s:%d Looks interesting: failed thresh 2 %d %d %f" % (snpchr, snppos, both, smaller, pval)
            print >>LOGFP,  "SNPS %s/%s, COUNTS a:%d c:%d g:%d t:%d" % (a1, a2, counts['a'], counts['c'], counts['g'], counts['t'])
            print >>LOGFP,  "Phasing P:%s M:%s D:%s" % (pat_allele, mat_allele, snprec)
            print >>LOGFP,  "\n"
            return (ASYMMETRIC, pval, a1, a2, counts, winningParent)
        else:
            return (SYMMETRIC, pval, a1, a2, counts, winningParent)
    else:
        return (HOMOZYGOUS, pval, a1, a2, counts, winningParent)


TABLE=maketrans('ACGTacgt', 'TGCAtgca')

def reverseComplement(seq):
    tmp=seq[::-1]
    return tmp.translate(TABLE)

def doChrom(chrm, readfile, maptmplt):

    hetSnps=0
    interestingSnps=0

    readlen=None

    g=_1000G.getAnnotationsGenerator(chrm) #'22'

    # FIX abs path
    mapper = Mapping2.Mapping(maptmplt % ('chr%s' % chrm))

    try:
        reads=parse_bowtie(readfile, mapper)
    except IOError:
        print >> sys.stderr, "Failed to open %s, skipping" % readfile
        return
    for snp in g: # skip non het?
        KSVals=[]
        counts={'a':0, 'c':0, 'g':0, 't':0, 'n':0}
        snpchrm, snppos, snprec = snp
        mat_genotype, pat_genotype, child_genotype, mat_allele, pat_allele, typ, ref, hetSNP = snprec

        if not hetSNP:
            print >>LOGFP,  "Position %s %d failed het test" % (chrm, snppos)
            continue # skip non-Het snps
        snppos=int(snppos)
        
        chrrec=reads['chr%s'%chrm] # 'chr22.fa'
        readset=chrrec.getOverlappingReads(snppos)
        if len(readset)<mindepth: 
            print >>LOGFP,  "Position %s %d failed depth test with %d" % (chrm, snppos, len(readset))
            continue
        passed=0

        for read in readset:
            end, start, hap_end, hap_start, nmm, seq, strand, ident, hap = read
            readlen=len(seq)
            start=int(start)
            end=int(end)
            hap_start=int(hap_start)
            hap_end=int(hap_end)

            hap_snppos=mapper.trans(mapper.mt["REF"], mapper.mt[hap_trans[hap]], snppos)
            if hap_snppos==0: # This hap occurs in a gap in this haplotype.
                continue
            # bowtie appears to hand this internally
            '''
            if strand=='-':
                seq=reverseComplement(seq)
                KSVals.append(end)
            else:
                KSVals.append(start)
                '''
            # changed to use % location of snppos within read
            p=float(hap_snppos-hap_start)/readlen
            assert(p>=0.0 and p<=1.0)
            KSVals.append(p)
            allele=seq[hap_snppos-hap_start]
            allele=allele.lower()
            counts[allele]=counts[allele]+1
            passed+=1

        if passed==0: continue # FIX not needed??
        np=1.0 
        if len(KSVals):
            nd, np = scipy.stats.kstest(KSVals, 'uniform', (0.0, 1.0))  
        t, pval, a1, a2, counts, winningParent = testCounts(counts, snp)
        if t==ASYMMETRIC or t==SYMMETRIC:
            hetSnps+=1
        if t==ASYMMETRIC:
            interestingSnps+=1
        if BShandler:
            inBS=1 if BShandler.check("chr%s"%chrm, snppos) else 0
        else:
            inBS=-1

        cnv=CNVhandler.getAnnotation("chr%s"%chrm, snppos)

        if cnv:
            cnv=cnv[2]
        else:
            cnv='1.0'

        print >>OUTFP, "%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%s\t%s\t%f\t%f\t%d\t%s\t%s" % (chrm, snppos, ref, mat_genotype, pat_genotype, child_genotype, typ, mat_allele, pat_allele, counts['a'], counts['c'], counts['g'], counts['t'], winningParent, t, pval, np, inBS, "DUMMY", cnv)
        print >>KSFP, "%s,%d,%d,%d,%s" % (chrm, snppos, len(KSVals), readlen, ",".join(map(str, KSVals)))
        OUTFP.flush()

chrms=[str(c) for c in range(1, 23)] + ['X'] # +['X', 'Y', 'M']
#chrms=['22']

USAGE="%s mindepth snpfile readfiletmplt maptmplt bindingsites cnvfile outfile logfile ksfile"

if __name__=='__main__':

    if len(sys.argv) != 10:
        print USAGE % sys.argv[0]
        sys.exit(-1)

    mindepth=int(sys.argv[1])
    snpfile=sys.argv[2]
    readfiletmplt=sys.argv[3]
    maptmplt=sys.argv[4]
    BindingSitefile=sys.argv[5]
    CNVFile=sys.argv[6]
    OUTFP = open(sys.argv[7], 'w')
    LOGFP = open(sys.argv[8], 'w')
    KSFP = open(sys.argv[9], 'w')

    gc.disable()

    print >>OUTFP, '\t'.join(('chrm', 'snppos   ', 'ref', 'mat_gtyp', 'pat_gtyp', 'c_gtyp', 'phase', 'mat_all', 'pat_all', 'cA', 'cC', 'cG', 'cT', 'winning', 'SymCls', 'SymPval',  'KSPval', 'BindingSite', 'SymQval', 'cnv'))
    _1000G=getNew1000GSNPAnnotations.Handler(snpfile, hasHeader=True)
    # absence of bs file signals that we don't do this test
    if os.access(BindingSitefile, os.R_OK):
        BShandler=InBindingSite.BSHandler(BindingSitefile)
    else:
        BShandler=None

    CNVhandler=GetCNVAnnotations.Handler(CNVFile)

    for chrm in chrms:
        try:
            readfile=readfiletmplt%chrm
        except TypeError:
            readfile=readfiletmplt
        doChrom(chrm, readfile, maptmplt)
