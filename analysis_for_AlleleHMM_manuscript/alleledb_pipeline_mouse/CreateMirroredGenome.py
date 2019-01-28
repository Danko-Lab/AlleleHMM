
import random
import getNew1000GSNPAnnotations
import singleFasta

splitTable={
'a':('a','a'),
'c':('c','c'),
'g':('g','g'),
't':('t','t'),
'r':('a','g'),
'y':('c','t'),
's':('c','g'),
'w':('a','t'),
'k':('g','t'),
'm':('a','c'),
'n':('n','n'),
}

def flip(p):
    if random.random() < 0.50:
        return (p[0].lower(), p[1].lower())
    else:
        return (p[1].lower(), p[0].lower())

def matchCase(case, base):
    if case.isupper():
        return base.upper()
    else:
        return base.lower()

def splitGenotype(annot):
    m_genotype, p_genotype, c_genotype, m_allele, p_allele, typ, ref, hetSNP = annot
    if m_allele:
        return m_allele.lower(), p_allele.lower()
    else:
        return flip(splitTable[c_genotype.lower()])

def handleChrom(snphandler, chrom, inputfile, m_output, p_output):
    rdr=singleFasta.fastaReader(inputfile)
    hdr=rdr.getHdr()
    mwriter=singleFasta.fastaWriter(m_output, hdr)
    fwriter=singleFasta.fastaWriter(p_output, hdr)
    
    for s, e, seq in rdr.generator():
        annots = snphandler.getAnnotationsAsDict(chrom, s, e)
        if annots:
            mbuf = list(seq)
            fbuf = list(seq)
            for p in xrange(s, e):
                if p in annots:
                    m_allele, p_allele = splitGenotype(annots[p])
                    if annots[p][7] and seq[p-s].lower() not in (m_allele, p_allele):
                        print "huh... pos %d, %s %s <> %s" % (p, m_allele, p_allele, seq[p-s])
                        
                    mbuf[p-s]=matchCase(seq[p-s], m_allele)
                    fbuf[p-s]=matchCase(seq[p-s], p_allele)
                    
            mstr="".join(mbuf)
            fstr="".join(fbuf)
        else:
            mstr=seq
            fstr=seq

        mwriter.writeLine(mstr)
        fwriter.writeLine(fstr)

if __name__=='__main__':

    chroms=[str(i) for i in range(1,23)]+['X','M'] # no Y
    snphandler=getNew1000GSNPAnnotations.Handler('../GM12878/Snps/snp.calls', hasHeader=True)
    for chrom in chroms:
        print 'doing %s' % chrom
        handleChrom(snphandler, chrom, '../GM12878/Ref/chr%s.fa'%chrom, '../GM12878/AltRefMother/chr%s_mother.fa'%chrom, '../GM12878/AltRefFather/chr%s_father.fa'%chrom)
    
