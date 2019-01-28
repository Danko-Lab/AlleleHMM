

import getAnnotations, re
import utils

iupac2code={
'A':1,
'C':2,
'G':4,
'T':8,
'R':1|4,
'Y':2|8,
'S':2|4,
'W':1|8,
'K':4|8,
'M':1|2
}

code2iupac={
1:'A',
2:'C',
4:'G',
8:'T',
1|4:'R',
2|8:'Y',
2|4:'S',
1|8:'W',
4|8:'K',
1|2:'M'
}

def all2code(s):
    a1, a2=s[0], s[1]
    return code2iupac[iupac2code[a1.upper()]|iupac2code[a2.upper()]]

class Handler(getAnnotations.annotationHandler):

    def __init__(self, snpfile, mappers, parent, hasHeader, onlyHets=False):
        self.mappers={}
        self.onlyHets=onlyHets
        if mappers:
            for chr, mapper in mappers.iteritems():
                self.mappers[chr]=mapper.transFactory(mapper.mt['REF'], mapper.mt[parent])
        else:
            self.mappers=None
        super(Handler, self).__init__(snpfile, hasHeader)
        
    def process(self, l):
        e=l.rstrip().split()
        chrom=e[0]
        ref_pos=int(e[1])
        if self.mappers:
            mapper=self.mappers[utils.formatChrom(chrom)]
            pos=mapper(ref_pos)
        else:
            pos=ref_pos
        ref=e[2]
        mat_genotype=all2code(e[3])
        pat_genotype=all2code(e[4])
        child_genotype=all2code(e[5])
        typ=e[6]

        if typ in ['HOMO', 'PHASED']:
            m_allele = e[5][0]
            p_allele = e[5][1]
        else:
            m_allele = p_allele = "None"

        # this is checking that the phasing from the file matches what I calculate

        hetSNP=child_genotype in 'RYSWKM'

        if self.onlyHets and not hetSNP: return (None, None, None) # This will signal to the superclass to skip this record

        return chrom, pos, (ref_pos, mat_genotype, pat_genotype, child_genotype, m_allele, p_allele, typ, ref, hetSNP)
                
    def toDict(self, v):
        fields=('mat_genotype', 'pat_genotype', 'child_genotype', 'mat_allele', 'pat_allele', 'typ', 'ref', 'hetSNP')
        assert len(v)==len(fields)
        d=dict(zip(fields, v))
        return d
