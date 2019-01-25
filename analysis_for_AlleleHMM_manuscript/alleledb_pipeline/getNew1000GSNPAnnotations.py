

import getAnnotations, re

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

def homo(ltr):
    n=iupac2code[ltr]
    return n==1 or n==2 or n==4 or n==8

def hetero(ltr):
    return not homo(ltr)

def comp(a,b):
    return iupac2code[a]^iupac2code[b]

def mutant(f,m,c):
    '''check if child has any alleles not present in the parents, indicating
    a mutation or sequencing error'''
    return iupac2code[c]&~(iupac2code[f]|iupac2code[m])

class Handler(getAnnotations.annotationHandler):
    def process(self, l):

        e=l.rstrip().split()
        chrom=e[0]
        pos=int(e[1])
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
        '''
        my_m_allele, my_p_allele, my_typ = self.phase(mat_genotype, pat_genotype, child_genotype)
        assert (my_typ == typ)
        assert (my_m_allele==None or my_m_allele==m_allele)
        assert (my_p_allele==None or my_p_allele==p_allele)
        '''

        hetSNP=child_genotype in 'RYSWKM'

        return chrom, pos, (mat_genotype, pat_genotype, child_genotype, m_allele, p_allele, typ, ref, hetSNP)
                
    def phase(self, m, p, c):
        if mutant(p,m,c):
        #print "child mutant, skipping"
            return None, None, "MUTANT"
        if homo(c):
        #print "child homo, skipping"
            return None, None, "HOMO"
        if hetero(m) and hetero(p):
        #print "all hetero, skipping"
            return None, None, "HETERO"
        if homo(m):
        # since we're here, we know child and father are het, so 
            m_allele=m
            p_allele=code2iupac[comp(c, m_allele)]
        else:
            # since we're here, we know child and mother are het, so 
            p_allele=p
            m_allele=code2iupac[comp(c, p_allele)]

        return m_allele, p_allele, "PHASED"


    def toDict(self, v):
        fields=('mat_genotype', 'pat_genotype', 'child_genotype', 'mat_allele', 'pat_allele', 'typ', 'ref', 'hetSNP')
        assert len(v)==len(fields)
        d=dict(zip(fields, v))
        return d
