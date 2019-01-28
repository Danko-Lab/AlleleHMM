
'''
Handles this sort of file:

chrm	snppos   	ref	mat_gtyp	pat_gtyp	c_gtyp	phase	mat_all	pat_all	cA	cC	cG	cT	winning	SymCls	SymPval	BindingSite	cnv
1	934427	T	Y	Y	Y	HETERO	None	None	0	5	0	2	?	Sym	0.453125	-1	0.989249
1	939471	G	R	R	R	HETERO	None	None	9617	15	9295	3	?	Asym	0.0195857475458	-1	0.661849
1	949032	G	S	S	S	HETERO	None	None	0	1	7	0	?	Sym	0.0703125	-1	0.69965

'''

import getAnnotations

class Handler(getAnnotations.annotationHandler):
    def process(self, l):

        e=l.rstrip().split()
        chrom=e[0]
        pos=int(e[1])
        return chrom, pos, tuple(e[2:])
                
    def toDict(self, v):
        fields=('ref', 'm_genotype', 'f_genotype', 'c_genotype', 'phase', 'm_allele', 'f_allele', 'cA', 'cC', 'cG', 'cT', 'winning', 'SymClas', 'SymPval', 'BindingSite', 'cnv')
        if not v:
            v=[None for f in fields]
        assert len(v)==len(fields)
        d=dict(zip(fields, v))
        return d
