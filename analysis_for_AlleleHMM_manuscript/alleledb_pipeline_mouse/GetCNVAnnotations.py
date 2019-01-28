
'''
This file parses a cnv file of the following form:

chrm	snppos	rd
1	52066	0.902113
1	695745	0.909802
1	742429	0.976435
1	742584	0.978998

'''

import getAnnotations

class Handler(getAnnotations.annotationHandler):
    def __init__(self, f=None, hasHeader=True):
        self.cnt=0
        super(Handler,self).__init__(f, hasHeader)

    def process(self, l):
        e=l.rstrip().split()
        chrom='chr'+e[0]
        pos=int(e[1])
        return chrom, pos, tuple(e)
                
    def toDict(self, v):
        fields=('chrom', 'pos', 'depth')
        if not v:
            v=[None for f in fields]
        assert len(v)==len(fields)
        d=dict(zip(fields, v))
        return d

