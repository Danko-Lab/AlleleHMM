

import bisect

class BSHandler(object):
    # simple test to make sure the sorted regions are non-overlapping
    def test(self, v):
        le=-1
        for s, e in v:
            assert(s > le)
            le = e
            
    def __init__(self, hitFile):
        '''construct regions from hits file'''
        self.d={}
        self.k={}
        for l in open(hitFile):
            e=l.rstrip().split('\t')
            c,s,e=e[0:3]
            self.d.setdefault(c,[]).append((int(s),int(e)))
        for c, v in self.d.iteritems():
            v=sorted(v)
            self.test(v)
            self.d[c]=v # a sorted list of pairs (s, e)
            self.k[c]=[e[1] for e in v] # a list of the end values for bisection

    def check(self, c, p):
        k=self.k[c]
        d=self.d[c]
        i=bisect.bisect(k, p)
        while i < len(k):
            s, e = d[i]
            if s <= p and p < e: return True
            if s > p: return False
            i+=1
        return False
