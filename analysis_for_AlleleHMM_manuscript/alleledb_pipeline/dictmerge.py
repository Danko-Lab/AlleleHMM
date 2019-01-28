

def accum(d, d1, initfunc, mergefunc):
    for k, v in d1.iteritems():
        if isinstance(v, dict):
            dd=d.setdefault(k, {})
            accum(dd, v, initfunc, mergefunc)
        else:
            d[k]=mergefunc(d.get(k, initfunc()), v)

if __name__=='__main__':
    d1={'chr1':{'3333': {'a':1, 'c':3, 'g':5, 't':7}, '5555': {'a':2, 'c':4, 'g':6, 't':8}, '7777': {'a':20, 'c':40, 'g':60, 't':80}}, }
    d2={'chr1':{'3333': {'a':1, 'c':3, 'g':5, 't':7}, '6666': {'a':2, 'c':4, 'g':6, 't':8}, '7777': {'a':20, 'c':40, 'g':60, 't':80}}, } 

    def initfunc():
        return 0

    def mergefunc(a, b):
        return a+b

    d={}
    accum(d, d1, initfunc, mergefunc)
    accum(d, d2, initfunc, mergefunc)
    print d
    
            
