import os
import Mapping2

def myFormat(sep, lst):
    return sep.join([str(e) for e in lst])

def makeMappers(maptmplt):
    mappers={}
    req_cs=[str(c) for c in range(1,23)]
    opt_cs=['X', 'Y', 'M']

    for c in req_cs+opt_cs:
        f=maptmplt % c
        if os.path.exists(f):
            mappers[formatChrom(c)] = Mapping2.Mapping(f)
        else:
            if c in req_cs:
                raise Exception("Required map file %s not found" % f)
    return mappers

# This function will be passed just the chrom number, i.e. 1, 2, ... X, Y, M.
# it will also be passed an optional context
# it should return the correct format for that context
# This is used in several places: 
def formatChrom(num, context=None):
    #return "chr"+num
    return num

def parseChrom(c):
    # handle: 4_maternal
    chrm, parent = c.split('_')
    return (chrm, parent)
    # handle: chr4_maternal
    #chrm, parent = c[3:].split('_')
    
