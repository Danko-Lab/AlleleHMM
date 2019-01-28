
import sys

def parse(fn):
    if fn=='-':
        fp=sys.stdin
    else:
        fp=open(fn)
    for l in fp:
        yield l.split('\t')[0:7]

