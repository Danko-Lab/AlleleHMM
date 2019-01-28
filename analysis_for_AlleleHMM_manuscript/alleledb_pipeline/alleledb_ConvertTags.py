
import sys
import re

if __name__=='__main__':
    cnt=0
    while True:
        tag=sys.stdin.readline()
        seq=sys.stdin.readline()
        if not tag: break
        assert tag[0]=='>'
        newtag1='>%d:%s' % (cnt, tag[1:])
	newtag = re.sub('[ \t]+', '_', newtag1)
        print >> sys.stdout, newtag,
        print >> sys.stdout, seq,
        cnt+=1

        
