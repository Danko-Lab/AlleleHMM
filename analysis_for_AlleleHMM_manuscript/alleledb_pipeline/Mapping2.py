
'''
Robert Bjornson, Yale University

This utility provides facilities to use .map files (Abyzov) to map between genomes efficiently.
Mapping objects are initialized using the .map file. Then, trans can be used to map between 
columns, or transFactory can be used to stamp out particular mapping functions.

The first row of the file are the column headers.

'''


import sys
import SortedCollection
from operator import itemgetter

class Mapping(object):
  
  def __init__(self, mapfile):

      fp=open(mapfile)
      header=fp.readline().rstrip()
      assert header[0]=='#'
      self.mt={}
      for i, name in enumerate(header[1:].split()):
        self.mt[name]=i 
      
      self.lists=[[],[],[]]
      tmp=[[],[],[]]

      i=0
      lastSig=[None, None, None]
      for l in fp:
          idxs=[int(v) for v in l.rstrip().split()]
          #sig=[idx==0 for idx in idxs]
          #if sig==lastSig:
          #    continue
          #lastSig=sig
          for v, l, t in zip(idxs, self.lists, tmp):
              l.append((v, i))
              if v != 0:
                  t.append((v, i))
          i+=1

      self.scs=[SortedCollection.SortedCollection(t, key=itemgetter(0)) for t in tmp]

  def transFactory(self, i, j):
      def tmp(p):
          idx=self.scs[i].find_le(p)[1]
          to_val=self.lists[j][idx][0]
          if to_val==0: 
              return 0
          else:
              from_val=self.lists[i][idx][0]
              delta=p-from_val
              return to_val+delta
      return tmp

  def trans(self, i, j, p):
      idx=self.scs[i].find_le(p)[1]
      to_val=self.lists[j][idx][0]
      if to_val==0: 
          return 0
      else:
          from_val=self.lists[i][idx][0]
          delta=p-from_val
          return to_val+delta

import random, time

# just for testing

if __name__=='__main__':
    m=Mapping(sys.argv[1])
    
    m2r=m.transFactory(0,2)
    r2m=m.transFactory(2,0)
    print "GO"
    t1=time.clock()
    for i in range(int(sys.argv[2])):
        fm=random.randint(1,10000000)
        to=m.trans(0,2,fm)
        ttoo=m2r(fm)
        assert(to==ttoo)
        if to!=0:
            cand_fm=m.trans(2,0,to)
            ccand_fm=r2m(to)
            assert(ccand_fm==cand_fm)
            if fm!=cand_fm:
                print fm, to, cand_fm
    
    print time.clock()-t1
    
