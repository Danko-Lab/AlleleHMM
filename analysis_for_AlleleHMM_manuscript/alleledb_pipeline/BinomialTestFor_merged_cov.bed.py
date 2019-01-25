# output of BinomialTestFor_merged_cov.bed.py: if p-value <= 0.05, remain what it got from hmm (can ne M,P, or S), otherwise S.

import numpy as np
from math import *
import scipy.stats
from sys import argv

f_int = argv[1]
f_out = argv[2]

#f_int = "counts_minus_hmm_regions.merged_cov.bed"
#f_out = "counts_minus_hmm_regions.merged_cov_binomtest.bed"

mat  = np.loadtxt(f_int, dtype=int ,delimiter='\t', usecols=[4], skiprows=0)
pat = np.loadtxt(f_int, dtype=int ,delimiter='\t', usecols=[5], skiprows=0)
data = np.loadtxt(f_int, dtype=str ,delimiter='\t', usecols=range(0,7), skiprows=0)
hmm_states = data[:,3]
total = mat + pat


binomtest_dic={}
def binomtest(x, n, p):
    if (x,n,p) not in binomtest_dic:
        binomtest_dic[(x,n,p)] = scipy.stats.binom_test(x, n, p)
    return binomtest_dic[(x,n,p)]

p_value_list=[]
binomtest_state=[]
for i in xrange(len(mat)):
    p_value = binomtest(mat[i], total[i], 0.5)
    p_value_list.append(p_value)
    if p_value <= 0.05:
        binomtest_state.append(hmm_states[i])
    else:
        binomtest_state.append('S')
        

with open(f_out, 'w') as out:
    out.write("\t".join(['#chrm','chrmStart', 'chrmEnd', 'hmm_state','hmm+BinomialTest','mat_allele_count','pat_allele_count','identical_reads_count','Binom_p_value']))
    out.write("\n")
    for i in xrange(len(p_value_list)):
        out.write("\t".join(list(data[i,0:4])+[binomtest_state[i]]+list(data[i,4:7])+[str(p_value_list[i])]))
        out.write("\n")

