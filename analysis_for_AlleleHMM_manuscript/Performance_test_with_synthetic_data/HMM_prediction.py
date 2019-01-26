import numpy as np
from math import *
import scipy.stats
from sys import argv
###structure of hmm

##intial prob
I_s=0.5
I_m=0.25
I_p=0.25

##transition
# 0,1,2 = M, S, P

# emmision
binomlogpmf_dic={}
def get_emission_log_prob(x,n,p):
    if (x,n,p) not in binomlogpmf_dic:
        binomlogpmf_dic[(x,n,p)] = binomlogpmf(x, n, p)
    return binomlogpmf_dic[(x,n,p)]


def binomlogpmf(k, n, p):
    return scipy.stats.binom.logpmf(k, n, p)


def get_max_argmax(l):
    m = max(l)
    a = l.index(m)
    return m, a


binomtest_dic={}
def binomtest(x, n, p):
    x = min (x, n-x)
    if (x,n,p) not in binomtest_dic:
        binomtest_dic[(x,n,p)] = scipy.stats.binom_test(x, n, p)
    return binomtest_dic[(x,n,p)]


def viterbi (p, T, x, n):
    v = np.full((3, len(x)), float('-inf'))
    b = np.full((3, len(x)), int(3), dtype=int)
    #Initialization
    v[0, 0] = log(I_m) + get_emission_log_prob(x[0],n[0],p[0])
    v[1, 0] = log(I_s) + get_emission_log_prob(x[0],n[0],p[1])
    v[2, 0] = log(I_p) + get_emission_log_prob(x[0],n[0],p[2])
    b[0 ,0] = 0
    b[1 ,0] = 1
    b[2 ,0] = 2
    #Iteration
    for i in xrange(1, len(x)):
        v[0, i], b[0, i] = get_max_argmax([v[0, i-1] + T[0,0],  v[1, i-1] + T[1,0], v[2, i-1] + T[2,0]])
        v[1, i], b[1, i] = get_max_argmax([v[0, i-1] + T[0,1],  v[1, i-1] + T[1,1], v[2, i-1] + T[2,1]])
        v[2, i], b[2, i] = get_max_argmax([v[0, i-1] + T[0,2],  v[1, i-1] + T[1,2], v[2, i-1] + T[2,2]])
        v[0, i] += get_emission_log_prob(x[i],n[i],p[0])
        v[1, i] += get_emission_log_prob(x[i],n[i],p[1])
        v[2, i] += get_emission_log_prob(x[i],n[i],p[2])
    #track back pointer
    viterbi_path_backward=[]
    i = len(x) - 1
    f = np.argmax(v[:,i])
    viterbi_path_backward.append(f)
    while (i > 0):
        f = b[f,i]
        viterbi_path_backward.append(f)
        i -= 1
    
    viterbi_path_forward = np.array(viterbi_path_backward[::-1])
    return viterbi_path_forward


### input data for viterbi
# seperate the autosome
# seperate plus and minus strand
#def hmm_prediction(f_v, strand, t,new_T, new_P):
import sys, getopt
f_v = "similated_Counts.txt"
t = 1e-05
try:
   opts, args = getopt.getopt(sys.argv[1:],"i:t:")
except getopt.GetoptError:
   print 'HMM_prediction.py -i <inputfile> -t <tuning parameter t>'
   sys.exit(2)

for opt, arg in opts:
    if opt in ("-i"):
        f_v = arg
    elif opt in ("-t"):
        t = float(arg)
#print 'Input file:\t', f_v
#print 't:\t', t


#if len(argv)==1:
#    f_v = "similated_Counts.txt"
#else:
#    try:
#        f_v=argv[1]
#    except:
#        f_v = "similated_Counts.txt"

data_v = np.loadtxt(f_v, dtype=str ,delimiter='\t', usecols=range(0,6), skiprows=1)
mat_v  = np.loadtxt(f_v, dtype=int ,delimiter='\t', usecols=[5], skiprows=1)
total_v = np.loadtxt(f_v, dtype=int ,delimiter='\t', usecols=[4], skiprows=1)
#t = 1e-05
t_mm, t_ms, t_mp = 1-t, t/2, t/2
#t_sm, t_ss, t_sp =  1.84668156e-03, 9.96232518e-01, 1.92080089e-03
t_sm, t_ss, t_sp = t/2, 1-t, t/2
t_pm, t_ps, t_pp = t/2, t/2, 1-t
T = np.log(np.array( [[t_mm, t_ms, t_mp],[t_sm, t_ss, t_sp],[t_pm, t_ps, t_pp]]))


p_m, p_s,p_p = 0.9, 0.5, 0.1
p = [p_m, p_s,p_p]
#p= [0.89955685995989998, 0.50242206595539318, 0.10013546869624733] 
v_path=viterbi (p, T,x=mat_v, n=total_v)


state_map={0:"M", 1:"S", 2:"P"}
with open(f_v[0:-4]+'_AlleleHMM.txt', 'w') as out:
    out.write(open(f_v).readlines()[0].strip())
    out.write('\thmm_state\n')
    for i in xrange(data_v.shape[0]):
        out.write('\t'.join(data_v[i]))
        out.write('\t'+state_map[v_path[i]])
        out.write('\n')
