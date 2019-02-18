#python AlleleHMM.py prefix counts_plus_hmm.txt counts_minus_hmm.txt
import numpy as np
from math import *
import scipy.stats
#from sys import argv
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import time
import multiprocessing

import sys, getopt

counts_hmm="-"
counts_plus_hmm = "-"
counts_minus_hmm = "-"
predict=False
tao = "default"
prefix="AlleleHMM_output"
ITER=30 # number of iteration

def help_message():
    print "python AlleleHMM.py [options]"
    print ""
    print "options:"
    print ""
    print "To get help:"
    print "-h, --help             Show this brief help menu."
    print ""
    print "Required options:"
    print "For non-strand-specific data such as ChIP-seq:"
    print "-i, --input_hmm=PATH   Path to the non-strnad-specific, allele-specific read counts file (counts_hmm.txt)"
    print ""
    print "For strand-specific data such as PRO-seq:"
    print "-p, --input_plus_hmm=PATH    Path to the plus-strand allele-specific read counts file (counts_plus_hmm.txt)"
    print "-m, --input_minus_hmm=PATH   Path to the minus-strand allele-specific read counts file (counts_minus_hmm.txt)"
    print ""
    
    print "Optional operations:"
    print "-o, --output_prefix=STR      prefix for the output file. default=AlleleHMM_output"
    print "-t, --tao=FLOAT   AlleleHMM identify allele-specific blocks using 9 values of t (1E-01, 1E-02, ...,1E-09) by default."
    print "                  User can assign a specific tao for the calculation."
    print ""
    print "examples:"
    print "For strand-specific data such as PRO-seq use:"
    print "python AlleleHMM.py -p counts_plus_hmm.txt -m counts_minus_hmm.txt"
    
    print "For non-strand-specific data such as ChIP-seq use: "
    print "python AlleleHMM.py -i counts_hmm.txt"

#h_message='For strand-specific data such as PRO-seq use: \npython AlleleHMM.py -p counts_plus_hmm.txt -m counts_minus_hmm.txt \n\nFor non-strand-specific data such as ChIP-seq use: \npython AlleleHMM.py -i counts_hmm.txt'
#h_message='python AlleleHMM.py -h'

try:
   opts, args = getopt.getopt(sys.argv[1:],"ht:p:m:i:o:",["input_hmm=","input_plus_hmm=","input_minus_hmm=", "predict=", "tao=", "output_prefix="])
   if len(opts)== 0:
      help_message()
      sys.exit(2)
except getopt.GetoptError:
   help_message()
   sys.exit(2)

#print opts
i,p,m=0,0,0
for opt, arg in opts:
    if opt in ("-h"):
      help_message()
      sys.exit()
    elif opt in ("-p","--input_plus_hmm"):
      counts_plus_hmm = arg
      p+=1
    elif opt in ("-m","--input_minus_hmm"):
      counts_minus_hmm = arg
      assert counts_minus_hmm.count("plus") <1, 'please double check file name is correct'
      m+=1
    elif opt in ("-i","--input_hmm"):
      counts_hmm = arg
      i+=1
    elif opt in ("--predict="):
      predict = arg
    elif opt in ("-t", "--tao="):
      tao = float(arg)
    elif opt in ("-o", "--output_prefix="):
      prefix = arg

print 'Input non-strand-specific file :\t', counts_hmm
print 'Input plus file :\t', counts_plus_hmm
print 'Input minus file:\t', counts_minus_hmm
print 'tao:\t', tao

if (i+p+m) <1 or (i+p > 1) or (i+p+m >2) or (p != m):
   print '\nInput file error, please check the number of input files!\n'
   help_message()
   sys.exit()

if p==0 and i==1:
      counts_plus_hmm = counts_hmm


print 'output prefix:\t', prefix

input_i, input_p, input_m= i,p,m

#data
### input data for training
# comnined chromosomes
# combine plus strand and minus strand info for user
if counts_minus_hmm == "-":
    mat  = np.loadtxt(counts_plus_hmm, dtype=int ,delimiter='\t', usecols=[2], skiprows=1)
    pat  = np.loadtxt(counts_plus_hmm, dtype=int ,delimiter='\t', usecols=[3], skiprows=1)
    total = mat + pat
    #total = np.loadtxt(counts_plus_hmm, dtype=int ,delimiter='\t', usecols=[4], skiprows=1)
    #state = np.loadtxt(counts_plus_hmm, dtype=str ,delimiter='\t', usecols=[5], skiprows=1)
else:
    mat_1  = np.loadtxt(counts_plus_hmm, dtype=int ,delimiter='\t', usecols=[2], skiprows=1)
    pat_1 = np.loadtxt(counts_plus_hmm, dtype=int ,delimiter='\t', usecols=[3], skiprows=1)
    #state_1 = np.loadtxt(counts_plus_hmm, dtype=str ,delimiter='\t', usecols=[5], skiprows=1)
    mat_2  = np.loadtxt(counts_minus_hmm, dtype=int ,delimiter='\t', usecols=[2], skiprows=1)
    pat_2 = np.loadtxt(counts_minus_hmm, dtype=int ,delimiter='\t', usecols=[3], skiprows=1)
    #state_2 = np.loadtxt(counts_minus_hmm, dtype=str ,delimiter='\t', usecols=[5], skiprows=1)
    mat=np.concatenate((mat_1, mat_2))
    pat=np.concatenate((pat_1, pat_2))
    total = mat + pat
    #state=np.concatenate((state_1, state_2))

#n_state = np.full(len(state), int(3), dtype=int)
#n_state[state=="M"] = 0
#n_state[state=="S"] = 1
#n_state[state=="P"] = 2

###structure of hmm

##intial prob
I_s=0.5
I_m=0.25
I_p=0.25

##transition
# 0,1,2 = M, S, P
#t = 1e-05
#t_mm, t_ms, t_mp = 1-t, t/2, t/2
#t_sm, t_ss, t_sp = t/2, 1-t, t/2
#t_pm, t_ps, t_pp = t/2, t/2, 1-t
#T = np.log(np.array( [[t_mm, t_ms, t_mp],[t_sm, t_ss, t_sp],[t_pm, t_ps, t_pp]]))

## emmision
p_m, p_s,p_p = 0.7, 0.5, 0.3
p = [p_m, p_s,p_p]


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


def sumLogProb(a, b):
    # a and b are log probabilities
    # return log(exp(a)+exp(b))
    if b is None:
        return a
    elif a is None:
        return b
    elif a > b:
        return a + np.log1p(exp(b - a))
    else:
        return b + np.log1p(exp(a - b))



def forward_probability_calculation(x, n, p, T):
    f_p_m = np.full((3, len(x)), float('-inf'))
    # Initialization:
    f_p_m[0, 0] = log(I_m) + get_emission_log_prob(x[0],n[0],p[0])
    f_p_m[1, 0] = log(I_s) + get_emission_log_prob(x[0],n[0],p[1])
    f_p_m[2, 0] = log(I_p) + get_emission_log_prob(x[0],n[0],p[2])
    #Iteration
    for i in xrange(1, len(x)):
        f_p_m[0, i] = get_emission_log_prob(x[i],n[i],p[0]) + sumLogProb(sumLogProb(f_p_m[0, i-1] + T[0,0], f_p_m[1, i-1] + T[1,0]), f_p_m[2, i-1] + T[2,0])
        f_p_m[1, i] = get_emission_log_prob(x[i],n[i],p[1]) + sumLogProb(sumLogProb(f_p_m[0, i-1] + T[0,1], f_p_m[1, i-1] + T[1,1]), f_p_m[2, i-1] + T[2,1])
        f_p_m[2, i] = get_emission_log_prob(x[i],n[i],p[2]) + sumLogProb(sumLogProb(f_p_m[0, i-1] + T[0,2], f_p_m[1, i-1] + T[1,2]), f_p_m[2, i-1] + T[2,2])
    #Final value: 
    p_Y_f= sumLogProb(sumLogProb(f_p_m[0, len(x)-1], f_p_m[1,len(x)-1]),f_p_m[2,len(x)-1])
    return f_p_m, p_Y_f



def backward_probability_calculation(x, n, p, T):
    b_p_m = np.full((3, len(x)), float('-inf'))
    # Initialization:
    b_p_m[0, len(x)-1] = log(1) 
    b_p_m[1, len(x)-1] = log(1)
    b_p_m[2, len(x)-1] = log(1)
    #Iteration
    i = len(x) - 2
    while i >= 0:
        b_p_m[0, i] = sumLogProb(sumLogProb(T[0,0] + b_p_m[0, i+1] + get_emission_log_prob(x[i+1],n[i+1],p[0]), T[0,1] + b_p_m[1, i+1] + get_emission_log_prob(x[i+1],n[i+1],p[1])),T[0,2] + b_p_m[2, i+1] + get_emission_log_prob(x[i+1],n[i+1],p[2]))
        b_p_m[1, i] = sumLogProb(sumLogProb(T[1,0] + b_p_m[0, i+1] + get_emission_log_prob(x[i+1],n[i+1],p[0]), T[1,1] + b_p_m[1, i+1] + get_emission_log_prob(x[i+1],n[i+1],p[1])),T[1,2] + b_p_m[2, i+1] + get_emission_log_prob(x[i+1],n[i+1],p[2]))
        b_p_m[2, i] = sumLogProb(sumLogProb(T[2,0] + b_p_m[0, i+1] + get_emission_log_prob(x[i+1],n[i+1],p[0]), T[2,1] + b_p_m[1, i+1] + get_emission_log_prob(x[i+1],n[i+1],p[1])),T[2,2] + b_p_m[2, i+1] + get_emission_log_prob(x[i+1],n[i+1],p[2]))
        i -=1
    #Final value:
    p_Y_b = sumLogProb(sumLogProb(log(I_m) + b_p_m[0, 0] + get_emission_log_prob(x[0],n[0],p[0]), log(I_s) + b_p_m[1, 0] + get_emission_log_prob(x[0],n[0],p[1])),log(I_p)+ b_p_m[2, 0] + get_emission_log_prob(x[0],n[0],p[2]))
    return b_p_m, p_Y_b


def em_interate(T, p, x=mat, n=total):
    #t = time.time()
    f_p_m, p_Y_f = forward_probability_calculation(x, n, p, T)
    #print "forward: ", t- time.time()
    b_p_m, p_Y_b = backward_probability_calculation(x, n, p, T)
    #print "backward: ", t- time.time()
    
    #local P(Y)
    p_Y_l = np.full((1, len(x)), float('-inf'))
    for i in xrange(len(x)):
        p_Y_l[0,i] = sumLogProb(sumLogProb(b_p_m[0,i]+f_p_m[0,i], b_p_m[1,i]+f_p_m[1, i]),b_p_m[2, i]+f_p_m[2, i])
    #print "p_Y_l ", t- time.time()
    A = np.zeros((3,3))
    new_P = [None, None, None] #P_m, P_s, P_p 
    
    # can add multiple sequence
    for i in xrange(len(x)-1):
        for k in range(3):
            for l in range(3):
                A[k,l] = A[k,l] + exp(f_p_m[k, i] + T[k,l] + get_emission_log_prob(x[i+1],n[i+1],p[l]) + b_p_m[l, i+1] - p_Y_l[0,i])
    #print "A : ", t- time.time()
    new_T = np.zeros((3,3))
    for k in range(3):
        new_P[k] = np.sum(np.exp(f_p_m[k,] + b_p_m[k,]-p_Y_l) * x ) / np.sum(np.exp(f_p_m[k,] + b_p_m[k,]-p_Y_l) * n ) 
        new_T[k] = A[k]/A[k].sum()
    # adjust so that p_m + p_p = 1
    p0=new_P[0]
    p2=new_P[2]
    new_P[0]= (p0+1-p2)/2.0
    new_P[2]= (p2+1-p0)/2.0
    #print "secs: ", t- time.time()
    print new_T, new_P, p_Y_f
    return np.log(new_T), new_P, p_Y_f


def em_interate_T_mp_fixed(T, p, x=mat, n=total, update_state=1):
    #t = time.time()
    f_p_m, p_Y_f = forward_probability_calculation(x, n, p, T)
    #print "forward: ", t- time.time()
    b_p_m, p_Y_b = backward_probability_calculation(x, n, p, T)
    #print "backward: ", t- time.time()
    
    #local P(Y)
    p_Y_l = np.full((1, len(x)), float('-inf'))
    for i in xrange(len(x)):
        p_Y_l[0,i] = sumLogProb(sumLogProb(b_p_m[0,i]+f_p_m[0,i], b_p_m[1,i]+f_p_m[1, i]),b_p_m[2, i]+f_p_m[2, i])
    #print "p_Y_l ", t- time.time()
    A = np.zeros((3,3))
    new_P = [None, None, None] #P_m, P_s, P_p 
    
    # T of m and p fixed, but S to S,M,P update
    for i in xrange(len(x)-1):
        for l in range(3):
            k=update_state
            A[k,l] = A[k,l] + exp(f_p_m[k, i] + T[k,l] + get_emission_log_prob(x[i+1],n[i+1],p[l]) + b_p_m[l, i+1] - p_Y_l[0,i])
    #print "A : ", t- time.time()
    for k in range(3):
        new_P[k] = np.sum(np.exp(f_p_m[k,] + b_p_m[k,]-p_Y_l) * x ) / np.sum(np.exp(f_p_m[k,] + b_p_m[k,]-p_Y_l) * n ) 
    # adjust so that p_m + p_p = 1
    p0=new_P[0]
    p2=new_P[2]
    new_P[0]= (p0+1-p2)/2.0
    new_P[2]= (p2+1-p0)/2.0
    
    new_T = np.exp(T)
    k = update_state
    new_T[k] = A[k]/A[k].sum()
    #print "secs: ", t- time.time()
    #print new_T, new_P, p_Y_f
    return np.log(new_T), new_P, p_Y_f


def make_em_plot(em_p_Y_f_list, t, file_name='em_p_Y_f_list_plot.pdf', i=0):
    # i is the number of iteration
    plt.plot(xrange(i, len(em_p_Y_f_list)),em_p_Y_f_list[i:])
    plt.xlabel('# of iteration')
    plt.ylabel('log likelihood')
    plt.title(t)
    plt.savefig(file_name)
    plt.close()


### input data for viterbi
# seperate the autosome
# seperate plus and minus strand
f_data_dic={}
def hmm_prediction(f_v, strand, t,new_T, new_P):
    #f_v = "counts_plus_hmm.txt"
    if f_v in f_data_dic:
        data_v, chrom_v, snppos_v, mat_v, total_v = f_data_dic[f_v]
    else:
        data_v = np.loadtxt(f_v, dtype=str ,delimiter='\t', usecols=range(0,4), skiprows=1)
        chrom_v=data_v[:,0] # np.loadtxt(f_v, dtype=str ,delimiter='\t', usecols=[0], skiprows=1)
        snppos_v = data_v[:,1].astype(int) #np.loadtxt(f_v, dtype=int ,delimiter='\t', usecols=[1], skiprows=1)
        mat_v = data_v[:,2].astype(int)   #np.loadtxt(f_v, dtype=int ,delimiter='\t', usecols=[2], skiprows=1)
        pat_v = data_v[:,3].astype(int) 
        total_v = mat_v + pat_v  #np.loadtxt(f_v, dtype=int ,delimiter='\t', usecols=[4], skiprows=1)
        #state_v = data_v[:,5] #np.loadtxt(f_v, dtype=str ,delimiter='\t', usecols=[5], skiprows=1)
        #n_state_v = np.full(len(state_v), int(3), dtype=int)
        #n_state_v[state_v=="M"] = 0
        #n_state_v[state_v=="S"] = 1
        #n_state_v[state_v=="P"] = 2
        f_data_dic[f_v]=[data_v, chrom_v, snppos_v, mat_v, total_v]
    
    v_path=[]
    c_path=[]
    chrom_list= list(set(chrom_v))
    chrom_list.sort()
    
    for i in chrom_list:
        t_c = total_v[chrom_v == i]
        x_c = mat_v[chrom_v == i]
        v_path += (list(viterbi (x=x_c, n=t_c, p=new_P, T=new_T)))
        c_path.extend([i]*len(x_c))
    
    c_path=np.array(c_path)
    # output regions with neighbor sharing the same states as a bed file
    state_map = {0:'M', 1:'S', 2:'P'}
    region_list=[]
    for c in chrom_list:
        snppos_c = snppos_v[chrom_v == c]
        v_path_c = np.array(v_path)[c_path == c]
        u = snppos_c[0]
        for l in xrange(1,len(v_path_c)):
            if v_path_c[l] != v_path_c[l-1]:
                v = snppos_c[l-1]
                region_list.append([str(c), str(u-1), str(v), state_map[v_path_c[l-1]]])
                u = snppos_c[l]
        region_list.append([str(c), str(u-1), str(snppos_c[-1]), state_map[v_path_c[-1]]])
    
    #with open(f_v[0:-4]+'_regions_t'+str('%.0E' %t)+'.bed', 'w') as out:
    name="both"
    if strand== "+":
        name="plus"
    elif strand== "-":
        name="minus"
        
    with open(prefix+"_"+name+'_regions_t'+str('%.0E' %t)+'.bed', 'w') as out:
        for r in region_list:
            out.write('\t'.join(r+['111',strand]))
            out.write('\n')


### run em with Tmp fixed, Ts update
def run_em_T_mp_fixed(t, max_iter = ITER):
    t_mm, t_ms, t_mp = 1-t, t/2, t/2
    t_sm, t_ss, t_sp = t/2, 1-t, t/2
    t_pm, t_ps, t_pp = t/2, t/2, 1-t
    T = np.log(np.array( [[t_mm, t_ms, t_mp],[t_sm, t_ss, t_sp],[t_pm, t_ps, t_pp]]))
    new_T, new_P, p_Y_f = em_interate_T_mp_fixed(T, p, x=mat, n=total)
    p_Y_f_list = [p_Y_f]
    new_T_list = [new_T]
    new_P_list = [new_P]
    max_iter = max_iter
    for i in xrange(max_iter):
        print "iteration",i, new_T, new_P, p_Y_f
        new_T, new_P, p_Y_f = em_interate_T_mp_fixed(new_T, new_P, x=mat, n=total)
        p_Y_f_list.append(p_Y_f)
        new_T_list.append(new_T)
        new_P_list.append(new_P)
    #make_em_plot(p_Y_f_list,"count_min=1 Tmx, Tpx fixed, t="+str('%.0E' %t)+", Tsx allow change for EM", prefix+"_em_p_Y_f_list_plot_count_min=1_Tmpfixed_t="+str('%.0E' %t)+".pdf")
    with open(prefix+"_t="+str('%.0E' %t)+'_parameters.txt', 'w') as out:
        out.write("T="+str(new_T_list[-1])+"\n")
        out.write("P="+str(new_P_list[-1])+"\n")
    return [t, new_T_list,new_P_list, p_Y_f_list]


def run_all():
    t_list=[]
    for i in range(1,10):
        t_list.append(10**(-i))
    try:
        pool = multiprocessing.Pool(processes=9)
        pool_output = pool.map(run_em_T_mp_fixed, t_list )
        # pool_output looks like [[t, new_T_list,new_P_list, p_Y_f_list],...]
        pool.close() # no more tasks
        pool.join()
        
    except: # in case pool doesn't work
        result=[]
        for i in range(1,10):
            result.append(run_em_T_mp_fixed(10**(-i)))


def prediction(t):
    print prefix+"_t="+str('%.0E' %t)+'_parameters.txt'
    with open(prefix+"_t="+str('%.0E' %t)+'_parameters.txt') as p_in:
        l=p_in.readlines()
    #print i
    T=[]
    for ll in l[0:3]:
        for lll in ll.strip().strip('T=').strip('[').strip(']').split():
            #print lll
            T.append(float(lll))
        #print T
    new_T=np.array(T).reshape(3,3)
    new_P=[float(ll) for ll in l[-1].strip().strip('P=').strip('[').strip(']').split(",")]
    print "T", new_T
    print "P", new_P
    if input_i+input_m+input_p==1 : #counts_minus_hmm == "-":
        hmm_prediction(counts_plus_hmm, ".", t,new_T, new_P)
    else:
        hmm_prediction(counts_plus_hmm, "+", t,new_T, new_P)
        hmm_prediction(counts_minus_hmm, "-", t,new_T, new_P)

def prediction_all():
    for i in range(1,10):
        prediction(10**(-i))



if __name__ == '__main__':
    t = time.time()
    if tao != "default": # specific tao
        try:
            if predict !=False:
                prediction(tao)
            else:
                run_em_T_mp_fixed(tao)
                prediction(tao)
        except:
            run_em_T_mp_fixed(tao)
            prediction(tao)
    else: # all tao
        try:
            if predict !=False:
                prediction_all()
            else:
                run_all()
                prediction_all()
        except:
            run_all()
            prediction_all()
    print "AlleleHMM Run finished:", time.time() -t , "secs"