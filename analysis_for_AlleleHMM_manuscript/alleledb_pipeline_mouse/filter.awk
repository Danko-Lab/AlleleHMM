
# 16 is pval
# 18 is Binding site test
# 19 is qval
# 20 is cnv level
NR==1 { print $0 }
NR>1 && $20 >= 0.5 && $20 <= 1.5 { print $0 }