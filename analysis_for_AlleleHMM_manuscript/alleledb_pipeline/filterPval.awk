
# 16 is pval
# 19 is qval
# 20 is cnv level
# thresh must be passed in
NR==1 { print $0 }
NR>1 && $16 <= thresh { print $0 }