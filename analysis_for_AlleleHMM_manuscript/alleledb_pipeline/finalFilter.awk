# 16 is pval
# 18 is cnv level
# thresh is passed in
NR==1 { print $0 }
NR>1 && ($18 >= 0.5 && $18 <= 1.5 ) && $16 <= thresh { print $0 }