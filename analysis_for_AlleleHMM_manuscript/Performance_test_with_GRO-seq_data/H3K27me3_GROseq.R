#Fig6
#setwd("~/Box Sync/Danko_lab_work/mouse_AlleleDB/GM12878/H3K27me3")
m_input_f="H3K27me3_minus_counts_hmm_regions_t1e-05.merged_cov_binomtest.bed"
p_input_f="H3K27me3_plus_counts_hmm_regions_t1e-05.merged_cov_binomtest.bed"

mStrand_count = read.table(m_input_f,header=T,sep = "\t",comment.char = "\\")
pStrand_count = read.table(p_input_f,header=T,sep = "\t",comment.char = "\\")

dummy = 0.1
mStrand_count$H_mat_pat_ratio = with(mStrand_count, log2((mat_allele_count+dummy)/(pat_allele_count+dummy)))
mStrandGROseq <- data.frame(do.call('rbind', strsplit(as.character(mStrand_count$GROseq_hmm_state),',',fixed=TRUE)))
names(mStrandGROseq) <- c("hmm_state","pat_count","mat_count")
mStrandGROseq$pat_count <- as.numeric(as.character(mStrandGROseq$pat_count))
mStrandGROseq$mat_count <- as.numeric(as.character(mStrandGROseq$mat_count))
mStrand_count$G_mat_pat_ratio = with(mStrandGROseq, log2((mat_count+dummy)/(pat_count+dummy)))
mStrand_count$H_mat_pat_sum = with(mStrand_count, mat_allele_count+pat_allele_count)
mStrand_count$G_mat_pat_sum = with(mStrandGROseq, mat_count+pat_count)

pStrand_count$H_mat_pat_ratio = with(pStrand_count, log2((mat_allele_count+dummy)/(pat_allele_count+dummy)))
pStrandGROseq <- data.frame(do.call('rbind', strsplit(as.character(pStrand_count$GROseq_hmm_state),',',fixed=TRUE)))
names(pStrandGROseq) <- c("hmm_state","pat_count","mat_count")
pStrandGROseq$pat_count <- as.numeric(as.character(pStrandGROseq$pat_count))
pStrandGROseq$mat_count <- as.numeric(as.character(pStrandGROseq$mat_count))
pStrand_count$G_mat_pat_ratio = with(pStrandGROseq, log2((mat_count+dummy)/(pat_count+dummy)))
pStrand_count$H_mat_pat_sum = with(pStrand_count, mat_allele_count+pat_allele_count)
pStrand_count$G_mat_pat_sum = with(pStrandGROseq, mat_count+pat_count)
fdr10_p=0.005  # GRO-seq all HMM regions

pdf("out.pdf")
par(mar=c(6.1, 7.1, 2.1, 2.1)) #d l u r 5.1, 4.1, 4.1, 2.1
par(mgp=c(3,1,0))
par(cex.lab=2.2, cex.axis=2.2)
with (pStrand_count, plot(G_mat_pat_ratio[Binom_p_value <= fdr10_p], H_mat_pat_ratio[Binom_p_value <= fdr10_p],
                          col="dark blue",
                          pch=20,
                          cex=2,
                          xlab="GRO-seq log2(mat/pat)",
                          ylab="H3K27me3 ChIP-seq log2(mat/pat)",
                          bty='n',
                          ylim = c(-10,10)
))
abline(v=0)
abline(h=0)
with (pStrand_count, points(G_mat_pat_ratio[Binom_p_value <= fdr10_p], H_mat_pat_ratio[Binom_p_value <= fdr10_p],
                            col="dark blue",
                            pch=20,
                            cex=2,
                            xlab="GRO-seq log2(mat/pat)",
                            ylab="H3K27me3 ChIP-seq log2(mat/pat)",
                            bty='n'
))
with (mStrand_count, points(G_mat_pat_ratio[Binom_p_value <= fdr10_p], H_mat_pat_ratio[Binom_p_value <= fdr10_p],
                            col="dark blue",
                            pch=20,
                            cex=2
))

newdatap=with(pStrand_count, pStrand_count[Binom_p_value <= fdr10_p,])
newdatam=with(mStrand_count, mStrand_count[Binom_p_value <= fdr10_p,])
newdata=rbind(newdatam, newdatap)
y=newdata$H_mat_pat_ratio
x=newdata$G_mat_pat_ratio

### TLS total least square
# https://stats.stackexchange.com/questions/13152/how-to-perform-orthogonal-regression-total-least-squares-via-pca
#https://stat.ethz.ch/pipermail/r-help/2016-March/437386.html

tls <- function(X,y){
  
  v <- prcomp(cbind(X,y))$rotation
  beta <- -v[-ncol(v),ncol(v)] / v[ncol(v),ncol(v)]
  return(beta)
  
}

tls_beta0 <- function(X,y, beta){
  beta0 <- mean(y)-sum(colMeans(X)*beta)
  return(beta0)
  
}

beta = tls(matrix(x),y)
beta0 = tls_beta0(matrix(x), y, beta)
a=-15:15
lines(a, beta*a+beta0 , col="blue", lwd=3)

with(newdata, cor.test(H_mat_pat_ratio, G_mat_pat_ratio,method = "pearson"))

dev.off()