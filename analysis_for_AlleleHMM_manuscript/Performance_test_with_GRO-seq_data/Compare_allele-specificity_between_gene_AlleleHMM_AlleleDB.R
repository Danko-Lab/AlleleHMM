#Fig4 D,E, Sup_fig7

#Figure4D
pdf("fig4D.pdf")
par(mfrow=c(1,3))
# Simple Pie Chart
lbls <- c("Concordant", "Discordant", "Symmetric")
H <- c( 7957, 399, 2250)
A <- c( 5274, 512, 1993)
D <- c( 15926, 7681, 34356)
#par(mfrow=c(1,3))
pie(H, labels = lbls, main="H")
pie(A, labels = lbls, main="A")
pie(D, labels = lbls, main="D")
dev.off()

#Sup_fig7
## read counts of SNPs in H (AlleleHMM, not AlleleDB), A (AlleleHMM and AlleleDB), and D (AlleleDB, not AlleleHMM)
pdf("sup_fig7.pdf")
H_Con=read.table("H_Concordant_counts.txt")
H_Dis=read.table("H_Discordant_counts.txt")
H_Sym=read.table("H_Symmetric_counts.txt")
A_Con=read.table("A_Concordant_counts.txt")
A_Dis=read.table("A_Discordant_counts.txt")
A_Sym=read.table("A_Symmetric_counts.txt")
D_Con=read.table("D_Concordant_counts.txt")
D_Dis=read.table("D_Discordant_counts.txt")
D_Sym=read.table("D_Symmetric_counts.txt")
H=c(H_Con$V1, H_Dis$V1, H_Sym$V1)
A=c(A_Con$V1, A_Dis$V1, A_Sym$V1)
D=c(D_Con$V1, D_Dis$V1, D_Sym$V1)

u=100
H[H>=u]=u
A[A>=u]=u
D[D>=u]=u

par(mfrow=c(3,1))

b=1
hist(H, freq=F, col="red", breaks =  seq(1,u+b,b), xlab="Read counts per SNP", ylab="fraction of SNPs" )
hist(A, freq=F, col="purple", breaks =  seq(1,u+b,b), xlab="Read counts per SNP", ylab="fraction of SNPs")
hist(D, freq=F, col="blue", breaks =  seq(1,u+b,b), xlab="Read counts per SNP", ylab="fraction of SNPs")

dev.off()

#Fig4E
pdf("Fig4E.pdf")
A=read.table("interestingHets_AlleleDB_in_AlleleHMM_MP_switch_counts.txt")
D=read.table("interestingHets_AlleleDB_out_AlleleHMM_MP_switch_counts.txt")
H=read.table("counts_noX_MinCount1_inAlleleHMM_t1e-05_interestingHets_outAlleleDB_switch_counts.txt")
xmax = max(c(A$V1, D$V1, H$V1))
par(mfrow=c(3,1))
par(cex.lab=2.2, cex.axis=2.2)
u=10
H$V1[H$V1>=u]=u
A$V1[A$V1>=u]=u
D$V1[D$V1>=u]=u
hist(H$V1-1, breaks = seq(-0.5,xmax,1), freq=F,col="red", xlim=c(0,10), ylim=c(0,1), main=NA, xlab="number of switches")
hist(A$V1-1, breaks = seq(-0.5,xmax,1), freq=F,col="purple", xlim=c(0,10), ylim=c(0,1), main=NA, xlab="number of switches") #,density=20,angle=180
hist(D$V1-1, breaks = seq(-0.5,xmax,1), freq=F,col="blue", xlim=c(0,10), ylim=c(0,1), main=NA, xlab="number of switches") #,density=20,angle=45
dev.off()


