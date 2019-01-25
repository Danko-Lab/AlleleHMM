# Sup_Fig10

par(mar=c(6.1, 7.1, 2.1, 2.1)) #d l u r 5.1, 4.1, 4.1, 2.1
par(mgp=c(3,1,0))
par(cex.lab=2.2, cex.axis=2.2)
par(bty = 'n') 
par(mfrow=c(2,1)) 

setwd("~/Box Sync/Danko_lab_work/mouse_AlleleDB/GM12878/SNPsInBlocks/")
F1_m=read.table(file = "SRR4041366_total_minus_regions_t1E-05.sorted_binomtest_SigBlocks_NoS.bedmappedSNP.txt", header = F)
F1_p=read.table(file = "SRR4041366_total_plus_regions_t1E-05.sorted_binomtest_SigBlocks_NoS.bedmappedSNP.txt", header = F)
F1=c(F1_m$V1, F1_p$V1)

hist(log10(F1),xlab="log10(number of SNPs per AlleleHMM block)", cex=1.5,lwd=2, las=1, main="129/cast F1", col="gray",
     border = "white", xlim = c(0,3))
abline(v=log10(mean(F1)), col="red", lwd=4, lty=2)
abline(v=log10(median(F1)), col="blue", lwd=4, lty=2)

summary(F1)


#F1_m=read.table(file = "SRR4041366_counts_minus_hmm_regions_t1e-05_interestingHets_mappedSNP.txt", header = F)
#F1_p=read.table(file = "SRR4041366_counts_plus_hmm_regions_t1e-05_interestingHets_mappedSNP.txt", header = F)
#F1=c(F1_m$V1, F1_p$V1)

#hist(log10(F1),xlab="log10(number of SNPs per AlleleHMM block)", cex=1.5,lwd=2, las=1, main="129/cast F1", col="gray",
#     border = "white", xlim = c(0,3))
#abline(v=log10(mean(F1)), col="red", lwd=4, lty=2)
#abline(v=log10(median(F1)), col="blue", lwd=4, lty=2)

#summary(F1)






GM_m=read.table("GM_counts_minus_hmm_regions_t1e-05_interestingHets_mappedSNP.txt")
GM_p=read.table("GM_counts_plus_hmm_regions_t1e-05_interestingHets_mappedSNP.txt")
GM=c(GM_m$V1, GM_p$V1)
hist(log10(GM),xlab="log10(number of SNPs per AlleleHMM block)", cex=1.5,lwd=2, las=1, 
     xlim = c(0,3),
     main="GM12878", col="gray",border = "white")
abline(v=log10(mean(GM)), col="red", lwd=4, lty=2)
abline(v=log10(median(GM)), col="blue", lwd=4, lty=2)
summary(GM)

#GM_m=read.table("SRR1552485_total_minus_regions_t1E-05.sorted_binomtest_SigBlocks_NoS_mappedSNP.txt")
#GM_p=read.table("SRR1552485_total_plus_regions_t1E-05.sorted_binomtest_SigBlocks_NoS_mappedSNP.txt")
#GM=c(GM_m$V1, GM_p$V1)
#hist(log10(GM),xlab="log10(number of SNPs per AlleleHMM block)", cex=1.5,lwd=2, las=1, 
#     xlim = c(0,3),
#     main="GM12878", col="gray",border = "white")
#abline(v=log10(mean(GM)), col="red", lwd=4, lty=2)
#abline(v=log10(median(GM)), col="blue", lwd=4, lty=2)
#summary(GM)
