#Sup_Fig2B
par(mar=c(6.1, 7.1, 2.1, 2.1)) #d l u r 5.1, 4.1, 4.1, 2.1
par(mgp=c(3,1,0))
par(cex.lab=2.2, cex.axis=2.2)
par(bty = 'n') 
par(mfrow=c(2,1)) 

setwd("~/Box Sync/Danko_lab_work/mouse_AlleleDB/GeneSize")
F1=read.table(file = "gencode.vM17.annotation_gene_SNPCount.txt")


hist(log10(F1$V1),xlab="log10(number of SNPs per gene)", cex=1.5,lwd=2, las=1, main="129/cast F1", col="gray",
     border = "white", xlim = c(0,5))
abline(v=log10(mean(F1$V1)), col="red", lwd=4, lty=2)
abline(v=log10(median(F1$V1)), col="blue", lwd=4, lty=2)
       

GM=read.table(file = "gencode.v28lift37.annotation_gene_SNPCount.txt")
hist(log10(GM$V1),xlab="log10(number of SNPs per gene)", cex=1.5,lwd=2, las=1, 
     xlim = c(0,5),
     main="GM12878", col="gray",border = "white")
abline(v=log10(mean(GM$V1)), col="red", lwd=4, lty=2)
abline(v=log10(median(GM$V1)), col="blue", lwd=4, lty=2)