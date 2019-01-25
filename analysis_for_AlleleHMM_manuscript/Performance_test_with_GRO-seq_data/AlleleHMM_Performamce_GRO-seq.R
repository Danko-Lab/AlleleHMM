setwd("~/Box Sync/Danko_lab_work/AlleleHMM_performance_test")
par(mar=c(6.1, 7.1, 2.1, 2.1)) #d l u r 5.1, 4.1, 4.1, 2.1
par(mgp=c(3,1,0))
par(cex.lab=2.2, cex.axis=2.2)
par(bty = 'n') 

library("RColorBrewer")
red_cols = brewer.pal(n = 7, name = "Reds")
blue_cols = brewer.pal(n = 7, name = "Blues")


AlleleHMM_sensitivity_Specificity <- function (input_fp) {
  df=read.table(input_fp)
  df=df[df$V12+df$V13 >=10,]
  colnames(df)[7]="GeneAS"
  colnames(df)[9]="AlleleHMMAS"
  colnames(df)[8]="AlleleDBAS"
  diff=0.05
  AlleleHMM_result = data.frame(mat=seq(diff,1,diff), S_S=0, S_P=0, S_M=0, P_S=0, P_P=0, P_M=0, M_S=0, M_P=0, M_M=0)
  AS=c("S","P","M")
  for (m in 1:(dim(AlleleHMM_result)[1])){
    sub_df = df[df$V10 <= AlleleHMM_result$mat[m] & df$V10 >= (AlleleHMM_result$mat[m] - diff) ,]
    #V10 = mat/mat+pat
    for (i in seq(3)){
      for (j in seq(3)){
        AlleleHMM_result[m,3*(i-1)+j+1]=sum(sub_df$GeneAS==AS[i] & sub_df$AlleleHMMAS==AS[j])
      }
    }
  }
  AlleleHMM_result$mat = AlleleHMM_result$mat - diff/2
  AlleleHMM_result$sen = with(AlleleHMM_result, (P_P+M_M)/(P_S+P_P+M_M+M_S))
  AlleleHMM_result$spec = with(AlleleHMM_result, S_S/(S_S+S_P+S_M))
  AlleleHMM_result$prec = with(AlleleHMM_result, (P_P+M_M)/(P_P+M_M+S_M+S_P))
  return (list(AlleleHMM_result$mat, AlleleHMM_result$sen, AlleleHMM_result$spec, AlleleHMM_result$prec))
}

AlleleDB_sensitivity_Specificity <- function (input_fp) {
  df=read.table(input_fp)
  df=df[df$V12+df$V13 >=10,]
  colnames(df)[7]="GeneAS"
  colnames(df)[9]="AlleleHMMAS"
  colnames(df)[8]="AlleleDBAS"
  diff=0.05
  AlleleDB_result = data.frame(mat=seq(diff,1,diff), S_S=0, S_P=0, S_M=0, P_S=0, P_P=0, P_M=0, M_S=0, M_P=0, M_M=0)
  AS=c("S","P","M")
  for (m in 1:(dim(AlleleDB_result)[1])){
    sub_df = df[df$V10 <= AlleleDB_result$mat[m] & df$V10 >= (AlleleDB_result$mat[m] - diff) ,]
    
    for (i in seq(3)){
      for (j in seq(3)){
        AlleleDB_result[m,3*(i-1)+j+1]=sum(sub_df$GeneAS==AS[i] & sub_df$AlleleDBAS==AS[j])
      }
    }
  }
  AlleleDB_result$mat = AlleleDB_result$mat - diff/2
  AlleleDB_result$sen = with(AlleleDB_result, (P_P+M_M)/(P_S+P_P+M_M+M_S))
  AlleleDB_result$spec = with(AlleleDB_result, S_S/(S_S+S_P+S_M))
  AlleleDB_result$prec = with(AlleleDB_result, (P_P+M_M)/(P_P+M_M+S_M+S_P))
  return (list(AlleleDB_result$mat, AlleleDB_result$sen, AlleleDB_result$spec, AlleleDB_result$prec))
}

# Fig4 A,B
# Sup_Fig4 and Sup_Fig 6
## GM sensitivity
i=5+2
a=AlleleHMM_sensitivity_Specificity("SRR1552485_total_SNP_Allele_Specificity_matrix_t1E-05.bed")
plot(a[[1]], a[[2]], type='b', lty=i, ylab="Sensitivity", col=red_cols[i], ylim=c(0,1), xlab="Maternal reads fraction", cex=1.5,lwd=2, las=1)
b=AlleleDB_sensitivity_Specificity("SRR1552485_total_SNP_Allele_Specificity_matrix_t1E-05.bed")
lines(b[[1]], b[[2]], type='b', lty=i, ylab="Sensitivity", col=blue_cols[i], ylim=c(0,1), xlab="Maternal reads fraction", cex=1.5,lwd=2)

legend("top", pch=1, lty=i, col=c(red_cols[i], blue_cols[i]),
       legend=c("AlleleHMM", "AlleleDB"), cex=1.5,lwd=2, bty="n")

## GM precision
a=AlleleHMM_sensitivity_Specificity("SRR1552485_total_SNP_Allele_Specificity_matrix_t1E-05.bed")
plot(a[[1]], a[[4]], type='b', lty=i, ylab="Precision", col=red_cols[i], ylim=c(0,1), xlab="Maternal reads fraction", cex=1.5,lwd=2, las=1)
b=AlleleDB_sensitivity_Specificity("SRR1552485_total_SNP_Allele_Specificity_matrix_t1E-05.bed")
lines(b[[1]], b[[4]], type='b', lty=i, ylab="Precision", col=blue_cols[i], ylim=c(0,1), xlab="Maternal reads fraction", cex=1.5,lwd=2)

legend("bottomright", pch=1, lty=i, col=c(red_cols[i], blue_cols[i]),
       legend=c("AlleleHMM", "AlleleDB"), cex=1.5,lwd=2, bty="n")


# GM specificity
i=5+2
a=AlleleHMM_sensitivity_Specificity("SRR1552485_total_SNP_Allele_Specificity_matrix_t1E-05.bed")
plot(a[[1]], a[[3]], type='b', lty=i, ylab="Specificity", col=red_cols[i], ylim=c(0,1), xlab="Maternal reads fraction", cex=1.5,lwd=2,las=1)
a=AlleleDB_sensitivity_Specificity("SRR1552485_total_SNP_Allele_Specificity_matrix_t1E-05.bed")
lines(a[[1]], a[[3]], type='b', lty=i, col=blue_cols[i], ylim=c(0,1), xlab="Maternal reads fraction", cex=1.5,lwd=2)

legend("bottom", pch=1, lty=i, col=c(red_cols[i], blue_cols[i]),
       legend=c("AlleleHMM", "AlleleDB"), cex=1.5,lwd=2, bty="n")



## mouse sensitivity
i=5+2
a=AlleleHMM_sensitivity_Specificity("SNP_Allele_Specificity_matrix_mouse_2FDR.bed")
plot(a[[1]], a[[2]], type='b', lty=i, ylab="Sensitivity", col=red_cols[i], ylim=c(0,1), xlab="Maternal reads fraction", cex=1.5,lwd=2, las=1)
b=AlleleDB_sensitivity_Specificity("SNP_Allele_Specificity_matrix_mouse_2FDR.bed")
lines(b[[1]], b[[2]], type='b', lty=i, ylab="Sensitivity", col=blue_cols[i], ylim=c(0,1), xlab="Maternal reads fraction", cex=1.5,lwd=2)

legend("top", pch=1, lty=i, col=c(red_cols[i], blue_cols[i]),
       legend=c("AlleleHMM", "AlleleDB"), cex=1.5,lwd=2, bty="n")

## mouse precision
a=AlleleHMM_sensitivity_Specificity("SNP_Allele_Specificity_matrix_mouse_2FDR.bed")
plot(a[[1]], a[[4]], type='b', lty=i, ylab="Precision", col=red_cols[i], ylim=c(0,1), xlab="Maternal reads fraction", cex=1.5,lwd=2, las=1)
b=AlleleDB_sensitivity_Specificity("SNP_Allele_Specificity_matrix_mouse_2FDR.bed")
lines(b[[1]], b[[4]], type='b', lty=i, ylab="Precision", col=blue_cols[i], ylim=c(0,1), xlab="Maternal reads fraction", cex=1.5,lwd=2)

legend("bottomleft", pch=1, lty=i, col=c(red_cols[i], blue_cols[i]),
       legend=c("AlleleHMM", "AlleleDB"), cex=1.5,lwd=2, bty="n")

## mouse Specificity
i=5+2
a=AlleleHMM_sensitivity_Specificity("SNP_Allele_Specificity_matrix_mouse_2FDR.bed")
plot(a[[1]], a[[3]], type='b', lty=i, ylab="Specificity", col=red_cols[i], ylim=c(0,1), xlab="Maternal reads fraction", cex=1.5,lwd=2,las=1)
b=AlleleDB_sensitivity_Specificity("SNP_Allele_Specificity_matrix_mouse_2FDR.bed")
lines(b[[1]], b[[3]], type='b', lty=i, col=blue_cols[i], ylim=c(0,1), xlab="Maternal reads fraction", cex=1.5,lwd=2)

legend("bottom", pch=1, lty=i, col=c(red_cols[i], blue_cols[i]),
       legend=c("AlleleHMM", "AlleleDB"), cex=1.5,lwd=2, bty="n")


# GM
# AlleleHMM read depth vs sensitivity
display.brewer.pal(n = 5, name = "OrRd")
red_cols = brewer.pal(n = 7, name = "OrRd")
type='b'
i=1+2
a=AlleleHMM_sensitivity_Specificity("SRR1552485_sub16_SNP_Allele_Specificity_matrix_t1E-02.bed")
plot(a[[1]], a[[2]], type=type, pch=20, lty=1, ylab="Sensitivity", col=red_cols[i], ylim=c(0,1), xlab="Maternal reads fraction", cex=1.5,lwd=2, las=1)
i=2+2
a=AlleleHMM_sensitivity_Specificity("SRR1552485_sub8_SNP_Allele_Specificity_matrix_t1E-02.bed")
lines(a[[1]], a[[2]],type=type, pch=20, lty=1, ylab="Sensitivity", col=red_cols[i], cex=1.5,lwd=2)
i=3+2
a=AlleleHMM_sensitivity_Specificity("SRR1552485_sub4_SNP_Allele_Specificity_matrix_t1E-03.bed")
lines(a[[1]], a[[2]], type=type, pch=20, lty=1, ylab="Sensitivity", col=red_cols[i], cex=1.5,lwd=2)
i=4+2
a=AlleleHMM_sensitivity_Specificity("SRR1552485_sub2_SNP_Allele_Specificity_matrix_t1E-04.bed")
lines(a[[1]], a[[2]], type=type, pch=20, lty=1, ylab="Sensitivity", col=red_cols[i], cex=1.5,lwd=2)
i=5+2
a=AlleleHMM_sensitivity_Specificity("SRR1552485_total_SNP_Allele_Specificity_matrix_t1E-05.bed")
lines(a[[1]], a[[2]], type=type, pch=1, lty=1, ylab="Sensitivity", col=red_cols[i], ylim=c(0,1), xlab="Maternal reads fraction", cex=1.5,lwd=2)

legend("top", lty =1,col=rev(red_cols[3:7]),pch=c(1,20,20,20,20),
       legend=c("100%", "  50%", "  25%", "12.5%", "6.25%"), cex=1.5,lwd=2, bty="n")


# AlleleHMM read depth vs precision
display.brewer.pal(n = 5, name = "OrRd")
red_cols = brewer.pal(n = 7, name = "OrRd")
type='b'
i=1+2
a=AlleleHMM_sensitivity_Specificity("SRR1552485_sub16_SNP_Allele_Specificity_matrix_t1E-02.bed")
plot(a[[1]], a[[4]], type=type, pch=20, lty=1, ylab="Precision", col=red_cols[i], ylim=c(0,1), xlab="Maternal reads fraction", cex=1.5,lwd=2, las=1)
i=2+2
a=AlleleHMM_sensitivity_Specificity("SRR1552485_sub8_SNP_Allele_Specificity_matrix_t1E-02.bed")
lines(a[[1]], a[[4]],type=type, pch=20, lty=1, col=red_cols[i], cex=1.5,lwd=2)
i=3+2
a=AlleleHMM_sensitivity_Specificity("SRR1552485_sub4_SNP_Allele_Specificity_matrix_t1E-03.bed")
lines(a[[1]], a[[4]], type=type, pch=20, lty=1, col=red_cols[i], cex=1.5,lwd=2)
i=4+2
a=AlleleHMM_sensitivity_Specificity("SRR1552485_sub2_SNP_Allele_Specificity_matrix_t1E-04.bed")
lines(a[[1]], a[[4]], type=type, pch=20, lty=1, col=red_cols[i], cex=1.5,lwd=2)
i=5+2
a=AlleleHMM_sensitivity_Specificity("SRR1552485_total_SNP_Allele_Specificity_matrix_t1E-05.bed")
lines(a[[1]], a[[4]], type=type, pch=1, lty=1, col=red_cols[i], ylim=c(0,1), xlab="Maternal reads fraction", cex=1.5,lwd=2)

legend("bottomright", lty =1,col=rev(red_cols[3:7]),pch=c(1,20,20,20,20),
       legend=c("100%", "  50%", "  25%", "12.5%", "6.25%"), cex=1.5,lwd=2, bty="n")

# AlleleHMM read depth vs specificity
i=1+2
a=AlleleHMM_sensitivity_Specificity("SRR1552485_sub16_SNP_Allele_Specificity_matrix_t1E-02.bed")
plot(a[[1]], a[[3]], type=type, pch=20, lty=1, ylab="Specificity", col=red_cols[i], ylim=c(0,1), xlab="Maternal reads fraction", cex=1.5,lwd=2, las=1)
i=2+2
a=AlleleHMM_sensitivity_Specificity("SRR1552485_sub8_SNP_Allele_Specificity_matrix_t1E-02.bed")
lines(a[[1]], a[[3]], type=type, pch=20, lty=1, ylab="Sensitivity", col=red_cols[i], cex=1.5,lwd=2)
i=3+2
a=AlleleHMM_sensitivity_Specificity("SRR1552485_sub4_SNP_Allele_Specificity_matrix_t1E-03.bed")
lines(a[[1]], a[[3]], type=type, pch=20, lty=1, ylab="Sensitivity", col=red_cols[i], cex=1.5,lwd=2)
i=4+2
a=AlleleHMM_sensitivity_Specificity("SRR1552485_sub2_SNP_Allele_Specificity_matrix_t1E-04.bed")
lines(a[[1]], a[[3]], type=type, pch=20, lty=1, ylab="Sensitivity", col=red_cols[i], cex=1.5,lwd=2)
i=5+2
a=AlleleHMM_sensitivity_Specificity("SRR1552485_total_SNP_Allele_Specificity_matrix_t1E-05.bed")
lines(a[[1]], a[[3]], type=type, pch=1, lty=1, ylab="Sensitivity", col=red_cols[i], ylim=c(0,1), xlab="Maternal reads fraction", cex=1.5,lwd=2)

legend("bottom", lty =1,col=rev(red_cols[3:7]),pch=c(1,20,20,20,20),
       legend=c("100%", "  50%", "  25%", "12.5%", "6.25%"), cex=1.5,lwd=2, bty="n")


# AlleleDB sensitivity vs readDepth
i=1+2
a=AlleleDB_sensitivity_Specificity("SRR1552485_sub16_SNP_Allele_Specificity_matrix_t1E-02.bed")
plot(a[[1]], a[[2]], type=type, pch=20, lty=1, ylab="Sensitivity", col=blue_cols[i], ylim=c(0,1), xlab="Maternal reads fraction", cex=1.5,lwd=2, las=1)
i=2+2
a=AlleleDB_sensitivity_Specificity("SRR1552485_sub8_SNP_Allele_Specificity_matrix_t1E-02.bed")
lines(a[[1]], a[[2]], type=type, pch=20, lty=1, ylab="Sensitivity", col=blue_cols[i], cex=1.5,lwd=2)
i=3+2
a=AlleleDB_sensitivity_Specificity("SRR1552485_sub4_SNP_Allele_Specificity_matrix_t1E-03.bed")
lines(a[[1]], a[[2]], type=type, pch=20, lty=1, ylab="Sensitivity", col=blue_cols[i], cex=1.5,lwd=2)
i=4+2
a=AlleleDB_sensitivity_Specificity("SRR1552485_sub2_SNP_Allele_Specificity_matrix_t1E-04.bed")
lines(a[[1]], a[[2]],type=type, pch=20, lty=1, ylab="Sensitivity", col=blue_cols[i], cex=1.5,lwd=2)
i=5+2
a=AlleleDB_sensitivity_Specificity("SRR1552485_total_SNP_Allele_Specificity_matrix_t1E-05.bed")
lines(a[[1]], a[[2]], type=type, pch=1, lty=1, ylab="Sensitivity", col=blue_cols[i], ylim=c(0,1), xlab="Maternal reads fraction", cex=1.5,lwd=2)

legend("top", lty =1,col=rev(blue_cols[3:7]),pch=c(1,20,20,20,20),
       legend=c("100%", "  50%", "  25%", "12.5%", "6.25%"), cex=1.5,lwd=2, bty="n")

# AlleleDB precision vs readDepth
i=1+2
a=AlleleDB_sensitivity_Specificity("SRR1552485_sub16_SNP_Allele_Specificity_matrix_t1E-02.bed")
plot(a[[1]], a[[4]], type=type, pch=20, lty=1, ylab="Precision", col=blue_cols[i], ylim=c(0,1), xlab="Maternal reads fraction", cex=1.5,lwd=2, las=1)
i=2+2
a=AlleleDB_sensitivity_Specificity("SRR1552485_sub8_SNP_Allele_Specificity_matrix_t1E-02.bed")
lines(a[[1]], a[[4]], type=type, pch=20, lty=1, col=blue_cols[i], cex=1.5,lwd=2)
i=3+2
a=AlleleDB_sensitivity_Specificity("SRR1552485_sub4_SNP_Allele_Specificity_matrix_t1E-03.bed")
lines(a[[1]], a[[4]], type=type, pch=20, lty=1, col=blue_cols[i], cex=1.5,lwd=2)
i=4+2
a=AlleleDB_sensitivity_Specificity("SRR1552485_sub2_SNP_Allele_Specificity_matrix_t1E-04.bed")
lines(a[[1]], a[[4]],type=type, pch=20, lty=1, col=blue_cols[i], cex=1.5,lwd=2)
i=5+2
a=AlleleDB_sensitivity_Specificity("SRR1552485_total_SNP_Allele_Specificity_matrix_t1E-05.bed")
lines(a[[1]], a[[4]], type=type, pch=1, lty=1, col=blue_cols[i], ylim=c(0,1), xlab="Maternal reads fraction", cex=1.5,lwd=2)

legend("bottomright", lty =1,col=rev(blue_cols[3:7]),pch=c(1,20,20,20,20),
       legend=c("100%", "  50%", "  25%", "12.5%", "6.25%"), cex=1.5,lwd=2, bty="n")


# AlleleDB precision vs Specificity
i=1+2
a=AlleleDB_sensitivity_Specificity("SRR1552485_sub16_SNP_Allele_Specificity_matrix_t1E-02.bed")
plot(a[[1]], a[[3]], type=type, pch=20, lty=1, col=blue_cols[i], ylim=c(0,1), xlab="Maternal reads fraction", cex=1.5,lwd=2, ylab="Specificity", las=1)
i=2+2
a=AlleleDB_sensitivity_Specificity("SRR1552485_sub8_SNP_Allele_Specificity_matrix_t1E-02.bed")
lines(a[[1]], a[[3]], type=type, pch=20, lty=1, col=blue_cols[i], cex=1.5,lwd=2)
i=3+2
a=AlleleDB_sensitivity_Specificity("SRR1552485_sub4_SNP_Allele_Specificity_matrix_t1E-03.bed")
lines(a[[1]], a[[3]], type=type, pch=20, lty=1, col=blue_cols[i], cex=1.5,lwd=2)
i=4+2
a=AlleleDB_sensitivity_Specificity("SRR1552485_sub2_SNP_Allele_Specificity_matrix_t1E-04.bed")
lines(a[[1]], a[[3]],type=type, pch=20, lty=1, col=blue_cols[i], cex=1.5,lwd=2)
i=5+2
a=AlleleDB_sensitivity_Specificity("SRR1552485_total_SNP_Allele_Specificity_matrix_t1E-05.bed")
lines(a[[1]], a[[3]], type=type, pch=1, lty=1,  col=blue_cols[i], ylim=c(0,1), xlab="Maternal reads fraction", cex=1.5,lwd=2)

legend("bottom", lty =1,col=rev(blue_cols[3:7]),pch=c(1,20,20,20,20),
       legend=c("100%", "  50%", "  25%", "12.5%", "6.25%"), cex=1.5,lwd=2, bty="n")


## examine how tao affect sen and spec
# Sup_Fig1 D,E,F
## positive set (geneBiased):combine the regions with mat <=0.2 or >= 0.8 and gene biased
## negative set (geneSym) :combine the regions with 0.45 < mat <0.55 and gene symmetric
AlleleHMM_REGION_sensitivity_Specificity <- function (input_fp) {
  df=read.table(input_fp)
  df=df[df$V12+df$V13 >=10,]  # reads count within gene >=10
  colnames(df)[7]="GeneAS"
  colnames(df)[9]="AlleleHMMAS"
  colnames(df)[8]="AlleleDBAS"
  diff=0.05
  geneBiased=df[df$GeneAS!="S" & (df$V10<=0.2 |df$V10 >= 0.8 ),]
  geneSym=df[df$GeneAS=="S"& (df$V10>0.45 & df$V10 < 0.55 ),]
  
  df=rbind.data.frame(geneBiased, geneSym)
  AlleleHMM_result = data.frame(mat=seq(diff,1,diff), S_S=0, S_P=0, S_M=0, P_S=0, P_P=0, P_M=0, M_S=0, M_P=0, M_M=0)
  AS=c("S","P","M")
  for (m in 1:(dim(AlleleHMM_result)[1])){
    sub_df = df[df$V10 <= AlleleHMM_result$mat[m] & df$V10 >= (AlleleHMM_result$mat[m] - diff) ,]
    #V10 = mat/mat+pat
    for (i in seq(3)){
      for (j in seq(3)){
        AlleleHMM_result[m,3*(i-1)+j+1]=sum(sub_df$GeneAS==AS[i] & sub_df$AlleleHMMAS==AS[j])
      }
    }
  }
  AlleleHMM_result$mat = AlleleHMM_result$mat - diff/2
  AlleleHMM_result$sen = with(AlleleHMM_result, (sum(P_P)+sum(M_M))/(sum(P_S)+sum(P_P)+sum(M_M)+sum(M_S)))
  AlleleHMM_result$spec = with(AlleleHMM_result, sum(S_S)/(sum(S_S)+sum(S_P)+sum(S_M)))
  AlleleHMM_result$prec = with(AlleleHMM_result, (sum(P_P)+sum(M_M))/(sum(P_P)+sum(M_M)+sum(S_P)+sum(S_M)))
  return (list(AlleleHMM_result$mat, AlleleHMM_result$sen[1], AlleleHMM_result$spec[1], AlleleHMM_result$prec[1]))
}

AlleleDB_REGION_sensitivity_Specificity <- function (input_fp) {
  df=read.table(input_fp)
  df=df[df$V12+df$V13 >=10,]
  colnames(df)[7]="GeneAS"
  colnames(df)[9]="AlleleHMMAS"
  colnames(df)[8]="AlleleDBAS"
  diff=0.05
  geneSym=df[df$GeneAS=="S",]
  geneBiased=df[df$GeneAS!="S" & (df$V10<=0.2 |df$V10 >= 0.8 ),]
  df=rbind.data.frame(geneBiased, geneSym)
  AlleleDB_result = data.frame(mat=seq(diff,1,diff), S_S=0, S_P=0, S_M=0, P_S=0, P_P=0, P_M=0, M_S=0, M_P=0, M_M=0)
  AS=c("S","P","M")
  for (m in 1:(dim(AlleleDB_result)[1])){
    sub_df = df[df$V10 <= AlleleDB_result$mat[m] & df$V10 >= (AlleleDB_result$mat[m] - diff) ,]
    
    for (i in seq(3)){
      for (j in seq(3)){
        AlleleDB_result[m,3*(i-1)+j+1]=sum(sub_df$GeneAS==AS[i] & sub_df$AlleleDBAS==AS[j])
      }
    }
  }
  AlleleDB_result$mat = AlleleDB_result$mat - diff/2
  AlleleDB_result$sen = with(AlleleDB_result, (sum(P_P)+sum(M_M))/(sum(P_S)+sum(P_P)+sum(M_M)+sum(M_S)))
  AlleleDB_result$spec = with(AlleleDB_result, sum(S_S)/(sum(S_S)+sum(S_P)+sum(S_M)))
  AlleleDB_result$prec = with(AlleleDB_result, (sum(P_P)+sum(M_M))/(sum(P_P)+sum(M_M)+sum(S_P)+sum(S_M)))
  return (list(AlleleDB_result$mat, AlleleDB_result$sen[1], AlleleDB_result$spec[1], AlleleDB_result$prec[1]))
}

for (PREFIX in c("SRR1552485_total","SRR1552485_sub4","SRR4041366_total")){
  pdf(paste("tao_performance_test/",PREFIX,"_tao_performance_geneSym_and_geneBiased28.pdf",sep=""))
  par(mar=c(6.1, 7.1, 2.1, 2.1)) #d l u r 5.1, 4.1, 4.1, 2.1
  par(mgp=c(3,1,0))
  par(cex.lab=2.2, cex.axis=2.2)
  par(bty = 'n') 
  i=9
  a=AlleleHMM_REGION_sensitivity_Specificity(paste(PREFIX,"_SNP_Allele_Specificity_matrix_t1E-0",i,".bed",sep=""))
  d=AlleleDB_REGION_sensitivity_Specificity(paste(PREFIX,"_SNP_Allele_Specificity_matrix_t1E-0",i,".bed",sep=""))
  tao=paste("1E-0",i,sep="")
  sen=a[[2]]
  spec=a[[3]]
  prec=a[[4]]
  for (i in 8:1){
    a=AlleleHMM_REGION_sensitivity_Specificity(paste(PREFIX,"_SNP_Allele_Specificity_matrix_t1E-0",i,".bed",sep=""))
    tao=c(tao,paste("1E-0",i,sep=""))
    sen=c(sen, a[[2]])
    spec=c(spec, a[[3]])
    prec=c(prec,a[[4]])
  }
  plot(tao, sen, log="x", type='b', pch=1, ylab="Value", col="red", ylim=c(0,1), xlab="Tuning parameter Tao", cex=1.5,lwd=2, las=1)
  lines(tao, spec, type='b', pch=0, lty=1, col="blue", ylim=c(0,1), main=PREFIX, xlab="Maternal reads fraction", cex=1.5,lwd=2, las=1)
  legend("bottomright", lty =1,col=c("blue", "red"),pch=c(0,1),
         legend=c("Specificity", "Sensitivity"), cex=1.5,lwd=2, bty="n", title = PREFIX)
  dev.off()
}


## examine how read depth affect sen, spec and precision
# Sup_Fig5
## positive set (geneBiased):combine the regions with mat <=0.2 or >= 0.8 and gene biased
## negative set (geneSym) :combine the regions with 0.45 < mat <0.55 and gene symmetric
#read depth

f_list=c("SRR1552485_sub16_SNP_Allele_Specificity_matrix_t1E-02.bed","SRR1552485_sub8_SNP_Allele_Specificity_matrix_t1E-02.bed",
         "SRR1552485_sub4_SNP_Allele_Specificity_matrix_t1E-03.bed","SRR1552485_sub2_SNP_Allele_Specificity_matrix_t1E-04.bed",
         "SRR1552485_total_SNP_Allele_Specificity_matrix_t1E-05.bed")

readDepth=c(1/16, 1/8, 1/4, 1/2,1)


sen=c()
spec=c()
prec=c()
for (f in f_list){
  a=AlleleHMM_REGION_sensitivity_Specificity(f)
  sen=c(sen, a[[2]])
  spec=c(spec, a[[3]])
  prec=c(prec,a[[4]])
}
plot(readDepth, sen, main = "Sym+Bias" ,type='b',  pch=1, lty=1, col="red", ylim=c(0,1), xlab="Subsampling", cex=1.5,lwd=2, las=1, ylab="value")
lines(readDepth, sen, type='b',  pch=16, lty=1, col="red", ylim=c(0,1), xlab="Subsampling", cex=1.5,lwd=2, las=1)
lines(readDepth, spec, type='b', pch=15, lty=1, col="blue", ylim=c(0,1), xlab="Maternal reads fraction", cex=1.5,lwd=2, las=1)
lines(readDepth, prec, type='b', pch=17, lty=1, col="dark green", ylim=c(0,1), xlab="Maternal reads fraction", cex=1.5,lwd=2, las=1)

AlleleDB_sen=c()
AlleleDB_spec=c()
AlleleDB_prec=c()
for (f in f_list){
  a=AlleleDB_REGION_sensitivity_Specificity(f)
  AlleleDB_sen=c(AlleleDB_sen, a[[2]])
  AlleleDB_spec=c(AlleleDB_spec, a[[3]])
  AlleleDB_prec=c(AlleleDB_prec,a[[4]])
}
lines(readDepth, AlleleDB_sen, type='b', pch=1 ,lty=2, ylab="Value", col="red", ylim=c(0,1), xlab="Subsampling", cex=1.5,lwd=2, las=1)
lines(readDepth, AlleleDB_spec, type='b', pch=0, lty=2, col="blue", ylim=c(0,1), xlab="Maternal reads fraction", cex=1.5,lwd=2, las=1)
lines(readDepth, AlleleDB_prec, type='b', pch=2, lty=2, col="dark green", ylim=c(0,1), xlab="Maternal reads fraction", cex=1.5,lwd=2, las=1)

legend("right", lty =2,col=c("dark green","blue", "red"),pch=c(16,15,17),
       legend=c("Precision","Specificity", "Sensitivity"), cex=1.5,lwd=2, bty="n", title = "SNP independt;y")

