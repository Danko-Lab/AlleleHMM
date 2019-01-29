# Figure 5B,C and Sup_Fig9
# figure patameter
# export  7.92 x 5.92 inches
e_col= "red" #t_col("red", perc = 20)
g_col="blue" #t_col("blue", perc = 20)
par(mar=c(6.1, 7.1, 2.1, 2.1)) #d l u r 5.1, 4.1, 4.1, 2.1
par(mgp=c(3,1,0))
par(cex.lab=2.2, cex.axis=2.2)

# Figure 5B
# AlleleHMM Block size distribution
#setwd("~/Box Sync/Danko_lab_work/mouse_AlleleDB/GM12878/TssInBlocks/")
hist_of_HMM_block_from_bed6_solid <- function(m_all, p_all, add=F, col="blue",freq=F){
  mp_all = rbind.data.frame(m_all, p_all)
  mp_all$blockSize = mp_all$V3 - mp_all$V2
  hist(log10(mp_all$blockSize),
       freq=freq,
       col = col,
       xlim=c(0,8),
       xlab="log10(AlleleHMM block size in bp)",
       main="",
       add=add,
       las=1,
       bty="l")
}


hist_of_HMM_block_from_bed6_strip <- function(m_all, p_all, add=F, col="blue", freq=F){
  mp_all = rbind.data.frame(m_all, p_all)
  mp_all$blockSize = mp_all$V3 - mp_all$V2
  hist(log10(mp_all$blockSize),
       freq=freq,
       col = col,
       density=25,angle=45,
       xlab="log10(HMM block size in bp)",
       las=1,
       add=add)
}

###color transparant
t_col <- function(color, percent = 50, name = NULL) {
  #      color = color name
  #       percent = % transparency
  #          name = an optional name for the color
  ## Get RGB values for named color
  rgb.val <- col2rgb(color)
  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100-percent)*255/100,
               names = name)
  ## Save the color
  invisible(t.col)
}


#Engreitz_F1
m_all_SRR4041366 = read.table("SRR4041366_dedup_2/counts_minus_hmm_regions_t1e-05_interestingHets_IGV.bed")
p_all_SRR4041366 = read.table("SRR4041366_dedup_2/counts_plus_hmm_regions_t1e-05_interestingHets_IGV.bed")
hist_of_HMM_block_from_bed6_solid(m_all_SRR4041366, p_all_SRR4041366, col= e_col)
mp_all = rbind.data.frame(m_all_SRR4041366, p_all_SRR4041366)
mp_all$blockSize = mp_all$V3 - mp_all$V2
summary(mp_all$blockSize)

F1_ma_all=mp_all


#GM12878
GM_m_all=read.table("SRR1552485_total/counts_minus_hmm_regions_t1e-05_interestingHets_IGV.bed")
GM_p_all=read.table("SRR1552485_total/counts_plus_hmm_regions_t1e-05_interestingHets_IGV.bed")

hist_of_HMM_block_from_bed6_strip(GM_m_all, GM_p_all, T, g_col)
mp_all = rbind.data.frame(GM_m_all, GM_p_all)
mp_all$blockSize = mp_all$V3 - mp_all$V2
summary(mp_all$blockSize)

legend("topright", 
       legend = c("F1", "GM"),
       #pch=c(15,15),
       cex=2, 
       lty=c(0,0),
       #bty="n",
       lwd=1.5, 
       density=c(10000,25),
       angle=c(180,45),
       #angle=45,
       fill=c(e_col,g_col)
       , bty = "n"
)



##sup_fig_9

GM_m_all=read.table("SRR1552485_total/counts_minus_hmm_regions_t1e-05_interestingHets_IGV.bed")
GM_p_all=read.table("SRR1552485_total/counts_plus_hmm_regions_t1e-05_interestingHets_IGV.bed")

GM_m_sub4=read.table("SRR1552485_sub4/SRR1552485_sub4_minus_regions_t1E-03.sorted_binomtest_SigBlocks.bed")
GM_p_sub4=read.table("SRR1552485_sub4/SRR1552485_sub4_plus_regions_t1E-03.sorted_binomtest_SigBlocks.bed")

freq=F
hist_of_HMM_block_from_bed6_strip(GM_m_all, GM_p_all, add=F, g_col, freq)
hist_of_HMM_block_from_bed6_solid(GM_m_sub4, GM_p_sub4, add=T, "yellow", freq)
hist_of_HMM_block_from_bed6_strip(GM_m_all, GM_p_all, add=T, g_col, freq)

legend("topright", 
       legend = c("100%", "25%"),
       #pch=c(15,15),
       cex=2, 
       lty=c(0,0),
       #bty="n",
       lwd=1.5, 
       density=c(25,10000),
       angle=c(45,180),
       #angle=45,
       fill=c("blue", "yellow")
       , bty = "n"
)

freq=T
hist_of_HMM_block_from_bed6_strip(GM_m_all, GM_p_all, add=F, g_col, freq)
hist_of_HMM_block_from_bed6_solid(GM_m_sub4, GM_p_sub4, add=T, "yellow", freq)
hist_of_HMM_block_from_bed6_strip(GM_m_all, GM_p_all, add=T, g_col, freq)

legend("topright", 
       legend = c("100%", "25%"),
       #pch=c(15,15),
       cex=2, 
       lty=c(0,0),
       #bty="n",
       lwd=1.5, 
       density=c(25,10000),
       angle=c(45,180),
       #angle=45,
       fill=c("blue", "yellow")
       , bty = "n"
)

# Figure 5C

## gene per block
#setwd("~/Box Sync/Danko_lab_work/mouse_AlleleDB/GM12878/GeneInBlocks/")
#GM12878
m_withG=read.table("counts_minus_hmm_regions_t1e-05_interestingHets_IGV_genePerBlock.txt")
p_withG=read.table("counts_plus_hmm_regions_t1e-05_interestingHets_IGV_genePerBlock.txt")
Block_with_no_gene=4028-2-(dim(m_withG)[1] + dim(p_withG)[1])

##SRR4041366_dedup_2
Em_withG=read.table("SRR4041366_dedup_2/counts_minus_hmm_regions_t1e-05_interestingHets_IGV_genePerBlock.txt")
Ep_withG=read.table("SRR4041366_dedup_2/counts_plus_hmm_regions_t1e-05_interestingHets_IGV_genePerBlock.txt")
EBlock_with_no_gene= 3485-2-(dim(Em_withG)[1] + dim(Ep_withG)[1])

hist(c(rep(0,Block_with_no_gene), m_withG$V1, p_withG$V1)
     ,col= g_col
     ,density=20,angle=45
     , breaks = seq(-0.0001,200,1)
     , freq = F
     , xlim=c(0,20)
     ,xlab="Number of genes in each HMM block"
     ,main=" "
     ,las=1
    # , add=T
)
 
 hist(c(rep(0,EBlock_with_no_gene), Em_withG$V1, Ep_withG$V1)
     , breaks = seq(-0.0001,200,1)
     #,density=20,angle=135
     , freq = F
     , xlim=c(0,20)
     #,ylim=c(0,0.6)
     ,xlab="Number of genes in each HMM block"
     ,main=""
     ,col= e_col
     ,add =T
    
)


hist(c(rep(0,Block_with_no_gene), m_withG$V1, p_withG$V1)
     ,col= g_col
     ,density=20,angle=45
     , breaks = seq(-0.0001,200,1)
     , freq = F
     , xlim=c(0,20)
     ,xlab="Number of genes in each HMM block"
     ,main=" "
     , add=T
)






