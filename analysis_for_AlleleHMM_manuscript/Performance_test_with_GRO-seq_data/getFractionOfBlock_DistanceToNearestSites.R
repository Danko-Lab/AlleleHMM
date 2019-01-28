#R --vanilla --slave --args $(pwd) "counts_hmm_regions_t*_interestingHets_5head_distance_toclosest-dReg_AccumulateCounts.txt" counts_hmm_regions_tX_interestingHets_5head_distance_toclosest-dReg_AccumulateCounts.pdf counts_hmm_regions_tX_interestingHets_5head_distance_toclosest-dReg_At5Kb.pdf < getFractionOfBlock_DistanceToNearestSites.R

#arguments here
args=(commandArgs(TRUE))
setwd(args[1])
f_pattern = glob2rx(args[2])


library("colorspace")

#f_pattern = glob2rx("counts_hmm_regions_t*_interestingHets_5head_distance_toclosest-dReg_AccumulateCounts.txt")
f_list = list.files(pattern=f_pattern) 
print(f_list)
i=1
tmp=strsplit(f_list[i], "_")[[1]]
t_list=as.numeric(substr(tmp[length(tmp)-5],2,nchar(tmp[length(tmp)-5])))


data = read.table(f_list[i])
pdf(args[3])
fraction_at5K=max((data$V2/data$V2[length(data$V2)])[data$V1 <= 5000])
fraction_at10K=max((data$V2/data$V2[length(data$V2)])[data$V1 <= 10000])
plot(main = "M|P Block 5 prime head to nearest dREG sites (Both up and downstream)",
     data$V1[data$V1 <= 20000],  (data$V2/data$V2[length(data$V2)])[data$V1 <= 20000], 
     #data$V1,  (data$V2/data$V2[length(data$V2)]), 
     type='l', pch=20,
     xlab='distance to closest dREG site', ylab='fraction of M|P blocks',
     col=rainbow_hcl(length(f_list))[1])

for (i in 2:length(f_list)){
  #print(f_list[i])
  data = read.table(f_list[i])
  #print(min(data$V1))
  #print(max(data$V1))
  tmp=strsplit(f_list[i], "_")[[1]]
  t_list=c(t_list, as.numeric(substr(tmp[length(tmp)-5],2,nchar(tmp[length(tmp)-5]))))
  fraction_at5K=c(fraction_at5K, max((data$V2/data$V2[length(data$V2)])[data$V1 <= 5000]))
  fraction_at10K=c(fraction_at10K,max((data$V2/data$V2[length(data$V2)])[data$V1 <= 10000]))
  lines(data$V1[data$V1 <= 20000],  (data$V2/data$V2[length(data$V2)])[data$V1 <= 20000],lty=i, pch=20, col=rainbow_hcl(length(f_list))[i])
  #lines(data$V1,  (data$V2/data$V2[length(data$V2)]),lty=i, pch=20, col=rainbow_hcl(length(f_list))[i])
}
legend("bottomright", pch=NA, lty = 1:9,col=rainbow_hcl(length(f_list)),
       legend=t_list)

dev.off()
###choose a distance, say 5kb, and plot the fraction close to a dREG site at each of the thresholds 
pdf(args[4])
plot(log10(t_list), fraction_at5K, xlab='log10(t)', ylab='fraction of M|P blocks at 5Kb', type="b")
#plot(log10(t_list), fraction_at10K, xlab='log10(t)', ylab='fraction of M|P blocks at 10Kb')
dev.off()


#text(75, data$V2[data$V1 ==75]/mappable_Reads, labels= round(data$V2[data$V1 ==75]/mappable_Reads, 2), cex=1, pos= 3, col='red')
