# R --vanilla --slave --args $(pwd) counts.txt Min_count Max_Pvalus< test.R 

args=(commandArgs(TRUE))

folder = args[1]
count_fp = args[2]
#folder = '/Users/shaopei/Box\ Sync/Danko_lab_work/Script/After_AlleleDB_pipeline/' #args[1]
#count_fp = 'interestingHets_plus.txt' #args[2]


if(length(args)==2){
  Min_count = 5
  Max_Pvalus = 1
}else if(length(args)==3){
  Min_count = as.numeric(args[3])
  Max_Pvalus = 1
}else{
  Min_count = as.numeric(args[3])
  Max_Pvalus = as.numeric(args[4])
}
cat('pwd=',folder,'\n')
cat('input_file =',count_fp, '\n')
cat('min reads count = ',Min_count, '\n')
cat('max p-value = ',Max_Pvalus, '\n')

####functions
get.mat_allele_count <- function(count){
  mat_allele_count = apply(count, 1,function(x){
    if(x['mat_all']=='T') return (as.numeric(x['cT']));
    if(x['mat_all']=='C') return (as.numeric(x['cC']));
    if(x['mat_all']=='G') return (as.numeric(x['cG']));
    if(x['mat_all']=='A') return (as.numeric(x['cA']))
    else return (NA)})
  return (mat_allele_count)
}
get.pat_allele_count <- function(count){
  pat_allele_count = apply(count, 1,function(x){
    if(x['pat_all']=='T') return (as.numeric(x['cT']));
    if(x['pat_all']=='C') return (as.numeric(x['cC']));
    if(x['pat_all']=='G') return (as.numeric(x['cG']));
    if(x['pat_all']=='A') return (as.numeric(x['cA']))
    else return (NA)})
  return (pat_allele_count)
}

###main

setwd(folder)
count <- read.table(count_fp,header=T,sep = "\t")
dim(count)
count <- count[count$SymCls !='Weird', ]
dim(count)
count$mat_allele_count <- get.mat_allele_count(count)
count$pat_allele_count <- get.pat_allele_count(count)
count$total_reads_count <-count$mat_allele_count + count$pat_allele_count
#count$total_reads_count <-count$cA + count$cC + count$cG + count$cT
dim(count)
count <- count[count$total_reads_count >= Min_count, ]
dim(count)
summary(count$SymPval)
count <- count[count$SymPval <= Max_Pvalus, ]
dim(count)

write.table(count,paste(substr(count_fp, 1, nchar(count_fp)-4),'_MinCount',Min_count,'_MaxPvalue',Max_Pvalus,'.txt',sep=''), row.names=FALSE, sep="\t", quote = FALSE)

