.libPaths(c("/home/shenlm/R/x86_64-pc-linux-gnu-library/4.1","/opt/R/4.1.2/lib/R/library"))
library(getopt)
spec = matrix(c('input','i',1,"character",'work_dir','w',1,"character"),byrow = TRUE,ncol = 4)
opt = getopt(spec)
setwd(opt$work_dir)

countmatrix_df <- read.table(opt$input,header = TRUE,sep=",")
out_fileName <- sapply(colnames(countmatrix_df),function(x){paste(x,".count", sep='')})

for(r in 2:ncol(countmatrix_df)){
  write.table(countmatrix_df[,c(1,r)], file=out_fileName[r],row.names = F,col.names=F,sep="\t",quote =FALSE)
}

sample_name <- colnames(countmatrix_df)[-1]
group <- rep("Group1",length(sample_name))
work_dir <- getwd()
filename <- sapply(colnames(countmatrix_df),function(x){paste(work_dir,"/",x,".count", sep='')})[-1]
count <- rep(2,length(sample_name))

df <- data.frame(sample_name,group,filename,count)
write.table(df, file="qualimeta.txt",row.names = F,col.names=F,sep="\t",quote =FALSE)
