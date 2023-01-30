.libPaths(c("/home/shenlm/R/x86_64-pc-linux-gnu-library/4.1","/opt/R/4.1.2/lib/R/library"))
library(getopt)
library(reshape2)

spec = matrix(c('dis','d',1,"character",'total','t',1,"character",'high','h',1,"character",'low','l',1,"character",'work_dir','w',1,"character",'function_file','f',1,"character",'sample','s',1,"character"),byrow = TRUE,ncol = 4)
opt = getopt(spec)
setwd(opt$work_dir)
source(opt$function_file)

quali_dis <- as.data.frame(t(na.omit(t(read.table(opt$dis,header=FALSE,sep=" ")))))
colnames(quali_dis) <- c("type","count","perc") 
quali_dis$perc <- paste(quali_dis$type,"(",quali_dis$count,",",quali_dis$perc,"%",")",sep = " ")
quali_dis$count <- as.numeric(gsub(" ","",quali_dis$count))
quali_region <- data.frame(counts=quali_dis$count,type=quali_dis$type)
quali_region$sample <- rep(opt$sample,nrow(quali_region))

quali_total <- read.delim(opt$total,header = TRUE,sep="\t")
quali_high <- read.delim(opt$high,header = TRUE,sep="\t")
quali_low <- read.delim(opt$low,header = TRUE,sep="\t")
colnames(quali_total) <- c("position","coverage_total")
colnames(quali_high) <- c("position","coverage_high")
colnames(quali_low) <- c("position","coverage_low")
quali <- merge(quali_total,quali_high,by="position")
quali<-merge(quali,quali_low,by="position")
quali_df<-melt(quali,id=c("position"))
quali_df$value <- log10(quali_df$value)
colnames(quali_df) <- c("position","type","log10_coverage")


pdf("tmp_readdis_pie.pdf",width=5,height=5)
pieplot(quali_dis)
dev.off()

pdf("tmp_readdis.pdf",width=10,height=5)
read_dis_lineplot(quali_df)
dev.off()

write.table(quali_region,"tmp_quali.txt",quote = FALSE,sep=" ",row.names = FALSE,col.names = FALSE)
