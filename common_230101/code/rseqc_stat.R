.libPaths(c("/home/shenlm/R/x86_64-pc-linux-gnu-library/4.1","/opt/R/4.1.2/lib/R/library"))
library(getopt)

spec = matrix(c('input','i',1,"character",'work_dir','w',1,"character",'function_file','f',1,"character",'sample','s',1,"character"),byrow = TRUE,ncol = 4)
opt = getopt(spec)
setwd(opt$work_dir)
source(opt$function_file)

rseqc_stat <- read.delim(opt$input,header = TRUE,sep=" ")
rseqc_tmp <- rseqc_stat
colnames(rseqc_tmp) <- c("type","tbase","tcount","tag_per_kb","other")

rseqc <- data.frame(counts=rseqc_tmp$tag_per_kb,type=rseqc_tmp$type)
rseqc$sample <- rep(opt$sample,nrow(rseqc))

pdf("tmp_readdis.pdf",width=10,height=5)
rseqc_barplot(rseqc_stat)
dev.off()

write.table(rseqc,"tmp_rseqc.txt",quote = FALSE,sep=" ",row.names = FALSE,col.names = FALSE)
