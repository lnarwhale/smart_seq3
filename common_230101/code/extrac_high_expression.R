suppressMessages(library(getopt))
suppressMessages(library(tidyverse))
suppressMessages(library(dplyr))
suppressMessages(library(edgeR))
suppressMessages(library(clusterProfiler))
suppressMessages(library(AnnotationDbi))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(org.Mm.eg.db))
suppressMessages(library(AnnotationHub))
suppressMessages(library(stringr))

spec = matrix(c('genefpkm','t',1,"character",
                'work_dir','w',1,"character",
                'control','c',1,"character",
                'sample','s',1,"character"),byrow = TRUE,ncol = 4)
opt=getopt(spec)
print(opt$genefpkm)
print(opt$work_dir)
print(opt$control)
print(opt$sample)
setwd(opt$work_dir)
raw_count <- read.csv(opt$genefpkm)
control <- opt$control
all_sample <- opt$sample


#---------------------------------test-----------------------------------------
# getwd()
# raw_count <- read.csv("gene_count_matrix.csv")
# control <- c("monmouMII")
# all_sample <- c("monmou3C","monmou1C","monmou2C")
#------------------------------------------------------------------------------
#---------------------------------tip-----------------------------------------
#--------it is begin for lots of samples,but it is hard to achieve it just by R
#so I achieve it by shell for , hope i could achieve it just by R in the future

for (i in 1:nrow(raw_count)) {
  raw_count$gene_id[i] <- str_sub(raw_count$gene_id[i],0,str_locate(raw_count$gene_id[i],"[|]")[1]-1)
}
raw_count_name <- colnames(raw_count)
for (i in 1:length(all_sample)) {
  name <- all_sample[i]
  con_loc <- which(raw_count_name==control)
  sample_loc <- which(raw_count_name==name)
  tmp <- data.frame(raw_count$gene_id,raw_count[,con_loc],raw_count[,sample_loc])
  colnames(tmp) <- c("gene_id",control,name)
  group <- c(name,control)
  y <- DGEList(counts = tmp[,2:3],genes = tmp[,1],group = group)
  keep <- rowSums(cpm(y)>1) >= 2
  y <- y[keep, , keep.lib.sizes=FALSE]
  y <- calcNormFactors(y)
  y_bcv <- y
  bcv <- 0.1
  et <- exactTest(y_bcv, dispersion = bcv ^ 2)
  gene1 <- decideTestsDGE(et, p.value = 0.05, lfc = 0)
  colnames(gene1) <- "Signifi"
  gene_res_all <- cbind(y$genes,y$counts,et$table,gene1)
  gene_res_all <- gene_res_all[order(gene_res_all$PValue),]
  gene_res_up <- gene_res_all[gene_res_all$Signifi==1,]
  gene_res_down <- gene_res_all[gene_res_all$Signifi==-1,]
  gene_res <- rbind(gene_res_up,gene_res_down)
  write.table(gene_res_up$genes,file=paste0(name,"_high.csv"),row.names = FALSE,col.names = FALSE)
}








