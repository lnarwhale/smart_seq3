library(dplyr)
library(tidyr)
setwd("/Users/mosh/Desktop/ssh/221106/result/")
in1 <- "/Users/mosh/Desktop/ssh/221106/result/"
gene_list <- read.csv("/Users/mosh/Desktop/ssh/221106/gtf/mm_DNA-binding_transcription_repressor_activity.csv",
                      header = FALSE)
mode <- read.csv("../mode_test/1zga.csv",header = TRUE)
colnames(gene_list) <- "gene"
path1 <- list.files(in1,full.names = TRUE)
path1 <- path1[grep(pattern = "*_cager_tss.csv",path1)]
names <- sub("_cager_tss.csv","",basename(path1))

gene_tpm <- list()
gene_loc <- list()
all <- list()
for (i in 1:length(names)) {
  gene_tpm[[i]] <- read.csv(paste0(names[[i]],"_cager_tss.csv"),header = TRUE)
  gene_tpm[[i]] <- select(gene_tpm[[i]],"seqnames","start","end","tpm")
  gene_tpm[[i]] <- unite(gene_tpm[[i]],"symbol",c("seqnames","start","end"),sep = "_",remove = TRUE)
  gene_loc[[i]] <- read.csv(paste0("tss/",names[i],"_cager_tss_x2_locgene.csv"),
                            header = FALSE)
  gene_loc[[i]] <- select(gene_loc[[i]],"V1","V2","V3","V5")
  gene_loc[[i]] <- unite(gene_loc[[i]],"symbol",c("V1","V2","V3"),sep = "_",remove = TRUE)
  all[[i]] <- merge(gene_tpm[[i]],gene_loc[[i]],by="symbol")
  all[[i]] <- select(all[[i]],"tpm","V5")
  colnames(all[[i]]) <- c("tpm","gene")
}
names(gene_tpm) <- names
names(gene_loc) <- names
names(all) <- names


gene_t <- merge(gene_list,all[[1]],by="gene",all = TRUE)
for (i in 2:length(names)) {
  gene_t <- merge(gene_t,all[[i]],by="gene",all = TRUE)
}


gene_t$gene <- make.unique(gene_t$gene)
rownames(gene_t) <- gene_t[,1]
gene_t <- gene_t[,-1]
colnames(gene_t) <- c(names,"x2","MII")
gene_t[is.na(gene_t)] <- 0
gene_t <- log10(gene_t+1)
gene_t <- scale(gene_t)

png("heatmap_ling.png",width = 1500,height = 1500)
p <- pheatmap::pheatmap(gene_t,cluster_rows = FALSE,cluster_col=TRUE ,border=FALSE)
dev.off()
