.libPaths(c("/home/shenlm/R/x86_64-pc-linux-gnu-library/4.1","/opt/R/4.1.2/lib/R/library"))
library(getopt)
library(ggplot2)
library(ggbiplot)
spec = matrix(c('matrix','m',1,"character",
                'gene_exp','g',1,"character",
                'out_dir','o',1,"character"),
              byrow=TRUE,ncol=4)
opt=getopt(spec)
setwd(opt$out_dir)

matrix <- read.csv(opt$matrix,header = TRUE)
rownames(matrix) <- matrix[,2]
matrix$Group <- NULL

fpkm <- read.csv(opt$gene_exp,header = T,row.names = 1)
fpkm <- t(fpkm)
pca_fp <- prcomp(fpkm)
a <- ggbiplot(pca_fp,
         var.axes = F,
         obs.scale = 1,
         groups = matrix[,1],
         ellipse = T,
         circle = F)+
  geom_text(
    aes(label=rownames(fpkm)),
    vjust=1.5,
    size=2
  )
png("pca_plot.png",width = 1000,height = 1000)
a
dev.off()
