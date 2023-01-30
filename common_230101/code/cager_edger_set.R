.libPaths(c("/home/shenlm/R/x86_64-pc-linux-gnu-library/4.1","/opt/R/4.1.2/lib/R/library"))
#library
library(getopt)
library(CAGEr)
library(ggplot2)
library(edgeR)
library(tidyr)
library(RColorBrewer)
library(UpSetR)
library(grid)
library(dplyr)
library(ChIPseeker)
library(biomaRt)
library(GenomicFeatures)
#species
#out_dir
#bam_dir
#pwm_dir
#mode(high-1 low-2)
#read
spec = matrix(c('species','s',1,"character",
                'out_dir','o',1,"character",
                'bam_dir','b',1,"character",
                'pwm_dir','p',1,"character",
                'mode','m',1,"character"),
              byrow = TRUE,ncol=4)
opt=getopt(spec)
#---
species <- opt$species
setwd(opt$out_dir)
inputDir1 <- opt$bam_dir
pwmdir <- opt$pwm_dir
select_sample <- opt$mode

#------------------------------------------------------------motif-scan-read-in-----------------------------
TBP <-read.csv(paste0(pwmdir,"/tbp.txt"),header = TRUE)
TBP <- t(TBP)
rownames(TBP) <- c('A','C','G','T')
tct <- read.csv(paste0(pwmdir,"/human_TCT.csv"),header = TRUE)
tct <- as.matrix(tct)


paths1 <- list.files(inputDir1,full.names = T)
pathsToInputFiles1 <- paths1[grep(pattern = "*_sorted.bam", paths1)]
pathsToInputFiles <- c(pathsToInputFiles1)
edger="TRUE"

samples <- sub("_sorted.bam","",basename(pathsToInputFiles))
names <- sub("_sorted.bam","",basename(pathsToInputFiles))
samples <- gsub("x1_","",samples)
samples <- gsub("x2_","",samples)
names <- gsub("x1_","",names)
names <- gsub("x2_","",names)

if (select_sample=="1") {
  name_s <- names[1]
}else if (select_sample=="2") {
  name_s <- names[2]
}else if (select_sample=="NULL") {
  name_s <- "NULL"
}

#-----------select---------------------------
group <- c(1,2)
merge_sam <- select_sample

#library
if (species == "human") {
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  seqnames(BSgenome.Hsapiens.UCSC.hg38) <- gsub("chr","",seqnames(BSgenome.Hsapiens.UCSC.hg38))
  BSname <- "BSgenome.Hsapiens.UCSC.hg38"
  BS <- BSgenome.Hsapiens.UCSC.hg38
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  org <- "org.Hs.eg.db"
}else if (species == "mouse") {
  library(BSgenome.Mmusculus.UCSC.mm10)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(org.Mm.eg.db)
  seqnames(BSgenome.Mmusculus.UCSC.mm10) <- gsub("chr","",seqnames(BSgenome.Mmusculus.UCSC.mm10))
  BSname <- "BSgenome.Mmusculus.UCSC.mm10"
  BS <- BSgenome.Mmusculus.UCSC.mm10
  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  org <- "org.Mm.eg.db"
}else if (species== "monkey") {
  library(BSgenome.Mfascicularis.NCBI.5.0)
  seqnames(BSgenome.Mfascicularis.NCBI.5.0) <- gsub("MFA","",seqnames(BSgenome.Mfascicularis.NCBI.5.0))
  BS <- BSgenome.Mfascicularis.NCBI.5.0
  BSname <- "BSgenome.Mfascicularis.NCBI.5.0"
  txdb <- makeTxDbFromGFF("/data1/shenluemou/reference/gtf/Macaca_fascicularis/MFA1912/macaca_MFA1912.gtf",format = "gtf")
}

#quality control
LICAGE <- new("CAGEset", genomeName = BSname, 
              inputFiles = pathsToInputFiles, inputFilesType = "bam", 
              sampleLabels = names)
getCTSS(LICAGE)


#----------------------------------------------------edger-------------------------------
if (edger) {
  raw_count <- data.frame(LICAGE@CTSScoordinates,LICAGE@tagCountMatrix)
  new_raw <- unite(raw_count,"symbol",c("chr","pos","strand"),sep = "_",remove = TRUE)
  test <- new_raw
  rownames(test) <- test[,1]
  test$symbol <- NULL
  y <- DGEList(counts = test[,1:length(group)],genes = rownames(test),group = group)
  keep <- rowSums(cpm(y)>1) >= 2
  y <- y[keep, , keep.lib.sizes=FALSE]
  y <- calcNormFactors(y)
  y_bcv <- y
  bcv <- 0.1
  et <- exactTest(y_bcv,dispersion = bcv^2)
  #gene1 <- decideTestsDGE(et,p.value = 0.05,lfc = 1)
  #colnames(gene1) <- "Signifi"
  results <- cbind(y$genes,y$counts,et$table)
  for (i in 1:nrow(results)) {
    if (results[i,6] < 0.05) {
      if (results[i,4]> 1) {
        results[i,7] <- 1
      }else if (results[i,4] < -1) {
        results[i,7] <- -1
      }else{
        results[i,7] <- 0
      }
    }else{
      results[i,7] <- "0"
    }
  }
  colnames(results)[7] <- "Signifi"
  if (select_sample=="2") {
    different <- results[which(results$Signifi==1),]
  }else if (select_sample=="1") {
    different <- results[which(results$Signifi==-1),]
  }else if (select_sample=="NULL") {
    different <- results[which(results$Signifi==0),]
  }
  count <- separate(data = different,col = genes,into = c("chr","pos","strand"),sep = "_")
  count_to_taq_count <- subset(count,select=samples)
  count_to_tss <- subset(count,select=c(chr,pos,strand))
  count_to_tss$pos <- as.numeric(count_to_tss$pos)
  rownames(count_to_taq_count) <- NULL
  rownames(count_to_tss) <- NULL    
  LICAGE@tagCountMatrix <- count_to_taq_count
  LICAGE@CTSScoordinates <- count_to_tss
  LICAGE@inputFiles <- LICAGE@sampleLabels
  LICAGE@sampleLabels <- LICAGE@sampleLabels
  LICAGE@librarySizes <- LICAGE@librarySizes
  normalizeTagCount(LICAGE, method = "powerLaw", fitInRange = c(100, 1000), alpha = 1.14, T = 1000000)
}

#exportCTSStoBedGraph(LICAGE, values = "normalized", format = "bedGraph", oneFile = FALSE)
clusterCTSS(object = LICAGE, threshold = 1, thresholdIsTpm = TRUE, nrPassThreshold = 1, 
            method = "distclu", maxDist = 20, removeSingletons = TRUE, keepSingletonsAbove = 5)
cumulativeCTSSdistribution(LICAGE, clusters = "tagClusters")
quantilePositions(LICAGE, clusters = "tagClusters", qLow = 0.1, qUp = 0.9)

# Extract tag clusters 
sampleLabels <- unname(sampleLabels(LICAGE))
tc.l <- lapply(sampleLabels, 
               function (x) tagClusters(LICAGE, sample = x, 
                                        returnInterquantileWidth = TRUE, qLow = 0.1, qUp = 0.9))
names(tc.l) <- sampleLabels

## Convert tag clusters list into a GRanges list
tc.grl <- GRangesList()
tc.grl <- lapply(tc.l, function(x) GRanges(seqnames = x$chr, 
                                           ranges = IRanges(start = x$start, end = x$end),
                                           strand = x$strand,
                                           nr_ctss = x$nr_ctss,
                                           dominant_ctss = x$dominant_ctss,
                                           tpm = x$tpm,
                                           tpm.dominant_ctss = x$tpm.dominant_ctss,
                                           q_0.1 = x$q_0.1,
                                           q_0.9 = x$q_0.9,
                                           interquantile_width = x$interquantile_width,
                                           seqlengths = seqlengths(BS)))
# add seqinfo information                
for (i in 1:length(tc.grl)) {
  seqinfo(tc.grl[[i]]) <- seqinfo(BS) 
}
## Dominant TSS position centered GRanges object
domTSS.grl <- list()
for (i in 1:length(tc.grl)) {
  domTSS.grl[[i]] <- GRanges(seqnames = seqnames(tc.grl[[i]]),
                             ranges = IRanges(start = tc.grl[[i]]$dominant_ctss, 
                                              end = tc.grl[[i]]$dominant_ctss),
                             strand = strand(tc.grl[[i]]),
                             nr_ctss = tc.grl[[i]]$nr_ctss,
                             dominant_ctss = tc.grl[[i]]$dominant_ctss,
                             tpm = tc.grl[[i]]$tpm, 
                             tpm.dominant_ctss = tc.grl[[i]]$tpm.dominant_ctss, 
                             q_0.1 = tc.grl[[i]]$q_0.1,
                             q_0.9 = tc.grl[[i]]$q_0.9,
                             interquantile_width = tc.grl[[i]]$interquantile_width,
                             tc_start = start(tc.grl[[i]]),
                             tc_end = end(tc.grl[[i]]))
  
  seqlevels(domTSS.grl[[i]]) <- seqlevels(tc.grl[[i]])
  seqlengths(domTSS.grl[[i]]) <- seqlengths(tc.grl[[i]])
  genome(domTSS.grl[[i]]) <- genome(BS)
}
names(domTSS.grl) <- names(tc.grl)

#------------------------------------------------------------------------------
#03_genomic_location_tc.R
# - rename chromosome names in txdb to match NCBI naming scheme - #
seqlevels(txdb) <- gsub("chr","",seqlevels(txdb))
tc_selected <- tc.grl[samples]
names(tc_selected) <- names
peakAnno_list <- lapply(tc_selected, function(x) annotatePeak(x, TxDb = txdb,  tssRegion = c(-1000, 1000), annoDb = org, sameStrand = TRUE, verbose = FALSE))
names(peakAnno_list) <- names


#select the promoter
if (select_sample!="NULL") {
  tmp <- as.data.frame(peakAnno_list[[as.numeric(merge_sam)]])
  tmp <- tmp[which(tmp$annotation=="Promoter"),]  
}else if (select_sample=="NULL") {
  tmp1 <- as.data.frame(peakAnno_list[[as.numeric(1)]])
  tmp1 <- tmp1[which(tmp1$annotation=="Promoter"),]
  tmp2 <- as.data.frame(peakAnno_list[[as.numeric(2)]])
  tmp2 <- tmp2[which(tmp2$annotation=="Promoter"),]
  tmp <- union(tmp1,tmp2)
}


a_dom <- subset(tmp,select=c(dominant_ctss))
a_all <- merge(domTSS.grl[[as.numeric(merge_sam)]],a_dom,by="dominant_ctss")
domTSS_slic_mergLandR.grl <- GRanges(seqnames = a_all$seqnames,
                                      ranges = IRanges(start = a_all$dominant_ctss, 
                                                      end = a_all$dominant_ctss),
                                      strand = a_all$strand,
                                      nr_ctss = a_all$nr_ctss,
                                      dominant_ctss = a_all$dominant_ctss,
                                      tpm = a_all$tpm, 
                                      tpm.dominant_ctss = a_all$tpm.dominant_ctss, 
                                      q_0.1 = a_all$q_0.1,
                                      q_0.9 = a_all$q_0.9,
                                      interquantile_width = a_all$interquantile_width,
                                      tc_start = a_all$start,
                                      tc_end = a_all$end)

#-----------------------------------
#Figures distribution of tag cluster/promoter interquantile widths
iq_all.l <- as.data.frame(domTSS_slic_mergLandR.grl$interquantile_width)
names(iq_all.l) <- "iq_width"
iq_all.l$sample <- rep(names[as.numeric(merge_sam)], nrow(iq_all.l))
p <- ggplot(iq_all.l, aes(iq_width)) +
  geom_histogram(aes(x = iq_width, y = ..density..), binwidth = 1, 
                 fill = "gray60", col = "black", size = 0.1) +
  theme_bw() +
  theme(text = element_text(size = 16, colour = "black"),
        legend.title = element_blank(),
        axis.title.x = element_text(colour = "black"),
        axis.title.y = element_text(colour = "black"),
        axis.text.x = element_text(colour = "black", size = 14),
        axis.text.y = element_text(colour = "black", size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(colour = "black", fill = "gray87")) +
  labs(x = "Interquantile width", y = "relative trequency") +
  #coord_equal(ratio = 1) +
  xlim(0, 50)
png(file = paste0(name_s,"_IQ_width_all_slic.png"), height = 500, width = 500)
p
dev.off()

tmp <- as.data.frame(domTSS_slic_mergLandR.grl)
write.csv(tmp,paste0(name_s,"_tss_count.csv"))

#---------------motif scan------------------------------------------------

up <- 40
down <- 0
range <- c(-up, down)
win <- up + down
domTSS_slic_win.grl <- promoters(domTSS_slic_mergLandR.grl,up=up,down=down)
domTSS_slic_win.grl <- domTSS_slic_win.grl[width(trim(domTSS_slic_win.grl))==win]
domTSS_slic_seq.l <- getSeq(BS,domTSS_slic_win.grl)
  
sum=0
simple=1
TBP_sy <- c()
all_sy <- c()
for (i in 1:length(domTSS_slic_seq.l)) {
  tmp_seq <- DNAString(as.character(domTSS_slic_seq.l[i]))
  tmp_count <- countPWM(TBP,tmp_seq,min.score = "80%")
  if (tmp_count > 0) {
    sum=sum+1
    TBP_sy <- c(TBP_sy,simple)
  }
  all_sy <- c(all_sy,simple)
  simple=simple+1
}
tbp <- sum/length(domTSS_slic_seq.l)
#plot
TBP_topie_name <- c('TATA-box','not TATA-box')
TBP_topie_percent <- c(tbp,1-tbp)
TBP_topie_percent <- round(TBP_topie_percent*100,2)
TBP_topie_name <- paste(TBP_topie_name," | ",TBP_topie_percent,"%",sep = "")
png(paste0(name_s,"_TATA-box_pie.png"),height = 800,width = 800)
pie(TBP_topie_percent,labels = TBP_topie_name,main=paste0(name_s,"_TATA-box_pie"),radius = 0.8,
    cex=1.0,clockwise = TRUE)
dev.off()

up <- 10
down <- 10
range <- c(-up, down)
win <- up + down
domTSS_slic_win.grl <- promoters(domTSS_slic_mergLandR.grl,up=up,down=down)
domTSS_slic_win.grl <- domTSS_slic_win.grl[width(trim(domTSS_slic_win.grl))==win]
domTSS_slic_seq.l <- getSeq(BS,domTSS_slic_win.grl)
  
TCT_sy <- c()
sum=0
simple=1
for (i in 1:length(domTSS_slic_seq.l)) {
  tmp_seq <- DNAString(as.character(domTSS_slic_seq.l[i]))
  tmp_count <- countPWM(tct,tmp_seq,min.score = "80%")
  if (tmp_count > 0) {
    sum=sum+1
    TCT_sy <- c(TCT_sy,simple)
  }
  simple=simple+1
}
TCT <- sum/length(domTSS_slic_seq.l)
#plot
TCT_topie_name <- c('TCT','not TCT')
TCT_topie_percent <- c(TCT,1-TCT)
TCT_topie_percent <- round(TCT_topie_percent*100,2)
TCT_topie_name <- paste(TCT_topie_name," | ",TCT_topie_percent,"%",sep = "")
png(paste0(name_s,"_TCT_pie.png"),height = 800,width = 800)
pie(TCT_topie_percent,labels = TCT_topie_name,main=paste0(name_s,"_TCT_pie"),radius = 0.8,
    cex=1.0,clockwise = TRUE)
dev.off()


TBP_a_TCT <- union(TBP_sy,TCT_sy)
others_sy <- setdiff(all_sy,TBP_a_TCT)
all <- list(TBP = TBP_sy,TCT=TCT_sy,others=others_sy)
p1 <- upset(fromList(all),text.scale = c(2,5,2,2,4,5))
png(paste0(name_s,"_upset.png"),width = 800,height = 700)
p1
dev.off()


