.libPaths(c("/home/shenlm/R/x86_64-pc-linux-gnu-library/4.1","/opt/R/4.1.2/lib/R/library"))
#output:1.width_picture. 2
#library
library(CAGEr)
library(CAGEr)
library(ggplot2)
library(edgeR)
library(tidyr)
library(RColorBrewer)
library(UpSetR)
library(ggplot2)
library(grid)
library(dplyr)
library(getopt)
library(ggbiplot)
library(ChIPseeker)
library(biomaRt)
library(GenomicFeatures)

spec = matrix(c('species','s',1,"character",
                'out_dir','o',1,"character",
                'bam_dir','b',1,"character",
                'pwm_dir','p',1,"character"),
              byrow = TRUE,ncol=4)
opt=getopt(spec)

setwd(opt$out_dir)
species <- opt$species
inputDir1 <- opt$bam_dir
pwm_dir <- opt$pwm_dir
TBM_score <- 10

paths1 <- list.files(inputDir1,full.names = T)
pathsToInputFiles1 <- paths1[grep(pattern = "*_sorted.bam", paths1)]
pathsToInputFiles <- c(pathsToInputFiles1)
promoter_c="TRUE"

samples <- sub("_sorted.bam","",basename(pathsToInputFiles))
names <- sub("_sorted.bam","",basename(pathsToInputFiles))
samples <- gsub("x1_","",samples)
samples <- gsub("x2_","",samples)
names <- gsub("x1_","",names)
names <- gsub("x2_","",names)

#----------------motif in-------------------------------------
TBP <-read.csv(paste0(pwm_dir,"/tbp.txt"),header = TRUE)
TBP <- t(TBP)
rownames(TBP) <- c('A','C','G','T')
tct <- read.csv(paste0(pwm_dir,"/human_TCT.csv"),header = TRUE)
tct <- as.matrix(tct)

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
seqlevels(txdb) <- gsub("chr","",seqlevels(txdb))
#quality control
LICAGE <- new("CAGEset", genomeName = BSname, 
              inputFiles = pathsToInputFiles, inputFilesType = "bam", 
              sampleLabels = names)
getCTSS(LICAGE)


normalizeTagCount(LICAGE, method = "powerLaw", fitInRange = c(100, 1000), alpha = 1.14, T = 1000000)
exportCTSStoBedGraph(LICAGE, values = "normalized", format = "bedGraph", oneFile = FALSE)
clusterCTSS(object = LICAGE, threshold = 1, thresholdIsTpm = TRUE, nrPassThreshold = 1, 
            method = "distclu", maxDist = 20, removeSingletons = TRUE, keepSingletonsAbove = 5)
cumulativeCTSSdistribution(LICAGE, clusters = "tagClusters")
quantilePositions(LICAGE, clusters = "tagClusters", qLow = 0.1, qUp = 0.9)
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
# save tc.grl as RDS object
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
tc_selected <- tc.grl[samples]
names(tc_selected) <- names
peakAnno_list <- lapply(tc_selected, function(x) annotatePeak(x, TxDb = txdb,  tssRegion = c(-1000, 1000), annoDb = org, sameStrand = TRUE, verbose = FALSE))
names(peakAnno_list) <- names

tmp <- list()
for (i in 1:length(domTSS.grl)) {
  ttmp <- as.data.frame(peakAnno_list[[i]])
  tmp[[i]] <- ttmp[which(ttmp$annotation=="Promoter"),]
}

#--------------------------------cager_genomic_picture--------------------------------------
feats.l <- lapply(peakAnno_list, function(x) x@annoStat)
names(feats.l) <- names
for (i in 1:length(names)) {
  feats.l[[i]]$sample <- rep(names[[i]], nrow(feats.l[[i]]))
}
feats.df <- do.call("rbind", feats.l)
feats.df$sample <- factor(feats.df$sample, levels = names[length(names):1])

domTSS_slic_mergLandR.grl <- list()
for (i in 1:length(domTSS.grl)) {
  a_dom <- subset(tmp[[i]],select=c(dominant_ctss))
  a_all <- merge(domTSS.grl[[i]],a_dom,by="dominant_ctss")
  domTSS_slic_mergLandR.grl[[i]] <- GRanges(seqnames = a_all$seqnames,
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
}
names(domTSS_slic_mergLandR.grl) <- names
#------------------------------------------------------------motif-scan-read-in-----------------------------
#chose the TPM > score
for (i in 1:length(domTSS_slic_mergLandR.grl)) {
  tpm_tmp <- domTSS_slic_mergLandR.grl[[i]]
  domTSS_slic_mergLandR.grl[[i]] <- tpm_tmp[which(tpm_tmp$tpm.dominant_ctss>TBM_score),]
}

#different TPM to choose
diff_tmp_1 <- as.data.frame(domTSS_slic_mergLandR.grl[[1]])
diff_tmp_2 <- as.data.frame(domTSS_slic_mergLandR.grl[[2]])
diff_tmp_1 <- unite(diff_tmp_1,"symbol",c("seqnames","start","strand"),sep = "_",remove = FALSE)
new_raw_1 <- diff_tmp_1$symbol
diff_tmp_2 <- unite(diff_tmp_2,"symbol",c("seqnames","start","strand"),sep = "_",remove = FALSE)
new_raw_2 <- diff_tmp_2$symbol
sy_1 <- setdiff(new_raw_1,new_raw_2)
sy_2 <- setdiff(new_raw_2,new_raw_1)
sy_1 <- as.data.frame(sy_1)
sy_2 <- as.data.frame(sy_2)
colnames(sy_1) <- "symbol"
colnames(sy_2) <- "symbol"
diff_table <- list()
diff_table[[1]] <- merge(diff_tmp_1,sy_1,by="symbol")
diff_table[[2]] <- merge(diff_tmp_2,sy_2,by="symbol")
for (merge_sam in 1:length(domTSS_slic_mergLandR.grl)) {
  domTSS_slic_mergLandR.grl[[merge_sam]] <- GRanges(seqnames = diff_table[[merge_sam]]$seqnames,
                                                    ranges = IRanges(start = diff_table[[merge_sam]]$dominant_ctss, 
                                                                     end = diff_table[[merge_sam]]$dominant_ctss),
                                                    strand = diff_table[[merge_sam]]$strand,
                                                    nr_ctss = diff_table[[merge_sam]]$nr_ctss,
                                                    dominant_ctss = diff_table[[merge_sam]]$dominant_ctss,
                                                    tpm = diff_table[[merge_sam]]$tpm, 
                                                    tpm.dominant_ctss = diff_table[[merge_sam]]$tpm.dominant_ctss, 
                                                    q_0.1 = diff_table[[merge_sam]]$q_0.1,
                                                    q_0.9 = diff_table[[merge_sam]]$q_0.9,
                                                    interquantile_width = diff_table[[merge_sam]]$interquantile_width,
                                                    tc_start = diff_table[[merge_sam]]$start,
                                                    tc_end = diff_table[[merge_sam]]$end)
}

#

#Figures 2B and F - distribution of tag cluster/promoter interquantile widths
iq.l <- lapply(domTSS_slic_mergLandR.grl, function(x) data.frame(iq_width = x$interquantile_width))
for (i in 1:length(names)) {
  iq.l[[i]]$sample <- rep(names[[i]], nrow(iq.l[[i]]))
}
iq.df <- do.call(rbind, iq.l)
iq.df$sample <- factor(iq.df$sample, levels = names)
iq_all.l <- lapply(tc.grl, function(x) data.frame(iq_width = x$interquantile_width))
names(iq_all.l) <- names
for (i in 1:length(names)) {
  iq_all.l[[i]]$sample <- rep(names[[i]], nrow(iq_all.l[[i]]))
}
iq_all.df <- do.call(rbind, iq_all.l)
iq_all.df$sample <- factor(iq_all.df$sample, levels = names)
# plotting
p <- ggplot(iq_all.df, aes(iq_width)) +
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
  labs(x = "Interquantile width", y = "percent") +
  xlim(0, 50)
pdf(paste0(names[1]," and ",names[2],"_after.pdf"), height = 4, width = 15)
p + facet_wrap(~ sample, ncol = 5)
dev.off()

#write tss
for (i in 1:length(domTSS_slic_mergLandR.grl)) {
  tmp <- as.data.frame(domTSS_slic_mergLandR.grl[[i]])
  write.csv(tmp,paste0(names[i],"_tss.csv"))
}

#TATA_TCT scan
domTSS_slic_win.grl <- list()
domTSS_slic_seq.l <- list()
for (i in 1:length(domTSS_slic_mergLandR.grl)) {
  up <- 40
  down <- 0
  range <- c(-up, down)
  win <- up + down
  domTSS_slic_win.grl[[i]] <- promoters(domTSS_slic_mergLandR.grl[[i]],up=up,down=down)
  scan_tmp <- domTSS_slic_win.grl[[i]]
  domTSS_slic_win.grl[[i]] <- scan_tmp[width(trim(domTSS_slic_win.grl[[i]]))==win]
  domTSS_slic_seq.l[[i]] <- getSeq(BS,domTSS_slic_win.grl[[i]])
  sum=0
  simple=1
  TBP_sy <- c()
  all_sy <- c()
  for (t in 1:length(domTSS_slic_seq.l[[i]])) {
    tmp_seq_file <- domTSS_slic_seq.l[[i]]
    tmp_seq <- DNAString(as.character(tmp_seq_file[t]))
    tmp_count <- countPWM(TBP,tmp_seq,min.score = "80%")
    if (tmp_count > 0) {
      sum=sum+1
      TBP_sy <- c(TBP_sy,simple)
    }
    simple=simple+1
    all_sy <- c(all_sy,simple)
  }
  tbp <- sum/length(domTSS_slic_seq.l[[i]])
  #plot
  TBP_topie_name <- c('TATA-box','not TATA-box')
  TBP_topie_percent <- c(tbp,1-tbp)
  TBP_topie_percent <- round(TBP_topie_percent*100,2)
  TBP_topie_name <- paste(TBP_topie_name," | ",TBP_topie_percent,"%",sep = "")
  png(paste0(names[i],"_TATA-box_pie.png"),height = 800,width = 800)
  pie(TBP_topie_percent,labels = TBP_topie_name,main=paste0(names[i],"_TATA-box_pie"),radius = 0.8,
      cex=1.0,clockwise = TRUE)
  dev.off()
  
  up <- 10
  down <- 10
  range <- c(-up, down)
  win <- up + down
  domTSS_slic_win.grl[[i]] <- promoters(domTSS_slic_mergLandR.grl[[i]],up=up,down=down)
  scan_tmp <- domTSS_slic_win.grl[[i]]
  domTSS_slic_win.grl[[i]] <- scan_tmp[width(trim(domTSS_slic_win.grl[[i]]))==win]
  domTSS_slic_seq.l[[i]] <- getSeq(BS,domTSS_slic_win.grl[[i]])
  
  TCT_sy <- c()
  sum=0
  simple=1
  for (t in 1:length(domTSS_slic_seq.l[[i]])) {
    tmp_seq_file <- domTSS_slic_seq.l[[i]]
    tmp_seq <- DNAString(as.character(tmp_seq_file[t]))
    tmp_count <- countPWM(tct,tmp_seq,min.score = "80%")
    if (tmp_count > 0) {
      sum=sum+1
      TCT_sy <- c(TCT_sy,simple)
    }
    simple=simple+1
  }
  TCT <- sum/length(domTSS_slic_seq.l[[i]])
  #plot
  TCT_topie_name <- c('TCT','not TCT')
  TCT_topie_percent <- c(TCT,1-TCT)
  TCT_topie_percent <- round(TCT_topie_percent*100,2)
  TCT_topie_name <- paste(TCT_topie_name," | ",TCT_topie_percent,"%",sep = "")
  png(paste0(names[i],"_TCT_pie.png"),height = 800,width = 800)
  pie(TCT_topie_percent,labels = TCT_topie_name,main=paste0(names[i],"_TCT_pie"),radius = 0.8,
      cex=1.0,clockwise = TRUE)
  dev.off()
  
  TBP_a_TCT <- union(TBP_sy,TCT_sy)
  others_sy <- setdiff(all_sy,TBP_a_TCT)
  all <- list(TBP = TBP_sy,TCT=TCT_sy,others=others_sy)
  if (i==1) {
    p1 <- upset(fromList(all),text.scale = c(2,5,2,2,4,5))
    dev.off()   
  }else if (i==2) {
    p2 <- upset(fromList(all),text.scale = c(2,5,2,2,4,5))
  }
}

png(paste0(names[1],"_upset.png"),width = 800,height = 700)
p1
dev.off()
png(paste0(names[2],"_upset.png"),width = 800,height = 700)
p2
dev.off()

domTSS_slic_win.grl <- list()
domTSS_slic_seq.l <- list()
for (i in 1:length(domTSS_slic_mergLandR.grl)) {
  up <- 100
  down <- 100
  range <- c(-up, down)
  win <- up + down
  domTSS_slic_win.grl[[i]] <- promoters(domTSS_slic_mergLandR.grl[[i]],up=up,down=down)
  scan_tmp <- domTSS_slic_win.grl[[i]]
  domTSS_slic_win.grl[[i]] <- scan_tmp[width(trim(domTSS_slic_win.grl[[i]]))==win]
  domTSS_slic_seq.l[[i]] <- getSeq(BS,domTSS_slic_win.grl[[i]])
  tomotif_tmp <- as.data.frame(domTSS_slic_win.grl[[i]])
  for (line in 1:nrow(tomotif_tmp)) {
    tomotif_tmp[line,15] <- paste0("peak_",line)
  }
  tomotif_tmp <- select_(tomotif_tmp,15,1,2,3,5)
  write.csv(tomotif_tmp,paste0(names[i],"_sx100.csv"))
}







