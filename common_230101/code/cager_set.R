#set
.libPaths(c("/home/shenlm/R/x86_64-pc-linux-gnu-library/4.1","/opt/R/4.1.2/lib/R/library"))
library(CAGEr)
setwd("/data1/shenluemou/ing_data/sk/220927_nss")
species <- "mouse"
inputDir1 <- "bam/"
paths1 <- list.files(inputDir1,full.names = T)
pathsToInputFiles1 <- paths1[grep(pattern = "*_sorted.bam", paths1)]
pathsToInputFiles <- c(pathsToInputFiles1)
#edger: 是否使用差异分析
#merge_c: 是否进行样本merge
#promoter_c : 是否取出注释为promoter的区域
#de_c: 是否删除某些特定的基因
edger="TRUE"
merge_c="TRUE"
promoter_c="TRUE"
de_c="FALSE"
tss_c="FALSE"
if (tss_c) {
  tss <- read.csv("")
}

samples <- sub("_sorted.bam","",basename(pathsToInputFiles))
names <- sub("_sorted.bam","",basename(pathsToInputFiles))
if (de_c) {
  qu_gene <- read.csv("qu_MIIgene.csv",header = TRUE)
  qu_gene$X <- NULL
}

#-----------select---------------------------
if (merge_c) {
  merge_index <- c(1,2)
  merge_sample <- c("1c","2c")
  merge_sam <- "1"
}

if (edger) {
  group <- merge_index#2-1
  select_sample="1"
  merge_sam <- select_sample
}

#library
if (species == "human") {
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  seqnames(BSgenome.Hsapiens.UCSC.hg38) <- gsub("chr","",seqnames(BSgenome.Hsapiens.UCSC.hg38))
  BSname <- "BSgenome.Hsapiens.UCSC.hg38"
  BS <- BSgenome.Hsapiens.UCSC.hg38
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  org <- "org.Mm.eg.db"
}else if (species == "mouse") {
  library(BSgenome.Mmusculus.UCSC.mm10)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(org.Mm.eg.db)
  seqnames(BSgenome.Mmusculus.UCSC.mm10) <- gsub("chr","",seqnames(BSgenome.Mmusculus.UCSC.mm10))
  BSname <- "BSgenome.Mmusculus.UCSC.mm10"
  BS <- BSgenome.Mmusculus.UCSC.mm10
  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  org <- "org.Hs.eg.db"
}
library(CAGEr)
library(ggplot2)
library(edgeR)
library(tidyr)
library(RColorBrewer)
library(UpSetR)
library(ggplot2)
library(grid)
library(dplyr)

#quality control
LICAGE <- new("CAGEset", genomeName = BSname, 
              inputFiles = pathsToInputFiles, inputFilesType = "bam", 
              sampleLabels = sub("_sorted.bam","",basename(pathsToInputFiles)))
getCTSS(LICAGE)


#----------------------------------------------------edger-------------------------------
if (edger) {
  raw_count <- data.frame(LICAGE@CTSScoordinates,LICAGE@tagCountMatrix)
  new_raw <- unite(raw_count,"symbol",c("chr","pos","strand"),sep = "_",remove = TRUE)
  test <- new_raw
  rownames(test) <- test[,1]
  test$symbol <- NULL
  y <- DGEList(counts = test[,1:length(group)],genes = rownames(test),group = group)
  keep <- rowSums(cpm(y)>1) >= 1
  y <- y[keep, , keep.lib.sizes=FALSE]
  y <- calcNormFactors(y)
  y_bcv <- y
  bcv <- 0.1
  et <- exactTest(y_bcv,dispersion = bcv^2)
  gene1 <- decideTestsDGE(et,p.value = 0.05,lfc = 0)
  colnames(gene1) <- "Signifi"
  results <- cbind(y$genes,y$counts,et$table,gene1)
  if (select_sample=="2") {
    different <- results[which(results$Signifi==1),]
  }else if (select_sample=="1") {
    different <- results[which(results$Signifi==-1),]
  }
  count <- separate(data = different,col = genes,into = c("chr","pos","strand"),sep = "_")
  if (de_c) {
    count_to_tss <- subset(count,select=c(chr,pos,strand))
    count_to_tss$pos <- as.numeric(count_to_tss$pos)
    count_to_tss <- setdiff(count_to_tss,qu_gene)
    tmp_count <- merge(count_to_tss,count,by=c("chr","pos","strand"))
    count_to_taq_count <- subset(tmp_count,select=samples)
    rownames(count_to_taq_count) <- NULL
    rownames(count_to_tss) <- NULL    
  }else{
    count_to_taq_count <- subset(count,select=samples)
    count_to_tss <- subset(count,select=c(chr,pos,strand))
    count_to_tss$pos <- as.numeric(count_to_tss$pos)
    rownames(count_to_taq_count) <- NULL
    rownames(count_to_tss) <- NULL    
  }
  LICAGE@tagCountMatrix <- count_to_taq_count
  LICAGE@CTSScoordinates <- count_to_tss
  LICAGE@inputFiles <- LICAGE@sampleLabels
  LICAGE@sampleLabels <- LICAGE@sampleLabels
  LICAGE@librarySizes <- LICAGE@librarySizes
  normalizeTagCount(LICAGE, method = "powerLaw", fitInRange = c(100, 1000), alpha = 1.14, T = 1000000)
}
if (!edger&&de_c) {
  count <- data.frame(LICAGE@CTSScoordinates,LICAGE@tagCountMatrix)
  count_to_tss <- as.data.frame(LICAGE@CTSScoordinates)
  count_to_tss <- setdiff(count_to_tss,qu_gene)
  tmp_count <- merge(count_to_tss,count,by=c("chr","pos","strand"))
  count_to_taq_count <- subset(tmp_count,select=samples)
  rownames(count_to_taq_count) <- NULL
  rownames(count_to_tss) <- NULL
  LICAGE@tagCountMatrix <- count_to_taq_count
  LICAGE@CTSScoordinates <- count_to_tss
  LICAGE@inputFiles <- LICAGE@sampleLabels
  LICAGE@sampleLabels <- LICAGE@sampleLabels
  LICAGE@librarySizes <- LICAGE@librarySizes
  normalizeTagCount(LICAGE, method = "powerLaw", fitInRange = c(100, 1000), alpha = 1.14, T = 1000000)
}

#merge
if (merge_c) {
  mergeSamples(LICAGE,mergeIndex = merge_index,mergedSampleLabels = merge_sample)
  normalizeTagCount(LICAGE, method = "powerLaw", fitInRange = c(100, 1000), alpha = 1.14, T = 1000000)
  #plotReverseCumulatives(LICAGE, fitInRange = c(100, 1000), onePlot = TRUE, values = "normalized")
  names <- merge_sample
  samples <- merge_sample
}

normalizeTagCount(LICAGE, method = "powerLaw", fitInRange = c(100, 1000), alpha = 1.14, T = 1000000)
exportCTSStoBedGraph(LICAGE, values = "normalized", format = "bedGraph", oneFile = FALSE)
clusterCTSS(object = LICAGE, threshold = 1, thresholdIsTpm = TRUE, nrPassThreshold = 1, 
            method = "distclu", maxDist = 20, removeSingletons = TRUE, keepSingletonsAbove = 5)
cumulativeCTSSdistribution(LICAGE, clusters = "tagClusters")
quantilePositions(LICAGE, clusters = "tagClusters", qLow = 0.1, qUp = 0.9)
saveRDS(LICAGE, file = "LICAGE.RDS")

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
# save tc.grl as RDS object
saveRDS(object = tc.grl, file = "tagClusters_list.RDS")
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

# save dominant TSS cluster GRanges object for later
saveRDS(domTSS.grl, file ="domTSS_grl.RDS")

#------------------------------------------------------------------------------
#03_genomic_location_tc.R
library(ChIPseeker)
library(biomaRt)
library(GenomicFeatures)

# - rename chromosome names in txdb to match NCBI naming scheme - #
seqlevels(txdb) <- gsub("chr","",seqlevels(txdb))
tc_selected <- tc.grl[samples]
names(tc_selected) <- names
peakAnno_list <- lapply(tc_selected, function(x) annotatePeak(x, TxDb = txdb,  tssRegion = c(-1000, 1000), annoDb = org, sameStrand = TRUE, verbose = FALSE))
names(peakAnno_list) <- names


#select the promoter
if (promoter_c) {
  tmp <- as.data.frame(peakAnno_list[[as.numeric(merge_sam)]])
  tmp <- tmp[which(tmp$annotation=="Promoter"),]
}


#--------------------------------cager_genomic_picture--------------------------------------
feats.l <- lapply(peakAnno_list, function(x) x@annoStat)
names(feats.l) <- names
for (i in 1:length(names)) {
  feats.l[[i]]$sample <- rep(names[[i]], nrow(feats.l[[i]]))
}
feats.df <- do.call("rbind", feats.l)
feats.df$sample <- factor(feats.df$sample, levels = names[length(names):1])
write.csv(feats.df,file = "features.csv",row.names = FALSE)

col1 <- brewer.pal(5, "Greys")[3:5]
col2 <- brewer.pal(8, "Blues")[8:6]
col <- c(col2, col1)
# new colour scheme
col <- c("#F8F2AB", "#B4D6A4", "#67BC9A", "#13B0A5", "#0071A7", "#3E606F", "#88F9D4", "#18C29C",
         "#0B877D", "#126872", "#031727")
# set feature factors
feature_factors <- c("Promoter", "Promoter (<=1kb)", "Promoter (1-3kb)", "5' UTR", "1st Exon", "Other Exon", 
                     "1st Intron", "Other Intron", "3' UTR", "Downstream (<=3kb)", "Distal Intergenic")
names(col) <- feature_factors
col_sel <- col[names(col) %in% unique(feats.df$Feature)]
# plot genomic features
library(ggplot2)
p <- ggplot(feats.df, aes(x = sample, y = Frequency, fill = Feature), alpha = 0.7) +
  geom_bar(stat = "identity", width = 0.75, colour = "black", lwd = 0.125) +
  coord_flip() +
  scale_fill_manual("Features", values = col_sel) +
  theme_bw() +
  theme(text = element_text(size = 14, colour = "black"),
        legend.title = element_blank(),
        axis.title.x = element_text(colour = "black"),
        axis.title.y = element_text(colour = "black"),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(y = "Percentage", x = NULL)

pdf(file = "genomicFeatures_sc_all.pdf", height = 3, width = 5)
print(p)
dev.off()

#-----------------------------------
#Figures 2B and F - distribution of tag cluster/promoter interquantile widths
# extract interquantile widths
iq.l <- lapply(tc_selected, function(x) data.frame(iq_width = x$interquantile_width))

# add names to as a column to each sample
for (i in 1:length(names)) {
  iq.l[[i]]$sample <- rep(names[[i]], nrow(iq.l[[i]]))
}

# combine to 1 dataframe
iq.df <- do.call(rbind, iq.l)

# set levels for plotting - sample levels
iq.df$sample <- factor(iq.df$sample, levels = names)

# plotting
library(ggplot2)
p <- ggplot(iq.df, aes(iq_width)) +
  geom_histogram(aes(x = iq_width, y = ..ncount..* 100), binwidth = 1, 
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
  labs(x = "Interquantile width", y = "Percentage") +
  coord_equal(ratio = 1) +
  xlim(0, 150)

pdf(file = "IQ_width_sel_sc.pdf", height = 6, width = 6)
p + facet_wrap(~ sample, ncol = 2)
dev.off()


# plot all interquantile widths for supplementary
iq_all.l <- lapply(tc.grl, function(x) data.frame(iq_width = x$interquantile_width))

# set names for plotting

names(iq_all.l) <- names

# add names to as a column to each sample
for (i in 1:length(names)) {
  iq_all.l[[i]]$sample <- rep(names[[i]], nrow(iq_all.l[[i]]))
}

# combine to 1 dataframe
iq_all.df <- do.call(rbind, iq_all.l)

# set levels for plotting - sample levels
iq_all.df$sample <- factor(iq_all.df$sample, levels = names)

# plot iq-widths for all samples
# plotting
library(ggplot2)
p <- ggplot(iq_all.df, aes(iq_width)) +
  geom_histogram(aes(x = iq_width, y = ..ncount..*100), binwidth = 1, 
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
  labs(x = "Interquantile width", y = "count") +
  #coord_equal(ratio = 1) +
  xlim(0, 50)

pdf(file = "IQ_width_all_slic.pdf", height = 4, width = 15)
p + facet_wrap(~ sample, ncol = 5)
dev.off()

#------------------------------------------------------------motif-scan-read-in-----------------------------
domTSS_slic_mergLandR.grl <- readRDS(file = "domTSS_grl.RDS")

if (promoter_c) {
  a_dom <- subset(tmp,select=c(dominant_ctss))
  a_all <- merge(domTSS_slic_mergLandR.grl[[as.numeric(merge_sam)]],a_dom,by="dominant_ctss")
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
}

#lapply(domTSS_slic_mergLandR.grl,function(x) sum(x$interquantile_width <= 3)/length(x$interquantile_width))


#----------------motif in-------------------------------------
TBP <-read.csv("pwm/tbp.txt",header = TRUE)
TBP <- t(TBP)
rownames(TBP) <- c('A','C','G','T')
tct <- read.csv("pwm/human_TCT.csv",header = TRUE)
tct <- as.matrix(tct)



#---------------motif scan------------------------------------------------
if (promoter_c) {
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
  for (i in 1:length(domTSS_slic_seq.l)) {
    tmp_seq <- DNAString(as.character(domTSS_slic_seq.l[i]))
    tmp_count <- countPWM(TBP,tmp_seq,min.score = "80%")
    if (tmp_count > 0) {
      sum=sum+1
      TBP_sy <- c(TBP_sy,simple)
    }
    simple=simple+1
  }
  tbp <- sum/length(domTSS_slic_seq.l)
  
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
  
  all <- list(TBP = TBP_sy,TCT=TCT_sy)
  upset(fromList(all),text.scale = c(2,5,2,2,4,5))
  
  up <- 100
  down <- 100
  range <- c(-up, down)
  win <- up + down
  domTSS_slic_win.grl <- promoters(domTSS_slic_mergLandR.grl,up=up,down=down)
  domTSS_slic_win.grl <- domTSS_slic_win.grl[width(trim(domTSS_slic_win.grl))==win]
  domTSS_slic_seq.l <- getSeq(BS,domTSS_slic_win.grl)
}


#---------------------------------------------------------------
if (!promoter_c) {
  up <-40
  down <- 0
  range <- c(-up, down)
  win <- up + down
  domTSS_slic_win.grl <- sapply(domTSS_slic_mergLandR.grl, function(x) promoters(x, up = up, down = down))
  domTSS_slic_win.grl <- sapply(domTSS_slic_win.grl, function(x) x[width(trim(x)) == win])
  sapply(domTSS_slic_win.grl, length)
  domTSS_slic_seq.l <- sapply(domTSS_slic_win.grl,function(x) getSeq(BS, x))    
  motif <- list(tct)
  tomo_name <- names(domTSS_seq.l)
  tbp <- list()
  for (mo in 1:length(motif)) {
    for (nam in 1:length(tomo_name)) {
      tmp_dom <- domTSS_slic_seq.l[[nam]]
      sum=0
      for (i in 1:length(tmp_dom)) {
        tmp_seq <- DNAString(as.character(tmp_dom[i]))
        tmp_count <- countPWM(motif[[mo]],tmp_seq,min.score = "80%")
        if (tmp_count > 0) {
          sum=sum+1
        }
      }
      TCT <- c(TCT,sum/length(tmp_dom))
    }
  }
  
  
  up <-10
  down <- 10
  range <- c(-up, down)
  win <- up + down
  domTSS_slic_win.grl <- sapply(domTSS_slic_mergLandR.grl, function(x) promoters(x, up = up, down = down))
  domTSS_slic_win.grl <- sapply(domTSS_slic_win.grl, function(x) x[width(trim(x)) == win])
  sapply(domTSS_slic_win.grl, length)
  domTSS_slic_seq.l <- sapply(domTSS_slic_win.grl,function(x) getSeq(BS, x))    
  motif <- list(TBP)
  tomo_name <- names(domTSS_seq.l)
  tbp <- list()
  for (mo in 1:length(motif)) {
    for (nam in 1:length(tomo_name)) {
      tmp_dom <- domTSS_slic_seq.l[[nam]]
      sum=0
      for (i in 1:length(tmp_dom)) {
        tmp_seq <- DNAString(as.character(tmp_dom[i]))
        tmp_count <- countPWM(motif[[mo]],tmp_seq,min.score = "80%")
        if (tmp_count > 0) {
          sum=sum+1
        }
      }
      tbp <- sum/length(domTSS_slic_seq.l)
    }
  }
}

