.libPaths(c("/home/shenlm/R/x86_64-pc-linux-gnu-library/4.1","/opt/R/4.1.2/lib/R/library"))

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

#read
spec = matrix(c('work_dir','w',1,"character",
                'species','s',1,"character",
                'bam_dir','b',1,"character"),
              byrow = TRUE,ncol=4)
opt=getopt(spec)


#in
setwd(opt$work_dir)
species <- opt$species 
merge_c <- FALSE
inputDir1 <- opt$bam_dir

paths1 <- list.files(inputDir1,full.names = T)
pathsToInputFiles1 <- paths1[grep(pattern = "*_sorted.bam", paths1)]
pathsToInputFiles <- c(pathsToInputFiles1)
samples <- sub("_sorted.bam","",basename(pathsToInputFiles))
names <- sub("_sorted.bam","",basename(pathsToInputFiles))


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
}

#quality control
LICAGE <- new("CAGEset", genomeName = BSname, 
              inputFiles = pathsToInputFiles, inputFilesType = "bam", 
              sampleLabels = sub("_sorted.bam","",basename(pathsToInputFiles)))
getCTSS(LICAGE)

normalizeTagCount(LICAGE, method = "powerLaw", fitInRange = c(100, 1000), alpha = 1.14, T = 1000000)
exportCTSStoBedGraph(LICAGE, values = "normalized", format = "bedGraph", oneFile = FALSE)
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
# - rename chromosome names in txdb to match NCBI naming scheme - #
seqlevels(txdb) <- gsub("chr","",seqlevels(txdb))
tc_selected <- tc.grl[samples]
names(tc_selected) <- names
peakAnno_list <- lapply(tc_selected, function(x) annotatePeak(x, TxDb = txdb,  tssRegion = c(-1000, 1000), annoDb = org, sameStrand = TRUE, verbose = FALSE))
names(peakAnno_list) <- names

tmp <- list()
for (i in 1:length(domTSS.grl)) {
  ttmp <- as.data.frame(peakAnno_list[[i]])
  tmp[[i]] <- ttmp[which(ttmp$annotation=="Promoter"),]
}
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
#----------------------csv_tss------------------------------------------------
for (i in 1:length(domTSS_slic_mergLandR.grl)) {
  tmp <- as.data.frame(domTSS_slic_mergLandR.grl[[i]]) 
  write.csv(tmp,paste0(names[i],"_cager_tss.csv"))
}
write.csv(names,"names.csv")
