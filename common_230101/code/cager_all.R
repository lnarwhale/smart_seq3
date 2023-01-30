
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

#-----------select---------------------------
if (merge_c) {
  merge_index <- c(5,1,2,4,3)
  merge_sample <- c("fourC","MII","onec","threeC","twoC")
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
saveRDS(LICAGE, file = "LICAGE.RDS")

#corr
png("correlation2.png",width = length(names)*200,height = length(names)*200)
plotCorrelation2( LICAGE, samples = "all"
                  , tagCountThreshold = 1, applyThresholdBoth = FALSE
                  , method = "pearson")
dev.off()

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

for (i in 1:length(domTSS.grl)) {
  tmp <- as.data.frame(domTSS.grl[[i]]) 
  write.csv(tmp,paste0(names[i],"_cager_tss.csv"))
}

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



#--------------------------------cager_genomic_picture--------------------------------------
feats.l <- lapply(peakAnno_list, function(x) x@annoStat)
names(feats.l) <- names
for (i in 1:length(names)) {
  feats.l[[i]]$sample <- rep(names[[i]], nrow(feats.l[[i]]))
}
feats.df <- do.call("rbind", feats.l)
feats.df$sample <- factor(feats.df$sample, levels = names[length(names):1])

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
  #coord_equal(ratio = 1) +
  xlim(0, 50)
                                            
pdf(file = "IQ_width_all_slic.pdf", height = 4, width = 15)
p + facet_wrap(~ sample, ncol = 5)
dev.off()
