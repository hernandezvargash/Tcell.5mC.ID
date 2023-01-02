
# 151222
# remora model:
# "dna_r9.4.1_450bps_fast.cfg"
# mapped to GRCm38.p6.genome.fa

# libraries ---------------------------------------------------------------

rm(list=ls())

suppressPackageStartupMessages({
  library(rtracklayer)
  library(vioplot)
  library(reshape)
  library(ggplot2)
  library(gridExtra)
  library(BSgenome.Mmusculus.UCSC.mm10)
  library(viridis)
  library(DSS)
  library(minfi)
  library(stringr)
  library(bumphunter)
  library(bsseq)
  library(SummarizedExperiment)
  library(xlsx)
  library(dmrseq)
  library(DMRcate)
  library(biomaRt)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(ggfortify)
  library(pals)
  library(scales)
  library(EnrichedHeatmap)
  library(circlize)
  library(viridis)

})

setwd("~/Dropbox/BioInfo/Lab/Tcell_ID")

set.seed(123)


# load bed files ----------------------------------------------------------

genome <- BSgenome.Mmusculus.UCSC.mm10

bed.files <- list.files(path = "./remora/mouse_brain/", pattern = ".bed", recursive = T, full.names = T)

bed.list <- list()

for(i in 1:length(bed.files)){
  bed.list[[i]] <- import.bed(bed.files[i], which = NULL) # used "regions" or NULL
  bed.list[[i]] <- keepStandardChromosomes(bed.list[[i]], pruning.mode = "coarse")
  bed.list[[i]] <- sort(bed.list[[i]])
  bed.list[[i]] <- resize(bed.list[[i]], width = 3, fix = "start")
  bed.list[[i]]$context <- getSeq(genome, bed.list[[i]])
  bed.list[[i]]$blockSizes <- as.numeric(as.character(bed.list[[i]]$blockSizes))
}

names(bed.list) <- c("brain.5hmC", "brain.5mC")
save(bed.list, file="manuscript/bed.mouse.brain.RData")



# inspection --------------------------------------------------------------

load("manuscript/bed.mouse.brain.RData")

par(mar = c(10,15,5,15), mfrow = c(2,1))
barplot(as.numeric(lapply(bed.list, length)), las = 2,
        names.arg = names(bed.list),
        main = "Total reads per sample",
        col = viridis(4))

barplot(as.numeric(lapply(bed.list, function(x) sum(x$blockCount))), las = 2,
        names.arg = names(bed.list),
        main = "Total counts per sample",
        col = viridis(4))


# methylation distribution

lapply(bed.list, length)
# 5551760

block.sizes <- lapply(bed.list, function(x) as.numeric(x$blockSizes))
names(block.sizes)

block.sizes.df <- data.frame(Reduce(cbind, block.sizes))
colnames(block.sizes.df) <- names(block.sizes)
block.sizes.df <- as.matrix(block.sizes.df)
block.sizes.df <- block.sizes.df/100 # to get beta values

head(block.sizes.df)
tail(block.sizes.df)

hist(block.sizes.df[,1], main = "mouse brain 5hmC", col = "black", xlab = "beta values")
hist(block.sizes.df[,2], main = "mouse brain 5mC", col = "black", xlab = "beta values")




# 5mC/5hmC near genes -------------------------------------------------------

load("manuscript/bed.mouse.brain.RData")

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene


# TSS

genes <- genes(TxDb.Mmusculus.UCSC.mm10.knownGene)

tss = promoters(genes, upstream = 0, downstream = 1)
tss[1:5]

tss.2 <- subsetByOverlaps(tss, bed.list[[1]])

mat.tss.5hmC = normalizeToMatrix(bed.list[[1]], tss.2, value_column = "blockSizes", mean_mode = "absolute", extend = 5000, w = 200, background = NA, smooth = F)
mat.tss.5mC = normalizeToMatrix(bed.list[[2]], tss.2, value_column = "blockSizes", mean_mode = "absolute", extend = 5000, w = 200, background = NA, smooth = F)

#meth_col_fun = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
#EnrichedHeatmap(mat.tss.5hmC, col = meth_col_fun, name = "methylation", column_title = "methylation near TSS")

ht_list1 <- EnrichedHeatmap(mat.tss.5hmC, name = "5hmC", col = c("whitesmoke",  plasma(5)[1]),
                           column_title = "5hmC in brain cells\nTSS",
                           axis_name_rot = 90,
                           top_annotation = HeatmapAnnotation(height = unit(5,"cm"),
                                                              enrich = anno_enriched(ylim = c(0, 70), 
                                                                                     gp = gpar(col =  plasma(5)[1], lty = 1, lwd = 4),
                                                                                     axis_param = list(side = "left")))) +
  EnrichedHeatmap(mat.tss.5mC, name = "5mC", col = c("whitesmoke",  plasma(5)[3]),
                  column_title = "5mC in brain cells\nTSS",
                  axis_name_rot = 90,
                  top_annotation = HeatmapAnnotation(height = unit(5,"cm"),
                                                     enrich = anno_enriched(ylim = c(0, 70), 
                                                                            gp = gpar(col = plasma(5)[3], lty = 1, lwd = 4),
                                                                            axis_param = list(side = "left"))))

draw(ht_list1, ht_gap = unit(c(10,10), "mm"))
#jpeg("enriched.heatmap.brain.jpeg", height = 960, width = 800, quality = 100)
#draw(ht_list, ht_gap = unit(c(10,10), "mm"))
#dev.off()


# gene bodies

genes.2 <- subsetByOverlaps(genes, bed.list[[1]])

mat.genes.5hmC = normalizeToMatrix(bed.list[[1]], genes.2, extend = c(5000,5000), mean_mode = "absolute", w = 200, smooth = F, target_ratio = 0.5, value_column = "blockSizes", background = NA)
mat.genes.5mC = normalizeToMatrix(bed.list[[2]], genes.2, extend = c(5000,5000), mean_mode = "absolute", w = 200, smooth = F, target_ratio = 0.5, value_column = "blockSizes", background = NA)

ht_list2 <- EnrichedHeatmap(mat.genes.5hmC, name = "5hmC", col = c("whitesmoke",  plasma(5)[1]),
                  column_title = "5hmC in brain cells\ngene bodies",
                  axis_name_rot = 90,
                  top_annotation = HeatmapAnnotation(height = unit(5,"cm"),
                                                     enrich = anno_enriched(ylim = c(0, 70), 
                                                                            gp = gpar(col =  plasma(5)[1], lty = 1, lwd = 4),
                                                                            axis_param = list(side = "left")))) +
  EnrichedHeatmap(mat.genes.5mC, name = "5mC", col = c("whitesmoke",  plasma(5)[3]),
                  column_title = "5mC in brain cells\ngene bodies",
                  axis_name_rot = 90,
                  top_annotation = HeatmapAnnotation(height = unit(5,"cm"),
                                                     enrich = anno_enriched(ylim = c(0, 70), 
                                                                            gp = gpar(col = plasma(5)[3], lty = 1, lwd = 4),
                                                                            axis_param = list(side = "left"))))

draw(ht_list2, ht_gap = unit(c(10,10), "mm"))



# CGIs

session <- browserSession("UCSC")
genome(session) <- "mm10"
query <- ucscTableQuery(session, "CpG Islands", GRangesForUCSCGenome("mm10"))
cpg_islands <- getTable(query)
cgi <- GRanges(cpg_islands)
cgi <- keepStandardChromosomes(cgi, pruning.mode="coarse") # 16009

cgi.2 <- subsetByOverlaps(cgi, bed.list[[1]])

mat.cgi.5hmC = normalizeToMatrix(bed.list[[1]], cgi.2, extend = 5000, mean_mode = "absolute", w = 200, smooth = F, target_ratio = 0.3, value_column = "blockSizes", background = NA)
mat.cgi.5mC = normalizeToMatrix(bed.list[[2]], cgi.2, extend = 5000, mean_mode = "absolute", w = 200, smooth = F, target_ratio = 0.3, value_column = "blockSizes", background = NA)

ht_list3 <- EnrichedHeatmap(mat.cgi.5hmC, name = "5hmC", col = c("whitesmoke",  plasma(5)[1]),
                           column_title = "5hmC in brain cells\nCGIs",
                           axis_name_rot = 90,
                           top_annotation = HeatmapAnnotation(height = unit(5,"cm"),
                                                              enrich = anno_enriched(ylim = c(0, 70), 
                                                                                     gp = gpar(col =  plasma(5)[1], lty = 1, lwd = 4),
                                                                                     axis_param = list(side = "left")))) +
  EnrichedHeatmap(mat.cgi.5mC, name = "5mC", col = c("whitesmoke",  plasma(5)[3]),
                  column_title = "5mC in brain cells\nCGIs",
                  axis_name_rot = 90,
                  top_annotation = HeatmapAnnotation(height = unit(5,"cm"),
                                                     enrich = anno_enriched(ylim = c(0, 70), 
                                                                            gp = gpar(col = plasma(5)[3], lty = 1, lwd = 4),
                                                                            axis_param = list(side = "left"))))

draw(ht_list3, ht_gap = unit(c(10,10), "mm"))



# brain enhancers

#enh <- read.csv("mouse.brain.enhancers.csv") # from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3660042/
#enh <- read.table("SE_11_0008_SE_mm10.bed", header = T) # from http://www.licpathway.net/sedb/index.php cortex
#enh <- read.table("SE_11_0008_SE_ele_mm10.bed", header = T) # from http://www.licpathway.net/sedb/index.php cortex
#enh <- read.delim("SE_11_0008_TE_mm10.bed", header = T) # from http://www.licpathway.net/sedb/index.php cortex
#enh <- read.delim("Mouse_Cerebral_cortex_mm10/SE_12_0294_SE_mm10.bed", header = T) # from http://www.licpathway.net/sedb/index.php cerebral cortex
#enh <- read.delim("Mouse_Cerebral_cortex_mm10_ele/SE_12_0294_SE_ele_mm10.bed", header = T) # from http://www.licpathway.net/sedb/index.php cerebral cortex
enh <- read.delim("Mouse_Forebrain_mm10/SE_12_0226_SE_mm10.bed", header = T) # from http://www.licpathway.net/sedb/index.php forebrain

colnames(enh)[1:3] <- c("chr","start","end")
head(enh)
enh <- GRanges(enh)

enh.2 <- subsetByOverlaps(enh, bed.list[[1]])

mat.enh.5hmC = normalizeToMatrix(bed.list[[1]], enh.2, extend = 5000, mean_mode = "absolute", w = 200, smooth = F, target_ratio = 0.3, value_column = "blockSizes", background = NA)
mat.enh.5mC = normalizeToMatrix(bed.list[[2]], enh.2, extend = 5000, mean_mode = "absolute", w = 200, smooth = F, target_ratio = 0.3, value_column = "blockSizes", background = NA)

ht_list4 <- EnrichedHeatmap(mat.enh.5hmC, name = "5hmC", col = c("whitesmoke",  plasma(5)[1]),
                            column_title = "5hmC in brain cells\nForebrain Enhancers",
                            axis_name_rot = 90,
                            top_annotation = HeatmapAnnotation(height = unit(5,"cm"),
                                                               enrich = anno_enriched(ylim = c(0, 70), 
                                                                                      gp = gpar(col =  plasma(5)[1], lty = 1, lwd = 4),
                                                                                      axis_param = list(side = "left")))) +
  EnrichedHeatmap(mat.enh.5mC, name = "5mC", col = c("whitesmoke",  plasma(5)[3]),
                  column_title = "5mC in brain cells\nForebrain Enhancers",
                  axis_name_rot = 90,
                  top_annotation = HeatmapAnnotation(height = unit(5,"cm"),
                                                     enrich = anno_enriched(ylim = c(0, 70), 
                                                                            gp = gpar(col = plasma(5)[3], lty = 1, lwd = 4),
                                                                            axis_param = list(side = "left"))))

draw(ht_list4, ht_gap = unit(c(10,10), "mm"))


# control enhancers

#enh <- read.table("SE_11_0004_SE_mm10.bed", header = T) # from http://www.licpathway.net/sedb/index.php adipose
#enh <- read.delim("SE_11_0004_TE_mm10.bed", header = T) # from http://www.licpathway.net/sedb/index.php adipose
enh <- read.delim("Mouse_Adipose_mm10_ele/SE_11_0004_SE_ele_mm10.bed", header = T) # from http://www.licpathway.net/sedb/index.php adipose
#enh <- read.table("SE_11_0002_SE_mm10.bed", header = T) # from http://www.licpathway.net/sedb/index.php blood
#enh <- read.table("SE_11_0002_SE_ele_mm10.bed", header = T) # from http://www.licpathway.net/sedb/index.php blood
#enh <- read.table("SE_11_0007_SE_mm10.bed", header = T) # from http://www.licpathway.net/sedb/index.php cerebellum

colnames(enh)[1:3] <- c("chr","start","end")
head(enh)
enh <- GRanges(enh)

enh.2 <- subsetByOverlaps(enh, bed.list[[1]])

mat.enh.5hmC = normalizeToMatrix(bed.list[[1]], enh.2, extend = 5000, mean_mode = "absolute", w = 200, smooth = F, target_ratio = 0.3, value_column = "blockSizes", background = NA)
mat.enh.5mC = normalizeToMatrix(bed.list[[2]], enh.2, extend = 5000, mean_mode = "absolute", w = 200, smooth = F, target_ratio = 0.3, value_column = "blockSizes", background = NA)

ht_list4 <- EnrichedHeatmap(mat.enh.5hmC, name = "5hmC", col = c("whitesmoke",  plasma(5)[1]),
                            column_title = "5hmC in brain cells\nAdipose Enhancers",
                            axis_name_rot = 90,
                            top_annotation = HeatmapAnnotation(height = unit(5,"cm"),
                                                               enrich = anno_enriched(ylim = c(0, 70), 
                                                                                      gp = gpar(col =  plasma(5)[1], lty = 1, lwd = 4),
                                                                                      axis_param = list(side = "left")))) +
  EnrichedHeatmap(mat.enh.5mC, name = "5mC", col = c("whitesmoke",  plasma(5)[3]),
                  column_title = "5mC in brain cells\nAdipose Enhancers",
                  axis_name_rot = 90,
                  top_annotation = HeatmapAnnotation(height = unit(5,"cm"),
                                                     enrich = anno_enriched(ylim = c(0, 70), 
                                                                            gp = gpar(col = plasma(5)[3], lty = 1, lwd = 4),
                                                                            axis_param = list(side = "left"))))

draw(ht_list4, ht_gap = unit(c(10,10), "mm"))



# end ---------------------------------------------------------------------

sessionInfo()

