
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


# preprocessing ----------------------------------------------------------

## load bed files ----

genome <- BSgenome.Mmusculus.UCSC.mm10

bed.files.full <- list.files(path = "./dorado/", pattern = "calls.bam.sorted_mouse.brain.bed", recursive = T, full.names = T)
colnames.bed <- c("chr","start","end","base","score","strand","tstart","tend","color","coverage","freq","mod","canon","other","del","fail","diff","nocall")

bedtemp <- read.delim(bed.files.full, header = F)
colnames(bedtemp) <- colnames.bed
bedtemp <- GRanges(bedtemp)
bedtemp <- keepStandardChromosomes(bedtemp, pruning.mode = "coarse")
bedtemp <- sort(bedtemp)
bedtemp <- resize(bedtemp, width = 3, fix = "start")
bedtemp$context <- getSeq(genome, bedtemp)


bedlist.m <- bedtemp[bedtemp$base == "m", ]
head(bedlist.m)
bedlist.h <- bedtemp[bedtemp$base == "h", ]
head(bedlist.h)


## inspection ----


densityPlot(as.numeric(bedlist.m$freq), main = "Density Plot 5mCpG")

par(mar=c(5,5,5,5), mfrow = c(2,1))
densityPlot(as.matrix(bedlist.m$freq), main = "Density Plot 5mCpG")
densityPlot(as.matrix(bedlist.h$freq), main = "Density Plot 5hmCpG")

sum(bedlist.m$coverage) # 5071488
sum(bedlist.m$mod) # 3551368
sum(bedlist.h$coverage) # 5071488
sum(bedlist.h$mod) # 217493

hist(as.matrix(bedlist.m$freq), main = "mouse brain 5mC", col = "black", xlab = "beta values")
hist(as.matrix(bedlist.h$freq), main = "mouse brain 5hmC", col = "black", xlab = "beta values")



saveRDS(bedtemp, file="manuscript/mouse_brain.bedtemp.rds")
saveRDS(bedlist.m, file="manuscript/mouse_brain.bedlist.m.rds") # input of MIRA
saveRDS(bedlist.h, file="manuscript/mouse_brain.bedlist.h.rds") # input of MIRA



# 5mC/5hmC near genes -------------------------------------------------------

rm(list=ls())

brain.bed.5mC <- readRDS("manuscript/mouse_brain.bedlist.m.rds")
brain.bed.5hmC <- readRDS("manuscript/mouse_brain.bedlist.h.rds")

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene


# TSS

genes <- genes(TxDb.Mmusculus.UCSC.mm10.knownGene)

tss = promoters(genes, upstream = 0, downstream = 1)
tss[1:5]

tss.2 <- subsetByOverlaps(tss, brain.bed.5mC)

mat.tss.5mC = normalizeToMatrix(brain.bed.5mC, tss.2, value_column = "freq", mean_mode = "absolute", extend = 5000, w = 200, background = NA, smooth = F)
mat.tss.5hmC = normalizeToMatrix(brain.bed.5hmC, tss.2, value_column = "freq", mean_mode = "absolute", extend = 5000, w = 200, background = NA, smooth = F)

#meth_col_fun = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
#EnrichedHeatmap(mat.tss.5hmC, col = meth_col_fun, name = "methylation", column_title = "methylation near TSS")

ht_list1 <- EnrichedHeatmap(mat.tss.5hmC, name = "5hmC", col = c("whitesmoke",  plasma(5)[1]),
                           column_title = "5hmC in brain cells\nTSS",
                           axis_name_rot = 90,
                           top_annotation = HeatmapAnnotation(height = unit(5,"cm"),
                                                              enrich = anno_enriched(ylim = c(0, 80), 
                                                                                     gp = gpar(col =  plasma(5)[1], lty = 1, lwd = 4),
                                                                                     axis_param = list(side = "left")))) +
  EnrichedHeatmap(mat.tss.5mC, name = "5mC", col = c("whitesmoke",  plasma(5)[3]),
                  column_title = "5mC in brain cells\nTSS",
                  axis_name_rot = 90,
                  top_annotation = HeatmapAnnotation(height = unit(5,"cm"),
                                                     enrich = anno_enriched(ylim = c(0, 80), 
                                                                            gp = gpar(col = plasma(5)[3], lty = 1, lwd = 4),
                                                                            axis_param = list(side = "left"))))

draw(ht_list1, ht_gap = unit(c(10,10), "mm"))
#jpeg("enriched.heatmap.brain.jpeg", height = 960, width = 800, quality = 100)
#draw(ht_list, ht_gap = unit(c(10,10), "mm"))
#dev.off()


# gene bodies

genes.2 <- subsetByOverlaps(genes, brain.bed.5mC)

mat.genes.5hmC = normalizeToMatrix(brain.bed.5hmC, genes.2, extend = c(5000,5000), mean_mode = "absolute", w = 200, smooth = F, target_ratio = 0.5, value_column = "freq", background = NA)
mat.genes.5mC = normalizeToMatrix(brain.bed.5mC, genes.2, extend = c(5000,5000), mean_mode = "absolute", w = 200, smooth = F, target_ratio = 0.5, value_column = "freq", background = NA)

ht_list2 <- EnrichedHeatmap(mat.genes.5hmC, name = "5hmC", col = c("whitesmoke",  plasma(5)[1]),
                  column_title = "5hmC in brain cells\ngene bodies",
                  axis_name_rot = 90,
                  top_annotation = HeatmapAnnotation(height = unit(5,"cm"),
                                                     enrich = anno_enriched(ylim = c(0, 80), 
                                                                            gp = gpar(col =  plasma(5)[1], lty = 1, lwd = 4),
                                                                            axis_param = list(side = "left")))) +
  EnrichedHeatmap(mat.genes.5mC, name = "5mC", col = c("whitesmoke",  plasma(5)[3]),
                  column_title = "5mC in brain cells\ngene bodies",
                  axis_name_rot = 90,
                  top_annotation = HeatmapAnnotation(height = unit(5,"cm"),
                                                     enrich = anno_enriched(ylim = c(0, 80), 
                                                                            gp = gpar(col = plasma(5)[3], lty = 1, lwd = 4),
                                                                            axis_param = list(side = "left"))))

draw(ht_list2, ht_gap = unit(c(10,10), "mm"))



# CGIs

session <- browserSession(url="https://genome.ucsc.edu/cgi-bin/")
#session <- browserSession("UCSC")
genome(session) <- "mm10"
ucscTables("mm10", "cpgIsland")
query <- ucscTableQuery(session, table = "cpgIslandExt", genome = "mm10")
cpg_islands <- getTable(query)
cgi <- GRanges(cpg_islands)
cgi <- keepStandardChromosomes(cgi, pruning.mode="coarse") # 16009

cgi.2 <- subsetByOverlaps(cgi, brain.bed.5mC)

mat.cgi.5hmC = normalizeToMatrix(brain.bed.5hmC, cgi.2, extend = 5000, mean_mode = "absolute", w = 200, smooth = F, target_ratio = 0.3, value_column = "freq", background = NA)
mat.cgi.5mC = normalizeToMatrix(brain.bed.5mC, cgi.2, extend = 5000, mean_mode = "absolute", w = 200, smooth = F, target_ratio = 0.3, value_column = "freq", background = NA)

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

enh.2 <- subsetByOverlaps(enh, brain.bed.5mC)

mat.enh.5hmC = normalizeToMatrix(brain.bed.5hmC, enh.2, extend = 5000, mean_mode = "absolute", w = 200, smooth = F, target_ratio = 0.3, value_column = "freq", background = NA)
mat.enh.5mC = normalizeToMatrix(brain.bed.5mC, enh.2, extend = 5000, mean_mode = "absolute", w = 200, smooth = F, target_ratio = 0.3, value_column = "freq", background = NA)

ht_list4 <- EnrichedHeatmap(mat.enh.5hmC, name = "5hmC", col = c("whitesmoke",  plasma(5)[3]),
                            column_title = "5hmC in brain cells\nForebrain Enhancers",
                            axis_name_rot = 90,
                            top_annotation = HeatmapAnnotation(height = unit(5,"cm"),
                                                               enrich = anno_enriched(ylim = c(0, 20), 
                                                                                      gp = gpar(col =  plasma(5)[3], lty = 1, lwd = 4),
                                                                                      axis_param = list(side = "left")))) +
  EnrichedHeatmap(mat.enh.5mC, name = "5mC", col = c("whitesmoke",  plasma(5)[3]),
                  column_title = "5mC in brain cells\nForebrain Enhancers",
                  axis_name_rot = 90,
                  top_annotation = HeatmapAnnotation(height = unit(5,"cm"),
                                                     enrich = anno_enriched(ylim = c(0, 80), 
                                                                            gp = gpar(col = plasma(5)[3], lty = 1, lwd = 4),
                                                                            axis_param = list(side = "left"))))

draw(ht_list4, ht_gap = unit(c(10,10), "mm"))


# control enhancers

#enh <- read.table("SE_11_0004_SE_mm10.bed", header = T) # from http://www.licpathway.net/sedb/index.php adipose
#enh <- read.delim("SE_11_0004_TE_mm10.bed", header = T) # from http://www.licpathway.net/sedb/index.php adipose
enh.ctrl <- read.delim("Mouse_Adipose_mm10_ele/SE_11_0004_SE_ele_mm10.bed", header = T) # from http://www.licpathway.net/sedb/index.php adipose
#enh <- read.table("SE_11_0002_SE_mm10.bed", header = T) # from http://www.licpathway.net/sedb/index.php blood
#enh <- read.table("SE_11_0002_SE_ele_mm10.bed", header = T) # from http://www.licpathway.net/sedb/index.php blood
#enh <- read.table("SE_11_0007_SE_mm10.bed", header = T) # from http://www.licpathway.net/sedb/index.php cerebellum

colnames(enh.ctrl)[1:3] <- c("chr","start","end")
head(enh.ctrl)
enh.ctrl <- GRanges(enh.ctrl)

enh.ctrl.2 <- subsetByOverlaps(enh.ctrl, brain.bed.5mC)

mat.enh.ctrl.5hmC = normalizeToMatrix(brain.bed.5hmC, enh.ctrl.2, extend = 5000, mean_mode = "absolute", w = 200, smooth = F, target_ratio = 0.3, value_column = "freq", background = NA)
mat.enh.ctrl.5mC = normalizeToMatrix(brain.bed.5mC, enh.ctrl.2, extend = 5000, mean_mode = "absolute", w = 200, smooth = F, target_ratio = 0.3, value_column = "freq", background = NA)

ht_list4b <- EnrichedHeatmap(mat.enh.ctrl.5hmC, name = "5hmC", col = c("whitesmoke",  plasma(5)[1]),
                            column_title = "5hmC in brain cells\nAdipose Enhancers",
                            axis_name_rot = 90,
                            top_annotation = HeatmapAnnotation(height = unit(5,"cm"),
                                                               enrich = anno_enriched(ylim = c(0, 20), 
                                                                                      gp = gpar(col =  plasma(5)[1], lty = 1, lwd = 4),
                                                                                      axis_param = list(side = "left")))) +
  EnrichedHeatmap(mat.enh.ctrl.5mC, name = "5mC", col = c("whitesmoke",  plasma(5)[1]),
                  column_title = "5mC in brain cells\nAdipose Enhancers",
                  axis_name_rot = 90,
                  top_annotation = HeatmapAnnotation(height = unit(5,"cm"),
                                                     enrich = anno_enriched(ylim = c(0, 80), 
                                                                            gp = gpar(col = plasma(5)[1], lty = 1, lwd = 4),
                                                                            axis_param = list(side = "left"))))

draw(ht_list4b, ht_gap = unit(c(10,10), "mm"))


#


## stats ----

rm(list=ls())

brain.bed.5hmC <- readRDS("manuscript/mouse_brain.bedlist.h.rds")

prop.table(table(brain.bed.5hmC$freq > 0)) # 3.2%
prop.table(table(brain.bed.5hmC$freq > 10)) # 3.2%
prop.table(table(brain.bed.5hmC$freq > 50)) # 2.5%

enh <- read.delim("Mouse_Forebrain_mm10/SE_12_0226_SE_mm10.bed", header = T) # from http://www.licpathway.net/sedb/index.php forebrain
colnames(enh)[1:3] <- c("chr","start","end")
head(enh)
enh <- GRanges(enh)

enh.ctrl <- read.delim("Mouse_Adipose_mm10_ele/SE_11_0004_SE_ele_mm10.bed", header = T) # from http://www.licpathway.net/sedb/index.php adipose
colnames(enh.ctrl)[1:3] <- c("chr","start","end")
head(enh.ctrl)
enh.ctrl <- GRanges(enh.ctrl)

enh.brain <- subsetByOverlaps(brain.bed.5hmC, enh)
enh.adipo <- subsetByOverlaps(brain.bed.5hmC, enh.ctrl)

hist(enh.brain$freq)
hist(enh.adipo$freq)

par(mar=c(5,5,5,5), mfrow = c(2,1))
densityPlot(as.matrix(enh.brain$freq), main = "5hmCpG Density Brain Enhancers")
densityPlot(as.matrix(enh.adipo$freq), main = "5hmCpG Density Adipose Enhancers")

table(enh.brain$freq > 10)
#FALSE   TRUE 
#217725   9605 = 4.225135 %
table(enh.adipo$freq > 10)
#FALSE  TRUE 
#30136  1081 = 3.462857 %
table(enh.brain$freq > 50)
#FALSE   TRUE 
#219791   7539 = 3.31 %
table(enh.adipo$freq > 50)
#FALSE  TRUE 
#30353  864 = 2.76 %

summary(enh.brain$freq) # Mean: 3.744
summary(enh.adipo$freq) # Mean: 3.102

enh.brain.pos <- enh.brain$freq[enh.brain$freq > 0]
enh.adipo.pos <- enh.adipo$freq[enh.adipo$freq > 0]

plot(density(enh.brain$freq))

par(mar=c(5,5,5,5), mfrow = c(2,1))
densityPlot(as.matrix(enh.brain.pos), main = "5hmCpG Density Brain Enhancers", ylim = c(0,0.1))
densityPlot(as.matrix(enh.adipo.pos), main = "5hmCpG Density Adipose Enhancers", ylim = c(0,0.1))

t.test(enh.brain$freq, enh.adipo$freq)
# p-value = 4.279e-10
# mean of x mean of y 
# 3.743609  3.101673 
wilcox.test(enh.brain$freq, enh.adipo$freq)
# p-value = 2.566e-10
# 

par(mar=c(5,5,5,5), mfrow = c(1,1))
plot(NA, xlim=range(0,110), ylim=range(0,0.12), 
     main= "5hmCpG frequency in enhancer regions", 
     ylab = "Density", xlab = "5hmCpG methylation")
#lines(density(enh.adipo.pos), col = "lightblue")
#lines(density(enh.brain.pos), col = "blue")
polygon(density(enh.adipo.pos), col = rgb(0.78, 0.89, 1, alpha = 0.2), lty = 2)
polygon(density(enh.brain.pos), col = rgb(1, 0, 0, alpha = 0.1))
legend("topleft",
       legend = c("Adipose Enhancers", "Brain Enhancers"),
       fill = c(rgb(0.78, 0.89, 1, alpha = 0.3), rgb(1, 0, 0, alpha = 0.2)),
       border = "black")
mtext("Wilcox p value = 2.6e-10")



# end ---------------------------------------------------------------------

sessionInfo()

