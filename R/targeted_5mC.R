
# main --------------------------------------------------------------------

# updated with dorado basecaller

# in vitro differentiation experiments:

# 20210113: second set of samples, replicates from the AK28 set
# 20210312: second set of samples sent, replicates from the AKXX set
# 20210317: second set of samples sent, replicates from the ARJP28 set

# Crispr targeted LIGATION sequencing SQK-LSK109
# EXP-NBD104 barcoding
# RB1=Th0, RB2=Th1, RB3=Th2, RB4=TH17

# 

# Treg experiments:

# 20201215: CAR157, RB8=ID29 (Treg), RB9=ID32 (DMK), RB10=ID40 (tTreg)
# 20201210: CAR158, RB5=ID43 (Treg), RB6=ID46 (DMK), RB7=ID49 (tTreg)

# Crispr targeted LIGATION sequencing SQK-LSK109
# EXP-NBD104


#

# basecalled with megalodon + modification called with remora
# remora model (fast): dna_r9.4.1_e8 fast 0.0.0 5hmc_5mc CG 0
# more info on remora: https://www.youtube.com/watch?v=y8nywBuSbpU


# Targeted loci: Foxp3, Il17a, Rorc, Ifng, Il10, Maf, Il4, Gata3, Zfp362, Tbx21, Itgb8



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
  
})


annoTrack <- getAnnot("mm10")

genome <- BSgenome.Mmusculus.UCSC.mm10

genes <- annotateTranscripts(TxDb.Mmusculus.UCSC.mm10.knownGene)


setwd("~/Dropbox/BioInfo/Lab/Tcell_ID")

set.seed(123)



# prepare phenotype data -----------------------------------------------------------

bed.files <- list.files(path = "./dorado/", pattern = ".bed", recursive = T, full.names = F)
# remove files not used in this analysis: uncalled Th17, DMK, and tTreg
bed.files <- bed.files[-c(2,3,5,6,21:23)]

pdata <- as.data.frame(str_split(bed.files, "/", simplify = T))
pdata <- cbind(pdata[, 1, drop = F], str_split(pdata$V2, "_", simplify = T))
colnames(pdata) <- c("run","kit","barcode","group")
pdata$barcode <- str_sub(pdata$barcode, 1, -12)
pdata$group <- str_sub(pdata$group, 1, -5)
pdata$bed.file <- bed.files
rownames(pdata) <- paste0(pdata$group, "_", pdata$run)
pdata$group[15] <- "Th17_noTgfb"
pdata$group[16] <- "Th17_noTgfb"
table(pdata$group)

group.colors <- alphabet(20)
show_col(group.colors)

pdata$col <- pdata$group
pdata$col <- gsub("Th0", group.colors[5], pdata$col)
pdata$col <- gsub("Th17_noTgfb", group.colors[19], pdata$col)
pdata$col <- gsub("Th2", group.colors[3], pdata$col)
pdata$col <- gsub("Treg", group.colors[7], pdata$col)
pdata$col <- gsub("Th17", group.colors[15], pdata$col)
pdata$col <- gsub("Th1", group.colors[18], pdata$col)

pdata$group <- as.factor(pdata$group)
pdata$run <- as.factor(paste0("Run_", pdata$run))

head(pdata)

#write.csv(pdata, file = "manuscript/pdata.5mC.csv")
#write.csv(pdata, file = "manuscript/pdata.5mC.full.csv")
write.csv(pdata, file = "manuscript/pdata.csv")


# preprocessing ----------------------------------------------------------

pdata <- read.csv("manuscript/pdata.csv", row.names = 1)

## load bed files----

bed.files.full <- list.files(path = "./dorado/", pattern = ".bed", recursive = T, full.names = T)
bed.files.full <- bed.files.full[-c(2,3,5,6,21:23)]

colnames.bed <- c("chr","start","end","base","score","strand","tstart","tend","color","coverage","freq","mod","canon","other","del","fail","diff","nocall")
bed.names <- rownames(pdata)

bedlist <- list()
for(i in 1:length(bed.files.full)){
  bedtemp <- read.delim(bed.files.full[i], header = F)
  colnames(bedtemp) <- colnames.bed
  bedtemp <- GRanges(bedtemp)
  bedtemp <- keepStandardChromosomes(bedtemp, pruning.mode = "coarse")
  bedtemp <- sort(bedtemp)
  bedtemp <- resize(bedtemp, width = 3, fix = "start")
  bedtemp$context <- getSeq(genome, bedtemp)
  bedlist[[i]] <- bedtemp
  names(bedlist)[i] <- bed.names[i]
  rm(bedtemp)
}

par(mar=c(15,5,5,5), mfrow = c(2,1))
barplot(unlist(lapply(bedlist, length)), las = 2, col = as.factor(pdata$group))
barplot(unlist(lapply(bedlist, length)), las = 2, col = as.factor(pdata$run))

bedlist.m <- lapply(bedlist, function(x) x[x$base == "m", ])
head(bedlist.m$Th0_20210113)
bedlist.h <- lapply(bedlist, function(x) x[x$base == "h", ])
head(bedlist.h$Th0_20210113)


## inspection ----

block.sizes <- lapply(bedlist.m, function(x) as.numeric(x$freq))
names(block.sizes)
block.sizes.df <- data.frame(Reduce(cbind, block.sizes))
colnames(block.sizes.df) <- names(block.sizes)
block.sizes.df <- as.matrix(block.sizes.df)
block.sizes.df <- block.sizes.df/100 # to get beta values
head(block.sizes.df)
tail(block.sizes.df)
dim(block.sizes.df)

block.sizes2 <- lapply(bedlist.h, function(x) as.numeric(x$freq))
block.sizes.df2 <- data.frame(Reduce(cbind, block.sizes2))
colnames(block.sizes.df2) <- names(block.sizes2)
block.sizes.df2 <- as.matrix(block.sizes.df2)
block.sizes.df2 <- block.sizes.df2/100 # to get beta values

par(mar=c(5,5,5,5), mfrow = c(2,1))
densityPlot(block.sizes.df, main = "Density Plot 5mCpG")
densityPlot(block.sizes.df2, main = "Density Plot 5hmCpG")

par(mar=c(15,5,5,5), mfrow = c(3,1))
barplot(as.numeric(lapply(bedlist.m, function(x) sum(x$coverage))), las = 2,
        names.arg = names(bedlist.m),
        main = "Total reads per sample",
        col = pdata$run)
barplot(as.numeric(lapply(bedlist.m, function(x) sum(x$mod))), las = 2,
        names.arg = names(bedlist.m),
        main = "Total 5mCpG counts per sample",
        col = pdata$run)
barplot(as.numeric(lapply(bedlist.h, function(x) sum(x$mod))), las = 2,
        names.arg = names(bedlist.h),
        main = "Total 5hmCpG counts per sample",
        col = pdata$run)


saveRDS(bedlist, file="manuscript/bedlist.rds")
saveRDS(bedlist.m, file="manuscript/bedlist.m.rds") # input of MIRA
saveRDS(bedlist.h, file="manuscript/bedlist.h.rds") # input of MIRA



## make BSseq objects ----

bedlist.df <- lapply(bedlist, function(x) as.data.frame(x))
bedlist.m <- lapply(bedlist.df, function(x) x[x$base == "m", ])
head(bedlist.m[[1]])
bedlist.h <- lapply(bedlist.df, function(x) x[x$base == "h", ])
head(bedlist.h[[1]])

pre.BSseq.m <- lapply(bedlist.m, function(x) x[,c("seqnames","start","coverage","mod")])
pre.BSseq.m <- lapply(pre.BSseq.m, function(x) {colnames(x)<-c("chr","pos","N","X");x})
BSobj.m <- makeBSseqData(pre.BSseq.m, sampleNames = names(bedlist.m))
pData(BSobj.m) <- pdata

pre.BSseq.h <- lapply(bedlist.h, function(x) x[,c("seqnames","start","coverage","mod")])
pre.BSseq.h <- lapply(pre.BSseq.h, function(x) {colnames(x)<-c("chr","pos","N","X");x})
BSobj.h <- makeBSseqData(pre.BSseq.h, sampleNames = names(bedlist.h))
pData(BSobj.h) <- pdata


saveRDS(BSobj.m, file= "manuscript/BSobj.m.rds")
saveRDS(BSobj.h, file= "manuscript/BSobj.h.rds")



## PCA----

meth.data <- bsseq::getMeth(BSseq = BSobj.m,
                               regions = NULL,
                               type = "raw",
                               what = "perBase")
head(meth.data)
dim(meth.data)

all(colnames(meth.data) == rownames(pdata))

meth.data <- na.omit(meth.data)

t.counts <- t(meth.data)
t.counts[1:5,1:10]

pca_res <- prcomp(t.counts, scale. = F)

a1 <- autoplot(pca_res, data = pdata, colour = 'group', size = 4) +
  theme_minimal() +
#  scale_color_manual(values = group.colors) +
  ggtitle("PCA of 5mCpG integrated data")

meth.data <- bsseq::getMeth(BSseq = BSobj.h,
                            regions = NULL,
                            type = "raw",
                            what = "perBase")
head(meth.data)
dim(meth.data)
all(colnames(meth.data) == rownames(pdata))
meth.data <- na.omit(meth.data)
t.counts <- t(meth.data)
t.counts[1:5,1:10]
pca_res <- prcomp(t.counts, scale. = F)

a2 <- autoplot(pca_res, data = pdata, colour = 'group', size = 4) +
  theme_minimal() +
  #  scale_color_manual(values = group.colors) +
  ggtitle("PCA of 5hmCpG integrated data")

grid.arrange(a1,a2, ncol=1)



# targeted regions ----

#rm(list=ls())

#pdata <- read.csv(file = "manuscript/pdata.csv", row.names = 1)

#bedlist <- readRDS(file="manuscript/bedlist.rds")

regions <- read.csv("targeted.regions.csv", row.names = 1)
head(regions)
regions <- GRanges(regions)

bedlist.t <- lapply(bedlist, function(x) subsetByOverlaps(x, regions))
lapply(bedlist.t, length)

bedlist.t <- lapply(bedlist.t, function(x) as.data.frame(x))

bedlist.m <- lapply(bedlist.t, function(x) x[x$base == "m", ])
head(bedlist.m[[1]])
bedlist.h <- lapply(bedlist.t, function(x) x[x$base == "h", ])
head(bedlist.h[[1]])

par(mar=c(15,5,5,5), mfrow = c(2,1))
barplot(unlist(lapply(bedlist.m, nrow)), las = 2, col = as.factor(pdata$group),
        main = "Informative on-target CpG sites per sample")
barplot(unlist(lapply(bedlist.m, nrow)), las = 2, col = as.factor(pdata$run),
        main = "Informative on-target CpG sites per sample")

barplot(unlist(lapply(bedlist.m, function(x) sum(x$coverage))), 
        las = 2, col = as.factor(pdata$group),
        main = "On-target CpG site coverage per sample")
barplot(unlist(lapply(bedlist.m, function(x) sum(x$coverage))), 
        las = 2, col = as.factor(pdata$run),
        main = "On-target CpG site coverage per sample")

pre.BSseq.m <- lapply(bedlist.m, function(x) x[,c("seqnames","start","coverage","mod")])
pre.BSseq.m <- lapply(pre.BSseq.m, function(x) {colnames(x)<-c("chr","pos","N","X");x})
BSobj.m <- makeBSseqData(pre.BSseq.m, sampleNames = names(bedlist.m))

pre.BSseq.h <- lapply(bedlist.h, function(x) x[,c("seqnames","start","coverage","mod")])
pre.BSseq.h <- lapply(pre.BSseq.h, function(x) {colnames(x)<-c("chr","pos","N","X");x})
BSobj.h <- makeBSseqData(pre.BSseq.h, sampleNames = names(bedlist.h))


saveRDS(BSobj.m, file= "manuscript/BSobj.m_regions.rds")
saveRDS(BSobj.h, file= "manuscript/BSobj.h_regions.rds")


meth.data <- bsseq::getMeth(BSseq = BSobj.m,
                            regions = NULL,
                            type = "raw",
                            what = "perBase")
head(meth.data)
dim(meth.data)
all(colnames(meth.data) == rownames(pdata))
meth.data <- na.omit(meth.data)
t.counts <- t(meth.data)
t.counts[1:5,1:10]
pca_res <- prcomp(t.counts, scale. = F)

a1 <- autoplot(pca_res, data = pdata, colour = 'group', size = 4) +
  theme_minimal() +
  #  scale_color_manual(values = group.colors) +
  ggtitle("PCA of 5mCpG integrated data")

meth.data <- bsseq::getMeth(BSseq = BSobj.h,
                            regions = NULL,
                            type = "raw",
                            what = "perBase")
head(meth.data)
dim(meth.data)
all(colnames(meth.data) == rownames(pdata))
meth.data <- na.omit(meth.data)
t.counts <- t(meth.data)
t.counts[1:5,1:10]
pca_res <- prcomp(t.counts, scale. = F)

a2 <- autoplot(pca_res, data = pdata, colour = 'group', size = 4) +
  theme_minimal() +
  #  scale_color_manual(values = group.colors) +
  ggtitle("PCA of 5hmCpG integrated data")

grid.arrange(a1,a2, ncol=1)





## DMRplots ----------------------------------------------------------------

rm(list=ls())

suppressPackageStartupMessages({
  library(bsseq)
  library(dmrseq)
  library(pals)
  library(scales)
  
})

annoTrack <- getAnnot("mm10")

pdata <- read.csv(file = "manuscript/pdata.csv", row.names = 1)

BSobj.m <- readRDS("manuscript/BSobj.m.rds")
BSobj.h <- readRDS("manuscript/BSobj.h.rds")

regions <- read.csv("targeted.regions.csv", row.names = 1)
head(regions)

# test
sel.num <- 1 # selected region
BSobj = BSobj.m
#BSobj = BSobj.h
plotDMRs(BSobj, regions=regions[sel.num,], testCovariate="group", addRegions = NULL, 
         main = paste0(rownames(regions)[sel.num], " Targeted Region_5mC"),
         extend = 1000, annoTrack = annoTrack, highlightMain = F, qval = F, stat = F, horizLegend = T)


for(i in 1:nrow(regions)){
  BSobj = BSobj.m
  jpeg(filename = paste0("5mCpG_", rownames(regions)[i], ".jpeg" ), height = 800, width = 1200, res = 150)
  plotDMRs(BSobj, regions=regions[i,], testCovariate="group", addRegions = NULL, 
           main = paste0(rownames(regions)[i], " Targeted Region_5mCpG"),
           extend = 0, annoTrack = annoTrack, highlightMain = F, qval = F, stat = F, horizLegend = T)
  dev.off()
  BSobj = BSobj.h
  jpeg(filename = paste0("5mhCpG_", rownames(regions)[i], ".jpeg" ), height = 800, width = 1200, res = 150)
  plotDMRs(BSobj, regions=regions[i,], testCovariate="group", addRegions = NULL, 
           main = paste0(rownames(regions)[i], " Targeted Region_5hmCpG"),
           extend = 0, annoTrack = annoTrack, highlightMain = F, qval = F, stat = F, horizLegend = T)
  dev.off()
  rm(BSobj)
}




## heatmap -----------------------------------------------------------------

rm(list=ls())

BSobj.m <- readRDS("manuscript/BSobj.m.rds")
BSobj.h <- readRDS("manuscript/BSobj.h.rds")

hasBeenSmoothed(BSobj.m)

regions <- read.csv("targeted.regions.csv", row.names = 1)
head(regions)

BSobj = BSobj.m
#BSobj = BSobj.h

# smoothing does not significantly change the clustering
#smoothed.BSoj <- BSmooth(BSobj)
#saveRDS(smoothed.BSoj, file= "manuscript/smoothed.BSobj.m.rds")
#heatmap.data <- bsseq::getMeth(BSseq = smoothed.BSoj, regions = GRanges(regions), type = "smooth", what = "perRegion")

heatmap.data <- bsseq::getMeth(BSseq = BSobj,
                               regions = GRanges(regions),
                               type = "raw",
                               what = "perRegion")
head(heatmap.data)
dim(heatmap.data)
rownames(heatmap.data) <- rownames(regions)

#heatmap.data %>%   na.omit() %>% as.matrix()

pheno <- as.data.frame(pData(BSobj))
pheno$group <- as.factor(pheno$group)
group.colors <- levels(as.factor(pheno$col))
names(group.colors) <- levels(pheno$group)
group.colors =  list(group = group.colors)

pheno <- pheno[, "group", drop = F]
levels(pheno$group)
head(pheno)

library(pheatmap)
pheatmap::pheatmap(heatmap.data,
                   scale = "row",
                   annotation_col = pheno,
                   annotation_colors = group.colors,
                   show_rownames = T, fontsize_row = 12,
                   show_colnames = T,
                   angle_col = 45,
                   border_color = "grey",
                   main = "Z-Scores of mean 5mC in all targeted regions",
#                   main = "Z-Scores of mean 5hmC in all targeted regions",
                   fontsize = 12,
                   cellwidth = 25,
                   cellheight = 20)


#write.csv(heatmap.data, file = "manuscript/heatmaps/5mC.heatmap.data.csv")



# differential methylation (5mC) --------------------------------------------------------

rm(list=ls())

genes <- annotateTranscripts(TxDb.Mmusculus.UCSC.mm10.knownGene)

BSobj <- readRDS("manuscript/BSobj.m.rds")

pheno <- as.data.frame(pData(BSobj))
pheno$group

# Treg vs Th0
dml.test <- DMLtest(BSobj, group1 = c(1,2), group2 = c(3,7,11), ncores = 12, smoothing = F) # 

# Th1 vs Th0
dml.test <- DMLtest(BSobj, group1 = c(4,8,12), group2 = c(3,7,11),  ncores = 12) # 

# Th2 vs Th0
dml.test <- DMLtest(BSobj, group1 = c(5,9,13), group2 = c(3,7,11),  ncores = 12) # 

# Th17 vs Th0
dml.test <- DMLtest(BSobj, group1 = c(6,10,14), group2 = c(3,7,11),  ncores = 12) # 

# Th17.no.TGFb vs Th17
dml.test <- DMLtest(BSobj, group1 = c(15:16), group2 = c(6,10,14),  ncores = 12) # 

# Th17.no.Tgfb vs Th0
dml.test <- DMLtest(BSobj, group1 = c(15:16), group2 = c(3,7,11),  ncores = 12) # 



dmls = callDML(dml.test, p.threshold=0.05, delta = 0.1)
table(dmls$fdr < 0.05)
head(dmls)
tail(dmls)

dmrs = callDMR(dml.test, p.threshold=0.05, delta = 0.1, minCG = 3, minlen = 50, dis.merge = 100, pct.sig = 0.5)
head(dmrs)
tail(dmrs)

dmls.gr <- dmls
colnames(dmls.gr)[2] <- "start"
dmls.gr$end <- dmls.gr$start
dmls.gr <- GRanges(dmls.gr)
match1 <- matchGenes(dmls.gr, genes, type = "fiveprime", promoterDist = 2500, skipExons = FALSE, verbose = TRUE)
dmls <- cbind(dmls, match1)
head(dmls)

dmrs.gr <- dmrs
colnames(dmrs.gr)[2] <- "start"
dmrs.gr$end <- dmrs.gr$start
dmrs.gr <- GRanges(dmrs.gr)
match1 <- matchGenes(dmrs.gr, genes, type = "fiveprime", promoterDist = 2500, skipExons = FALSE, verbose = TRUE)
dmrs <- cbind(dmrs, match1)
head(dmrs)
dmrs$name

# Treg vs Th0
write.csv(dmls, file = "manuscript/tables/DMLs_5mC_Treg.vs.Th0.csv")
write.csv(dmrs, file = "manuscript/tables/DMRs_5mC_Treg.vs.Th0.csv")
# Th1 vs Th0
write.csv(dmls, file = "manuscript/tables/DMLs_5mC_Th1.vs.Th0.csv")
write.csv(dmrs, file = "manuscript/tables/DMRs_5mC_Th1.vs.Th0.csv")
# Th2 vs Th0
write.csv(dmls, file = "manuscript/tables/DMLs_5mC_Th2.vs.Th0.csv")
write.csv(dmrs, file = "manuscript/tables/DMRs_5mC_Th2.vs.Th0.csv")
# Th17 vs Th0
write.csv(dmls, file = "manuscript/tables/DMLs_5mC_Th17.vs.Th0.csv")
write.csv(dmrs, file = "manuscript/tables/DMRs_5mC_Th17.vs.Th0.csv")
# Th17.no.TGFb vs Th17
write.csv(dmls, file = "manuscript/tables/DMLs_5mC_Th17.no.TGFb.vs.Th17.csv")
write.csv(dmrs, file = "manuscript/tables/DMRs_5mC_Th17.no.TGFb.vs.Th17.csv")
# Th17.no.Tgfb vs Th0
write.csv(dmls, file = "manuscript/tables/DMLs_5mC_Th17.no.TGFb.vs.Th0.csv")
write.csv(dmrs, file = "manuscript/tables/DMRs_5mC_Th17.no.Tgfb.vs.Th0.csv")




## Targeted DMRplots ----------------------------------------------------------------

# customized DMR plots
# edited dmrseq::plotDMRs functions to include chromatin status
# and minor format changes

rm(list=ls())

detach("package:dmrseq", unload = TRUE)

source("plotDMRs.R")
source("helper_functions.R")

library(bsseq)
library(RColorBrewer)
library(scales)

load("annoTrack2_mm10.RData")

BSobj <- readRDS("manuscript/BSobj.m.rds")


# Treg

sampleNames(BSobj)
BS.subset <- BSobj[, c(1:2, 3,7,11)]
sampleNames(BS.subset)

targets <- read.csv("manuscript/tables/DMRs/DMRs_5mC_Treg.vs.Th0.csv", row.names = 1)

region.cols <- c(rep(alpha("orange", 0.3), nrow(targets)))
sel.num <- 3 # selected region
plotDMRs(BS.subset, regions=targets[sel.num,], regionCol = region.cols, 
         testCovariate="group", addRegions = targets, 
         main = paste0(targets$name[sel.num], " Targeted Region_5mC"),
         extend = 1000, annoTrack = annoTrack2, highlightMain = F, qval = F, stat = F, horizLegend = T)

# Th1

sampleNames(BSobj)
BS.subset <- BSobj[, c(4,8,12, 3,7,11)]
sampleNames(BS.subset)

targets <- read.csv("manuscript/tables/DMRs/DMRs_5mC_Th1.vs.Th0.csv", row.names = 1)

sel.num <- 1 # selected region
plotDMRs(BS.subset, regions=targets[sel.num,], regionCol = region.cols, testCovariate="group", addRegions = targets, 
         main = paste0(targets$name[sel.num], " Targeted Region_5mC"),
         extend = 1000, annoTrack = annoTrack2, highlightMain = F, qval = F, stat = F, horizLegend = T)

# Th2

sampleNames(BSobj)
BS.subset <- BSobj[, c(5,9,13, 3,7,11)]
sampleNames(BS.subset)

targets <- read.csv("manuscript/tables/DMRs/DMRs_5mC_Th2.vs.Th0.csv", row.names = 1)

sel.num <- 1 # selected region
extended.region <- GRanges(targets[sel.num, ])
start(extended.region) <- start(extended.region)-2000
plotDMRs(BS.subset, regions=extended.region, testCovariate="group", regionCol = region.cols, addRegions = targets, 
         main = paste0(targets$name[sel.num], " Targeted Region_5mC"),
         extend = 3000, annoTrack = annoTrack2, highlightMain = F, qval = F, stat = F, horizLegend = T)

sel.num <- 5 # selected region
plotDMRs(BS.subset, regions=targets[sel.num,], testCovariate="group", regionCol = region.cols, addRegions = targets, 
         main = paste0(targets$name[sel.num], " Targeted Region_5mC"),
         extend = 8000, annoTrack = annoTrack2, highlightMain = F, qval = F, stat = F, horizLegend = T)


sel.num <- 8 # selected region
plotDMRs(BS.subset, regions=targets[sel.num,], testCovariate="group", regionCol = region.cols, addRegions = targets, 
         main = paste0(targets$name[sel.num], " Targeted Region_5mC"),
         extend = 9000, annoTrack = annoTrack2, highlightMain = F, qval = F, stat = F, horizLegend = T)


# Th17

sampleNames(BSobj)
BS.subset <- BSobj[, c(6,10,14, 3,7,11)]
sampleNames(BS.subset)

targets <- read.csv("manuscript/tables/DMRs/DMRs_5mC_Th17.vs.Th0.csv", row.names = 1)

sel.num <- 1 # selected region
plotDMRs(BS.subset, regions=targets[sel.num,], testCovariate="group", regionCol = region.cols, addRegions = targets, 
         main = paste0(targets$name[sel.num], " Targeted Region_5mC"),
         extend = 1000, annoTrack = annoTrack2, highlightMain = F, qval = F, stat = F, horizLegend = T)

sel.num <- 3 # selected region
plotDMRs(BS.subset, regions=targets[sel.num,], testCovariate="group", regionCol = region.cols, addRegions = targets, 
         main = paste0(targets$name[sel.num], " Targeted Region_5mC"),
         extend = 1000, annoTrack = annoTrack2, highlightMain = F, qval = F, stat = F, horizLegend = T)

sel.num <- 4 # selected region
extended.region <- GRanges(targets[sel.num, ])
end(extended.region) <- end(extended.region)+15000
plotDMRs(BS.subset, regions=extended.region, testCovariate="group", regionCol = region.cols, addRegions = targets, 
         main = paste0(targets$name[sel.num], " Targeted Region_5mC"),
         extend = 2000, annoTrack = annoTrack2, highlightMain = F, qval = F, stat = F, horizLegend = T)


# Th17 in absence of TGFb

sampleNames(BSobj)
BS.subset <- BSobj[, c(15,16,6,10,14)]
sampleNames(BS.subset)
pData(BS.subset)$group <- gsub("Th17_noTgfb", "Th17.noTgfb", pData(BS.subset)$group)

targets <- read.csv("manuscript/tables/DMRs/DMRs_5mC_Th17.no.TGFb.vs.Th17.csv", row.names = 1)

sel.num <- 1 # selected region
plotDMRs(BS.subset, regions=targets[sel.num,], testCovariate="group", regionCol = region.cols, addRegions = targets, 
         main = paste0(targets$name[sel.num], " Targeted Region_5mC"),
         extend = 1000, annoTrack = annoTrack2, highlightMain = F, qval = F, stat = F, horizLegend = T)



# Th17 in absence of TGFb vs Th0

sampleNames(BSobj)
BS.subset <- BSobj[, c(15,16, 3,7,11)]
sampleNames(BS.subset)
pData(BS.subset)$group
pData(BS.subset)$group <- gsub("Th17_noTgfb", "Th17.noTgfb", pData(BS.subset)$group)

targets <- read.csv("manuscript/tables/DMRs/DMRs_5mC_Th17.no.Tgfb.vs.Th0.csv", row.names = 1)

sel.num <- 1 # selected region
plotDMRs(BS.subset, regions=targets[sel.num,], testCovariate="group", regionCol = region.cols, addRegions = targets, 
         main = paste0(targets$name[sel.num], " Targeted Region_5mC"),
         extend = 1000, annoTrack = annoTrack2, highlightMain = F, qval = F, stat = F, horizLegend = T)




## Targeted heatmap -----------------------------------------------------------------

rm(list=ls())

BSobj <- readRDS("manuscript/BSobj.m.rds")
#BSobj <- readRDS("manuscript/smoothed.BSobj.m.rds")

# cell type specific DMRs

targets.1 <- read.csv("manuscript/tables/DMRs/DMRs_5mC_Treg.vs.Th0.csv", row.names = 1)
targets.2 <- read.csv("manuscript/tables/DMRs/DMRs_5mC_Th1.vs.Th0.csv", row.names = 1)
targets.3 <- read.csv("manuscript/tables/DMRs/DMRs_5mC_Th2.vs.Th0.csv", row.names = 1)
targets.4 <- read.csv("manuscript/tables/DMRs/DMRs_5mC_Th17.vs.Th0.csv", row.names = 1)
targets.5 <- read.csv("manuscript/tables/DMRs/DMRs_5mC_Th17.no.Tgfb.vs.Th0.csv", row.names = 1)
targets.6 <- read.csv("manuscript/tables/DMRs/DMRs_5mC_Th17.no.TGFb.vs.Th17.csv", row.names = 1)

targets <- rbind(targets.1, targets.2, targets.3, targets.4, targets.5, targets.6)
nrow(targets) # 22

# smoothing doesn't work

heatmap.data <- bsseq::getMeth(BSseq = BSobj,
                               regions = GRanges(targets),
                               type = "raw",
#                               type = "smooth",
                               what = "perRegion")
head(heatmap.data)
#rownames(heatmap.data) <- make.names(targets$name, unique = T)
rownames(heatmap.data) <- targets$name

heatmap.data %>%   na.omit() %>% as.matrix()

pheno <- as.data.frame(pData(BSobj))
pheno$group <- as.factor(pheno$group)
group.colors <- levels(as.factor(pheno$col))
names(group.colors) <- levels(pheno$group)
group.colors =  list(group = group.colors)

pheno <- pheno[, "group", drop = F]
levels(pheno$group)
head(pheno)

pheatmap::pheatmap(heatmap.data,
                   scale = "row",
                   annotation_col = pheno,
                   annotation_colors = group.colors,
                   show_rownames = T, fontsize_row = 12,
                   show_colnames = T,
                   angle_col = 45,
                   border_color = "grey",
                   main = "Cell type specific DMRs (5mC)",
                   fontsize = 14,
                   cellwidth = 20,
                   cellheight = 20)

write.csv(heatmap.data, file = "manuscript/heatmaps/5mC.targeted.heatmap.data.csv")

# no better clustering with smoothed data


# differential methylation (5hmC) --------------------------------------------------------

rm(list=ls())

genes <- annotateTranscripts(TxDb.Mmusculus.UCSC.mm10.knownGene)

#BSobj <- readRDS("manuscript/BSobj.h.rds")
BSobj <- readRDS("manuscript/BSobj.h_regions.rds")

pheno <- as.data.frame(pData(BSobj))
pheno$group

# Treg vs Th0
dml.test <- DMLtest(BSobj, group1 = c(1,2), group2 = c(3,7,11), ncores = 12, smoothing = T) # 

# Th1 vs Th0
dml.test <- DMLtest(BSobj, group1 = c(4,8,12), group2 = c(3,7,11),  ncores = 12, smoothing = T) # 

# Th2 vs Th0
dml.test <- DMLtest(BSobj, group1 = c(5,9,13), group2 = c(3,7,11),  ncores = 12, smoothing = T) # 

# Th17 vs Th0
dml.test <- DMLtest(BSobj, group1 = c(6,10,14), group2 = c(3,7,11),  ncores = 12, smoothing = T) # 

# Th17.no.TGFb vs Th17
dml.test <- DMLtest(BSobj, group1 = c(15:16), group2 = c(6,10,14),  ncores = 12, smoothing = T) # 

# Th17.no.Tgfb vs Th0
dml.test <- DMLtest(BSobj, group1 = c(15:16), group2 = c(3,7,11),  ncores = 12, smoothing = T) # 


dmls = callDML(dml.test, p.threshold=0.05, delta = 0) # no delta threshold used for 5hmC
table(dmls$fdr < 0.05)
head(dmls)
tail(dmls)

dmrs = callDMR(dml.test, p.threshold=0.05, delta = 0, minCG = 3, minlen = 50, dis.merge = 100, pct.sig = 0.5)  # no delta threshold used for 5hmC
head(dmrs)
tail(dmrs)

dmls.gr <- dmls
colnames(dmls.gr)[2] <- "start"
dmls.gr$end <- dmls.gr$start
dmls.gr <- GRanges(as.data.frame(dmls.gr))
match1 <- matchGenes(dmls.gr, genes, type = "fiveprime", promoterDist = 2500, skipExons = FALSE, verbose = TRUE)
dmls <- cbind(dmls, match1)
head(dmls)
dmls$name

dmrs.gr <- dmrs
colnames(dmrs.gr)[2] <- "start"
dmrs.gr$end <- dmrs.gr$start
dmrs.gr <- GRanges(dmrs.gr)
match1 <- matchGenes(dmrs.gr, genes, type = "fiveprime", promoterDist = 2500, skipExons = FALSE, verbose = TRUE)
dmrs <- cbind(dmrs, match1)
head(dmrs)
dmrs$name

# Treg vs Th0
write.csv(dmls, file = "manuscript/tables/DMLs/DMLs_5hmC_Treg.vs.Th0.csv")
write.csv(dmrs, file = "manuscript/tables/DMRs/DMRs_5hmC_Treg.vs.Th0.csv")
# Th1 vs Th0
write.csv(dmls, file = "manuscript/tables/DMLs/DMLs_5hmC_Th1.vs.Th0.csv")
write.csv(dmrs, file = "manuscript/tables/DMRs/DMRs_5hmC_Th1.vs.Th0.csv")
# Th2 vs Th0
write.csv(dmls, file = "manuscript/tables/DMLs/DMLs_5hmC_Th2.vs.Th0.csv")
write.csv(dmrs, file = "manuscript/tables/DMRs/DMRs_5hmC_Th2.vs.Th0.csv")
# Th17 vs Th0
write.csv(dmls, file = "manuscript/tables/DMLs/DMLs_5hmC_Th17.vs.Th0.csv")
write.csv(dmrs, file = "manuscript/tables/DMRs/DMRs_5hmC_Th17.vs.Th0.csv")
# Th17.no.TGFb vs Th17
write.csv(dmls, file = "manuscript/tables/DMLs/DMLs_5hmC_Th17.no.TGFb.vs.Th17.csv")
write.csv(dmrs, file = "manuscript/tables/DMRs/DMRs_5hmC_Th17.no.TGFb.vs.Th17.csv")
# Th17.no.Tgfb vs Th0
write.csv(dmls, file = "manuscript/tables/DMLs/DMLs_5hmC_Th17.no.TGFb.vs.Th0.csv")
write.csv(dmrs, file = "manuscript/tables/DMRs/DMRs_5hmC_Th17.no.Tgfb.vs.Th0.csv")




## Targeted heatmap -----------------------------------------------------------------

rm(list=ls())

BSobj <- readRDS("manuscript/BSobj.h.rds")

# cell type specific DMRs

#targets.1 <- read.csv("manuscript/tables/DMLs/DMLs_5hmC_Treg.vs.Th0.csv", row.names = 1)
#targets.2 <- read.csv("manuscript/tables/DMLs/DMLs_5hmC_Th1.vs.Th0.csv", row.names = 1)
#targets.3 <- read.csv("manuscript/tables/DMLs/DMLs_5hmC_Th2.vs.Th0.csv", row.names = 1)
#targets.4 <- read.csv("manuscript/tables/DMLs/DMLs_5hmC_Th17.vs.Th0.csv", row.names = 1)
#targets.5 <- read.csv("manuscript/tables/DMLs/DMLs_5hmC_Th17.no.TGFb.vs.Th17.csv", row.names = 1)
#targets.6 <- read.csv("manuscript/tables/DMLs/DMLs_5hmC_Th17.no.TGFb.vs.Th0.csv", row.names = 1)

targets.1 <- read.csv("manuscript/tables/DMRs/DMRs_5hmC_Treg.vs.Th0.csv", row.names = 1)
targets.2 <- read.csv("manuscript/tables/DMRs/DMRs_5hmC_Th1.vs.Th0.csv", row.names = 1)
targets.3 <- read.csv("manuscript/tables/DMRs/DMRs_5hmC_Th2.vs.Th0.csv", row.names = 1)
targets.4 <- read.csv("manuscript/tables/DMRs/DMRs_5hmC_Th17.vs.Th0.csv", row.names = 1)
targets.5 <- read.csv("manuscript/tables/DMRs/DMRs_5hmC_Th17.no.Tgfb.vs.Th0.csv", row.names = 1)
targets.6 <- read.csv("manuscript/tables/DMRs/DMRs_5hmC_Th17.no.TGFb.vs.Th17.csv", row.names = 1)

targets <- rbind(targets.1, targets.2, targets.3, targets.4, targets.5, targets.6)

targets <- targets[order(targets$nCG, decreasing = T), ]
dup <- duplicated(targets$start)
targets <- targets[!dup,]
head(targets)
nrow(targets) # 35
#colnames(targets)[2] <- "start"
#targets$end <- targets$start

# smoothing doesn't work

heatmap.data <- bsseq::getMeth(BSseq = BSobj,
                               regions = GRanges(targets),
                               type = "raw",
                               what = "perRegion")
head(heatmap.data)
#rownames(heatmap.data) <- make.names(targets$name, unique = T)
rownames(heatmap.data) <- targets$name

#heatmap.data <-  na.omit(heatmap.data) %>% as.matrix()
#heatmap.data <- heatmap.data + 0.001
#dim(heatmap.data) # 303 16

pheno <- as.data.frame(pData(BSobj))
pheno$group <- as.factor(pheno$group)
group.colors <- levels(as.factor(pheno$col))
names(group.colors) <- levels(pheno$group)
group.colors =  list(group = group.colors)

pheno <- pheno[, "group", drop = F]
levels(pheno$group)
head(pheno)

pheatmap::pheatmap(heatmap.data,
                   scale = "row",
                   annotation_col = pheno,
                   annotation_colors = group.colors,
                   show_rownames = T, fontsize_row = 12,
                   show_colnames = T,
                   angle_col = 45,
                   border_color = "grey",
                   main = "Cell type specific DMRs (5hmC)",
                   fontsize = 14,
                   cellwidth = 20,
                   cellheight = 15)

write.csv(heatmap.data, file = "manuscript/heatmaps/5hmC.targeted.heatmap.data.csv")




## Targeted DMRplots ----------------------------------------------------------------

# customized DMR plots
# edited dmrseq::plotDMRs functions to include chromatin status
# and minor format changes

rm(list=ls())

detach("package:dmrseq", unload = TRUE)

source("plotDMRs.R")
source("helper_functions_5hmC.R")

library(bsseq)
library(RColorBrewer)
library(scales)

load("annoTrack2_mm10.RData")

BSobj <- readRDS("manuscript/BSobj.h.rds")


# Treg

sampleNames(BSobj)
BS.subset <- BSobj[, c(1:2, 3,7,11)]
sampleNames(BS.subset)

targets <- read.csv("manuscript/tables/DMRs/DMRs_5hmC_Treg.vs.Th0.csv", row.names = 1)
targets <- targets[order(targets$nCG, decreasing = T), ]
head(targets)

region.cols <- c(rep(alpha("orange", 0.3), nrow(targets)))
sel.num <- 4 # selected region
plotDMRs(BS.subset, regions=targets[sel.num,], regionCol = region.cols, 
         testCovariate="group", addRegions = targets, 
         main = paste0(targets$name[sel.num], " Targeted Region_5hmC"),
         extend = 1000, annoTrack = annoTrack2, highlightMain = F, qval = F, stat = F, horizLegend = T)

# Th1

sampleNames(BSobj)
BS.subset <- BSobj[, c(4,8,12, 3,7,11)]
sampleNames(BS.subset)

targets <- read.csv("manuscript/tables/DMRs/DMRs_5hmC_Th1.vs.Th0.csv", row.names = 1)
targets <- targets[order(targets$nCG, decreasing = T), ]
head(targets)

sel.num <- 2 # selected region
plotDMRs(BS.subset, regions=targets[sel.num,], regionCol = region.cols, testCovariate="group", addRegions = targets, 
         main = paste0(targets$name[sel.num], " Targeted Region_5hmC"),
         extend = 1000, annoTrack = annoTrack2, highlightMain = F, qval = F, stat = F, horizLegend = T)

# Th2

sampleNames(BSobj)
BS.subset <- BSobj[, c(5,9,13, 3,7,11)]
sampleNames(BS.subset)

targets <- read.csv("manuscript/tables/DMRs/DMRs_5hmC_Th2.vs.Th0.csv", row.names = 1)
targets <- targets[order(targets$nCG, decreasing = T), ]
head(targets)

sel.num <- 2 # selected region
plotDMRs(BS.subset, regions=targets[sel.num,], testCovariate="group", regionCol = region.cols, addRegions = targets, 
         main = paste0(targets$name[sel.num], " Targeted Region_5hmC"),
         extend = 3000, annoTrack = annoTrack2, highlightMain = F, qval = F, stat = F, horizLegend = T)

# Th17

sampleNames(BSobj)
BS.subset <- BSobj[, c(6,10,14, 3,7,11)]
sampleNames(BS.subset)

targets <- read.csv("manuscript/tables/DMRs/DMRs_5hmC_Th17.vs.Th0.csv", row.names = 1)
targets <- targets[order(targets$nCG, decreasing = T), ]
head(targets)

sel.num <- 1 # selected region
plotDMRs(BS.subset, regions=targets[sel.num,], testCovariate="group", regionCol = region.cols, addRegions = targets, 
         main = paste0(targets$name[sel.num], " Targeted Region_5hmC"),
         extend = 1000, annoTrack = annoTrack2, highlightMain = F, qval = F, stat = F, horizLegend = T)

sel.num <- 2 # selected region
plotDMRs(BS.subset, regions=targets[sel.num,], testCovariate="group", regionCol = region.cols, addRegions = targets, 
         main = paste0(targets$name[sel.num], " Targeted Region_5hmC"),
         extend = 1000, annoTrack = annoTrack2, highlightMain = F, qval = F, stat = F, horizLegend = T)

sel.num <- 2 # selected region
plotDMRs(BS.subset, regions=targets[sel.num,], testCovariate="group", regionCol = region.cols, addRegions = targets, 
         main = paste0(targets$name[sel.num], " Targeted Region_5hmC"),
         extend = 1000, annoTrack = annoTrack2, highlightMain = F, qval = F, stat = F, horizLegend = T)


# Th17 in absence of TGFb

sampleNames(BSobj)
BS.subset <- BSobj[, c(15,16,6,10,14)]
sampleNames(BS.subset)
pData(BS.subset)$group <- gsub("Th17_noTgfb", "Th17.noTgfb", pData(BS.subset)$group)

targets <- read.csv("manuscript/tables/DMRs/DMRs_5hmC_Th17.no.TGFb.vs.Th17.csv", row.names = 1)
targets <- targets[order(targets$nCG, decreasing = T), ]
head(targets)

sel.num <- 1 # selected region
plotDMRs(BS.subset, regions=targets[sel.num,], testCovariate="group", regionCol = region.cols, addRegions = targets, 
         main = paste0(targets$name[sel.num], " Targeted Region_5hmC"),
         extend = 1000, annoTrack = annoTrack2, highlightMain = F, qval = F, stat = F, horizLegend = T)

sel.num <- 2 # selected region
plotDMRs(BS.subset, regions=targets[sel.num,], testCovariate="group", regionCol = region.cols, addRegions = targets, 
         main = paste0(targets$name[sel.num], " Targeted Region_5hmC"),
         extend = 1000, annoTrack = annoTrack2, highlightMain = F, qval = F, stat = F, horizLegend = T)

sel.num <- 3 # selected region
plotDMRs(BS.subset, regions=targets[sel.num,], testCovariate="group", regionCol = region.cols, addRegions = targets, 
         main = paste0(targets$name[sel.num], " Targeted Region_5hmC"),
         extend = 1000, annoTrack = annoTrack2, highlightMain = F, qval = F, stat = F, horizLegend = T)

# Th17 in absence of TGFb vs Th0

sampleNames(BSobj)
BS.subset <- BSobj[, c(15,16, 3,7,11)]
sampleNames(BS.subset)
pData(BS.subset)$group
pData(BS.subset)$group <- gsub("Th17_noTgfb", "Th17.noTgfb", pData(BS.subset)$group)

targets <- read.csv("manuscript/tables/DMRs/DMRs_5hmC_Th17.no.Tgfb.vs.Th0.csv", row.names = 1)
targets <- targets[order(targets$nCG, decreasing = T), ]
head(targets)

sel.num <- 2 # selected region
plotDMRs(BS.subset, regions=targets[sel.num,], testCovariate="group", regionCol = region.cols, addRegions = targets, 
         main = paste0(targets$name[sel.num], " Targeted Region_5hmC"),
         extend = 3000, annoTrack = annoTrack2, highlightMain = F, qval = F, stat = F, horizLegend = T)





# end ---------------------------------------------------------------------
sessionInfo()

