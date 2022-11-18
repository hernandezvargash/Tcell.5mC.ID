

# main --------------------------------------------------------------------

# in vitro differentiation experiments:

# 20210113: second set of samples, replicates from the AK28 set
# 20210312: second set of samples sent, replicates from the AKXX set
# 20210317: second set of samples sent, replicates from the ARJP28 set

# Crispr targeted LIGATION sequencing SQK-LSK109
# EXP-NBD104 barcoding
# RB1=Th0, RB2=Th1, RB3=Th2, RB4=TH17

# 

# Treg experiments (not using CAR157 nightmare run this time):

# 20201215: CAR157, RB8=ID29 (Treg), RB9=ID32 (DMK), RB10=ID40 (tTreg)
# 20201210: CAR158, RB5=ID43 (Treg), RB6=ID46 (DMK), RB7=ID49 (tTreg)

# Crispr targeted LIGATION sequencing SQK-LSK109
# EXP-NBD104

#

# Th17-specific experiments:

# 20210402: Th17 uncalled KO mice
# 20210419: Th17 uncalled WT mice

# 20210324: Crispr targeted LIGATION sequencing SQK-LSK109, EXP-NBD104, RB11: Th17nTGFB RB12: Th17nTGFB 

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

bed.files <- list.files(path = "./remora/", pattern = "5mC.bed", recursive = T, full.names = F)
# remove files not used in manuscript: uncalled Th17, DMK, and tTreg
bed.files <- bed.files[-c(1,3,4,6,21,22)]

pdata <- as.data.frame(str_split(bed.files, "/", simplify = T))
colnames(pdata) <- c("run","group","basename")
pdata$bed.file <- bed.files

rownames(pdata) <- paste0(pdata$group, "_", pdata$run)
pdata$group[15] <- "Th17_noTgfb"
pdata$group[16] <- "Th17_noTgfb"

head(pdata)

write.csv(pdata, file = "manuscript/pdata.5mC.csv")


# load bed files ----------------------------------------------------------

bed.files.full <- list.files(path = "./remora/", pattern = "5mC.bed", recursive = T, full.names = T)
bed.files.full <- bed.files.full[-c(1,3,4,6,21,22)]

# select targeted chromosomal locations to import

regions <- read.csv("targeted.regions.csv", row.names = 1)
head(regions)
regions <- GRanges(regions)


# import to a list

bed.list <- list()

for(i in 1:length(bed.files)){
  bed.list[[i]] <- import.bed(bed.files.full[i], which = regions) # used both "regions" and NULL
  bed.list[[i]] <- keepStandardChromosomes(bed.list[[i]], pruning.mode = "coarse")
  bed.list[[i]] <- sort(bed.list[[i]])
  bed.list[[i]] <- resize(bed.list[[i]], width = 3, fix = "start")
  bed.list[[i]]$context <- getSeq(genome, bed.list[[i]])
  bed.list[[i]]$blockSizes <- as.numeric(as.character(bed.list[[i]]$blockSizes))
}

names(bed.list) <- rownames(pdata)
save(bed.list, file="manuscript/bed.list.5mC.RData")



# inspection --------------------------------------------------------------

load("manuscript/bed.list.5mC.RData")
pdata <- read.csv("manuscript/pdata.5mC.csv", row.names = 1)

par(mar = c(15,10,8,10), mfrow = c(2,1))
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

block.sizes <- lapply(bed.list, function(x) as.numeric(x$blockSizes))
names(block.sizes)

block.sizes.df <- data.frame(Reduce(cbind, block.sizes))
colnames(block.sizes.df) <- names(block.sizes)
block.sizes.df <- as.matrix(block.sizes.df)
block.sizes.df <- block.sizes.df/100 # to get beta values
head(block.sizes.df)
tail(block.sizes.df)

# re-sort by group
all(rownames(pdata)==colnames(block.sizes.df))
pdata <- pdata[order(pdata$group), ]
block.sizes.df <- block.sizes.df[, rownames(pdata)]

#pdata <- pdata[order(pdata$group), ]
#block.sizes.df <- block.sizes.df[, rownames(pdata)]

par(mar = c(10,15,8,10), mfrow = c(1,1))
densityBeanPlot(block.sizes.df, 
                sampGroups = pdata$group,
#                sampNames = rownames(pdata),
                numPositions = 1000)




# prepare BSseq object ----------------------------------------------------

load("manuscript/bed.list.5mC.RData")
pdata <- read.csv("manuscript/pdata.5mC.csv", row.names = 1)

lapply(bed.list, length)

# reorder samples
all(rownames(pdata)==names(bed.list))
pdata <- pdata[order(pdata$group), ]
bed.list <- bed.list[rownames(pdata)]

length(bed.list) # 16
pre.BSseq <- lapply(bed.list, function(x) as.data.frame(x))
head(pre.BSseq[[1]])
pre.BSseq <- lapply(pre.BSseq, function(x) x[,c("seqnames","start","blockCount","blockSizes")])
pre.BSseq <- lapply(pre.BSseq, function(x) {colnames(x)<-c("chr", "pos", "N", "X");x})
pre.BSseq <- lapply(pre.BSseq, function(x)  transform(x, X = X*N/100))
lapply(pre.BSseq, dim)

head(pre.BSseq[[1]])

## make BSseq objects
BSobj <- makeBSseqData(pre.BSseq, sampleNames = rownames(pdata))
sampleNames(BSobj)
pData(BSobj) <- pdata

save(BSobj, file= "manuscript/BSobj_5mCpG.RData")



# PCA ---------------------------------------------------------------------

meth.data <- bsseq::getMeth(BSseq = BSobj,
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

autoplot(pca_res, data = pdata, colour = 'group', size = 4) +
  theme_minimal() +
#  scale_color_manual(values = group.colors) +
  ggtitle("PCA integrated data")

# one of the Treg samples cluster far appart from its replicate

sel.pdata <- pdata[-16, ]
sel.t.counts <- t.counts[-16, ]

pca_res2 <- prcomp(sel.t.counts, scale. = F)

autoplot(pca_res2, data = sel.pdata, colour = 'group', size = 4) +
  theme_minimal() +
  #  scale_color_manual(values = group.colors) +
  ggtitle("PCA integrated data")




# DMRplots ----------------------------------------------------------------

rm(list=ls())

suppressPackageStartupMessages({
  library(bsseq)
  library(dmrseq)
  library(pals)
  library(scales)
  
})

annoTrack <- getAnnot("mm10")

load("manuscript/BSobj_5mCpG.RData")

regions <- read.csv("targeted.regions.csv", row.names = 1)
head(regions)

pheno <- pData(BSobj)
table(pheno$group)

group.colors <- alphabet(20)
show_col(group.colors)

pheno$col <- pheno$group
pheno$col <- gsub("Th0", group.colors[5], pheno$col)
pheno$col <- gsub("Th17_noTgfb", group.colors[19], pheno$col)
pheno$col <- gsub("Th2", group.colors[3], pheno$col)
pheno$col <- gsub("Treg", group.colors[7], pheno$col)
pheno$col <- gsub("Th17", group.colors[15], pheno$col)
pheno$col <- gsub("Th1", group.colors[18], pheno$col)

pData(BSobj) <- pheno

# test
sel.num <- 1 # selected region
plotDMRs(BSobj, regions=regions[sel.num,], testCovariate="group", addRegions = NULL, 
         main = paste0(rownames(regions)[sel.num], " Targeted Region_5mC"),
         extend = 0, annoTrack = annoTrack, highlightMain = F, qval = F, stat = F, horizLegend = T)

for(i in 1:nrow(regions)){
  
  plotDMRs(BSobj, regions=regions[i,], testCovariate="group", addRegions = NULL, 
           main = paste0(rownames(regions)[i], " Targeted Region_5mC"),
           extend = 0, annoTrack = annoTrack, highlightMain = F, qval = F, stat = F, horizLegend = T)
  
}


save(BSobj, file= "manuscript/BSobj_5mCpG.RData")



# heatmap -----------------------------------------------------------------

rm(list=ls())

load("manuscript/BSobj_5mCpG.RData")
hasBeenSmoothed(BSobj)

regions <- read.csv("targeted.regions.csv", row.names = 1)
head(regions)


# smoothed.BSoj <- BSmooth(BSobj)
# smoothing this object gives an error

heatmap.data <- bsseq::getMeth(BSseq = BSobj,
                               regions = GRanges(regions),
                               type = "raw",
                               what = "perRegion")
head(heatmap.data)
rownames(heatmap.data) <- rownames(regions)

#heatmap.data %>%   na.omit() %>% as.matrix()

pheno <- as.data.frame(pData(BSobj))
pheno$group <- as.factor(pheno$group)
pheno <- pheno[, 2, drop = F]
levels(pheno$group)
head(pheno)

group.colors <- alphabet(20)
group.colors <- group.colors[c(5,18,15,19,3,7)]
names(group.colors) <- levels(pheno$group)
group.colors =  list(group = group.colors)

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
                   fontsize = 12,
                   cellwidth = 25,
                   cellheight = 20)



# differential methylation (5mC) --------------------------------------------------------

load("manuscript/BSobj_5mCpG.RData")

pheno <- as.data.frame(pData(BSobj))
pheno$group

# Treg vs Th0
dml.test <- DMLtest(BSobj, group1 = c(15,16), group2 = c(1:3), ncores = 12) # 
# Th1 vs Th0
dml.test <- DMLtest(BSobj, group1 = c(4:6), group2 = c(1:3),  ncores = 12) # 
# Th2 vs Th0
dml.test <- DMLtest(BSobj, group1 = c(12:14), group2 = c(1:3),  ncores = 12) # 
# Th17 vs Th0
dml.test <- DMLtest(BSobj, group1 = c(7:9), group2 = c(1:3),  ncores = 12) # 
# Th17.no.TGFb vs Th17
dml.test <- DMLtest(BSobj, group1 = c(10:11), group2 = c(7:9),  ncores = 12) # 
# Th17.no.Tgfb vs Th0
dml.test <- DMLtest(BSobj, group1 = c(10:11), group2 = c(1:3),  ncores = 12) # 

dmls = callDML(dml.test, p.threshold=0.05, delta = 0.1)
head(dmls)

dmrs = callDMR(dml.test, p.threshold=0.05, delta = 0.1)
# dmrs = callDMR(dml.test, p.threshold=0.05, delta = 0, minCG = 2) # criteria for Treg vs Th0
head(dmrs)

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




# Targeted DMRplots ----------------------------------------------------------------

rm(list=ls())

annoTrack <- getAnnot("mm10")

load("manuscript/BSobj_5mCpG.RData")

regions <- read.csv("targeted.regions.csv", row.names = 1)
head(regions)


# Treg

sampleNames(BSobj)
BS.subset <- BSobj[, c(1:3, 15:16)]
sampleNames(BS.subset)

targets <- read.csv("manuscript/tables/DMRs/DMRs_5mC_Treg.vs.Th0.csv", row.names = 1)

sel.num <- 1 # selected region
plotDMRs(BS.subset, regions=targets[sel.num,], testCovariate="group", addRegions = targets, 
         main = paste0(targets$name[sel.num], " Targeted Region_5mC"),
         extend = 3000, annoTrack = annoTrack, highlightMain = F, qval = F, stat = F, horizLegend = T)

# Th1

sampleNames(BSobj)
BS.subset <- BSobj[, c(1:3, 4:6)]
sampleNames(BS.subset)

targets <- read.csv("manuscript/tables/DMRs/DMRs_5mC_Th1.vs.Th0.csv", row.names = 1)

sel.num <- 1 # selected region
plotDMRs(BS.subset, regions=targets[sel.num,], testCovariate="group", addRegions = targets, 
         main = paste0(targets$name[sel.num], " Targeted Region_5mC"),
         extend = 1000, annoTrack = annoTrack, highlightMain = F, qval = F, stat = F, horizLegend = T)

# Th2

sampleNames(BSobj)
BS.subset <- BSobj[, c(1:3, 12:14)]
sampleNames(BS.subset)

targets <- read.csv("manuscript/tables/DMRs/DMRs_5mC_Th2.vs.Th0.csv", row.names = 1)

sel.num <- 1 # selected region
plotDMRs(BS.subset, regions=targets[sel.num,], testCovariate="group", addRegions = targets, 
         main = paste0(targets$name[sel.num], " Targeted Region_5mC"),
         extend = 2000, annoTrack = annoTrack, highlightMain = F, qval = F, stat = F, horizLegend = T)

sel.num <- 5 # selected region
plotDMRs(BS.subset, regions=targets[sel.num,], testCovariate="group", addRegions = targets, 
         main = paste0(targets$name[sel.num], " Targeted Region_5mC"),
         extend = 8000, annoTrack = annoTrack, highlightMain = F, qval = F, stat = F, horizLegend = T)


# Th17

sampleNames(BSobj)
BS.subset <- BSobj[, c(1:3, 7:9)]
sampleNames(BS.subset)

targets <- read.csv("manuscript/tables/DMRs/DMRs_5mC_Th17.vs.Th0.csv", row.names = 1)

sel.num <- 1 # selected region
plotDMRs(BS.subset, regions=targets[sel.num,], testCovariate="group", addRegions = targets, 
         main = paste0(targets$name[sel.num], " Targeted Region_5mC"),
         extend = 1000, annoTrack = annoTrack, highlightMain = F, qval = F, stat = F, horizLegend = T)

sel.num <- 3 # selected region
plotDMRs(BS.subset, regions=targets[sel.num,], testCovariate="group", addRegions = targets, 
         main = paste0(targets$name[sel.num], " Targeted Region_5mC"),
         extend = 1000, annoTrack = annoTrack, highlightMain = F, qval = F, stat = F, horizLegend = T)

sel.num <- 4 # selected region
plotDMRs(BS.subset, regions=targets[sel.num,], testCovariate="group", addRegions = targets, 
         main = paste0(targets$name[sel.num], " Targeted Region_5mC"),
         extend = 10000, annoTrack = annoTrack, highlightMain = F, qval = F, stat = F, horizLegend = T)


# Th17 in absence of TGFb

sampleNames(BSobj)
BS.subset <- BSobj[, c(10:11, 7:9)]
sampleNames(BS.subset)

targets <- read.csv("manuscript/tables/DMRs/DMRs_5mC_Th17.no.TGFb.vs.Th17.csv", row.names = 1)

sel.num <- 1 # selected region
plotDMRs(BS.subset, regions=targets[sel.num,], testCovariate="group", addRegions = targets, 
         main = paste0(targets$name[sel.num], " Targeted Region_5mC"),
         extend = 1000, annoTrack = annoTrack, highlightMain = F, qval = F, stat = F, horizLegend = T)

sel.num <- 2 # selected region
plotDMRs(BS.subset, regions=targets[sel.num,], testCovariate="group", addRegions = targets, 
         main = paste0(targets$name[sel.num], " Targeted Region_5mC"),
         extend = 4000, annoTrack = annoTrack, highlightMain = F, qval = F, stat = F, horizLegend = T)



# Th17 in absence of TGFb vs Th0

sampleNames(BSobj)
BS.subset <- BSobj[, c(10:11, 1:3)]
sampleNames(BS.subset)

targets <- read.csv("manuscript/tables/DMRs/DMRs_5mC_Th17.no.Tgfb.vs.Th0.csv", row.names = 1)

sel.num <- 1 # selected region
plotDMRs(BS.subset, regions=targets[sel.num,], testCovariate="group", addRegions = targets, 
         main = paste0(targets$name[sel.num], " Targeted Region_5mC"),
         extend = 3000, annoTrack = annoTrack, highlightMain = F, qval = F, stat = F, horizLegend = T)

sel.num <- 2 # selected region
plotDMRs(BS.subset, regions=targets[sel.num,], testCovariate="group", addRegions = targets, 
         main = paste0(targets$name[sel.num], " Targeted Region_5mC"),
         extend = 1000, annoTrack = annoTrack, highlightMain = F, qval = F, stat = F, horizLegend = T)




# Targeted heatmap -----------------------------------------------------------------

rm(list=ls())

load("manuscript/BSobj_5mCpG.RData")

# cell type specific DMRs

targets.1 <- read.csv("manuscript/tables/DMRs/DMRs_5mC_Treg.vs.Th0.csv", row.names = 1)
targets.2 <- read.csv("manuscript/tables/DMRs/DMRs_5mC_Th1.vs.Th0.csv", row.names = 1)
targets.3 <- read.csv("manuscript/tables/DMRs/DMRs_5mC_Th2.vs.Th0.csv", row.names = 1)
targets.4 <- read.csv("manuscript/tables/DMRs/DMRs_5mC_Th17.vs.Th0.csv", row.names = 1)
targets.5 <- read.csv("manuscript/tables/DMRs/DMRs_5mC_Th17.no.Tgfb.vs.Th0.csv", row.names = 1)
targets.6 <- read.csv("manuscript/tables/DMRs/DMRs_5mC_Th17.no.TGFb.vs.Th17.csv", row.names = 1)

targets <- rbind(targets.1, targets.2, targets.3, targets.4, targets.5, targets.6)
nrow(targets) # 23

# smoothing doesn't work

heatmap.data <- bsseq::getMeth(BSseq = BSobj,
                               regions = GRanges(targets),
                               type = "raw",
                               what = "perRegion")
head(heatmap.data)
#rownames(heatmap.data) <- make.names(targets$name, unique = T)
rownames(heatmap.data) <- targets$name

heatmap.data %>%   na.omit() %>% as.matrix()

pheno <- as.data.frame(pData(BSobj))
pheno$group <- as.factor(pheno$group)
pheno <- pheno[, 2, drop = F]
levels(pheno$group)
head(pheno)

group.colors <- alphabet(20)
group.colors <- group.colors[c(5,18,15,19,3,7)]
names(group.colors) <- levels(pheno$group)
group.colors =  list(group = group.colors)

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




# end ---------------------------------------------------------------------
sessionInfo()

