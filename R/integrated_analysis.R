

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
})


annoTrack <- getAnnot("mm10")

genome <- BSgenome.Mmusculus.UCSC.mm10

genes <- annotateTranscripts(TxDb.Mmusculus.UCSC.mm10.knownGene)


setwd("~/Dropbox/BioInfo/Lab/Tcell_ID")



# prepare phenotype data -----------------------------------------------------------

bed.files <- list.files(path = "./remora/", pattern = "5mC.bed", recursive = T, full.names = F)

pdata <- as.data.frame(str_split(bed.files, "/", simplify = T))
colnames(pdata) <- c("run","group","basename")
pdata$bed.file <- bed.files
pdata$origin <- c(rep("in_vitro", 2), "ex_vivo", rep("in_vitro", 2), "ex_vivo", rep("in_vitro", 14), rep("ex_vivo", 2))
pdata$protocol <- c(rep("Cas9", 20), rep("uncalled", 2))
rownames(pdata) <- paste0(pdata$group, "_", pdata$run)

head(pdata)

write.csv(pdata, file = "pdata.5mC.integrated.csv")


# load bed files ----------------------------------------------------------

bed.files.full <- list.files(path = "./remora/", pattern = "5mC.bed", recursive = T, full.names = T)


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
save(bed.list, file="bed.list.5mC.integrated.RData")
save(bed.list, file="bed.list.5mC.integrated.all.genome.RData")



# inspection --------------------------------------------------------------

load("bed.list.5mC.integrated.RData")
pdata <- read.csv("pdata.5mC.integrated.csv", row.names = 1)

par(mar = c(10,10,8,10), mfrow = c(2,1))
barplot(as.numeric(lapply(bed.list, length)), las = 2,
        names.arg = names(bed.list),
        main = "Total reads per sample",
        col = viridis(4))

barplot(as.numeric(lapply(bed.list, function(x) sum(x$blockCount))), las = 2,
        names.arg = names(bed.list),
        main = "Total counts per sample",
        col = viridis(4))


# methylation distribution

block.sizes <- lapply(bed.list, function(x) as.numeric(x$blockSizes))
names(block.sizes)

block.sizes.df <- data.frame(Reduce(cbind, block.sizes))
colnames(block.sizes.df) <- names(block.sizes)
block.sizes.df <- as.matrix(block.sizes.df)
block.sizes.df <- block.sizes.df/100 # to get beta values
head(block.sizes.df)

#pdata <- pdata[order(pdata$group), ]
#block.sizes.df <- block.sizes.df[, rownames(pdata)]

par(mar = c(10,10,8,10), mfrow = c(1,1))
densityBeanPlot(block.sizes.df, 
#                sampGroups = pdata$group,
#                sampNames = rownames(pdata),
                numPositions = 1000)




# prepare BSseq object ----------------------------------------------------

load("bed.list.5mC.integrated.RData")
pdata <- read.csv("pdata.5mC.integrated.csv", row.names = 1)

length(bed.list)
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

save(BSobj, file= "BSobj_5mCpG_integrated.RData")
save(BSobj, file= "BSobj_5mCpG_integrated.all.genome.RData")



# PCA ---------------------------------------------------------------------

# uncalled Th17 cluster way too far in preliminary tests
# DMK samples are next in clustering apart

meth.data <- bsseq::getMeth(BSseq = BSobj,
                               regions = NULL,
                               type = "raw",
                               what = "perBase")
head(meth.data)
dim(meth.data)

all(colnames(meth.data) == rownames(pdata))

sel.pdata <- pdata[-c(21:22), ]
sel.pdata <- sel.pdata[!sel.pdata$group == "DMK", ]

sel.meth.data <- meth.data[, rownames(sel.pdata)]
sel.meth.data <- na.omit(sel.meth.data)
head(sel.meth.data)
dim(sel.meth.data)

t.counts <- t(sel.meth.data)
t.counts[1:5,1:10]

pca_res <- prcomp(t.counts, scale. = F)


library(ggfortify)
autoplot(pca_res, data = sel.pdata, colour = 'group', size = 3) +
  theme_minimal() +
#  scale_color_manual(values = group.colors) +
  ggtitle("PCA integrated data")




# DMRplots ----------------------------------------------------------------

load("BSobj_5mCpG_integrated.RData")

regions <- read.csv("targeted.regions.csv", row.names = 1)
head(regions)

pheno <- pData(BSobj)
pheno$group[19] <- "Th17_noTGFb"
pheno$group[20] <- "Th17_noTGFb"

table(pheno$group)

library(pals)
group.colors <- alphabet(10)
group.colors <- rainbow(10)
group.colors <- c(rainbow(10), alphabet(10))
library(scales)
show_col(group.colors)
show_col(group.colors[c(1,2,3,6,8,9,15,13,19,20)])

pheno$col <- pheno$group
pheno$col <- gsub("DMK", group.colors[1], pheno$col)
pheno$col <- gsub("Th0", group.colors[2], pheno$col)
pheno$col <- gsub("Th17_noTGFb", group.colors[3], pheno$col)
pheno$col <- gsub("Th17_TGFb.KO", group.colors[6], pheno$col)
pheno$col <- gsub("Th17_TGFb.WT", group.colors[8], pheno$col)
pheno$col <- gsub("Th2", group.colors[9], pheno$col)
pheno$col <- gsub("tTreg", group.colors[15], pheno$col)
pheno$col <- gsub("Treg", group.colors[13], pheno$col)
pheno$col <- gsub("Th17", group.colors[19], pheno$col)
pheno$col <- gsub("Th1", group.colors[20], pheno$col)

pData(BSobj) <- pheno

# test
sel.num <- 11 # selected region
plotDMRs(BSobj, regions=regions[sel.num,], testCovariate="group", addRegions = NULL, 
         main = paste0(rownames(regions)[sel.num], " Targeted Region_5mC"),
         extend = 0, annoTrack = annoTrack, highlightMain = F, qval = F, stat = F, horizLegend = T)






# heatmap -----------------------------------------------------------------

load("BSobj_5mCpG_integrated.RData")

hasBeenSmoothed(BSobj)


# smoothed.BSoj <- BSmooth(BSobj)
# smoothing this object gives an error

heatmap.data <- bsseq::getMeth(BSseq = BSobj,
                               regions = GRanges(regions),
                               type = "raw",
                               what = "perRegion")
head(heatmap.data)
rownames(heatmap.data) <- rownames(regions)

#heatmap.data %>%   na.omit() %>% as.matrix()

pheatmap::pheatmap(heatmap.data,
                   scale = "row",
                   annotation_col =  as.data.frame(pData(BSobj)[, c(2,5)]),
                   show_rownames = T, fontsize_row = 12,
                   show_colnames = T,
                   angle_col = 45,
                   border_color = "grey",
                   main = "Z-Scores of 5mC in all targeted regions",
                   fontsize = 12,
                   cellwidth = 25,
                   cellheight = 20)



# differential methylation (5mC) --------------------------------------------------------

load("BSobj_5mCpG_integrated.RData")
load("BSobj_5mCpG_integrated.all.genome.RData")

pheno <- as.data.frame(pData(BSobj))
pheno$group

mParam = MulticoreParam(workers= 8, progressbar=TRUE)

# DMK vs Treg
dml.test <- DMLtest(BSobj, group1 = c(1,4), group2 = c(2,5), BPPARAM = mParam) # as a control for known differences (Valerie's paper), although using remora instead of rerio
# Treg vs tTreg
dml.test <- DMLtest(BSobj, group1 = c(2,5), group2 = c(3,6), BPPARAM = mParam) # as a control for known differences (Valerie's paper), although using remora instead of rerio
# Treg vs Th0
dml.test <- DMLtest(BSobj, group1 = c(2,5), group2 = c(7,11,15), BPPARAM = mParam) # 
# Th1 vs Th0
dml.test <- DMLtest(BSobj, group1 = c(8,12,16), group2 = c(7,11,15), BPPARAM = mParam) # 
# Th2 vs Th0
dml.test <- DMLtest(BSobj, group1 = c(10,14,18), group2 = c(7,11,15), BPPARAM = mParam) # 
# Th17 vs Th0
dml.test <- DMLtest(BSobj, group1 = c(9,13,17), group2 = c(7,11,15), BPPARAM = mParam) # 
# Th17.no.TGFb vs Th17
dml.test <- DMLtest(BSobj, group1 = c(19,20), group2 = c(9,13,17), BPPARAM = mParam) # 
# Th17.TGFb.KO vs Th17.TGFb.WT
dml.test <- DMLtest(BSobj, group1 = c(21), group2 = c(22), BPPARAM = mParam, smoothing = T) # when no replicates (uncalled samples)

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


# DMK vs Treg
write.csv(dmls, file = "DMLs_5mC_DMK.vs.Treg.csv")
write.csv(dmrs, file = "DMRs_5mC_DMK.vs.Treg.csv")
# Treg vs tTreg
write.csv(dmls, file = "DMLs_5mC_Treg.vs.tTreg.csv")
write.csv(dmrs, file = "DMRs_5mC_Treg.vs.tTreg.csv")
# Treg vs Th0
write.csv(dmls, file = "DMLs_5mC_Treg.vs.Th0.csv")
write.csv(dmrs, file = "DMRs_5mC_Treg.vs.Th0.csv")
# Th1 vs Th0
write.csv(dmls, file = "DMLs_5mC_Th1.vs.Th0.csv")
write.csv(dmrs, file = "DMRs_5mC_Th1.vs.Th0.csv")
# Th2 vs Th0
write.csv(dmls, file = "DMLs_5mC_Th2.vs.Th0.csv")
write.csv(dmrs, file = "DMRs_5mC_Th2.vs.Th0.csv")
# Th17 vs Th0
write.csv(dmls, file = "DMLs_5mC_Th17.vs.Th0.csv")
write.csv(dmrs, file = "DMRs_5mC_Th17.vs.Th0.csv")
# Th17.no.TGFb vs Th17
write.csv(dmls, file = "DMLs_5mC_Th17.no.TGFb.vs.Th17.csv")
write.csv(dmrs, file = "DMRs_5mC_Th17.no.TGFb.vs.Th17.csv")
# Th17.TGFb.KO vs Th17.TGFb.WT
write.csv(dmls, file = "DMLs_5mC_Th17.KO.vs.WT.csv")
write.csv(dmrs, file = "DMRs_5mC_Th17.KO.vs.WT.csv")




# Targeted DMRplots ----------------------------------------------------------------

load("BSobj_5mCpG_integrated.RData")

regions <- read.csv("targeted.regions.csv", row.names = 1)
head(regions)

pheno <- pData(BSobj)
pheno$group[19] <- "Th17_noTGFb"
pheno$group[20] <- "Th17_noTGFb"


# effect of TGFb during and after Th17 differentiation

BS.subset <- BSobj[, c(9, 13, 17, 19:22)]

targets.1 <- read.csv("DMRs_5mC_Th17.no.TGFb.vs.Th17.csv", row.names = 1)
targets.2 <- read.csv("DMRs_5mC_Th17.KO.vs.WT.csv", row.names = 1)

targets <- rbind(targets.1, targets.2)

sel.num <- 13 # selected region
plotDMRs(BS.subset, regions=targets[sel.num,], testCovariate="group", addRegions = targets, 
         main = paste0(targets$name[sel.num], " Targeted Region_5mC"),
         extend = 1000, annoTrack = annoTrack, highlightMain = F, qval = F, stat = F, horizLegend = T)



# DMK experiment

BS.subset <- BSobj[, c(1:6)]

sel.num <- 1 # selected region
plotDMRs(BS.subset, regions=regions[sel.num,], testCovariate="group", addRegions = targets, 
         main = paste0(rownames(regions)[sel.num], " Targeted Region_5mC"),
         extend = 1000, annoTrack = annoTrack, highlightMain = F, qval = F, stat = F, horizLegend = T)

all.regions <- data.frame(chr="chrX",start=7565986, end=7593799)
plotDMRs(BS.subset, regions=all.regions, testCovariate="group", addRegions = NULL, regionCol = alpha("orange", 0.2),
         extend = 1000, annoTrack = annoTrack, highlightMain = F, qval = F, stat = F, horizLegend = T, addLines = T)





# Targeted heatmap -----------------------------------------------------------------

load("BSobj_5mCpG_integrated.RData")

# cell type specific DMRs

targets.1 <- read.csv("DMRs_5mC_Treg.vs.Th0.csv", row.names = 1)
targets.2 <- read.csv("DMRs_5mC_Th1.vs.Th0.csv", row.names = 1)
targets.3 <- read.csv("DMRs_5mC_Th2.vs.Th0.csv", row.names = 1)
targets.4 <- read.csv("DMRs_5mC_Th17.vs.Th0.csv", row.names = 1)

targets <- rbind(targets.1, targets.2, targets.3, targets.4)

# smoothing doesn't work

heatmap.data <- bsseq::getMeth(BSseq = BSobj,
                               regions = GRanges(targets),
                               type = "raw",
                               what = "perRegion")
head(heatmap.data)
rownames(heatmap.data) <- targets$name

heatmap.data %>%   na.omit() %>% as.matrix()

pheatmap::pheatmap(heatmap.data,
                   scale = "row",
                   annotation_col =  as.data.frame(pData(BSobj)[, c(2,5)]),
                   show_rownames = T, fontsize_row = 12,
                   show_colnames = T,
                   angle_col = 45,
                   border_color = "grey",
                   main = "Mean methylation of cell type specific DMRs (5mC)\n",
                   fontsize = 14,
                   cellwidth = 20,
                   cellheight = 20)





# end ---------------------------------------------------------------------
sessionInfo()

