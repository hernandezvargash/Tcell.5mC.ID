

# main --------------------------------------------------------------------

# in vitro differentiation experiments:

# 20210113: second set of samples, replicates from the AK28 set
# 20210312: second set of samples sent, replicates from the AKXX set
# 20210317: second set of samples sent, replicates from the ARJP28 set

# Crispr targeted LIGATION sequencing SQK-LSK109
# EXP-NBD104 barcoding
# RB1=Th0, RB2=Th1, RB3=Th2, RB4=TH17


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


setwd("~/Dropbox/BioInfo/Lab/Tcell_ID")



# prepare phenotype data -----------------------------------------------------------

bed.files <- list.files(path = "./remora/", pattern = "5hmC.bed", recursive = T, full.names = F)

pdata <- as.data.frame(str_split(bed.files, "/", simplify = T))
colnames(pdata) <- c("run","group","basename")
pdata$bed.file <- bed.files
pdata$assay <- rep(paste0("Assay.", 1:3), each = 4)
rownames(pdata) <- paste0(pdata$group, "_", pdata$assay)

head(pdata)

write.csv(pdata, file = "pdata.5hmC.csv")


# load bed files ----------------------------------------------------------

bed.files.full <- list.files(path = "./remora/", pattern = "5hmC.bed", recursive = T, full.names = T)


# inspect 1st file

test <- import.bed(bed.files.full[1], which = NULL)
test <- keepStandardChromosomes(test, pruning.mode = "coarse")
test <- sort(test)
test <- resize(test, width = 3, fix = "start")
test$context <- getSeq(genome, test)
test$blockSizes <- as.numeric(as.character(test$blockSizes))
hist(test$blockCount, breaks = 100) # coverage
hist(test$blockSizes, breaks = 100) # methylation


# import to a list

bed.list <- list()

#region = GRanges("chr1:1-248956422")
#region =  GRanges("chrM:1-16569")

for(i in 1:length(bed.files)){
  bed.list[[i]] <- import.bed(bed.files.full[i], which = NULL)
  bed.list[[i]] <- keepStandardChromosomes(bed.list[[i]], pruning.mode = "coarse")
  bed.list[[i]] <- sort(bed.list[[i]])
  bed.list[[i]] <- resize(bed.list[[i]], width = 3, fix = "start")
  bed.list[[i]]$context <- getSeq(genome, bed.list[[i]])
  bed.list[[i]]$blockSizes <- as.numeric(as.character(bed.list[[i]]$blockSizes))
}

names(bed.list) <- rownames(pdata)
save(bed.list, file="bed.list.5hmC.RData")



# inspection --------------------------------------------------------------

load("bed.list.5hmC.RData")
pdata <- read.csv("pdata.5hmC.csv", row.names = 1)

par(mar = c(10,10,8,10), mfrow = c(2,1))
barplot(as.numeric(lapply(bed.list, length)), las = 2,
        names.arg = names(bed.list),
        main = "Total reads per sample",
        col = viridis(4))

barplot(as.numeric(lapply(bed.list, function(x) sum(x$blockCount))), las = 2,
        names.arg = names(bed.list),
        main = "Total counts per sample",
        col = viridis(4))


# global methylation

block.sizes <- lapply(bed.list, function(x) as.numeric(x$blockSizes))
names(block.sizes)

means.Th0 <- unlist(lapply(block.sizes, function(x) mean(x)))[c(1,5,9)]
means.Th1 <- unlist(lapply(block.sizes, function(x) mean(x)))[c(2,6,10)]
means.Th17 <- unlist(lapply(block.sizes, function(x) mean(x)))[c(3,7,11)]
means.Th2 <- unlist(lapply(block.sizes, function(x) mean(x)))[c(4,8,12)]

boxplot(means.Th0, means.Th1, means.Th2,means.Th17,
        col = viridis(4),
        names = c("Th0","Th1","Th2","Th17"),
        ylim = c(1,2),
        main = "mean methylation per condition",
        ylab = "5mCpG methylation [%]")


t.test(means.Th0, means.Th1, paired = F)
# p-value = 0.5551



# prepare BSseq object ----------------------------------------------------

load("bed.list.5hmC.RData")
pdata <- read.csv("pdata.5hmC.csv", row.names = 1)

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
save(BSobj, file= "BSobj_Th_5hmCpG.RData")



# fit DSS multi-factor --------------------------------------------------------

load("BSobj_Th_5hmCpG.RData")

genes <- annotateTranscripts(TxDb.Mmusculus.UCSC.mm10.knownGene)

pheno <- as.data.frame(pData(BSobj))
formula=~0+pheno$group+pheno$assay
model.matrix(formula)

DMLfit = DMLfit.multiFactor(BSobj, design= pheno, formula=formula, smoothing=F)

save(DMLfit, file="DMLfit_5hmC.RData")



# Th1 vs Th0 --------------------------------------------------------------

load("DMLfit_5hmC.RData")

head(DMLfit$X)

Contrast = matrix(c(-1,1,0,0,0,0), ncol=1)
DMLtest.Th1.vs.Th0 = DMLtest.multiFactor(DMLfit, Contrast=Contrast)
head(DMLtest.Th1.vs.Th0[order(DMLtest.Th1.vs.Th0$fdrs, decreasing = F), ])
table(DMLtest.Th1.vs.Th0$fdrs < 0.05) # 0


# Th2 vs Th0 --------------------------------------------------------------

load("DMLfit_5hmC.RData")

head(DMLfit$X)

Contrast = matrix(c(-1,0,0,1,0,0), ncol=1)
DMLtest.Th1.vs.Th0 = DMLtest.multiFactor(DMLfit, Contrast=Contrast)
head(DMLtest.Th1.vs.Th0[order(DMLtest.Th1.vs.Th0$fdrs, decreasing = F), ])
table(DMLtest.Th1.vs.Th0$fdrs < 0.05) # 1

DMLtest.Th1.vs.Th0.sites <- DMLtest.Th1.vs.Th0[DMLtest.Th1.vs.Th0$fdrs<0.05, ]
DMLtest.Th1.vs.Th0.sites <- na.omit(DMLtest.Th1.vs.Th0.sites)
gr1 <- DMLtest.Th1.vs.Th0.sites
colnames(gr1)[2] <- "start"
gr1$end <- gr1$start
gr1 <- GRanges(gr1)
match1 <- matchGenes(gr1, genes, type = "fiveprime", promoterDist = 2500, skipExons = FALSE, verbose = TRUE)
DMLtest.Th1.vs.Th0.sites <- cbind(DMLtest.Th1.vs.Th0.sites, match1)
head(DMLtest.Th1.vs.Th0.sites) # Il4

DMLtest.Th1.vs.Th0.regions <- callDMR(DMLtest.Th1.vs.Th0, p.threshold = 0.05, minCG = 3, dis.merge = 50, minlen = 10, pct.sig = 0.5)
match2 <- matchGenes(DMLtest.Th1.vs.Th0.regions, genes, type = "fiveprime", promoterDist = 2500, skipExons = FALSE, verbose = TRUE)
DMLtest.Th1.vs.Th0.regions <- cbind(DMLtest.Th1.vs.Th0.regions, match2)
head(DMLtest.Th1.vs.Th0.regions) # 3 (Il4)

save(DMLtest.Th1.vs.Th0, file="DMLtest_5hmC.Th2.vs.Th0.rData")
write.csv(DMLtest.Th1.vs.Th0.sites, file="DMLtest_5hmC.Th2.vs.Th0.sites.csv")
write.csv(DMLtest.Th1.vs.Th0.regions, file="DMLtest_5hmC.Th2.vs.Th0.regions.csv")


# Th17 vs Th0 --------------------------------------------------------------

load("DMLfit_5hmC.RData")

head(DMLfit$X)

Contrast = matrix(c(-1,0,1,0,0,0), ncol=1)
DMLtest.Th1.vs.Th0 = DMLtest.multiFactor(DMLfit, Contrast=Contrast)
head(DMLtest.Th1.vs.Th0[order(DMLtest.Th1.vs.Th0$fdrs, decreasing = F), ])
table(DMLtest.Th1.vs.Th0$fdrs < 0.05) # 0




# DMRplots ----------------------------------------------------------------

load("BSobj_Th_5hmCpG.RData")

# only significant 5hmC loci was Il4 in Th2 vs Th0 comparison

regions <- read.csv("targeted.regions.csv", row.names = 1)
regions <- regions[-11, ] # remove Itgb8
head(regions)

# test
plotDMRs(BSobj, regions=regions[1,], testCovariate="group", addRegions = regions, main = "test",
         extend = 5000, annoTrack = annoTrack, highlightMain = F, qval = F, stat = F, horizLegend = T)

for(i in 1:nrow(regions)){
  jpeg(filename = paste0("Targeted_Region_5hmC_", rownames(regions)[i], ".jpeg") , width = 960, height = 960, quality = 100)
  plotDMRs(BSobj, regions=regions[i,], testCovariate="group", addRegions = regions,
           main = paste0("5hmC: ", rownames(regions)[i], " locus"),
           extend = 5000, annoTrack = annoTrack, highlightMain = F, qval = F, stat = F, horizLegend = T)
  dev.off()
}


gr <- GRanges(regions[7, ])

plot(plotDMRs(BSobj, regions=regions[7,], testCovariate="group", addRegions = regions, main = "test",
         extend = 5000, annoTrack = annoTrack, highlightMain = F, qval = F, stat = F, horizLegend = T), 
     ylim = c(0,0.1), xlim = c(1, 10))




# heatmap -----------------------------------------------------------------


hasBeenSmoothed(BSobj)
smoothed.BSoj <- BSmooth(BSobj)

# save(smoothed.BSoj, file = "smoothed.5hmC.BSobj.RData")

# DMRichR::smoothPheatmap not working:
# smoothPheatmap(smoothed.BSoj, GRanges(all.DMRs), testCovariate = "group", filename = NULL)


load("smoothed.5hmC.BSobj.RData")

heatmap.data <- bsseq::getMeth(BSseq = smoothed.BSoj,
                               regions = GRanges(regions),
                               type = "smooth",
                               what = "perRegion")
head(heatmap.data)
rownames(heatmap.data) <- rownames(regions)

heatmap.data %>%   na.omit() %>% as.matrix()

pheatmap::pheatmap(heatmap.data,
                   scale = "row",
                   annotation_col =  as.data.frame(pData(smoothed.BSoj)[, c(2,5)]),
                   show_rownames = T, fontsize_row = 12,
                   show_colnames = T,
                   angle_col = 45,
                   border_color = "grey",
                   main = "Z-Scores of 5hmC in all targeted regions",
                   fontsize = 12,
                   cellwidth = 25,
                   cellheight = 20)



# combined 5mC/5hmC heatmap -----------------------------------------------------------------

regions <- read.csv("targeted.regions.csv", row.names = 1)
regions <- regions[-11, ] # remove Itgb8
head(regions)

load(file = "smoothed.BSobj.RData")

heatmap.data.5mC <- bsseq::getMeth(BSseq = smoothed.BSoj,
                               regions = GRanges(regions),
                               type = "smooth",
                               what = "perRegion")

load("smoothed.5hmC.BSobj.RData")

heatmap.data.5hmC <- bsseq::getMeth(BSseq = smoothed.BSoj,
                               regions = GRanges(regions),
                               type = "smooth",
                               what = "perRegion")

rownames(heatmap.data.5mC) <- paste0(rownames(regions),".5mC")
head(heatmap.data.5mC)
rownames(heatmap.data.5hmC) <- paste0(rownames(regions),".5hmC")
head(heatmap.data.5hmC)

# 5mC/5hmC aggregation

heatmap.data <- rbind(heatmap.data.5mC, heatmap.data.5hmC)

heatmap.data %>%   na.omit() %>% as.matrix()

pheatmap::pheatmap(heatmap.data,
                   scale = "row",
                   annotation_col =  as.data.frame(pData(smoothed.BSoj)[, c(2,5)]),
                   show_rownames = T, fontsize_row = 12,
                   show_colnames = T,
                   angle_col = 45,
                   border_color = "grey",
                   main = "Z-Scores of 5mC & 5hmC in all targeted regions",
                   fontsize = 12,
                   cellwidth = 25,
                   cellheight = 20)

# 5hmC/5mC "activity" ratio

heatmap.data <- heatmap.data.5hmC/ heatmap.data.5mC
rownames(heatmap.data) <- rownames(regions)
head(heatmap.data)

heatmap.data %>%   na.omit() %>% as.matrix()

pheatmap::pheatmap(heatmap.data,
                   scale = "row",
                   annotation_col =  as.data.frame(pData(smoothed.BSoj)[, c(2,5)]),
                   show_rownames = T, fontsize_row = 12,
                   show_colnames = T,
                   angle_col = 45,
                   border_color = "grey",
                   main = "'Activity' 5hmC/5mC ratio in all targeted regions",
                   fontsize = 12,
                   cellwidth = 25,
                   cellheight = 20)




# end ---------------------------------------------------------------------
sessionInfo()


