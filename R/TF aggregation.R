

# main --------------------------------------------------------------------

# TF aggregation

# The goal is to select regions (from off-target data) based on TF binding motifs and aggregate 5mC and 5hmC values for each TF
# Cas9-targeted loci: Foxp3, Il17a, Rorc, Ifng, Il10, Maf, Il4, Gata3, Zfp362, Tbx21, Itgb8
# TFs of interest: Foxp3, Rorc, Stat3, Maf, Gata3, Zfp362, Tbx21, Bcl6, Runx3, Hif1alpha
# including known lineage-specific master TFs (T-BET, GATA3, and ROR-γt) and STATs (STAT1 and STAT4, STAT6, and STAT3)

# cistromes from: https://bmcresnotes.biomedcentral.com/articles/10.1186/s13104-018-3856-x
# https://doi.org/10.6084/m9.figshare.7087697
# cistromes also used in this paper about pioneer TFs: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8184831/ , from that reference:
# Protective Pioneer Factors (PPFs): CTCF, REST, KLF4, KLF7, SOX2, SOX9, N-MYC, NRF1, OTX2, E2F1
# Super Pioneer Factors (SPFs): CTCF, REST, KLF4, SOX2, SOX9, SOX17, CREB, FOXA1, FOXD3, E2F1, N-MYC, GR


# libraries ---------------------------------------------------------------

rm(list=ls())

setwd("~/Dropbox/BioInfo/Lab/Tcell_ID")

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
  library(data.table)
  library(MIRA)
  
})

# The data set contains the cistrome genomic map (human, mm10 genome assembly) as standard genomic regions files (.bed). 
# Inside the zip archive, there are two folders (mm10_cistrome , mm10_cismotifs) corresponding to the base cistrome and the motif-annotated subset. 
# Inside each folder there are *.bed files names as <TF>.<reliability>.bed 
# where <TF> is the UniProt ID of a particular transcription factor and <reliability> is the reliability category from A (the highest) to D (limited). 
# The bed-files within mm10_cismotifs contain an additional column which lists the maximal HOCOMOCO motif score reachable for a particular cistrome segment 
# (-log10 scale, the filtering was performed using the score of 4 which corresponds to the motif P-value of 0.0001).

regions <- import.bed("cistrome_mm10/mm10_cistrome/GATA3_MOUSE.A.bed")
regions <- import.bed("cistrome_mm10/mm10_cistrome/TBX21_MOUSE.A.bed")
regions <- import.bed("cistrome_mm10/mm10_cistrome/FOXP3_MOUSE.A.bed")
regions <- import.bed("cistrome_mm10/mm10_cistrome/STAT3_MOUSE.A.bed")
regions <- import.bed("cistrome_mm10/mm10_cistrome/BCL6_MOUSE.A.bed")
regions <- import.bed("cistrome_mm10/mm10_cistrome/RUNX3_MOUSE.A.bed")
regions <- import.bed("cistrome_mm10/mm10_cistrome/STAT1_MOUSE.A.bed")
regions <- import.bed("cistrome_mm10/mm10_cistrome/KLF6_MOUSE.C.bed")
regions
mean(width(regions))
expanded.regions <- resize(regions, 2000, fix="center")
mean(width(expanded.regions))


# prepare phenotype data -----------------------------------------------------------

bed.files <- list.files(path = "./remora/", pattern = "5mC.bed", recursive = T, full.names = F)
bed.files <- bed.files[1:18] # omitting uncalled experiment

pdata <- as.data.frame(str_split(bed.files, "/", simplify = T))
colnames(pdata) <- c("run","group","basename")
pdata$bed.file <- bed.files
pdata$origin <- c(rep("in_vitro", 2), "ex_vivo", rep("in_vitro", 2), "ex_vivo", rep("in_vitro", 12))
pdata$protocol <- c(rep("Cas9", 18))
rownames(pdata) <- paste0(pdata$group, "_", pdata$run)

head(pdata)
colnames(pdata)[2] <- "sampleType"
table(pdata$sampleType)


# load bed files ----------------------------------------------------------

bed.files.full <- list.files(path = "./remora/", pattern = "5mC.bed", recursive = T, full.names = T)


# import to a list

bed.list <- list()

for(i in 1:18){ # omitting uncalled experiment
  bed.list[[i]] <- import.bed(bed.files.full[i], which = expanded.regions)
  bed.list[[i]] <- keepStandardChromosomes(bed.list[[i]], pruning.mode = "coarse")
  bed.list[[i]] <- sort(bed.list[[i]])
#  bed.list[[i]] <- resize(bed.list[[i]], width = 3, fix = "start")
#  bed.list[[i]]$context <- getSeq(genome, bed.list[[i]])
  bed.list[[i]]$blockSizes <- as.numeric(as.character(bed.list[[i]]$blockSizes))
}

names(bed.list) <- rownames(pdata)

lapply(bed.list, length)

#save(bed.list, file="bed.list.5mC.integrated.RData")

head(bed.list[[1]])
hist(bed.list[[1]]$blockCount)
hist(bed.list[[1]]$blockSizes)


# prepare data tables

pre.DT <- lapply(bed.list, function(x) as.data.frame(x))
head(pre.DT[[1]])
table(pre.DT[[1]]$coverage > 0)

pre.DT <- lapply(pre.DT, function(x) x[,c("seqnames","start","blockSizes","blockCount")])
pre.DT <- lapply(pre.DT, function(x) {colnames(x)<-c("chr", "start", "methylProp", "coverage");x})
pre.DT <- lapply(pre.DT, function(x)  transform(x, methylProp = methylProp/100))

BSDTList <- lapply(pre.DT, data.table)

# save(BSDTList, file = "BSDTList.RData")



# aggregate methylation data ----------------------------------------------------------

# load("BSDTList.RData")

bigBin <- lapply(X=BSDTList, FUN=aggregateMethyl, GRList=expanded.regions, binNum=21, minBaseCovPerBin = 0)
bigBin <- lapply(bigBin, function(x) transform(x, featureID = "GATA3") )

bigBinDT <- rbindNamedList(bigBin)
setkey(bigBinDT, sampleName)

annotDT <- data.table(sampleName = names(BSDTList), pdata)
setkey(annotDT, sampleName)

bigBinDT2 <- merge(bigBinDT, annotDT, all.x=TRUE)

plotMIRAProfiles(binnedRegDT=bigBinDT2, plotType = "jitter", colBlindOption = T)
#plotMIRAProfiles(binnedRegDT=bigBinDT2)
#bigBinDT2[, methylProp := methylProp - min(methylProp) + .05, by=.(featureID, sampleName)]
sampleScores <- calcMIRAScore(bigBinDT2,
                              shoulderShift="auto",
                              regionSetIDColName="featureID",
                              sampleIDColName="sampleName")
sampleScores
sampleScores$sampleType <- annotDT$sampleType
plotMIRAScores(sampleScores, colBlindOption = T)



# stats -------------------------------------------------------------------


score.1 <- sampleScores[sampleScores$sampleType == "Th0", ]
score.2 <- sampleScores[sampleScores$sampleType == "Th1", ]

wilcox.test(score.1$score, score.2$score)
t.test(score.1$score, score.2$score)



# loop --------------------------------------------------------------------

rm(list = ls())

# phenotype data

bed.files.full <- list.files(path = "./remora/", pattern = "5mC.bed", recursive = T, full.names = T)
bed.files.full <- bed.files.full[1:18] # omitting uncalled experiment
bed.files <- list.files(path = "./remora/", pattern = "5mC.bed", recursive = T, full.names = F)
bed.files <- bed.files[1:18] # omitting uncalled experiment
pdata <- as.data.frame(str_split(bed.files, "/", simplify = T))
colnames(pdata) <- c("run","group","basename")
pdata$bed.file <- bed.files
pdata$origin <- c(rep("in_vitro", 2), "ex_vivo", rep("in_vitro", 2), "ex_vivo", rep("in_vitro", 12))
pdata$protocol <- c(rep("Cas9", 18))
rownames(pdata) <- paste0(pdata$group, "_", pdata$run)
colnames(pdata)[2] <- "sampleType"
head(pdata)


# TF data

TFs.of.interest <- c("GATA3","TBX21","FOXP3","STAT3","BCL6","RUNX3","STAT1")

# list all available mm10 TFs
all.TFs <- list.files("cistrome_mm10/mm10_cistrome/")
# select motifs with highest confidence
all.TFs <- all.TFs[grep("MOUSE.A", all.TFs)]
all.TFs <- str_sub(all.TFs, 1, -13)

# check for available TFs of interest
sel.TFs <- intersect(TFs.of.interest, all.TFs)


for(t in 1:length(sel.TFs)){
  
  regions <- import.bed(paste0("cistrome_mm10/mm10_cistrome/", sel.TFs[t], "_MOUSE.A.bed"))
  expanded.regions <- resize(regions, 2000, fix="center")
  
  bed.list <- list()
  
  for(i in 1:18){ # omitting uncalled experiment
    bed.list[[i]] <- import.bed(bed.files.full[i], which = expanded.regions)
    bed.list[[i]] <- keepStandardChromosomes(bed.list[[i]], pruning.mode = "coarse")
    bed.list[[i]] <- sort(bed.list[[i]])
    bed.list[[i]]$blockSizes <- as.numeric(as.character(bed.list[[i]]$blockSizes))

  }
  
  names(bed.list) <- rownames(pdata)
  
  pre.DT <- lapply(bed.list, function(x) as.data.frame(x))
  pre.DT <- lapply(pre.DT, function(x) x[,c("seqnames","start","blockSizes","blockCount")])
  pre.DT <- lapply(pre.DT, function(x) {colnames(x)<-c("chr", "start", "methylProp", "coverage");x})
  pre.DT <- lapply(pre.DT, function(x)  transform(x, methylProp = methylProp/100))

  BSDTList <- lapply(pre.DT, data.table)
  
  bigBin <- lapply(X=BSDTList, FUN=aggregateMethyl, GRList=expanded.regions, binNum=21, minBaseCovPerBin = 1)
  bigBin <- lapply(bigBin, function(x) transform(x, featureID = sel.TFs[t]) )
  bigBinDT <- rbindNamedList(bigBin)
  setkey(bigBinDT, sampleName)
  annotDT <- data.table(sampleName = names(BSDTList), pdata)
  setkey(annotDT, sampleName)
  
  bigBinDT2 <- merge(bigBinDT, annotDT, all.x=TRUE)
  
  jpeg(filename = paste0("MIRA_Profile_", sel.TFs[t], ".jpeg"), height = 800, width = 800, quality = 100)
  p1 <- plotMIRAProfiles(binnedRegDT=bigBinDT2, plotType = "jitter", colBlindOption = T)
  plot(p1)
  dev.off()

  sampleScores <- calcMIRAScore(bigBinDT2,
                                shoulderShift="auto",
                                regionSetIDColName="featureID",
                                sampleIDColName="sampleName")
  sampleScores
  sampleScores$sampleType <- annotDT$sampleType
  
  
  jpeg(filename = paste0("MIRA_Scores_", sel.TFs[t], ".jpeg"), height = 800, width = 800, quality = 100)
  p2 <- plotMIRAScores(sampleScores, colBlindOption = T)
  plot(p2)
  dev.off()
  
  rm(regions, expanded.regions, bigBinDT2, sampleScores, p1, p2)
    
}



# 5hmC aggregation --------------------------------------------------------------------

rm(list = ls())

# phenotype data

bed.files.full <- list.files(path = "./remora/", pattern = "5hmC.bed", recursive = T, full.names = T)
bed.files.full <- bed.files.full[1:18] # omitting uncalled experiment
bed.files <- list.files(path = "./remora/", pattern = "5hmC.bed", recursive = T, full.names = F)
bed.files <- bed.files[1:18] # omitting uncalled experiment
pdata <- as.data.frame(str_split(bed.files, "/", simplify = T))
colnames(pdata) <- c("run","group","basename")
pdata$bed.file <- bed.files
pdata$origin <- c(rep("in_vitro", 2), "ex_vivo", rep("in_vitro", 2), "ex_vivo", rep("in_vitro", 12))
pdata$protocol <- c(rep("Cas9", 18))
rownames(pdata) <- paste0(pdata$group, "_", pdata$run)
colnames(pdata)[2] <- "sampleType"
head(pdata)


# TF data

TFs.of.interest <- c("GATA3","TBX21","FOXP3","STAT3","BCL6","RUNX3","STAT1")

# list all available mm10 TFs
all.TFs <- list.files("cistrome_mm10/mm10_cistrome/")
# select motifs with highest confidence
all.TFs <- all.TFs[grep("MOUSE.A", all.TFs)]
all.TFs <- str_sub(all.TFs, 1, -13)

# check for available TFs of interest
sel.TFs <- intersect(TFs.of.interest, all.TFs)


for(t in 1:length(sel.TFs)){
  
  regions <- import.bed(paste0("cistrome_mm10/mm10_cistrome/", sel.TFs[t], "_MOUSE.A.bed"))
  expanded.regions <- resize(regions, 2000, fix="center")
  
  bed.list <- list()
  
  for(i in 1:18){ # omitting uncalled experiment
    bed.list[[i]] <- import.bed(bed.files.full[i], which = expanded.regions)
    bed.list[[i]] <- keepStandardChromosomes(bed.list[[i]], pruning.mode = "coarse")
    bed.list[[i]] <- sort(bed.list[[i]])
    bed.list[[i]]$blockSizes <- as.numeric(as.character(bed.list[[i]]$blockSizes))
    
  }
  
  names(bed.list) <- rownames(pdata)
  
  pre.DT <- lapply(bed.list, function(x) as.data.frame(x))
  pre.DT <- lapply(pre.DT, function(x) x[,c("seqnames","start","blockSizes","blockCount")])
  pre.DT <- lapply(pre.DT, function(x) {colnames(x)<-c("chr", "start", "methylProp", "coverage");x})
  pre.DT <- lapply(pre.DT, function(x)  transform(x, methylProp = methylProp/100))
  
  BSDTList <- lapply(pre.DT, data.table)
  
  bigBin <- lapply(X=BSDTList, FUN=aggregateMethyl, GRList=expanded.regions, binNum=21, minBaseCovPerBin = 1)
  bigBin <- lapply(bigBin, function(x) transform(x, featureID = sel.TFs[t]) )
  bigBinDT <- rbindNamedList(bigBin)
  setkey(bigBinDT, sampleName)
  annotDT <- data.table(sampleName = names(BSDTList), pdata)
  setkey(annotDT, sampleName)
  
  bigBinDT2 <- merge(bigBinDT, annotDT, all.x=TRUE)
  
  jpeg(filename = paste0("MIRA_Profile_5hmC_", sel.TFs[t], ".jpeg"), height = 800, width = 800, quality = 100)
  p1 <- plotMIRAProfiles(binnedRegDT=bigBinDT2, plotType = "jitter", colBlindOption = T)
  plot(p1)
  dev.off()
  
  sampleScores <- calcMIRAScore(bigBinDT2,
                                shoulderShift="auto",
                                regionSetIDColName="featureID",
                                sampleIDColName="sampleName")
  sampleScores
  sampleScores$sampleType <- annotDT$sampleType
  
  
  jpeg(filename = paste0("MIRA_Scores_5hmC_", sel.TFs[t], ".jpeg"), height = 800, width = 800, quality = 100)
  p2 <- plotMIRAScores(sampleScores, colBlindOption = T)
  plot(p2)
  dev.off()
  
  rm(regions, expanded.regions, bigBinDT2, sampleScores, p1, p2)
  
}




# end ---------------------------------------------------------------------
sessionInfo()




