

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
  library(bsseq)
  library(SummarizedExperiment)
  library(stringr)
  library(xlsx)
  library(bsseq)
  library(SummarizedExperiment)
  library(DSS)
  library(dmrseq)
  library(DMRcate)
  library(biomaRt)
  library(bumphunter)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
})

annoTrack <- getAnnot("mm10")


setwd("~/Dropbox/BioInfo/Lab/Tcell_ID")



# DMRplots ----------------------------------------------------------------

load("BSobj_Th_5mCpG.RData")

Th1.vs.Th0.regions <- read.csv(file="DMLtest.Th1.vs.Th0.regions.csv", row.names = 1)
Th2.vs.Th0.regions <- read.csv(file="DMLtest.Th2.vs.Th0.regions.csv", row.names = 1)
Th17.vs.Th0.regions <- read.csv(file="DMLtest.Th17.vs.Th0.regions.csv", row.names = 1)
Th17.vs.Th1.regions <- read.csv(file="DMLtest.Th17.vs.Th1.regions.csv", row.names = 1)
Th17.vs.Th2.regions <- read.csv(file="DMLtest.Th17.vs.Th2.regions.csv", row.names = 1)
Th1.vs.Th2.regions <- read.csv(file="DMLtest.Th1.vs.Th2.regions.csv", row.names = 1)

all.DMRs <- rbind(Th1.vs.Th0.regions,
                  Th2.vs.Th0.regions,
                  Th17.vs.Th0.regions,
                  Th17.vs.Th1.regions,
                  Th17.vs.Th2.regions,
                  Th1.vs.Th2.regions)
all.DMRs$comparison <- c(rep("Th1.vs.Th0", nrow(Th1.vs.Th0.regions)),
                         rep("Th2.vs.Th0", nrow(Th2.vs.Th0.regions)),
                         rep("Th17.vs.Th0", nrow(Th17.vs.Th0.regions)),
                         rep("Th17.vs.Th1", nrow(Th17.vs.Th1.regions)),
                         rep("Th17.vs.Th2", nrow(Th17.vs.Th2.regions)),
                         rep("Th1.vs.Th2", nrow(Th1.vs.Th2.regions)))
head(all.DMRs) # 95
length(unique(all.DMRs$name)) # 15

write.csv(all.DMRs, file = "all.Th.DMRs.csv")

# test
plotDMRs(BSobj, regions=all.DMRs[37,], testCovariate="group", addRegions = all.DMRs, main = "test",
         extend = 5000, annoTrack = annoTrack, highlightMain = F, qval = F, stat = F, horizLegend = T)

# if customizing colors:
#pheno <- pData(BSobj)
#pheno$col <- pheno$group
#pheno$col <- gsub("nTreg", "green", pheno$col)
#pheno$col <- gsub("Treg", "blue", pheno$col)
#pheno$col <- gsub("DMK", "red", pheno$col)
#pData(BSobj) <- pheno


for(i in 1:nrow(all.DMRs)){
  jpeg(filename = paste0("DMR_", i, "_", all.DMRs[i, "comparison"], ".jpeg") , width = 960, height = 960, quality = 100)
  plotDMRs(BSobj, regions=all.DMRs[i,], testCovariate="group", addRegions = all.DMRs,
           main = paste0("pairwise comparison: ", all.DMRs[i, "comparison"]),
           extend = 5000, annoTrack = annoTrack, highlightMain = F, qval = F, stat = F, horizLegend = T)
  dev.off()
}

# error with DMRs # c(31,37.56)

all.DMRs <- all.DMRs[-c(31,37,56), ]

for(i in 1:nrow(all.DMRs)){
  jpeg(filename = paste0("DMR_", i, "_", all.DMRs[i, "comparison"], ".jpeg") , width = 960, height = 960, quality = 100)
  plotDMRs(BSobj, regions=all.DMRs[i,], testCovariate="group", addRegions = all.DMRs,
           main = paste0("pairwise comparison: ", all.DMRs[i, "comparison"]),
           extend = 5000, annoTrack = annoTrack, highlightMain = F, qval = F, stat = F, horizLegend = T)
  dev.off()
}



# heatmap -----------------------------------------------------------------


load("BSobj_Th_5mCpG.RData")

hasBeenSmoothed(BSobj)
smoothed.BSoj <- BSmooth(BSobj)

# DMRichR::smoothPheatmap not working:
# smoothPheatmap(smoothed.BSoj, GRanges(all.DMRs), testCovariate = "group", filename = NULL)

heatmap.data <- bsseq::getMeth(BSseq = smoothed.BSoj,
        regions = GRanges(all.DMRs),
        type = "smooth",
        what = "perRegion")
head(heatmap.data)
rownames(heatmap.data) <- all.DMRs$name

heatmap.data %>%   na.omit() %>% as.matrix()

pheatmap::pheatmap(heatmap.data,
                   scale = "row",
                   annotation_col =  as.data.frame(pData(smoothed.BSoj)[, c(2,5)]),
                   show_rownames = T, fontsize_row = 8,
                   show_colnames = T,
                   angle_col = 45,
                   border_color = "grey",
                   main = "Z-Scores of all pairwise DMRs",
                   fontsize = 14,
                   cellwidth = 25,
                   cellheight = 7)

save(smoothed.BSoj, file = "smoothed.BSobj.RData")


# heatmap with targeted regions

load(file = "smoothed.BSobj.RData")

regions <- read.csv("targeted.regions.csv", row.names = 1)
regions <- regions[-11, ] # remove Itgb8
head(regions)

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
                   main = "Z-Scores of 5mC in all targeted regions",
                   fontsize = 12,
                   cellwidth = 25,
                   cellheight = 20)



# end ---------------------------------------------------------------------
sessionInfo()


