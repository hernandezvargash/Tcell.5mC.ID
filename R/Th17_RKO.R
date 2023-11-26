
# main --------------------------------------------------------------------

# Th17 RKO model
# 


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



# preprocessing ----------------------------------------------------------

# keep only uncalled experiment


## load bed files----

bed.files.full <- list.files(path = "./dorado/", pattern = "Th17.TGFb", recursive = T, full.names = T)
colnames.bed <- c("chr","start","end","base","score","strand","tstart","tend","color","coverage","freq","mod","canon","other","del","fail","diff","nocall")
bed.names <- str_sub(bed.files.full, -6, -5)

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

par(mar=c(15,5,5,5), mfrow = c(1,1))
barplot(unlist(lapply(bedlist, length)))

bedlist.m <- lapply(bedlist, function(x) x[x$base == "m", ])
head(bedlist.m$WT)
bedlist.h <- lapply(bedlist, function(x) x[x$base == "h", ])
head(bedlist.h$KO)


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

par(mar=c(5,5,5,5), mfrow = c(1,3))
barplot(as.numeric(lapply(bedlist.m, function(x) sum(x$coverage))), las = 1,
        names.arg = names(bedlist.m),
        main = "Total reads per sample", col = c("orange","lightblue"))
barplot(as.numeric(lapply(bedlist.m, function(x) sum(x$mod))), las = 1,
        names.arg = names(bedlist.m),
        main = "Total 5mCpG counts per sample", col = c("orange","lightblue"))
barplot(as.numeric(lapply(bedlist.h, function(x) sum(x$mod))), las = 1,
        names.arg = names(bedlist.h),
        main = "Total 5hmCpG counts per sample",
        col = c("orange","lightblue"))


saveRDS(bedlist, file="manuscript/bedlist_uncalled.rds")
saveRDS(bedlist.m, file="manuscript/bedlist.m_uncalled.rds") # input of MIRA
saveRDS(bedlist.h, file="manuscript/bedlist.h_uncalled.rds") # input of MIRA



## make BSseq objects ----

bedlist.df <- lapply(bedlist, function(x) as.data.frame(x))
bedlist.m <- lapply(bedlist.df, function(x) x[x$base == "m", ])
head(bedlist.m[[1]])
bedlist.h <- lapply(bedlist.df, function(x) x[x$base == "h", ])
head(bedlist.h[[1]])

pre.BSseq.m <- lapply(bedlist.m, function(x) x[,c("seqnames","start","coverage","mod")])
pre.BSseq.m <- lapply(pre.BSseq.m, function(x) {colnames(x)<-c("chr","pos","N","X");x})
BSobj.m <- makeBSseqData(pre.BSseq.m, sampleNames = names(bedlist.m))

pre.BSseq.h <- lapply(bedlist.h, function(x) x[,c("seqnames","start","coverage","mod")])
pre.BSseq.h <- lapply(pre.BSseq.h, function(x) {colnames(x)<-c("chr","pos","N","X");x})
BSobj.h <- makeBSseqData(pre.BSseq.h, sampleNames = names(bedlist.h))


saveRDS(BSobj.m, file= "manuscript/BSobj.m_uncalled.rds")
saveRDS(BSobj.h, file= "manuscript/BSobj.h_uncalled.rds")




# targeted regions (n=11) ----

#rm(list=ls())

#bedlist <- readRDS(file="manuscript/bedlist.rds")

regions <- read.csv("targeted.regions.csv", row.names = 1)
head(regions)
regions <- GRanges(regions)

bedlist.t <- lapply(bedlist, function(x) subsetByOverlaps(x, regions))
lapply(bedlist.t, length)
#$KO 6286
#$WT 5278

bedlist.t <- lapply(bedlist.t, function(x) as.data.frame(x))

bedlist.m <- lapply(bedlist.t, function(x) x[x$base == "m", ])
head(bedlist.m[[1]])
bedlist.h <- lapply(bedlist.t, function(x) x[x$base == "h", ])
head(bedlist.h[[1]])

par(mar=c(5,5,5,5), mfrow = c(1,2))
barplot(unlist(lapply(bedlist.m, nrow)), las = 2, col = c("orange","lightblue"),
        main = "Informative on-target CpG sites per sample")
barplot(unlist(lapply(bedlist.m, function(x) sum(x$coverage))), 
        las = 2, col = c("orange","lightblue"),
        main = "On-target CpG site coverage per sample")

pre.BSseq.m <- lapply(bedlist.m, function(x) x[,c("seqnames","start","coverage","mod")])
pre.BSseq.m <- lapply(pre.BSseq.m, function(x) {colnames(x)<-c("chr","pos","N","X");x})
BSobj.m <- makeBSseqData(pre.BSseq.m, sampleNames = names(bedlist.m))

pre.BSseq.h <- lapply(bedlist.h, function(x) x[,c("seqnames","start","coverage","mod")])
pre.BSseq.h <- lapply(pre.BSseq.h, function(x) {colnames(x)<-c("chr","pos","N","X");x})
BSobj.h <- makeBSseqData(pre.BSseq.h, sampleNames = names(bedlist.h))


saveRDS(BSobj.m, file= "manuscript/BSobj.m_regions_uncalled.rds")
saveRDS(BSobj.h, file= "manuscript/BSobj.h_regions_uncalled.rds")



# targeted regions (n=139) ----

rm(list=ls())

bedlist <- readRDS(file="manuscript/bedlist_uncalled.rds")

regions <- read.csv("uncalled.regions.csv", row.names = 1)
head(regions)
regions <- GRanges(regions)

bedlist.t <- lapply(bedlist, function(x) subsetByOverlaps(x, regions))
lapply(bedlist.t, length)
#$KO 199746
#$WT 168710

bedlist.t <- lapply(bedlist.t, function(x) as.data.frame(x))

bedlist.m <- lapply(bedlist.t, function(x) x[x$base == "m", ])
head(bedlist.m[[1]])
bedlist.h <- lapply(bedlist.t, function(x) x[x$base == "h", ])
head(bedlist.h[[1]])

par(mar=c(5,5,5,5), mfrow = c(1,2))
barplot(unlist(lapply(bedlist.m, nrow)), las = 1, col = c("orange","lightblue"),
        main = "Informative on-target CpG sites per sample")
barplot(unlist(lapply(bedlist.m, function(x) sum(x$coverage))), 
        las = 1, col = c("orange","lightblue"),
        main = "On-target CpG site coverage per sample")

pre.BSseq.m <- lapply(bedlist.m, function(x) x[,c("seqnames","start","coverage","mod")])
pre.BSseq.m <- lapply(pre.BSseq.m, function(x) {colnames(x)<-c("chr","pos","N","X");x})
BSobj.m <- makeBSseqData(pre.BSseq.m, sampleNames = names(bedlist.m))

pre.BSseq.h <- lapply(bedlist.h, function(x) x[,c("seqnames","start","coverage","mod")])
pre.BSseq.h <- lapply(pre.BSseq.h, function(x) {colnames(x)<-c("chr","pos","N","X");x})
BSobj.h <- makeBSseqData(pre.BSseq.h, sampleNames = names(bedlist.h))


saveRDS(BSobj.m, file= "manuscript/BSobj.m_regions.139_uncalled.rds")
saveRDS(BSobj.h, file= "manuscript/BSobj.h_regions.139_uncalled.rds")



## DMRplots ----------------------------------------------------------------

rm(list=ls())

suppressPackageStartupMessages({
  library(bsseq)
  library(dmrseq)
  library(pals)
  library(scales)
  
})

annoTrack <- getAnnot("mm10")

BSobj.m <- readRDS("manuscript/BSobj.m_uncalled.rds")
BSobj.h <- readRDS("manuscript/BSobj.h_uncalled.rds")

regions <- read.csv("targeted.regions.csv", row.names = 1)
head(regions)

# test
sel.num <- 1 # selected region
BSobj = BSobj.m
#BSobj = BSobj.h
sampleNames(BSobj)
pData(BSobj) <- data.frame(group = sampleNames(BSobj))
plotDMRs(BSobj, regions=regions[sel.num,], testCovariate= "group", addRegions = NULL, 
         main = paste0(rownames(regions)[sel.num], " Targeted Region_5mC"),
         extend = 0, annoTrack = annoTrack, highlightMain = F, qval = F, stat = F, horizLegend = T)

output.dir <- "manuscript/uncalled/"

for(i in 1:nrow(regions)){
  BSobj = BSobj.m
  pData(BSobj) <- data.frame(group = sampleNames(BSobj))
  jpeg(filename = paste0(output.dir, "5mCpG_", rownames(regions)[i], ".jpeg" ), height = 800, width = 1200, res = 150)
  plotDMRs(BSobj, regions=regions[i,], testCovariate="group", addRegions = NULL, 
           main = paste0(rownames(regions)[i], " Targeted Region_5mCpG"),
           extend = 0, annoTrack = annoTrack, highlightMain = F, qval = F, stat = F, horizLegend = T)
  dev.off()
  BSobj = BSobj.h
  pData(BSobj) <- data.frame(group = sampleNames(BSobj))
  jpeg(filename = paste0(output.dir, "5mhCpG_", rownames(regions)[i], ".jpeg" ), height = 800, width = 1200, res = 150)
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
rownames(heatmap.data) <- rownames(regions)

#heatmap.data %>%   na.omit() %>% as.matrix()

pheno <- data.frame(group = sampleNames(BSobj))
pheno$group <- as.factor(pheno$group)
group.colors <- c("orange","lightblue")
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
                   #                   main = "Z-Scores of mean 5hmC in all targeted regions",
                   fontsize = 12,
                   cellwidth = 25,
                   cellheight = 20)


#write.csv(heatmap.data, file = "5mC.heatmap.data.csv")



# genome-wide differential methylation --------------------------------------------------------

## 5mCpG ----

rm(list=ls())

genes <- annotateTranscripts(TxDb.Mmusculus.UCSC.mm10.knownGene)

BSobj <- readRDS("manuscript/BSobj.m_uncalled.rds")

pheno <- data.frame(group = sampleNames(BSobj))
pheno$group

# Th17.KO vs Th17.WT
dml.test <- DMLtest(BSobj, group1 = 1, 
                    group2 = 2,  
                    equal.disp = T,
                    smoothing = T, # as there are no replicates
                    ncores = 12) # 

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

write.csv(dmls, file = "manuscript/tables/DMLs_5mC.uncalled_Th17.KO.vs.WT.csv")
write.csv(dmrs, file = "manuscript/tables/DMRs_5mC.uncalled_Th17.KO.vs.WT.csv")




## 5hmCpG ----

rm(list=ls())

genes <- annotateTranscripts(TxDb.Mmusculus.UCSC.mm10.knownGene)

BSobj <- readRDS("manuscript/BSobj.h_uncalled.rds")

pheno <- data.frame(group = sampleNames(BSobj))
pheno$group

# Th17.KO vs Th17.WT
dml.test <- DMLtest(BSobj, group1 = 1, 
                    group2 = 2,  
                    equal.disp = T,
                    smoothing = T, # as there are no replicates
                    ncores = 12) # 

dmls = callDML(dml.test, p.threshold=0.05, delta = 0.1)
table(dmls$fdr < 0.05)
head(dmls) # 8380
tail(dmls)

dmrs = callDMR(dml.test, p.threshold=0.05, delta = 0.1, minCG = 3, minlen = 50, dis.merge = 100, pct.sig = 0.5)
head(dmrs) # 463
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

write.csv(dmls, file = "manuscript/tables/DMLs_5hmC.uncalled_Th17.KO.vs.WT.csv")
write.csv(dmrs, file = "manuscript/tables/DMRs_5hmC.uncalled_Th17.KO.vs.WT.csv")



# 463 DMRs
unique(dmrs$name)
# corresponding to 441 unique gene symbols

# how many uncalled dmrs overlap with Cas9-targeted regions?
regions <- read.csv("targeted.regions.csv", row.names = 1)
head(regions)

intersect(dmrs$name, rownames(regions))
# 0




# regions (139) differential methylation --------------------------------------------------------


## 5mCpG ----

rm(list=ls())

genes <- annotateTranscripts(TxDb.Mmusculus.UCSC.mm10.knownGene)

BSobj <- readRDS("manuscript/BSobj.m_regions.139_uncalled.rds")

pheno <- data.frame(group = sampleNames(BSobj))
pheno$group

# Th17.KO vs Th17.WT
dml.test <- DMLtest(BSobj, group1 = 1, 
                    group2 = 2,  
                    equal.disp = F,
                    smoothing = T, # as there are no replicates
                    ncores = 12) # 

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
sort(dmrs$name)

write.csv(dmls, file = "manuscript/tables/DMLs_5mC.uncalled_regions.139_Th17.KO.vs.WT.csv")
write.csv(dmrs, file = "manuscript/tables/DMRs_5mC.uncalled_regions.139_Th17.KO.vs.WT.csv")



# 330 DMRs
unique(dmrs$name)
# corresponding to 136 unique gene symbols

# how many uncalled dmrs overlap with Cas9-targeted regions?
regions <- read.csv("targeted.regions.csv", row.names = 1)
head(regions)

intersect(dmrs$name, rownames(regions))
# "Zfp362"(3) "Tbx21"(8)  "Ifng"(2)   "Il10"(2)   "Rorc"(2)   "Gata3"(2)  "Itgb8"(1)  "Il4"(1)



## 5hmCpG ----

rm(list=ls())

genes <- annotateTranscripts(TxDb.Mmusculus.UCSC.mm10.knownGene)

BSobj <- readRDS("manuscript/BSobj.h_regions.139_uncalled.rds")

pheno <- data.frame(group = sampleNames(BSobj))
pheno$group

# Th17.KO vs Th17.WT
dml.test <- DMLtest(BSobj, group1 = 1, 
                    group2 = 2,  
                    equal.disp = F,
                    smoothing = T, # as there are no replicates
                    ncores = 12) # 

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
sort(dmls$name)

dmrs.gr <- dmrs
colnames(dmrs.gr)[2] <- "start"
dmrs.gr$end <- dmrs.gr$start
dmrs.gr <- GRanges(dmrs.gr)
match1 <- matchGenes(dmrs.gr, genes, type = "fiveprime", promoterDist = 2500, skipExons = FALSE, verbose = TRUE)
dmrs <- cbind(dmrs, match1)
head(dmrs)
sort(dmrs$name)
# "Pparg" "Zeb2"  "Zeb2" 

write.csv(dmls, file = "manuscript/tables/DMLs_5hmC.uncalled_regions.139_Th17.KO.vs.WT.csv")
write.csv(dmrs, file = "manuscript/tables/DMRs_5hmC.uncalled_regions.139_Th17.KO.vs.WT.csv")



# how many uncalled dmls overlap with Cas9-targeted regions?
regions <- read.csv("targeted.regions.csv", row.names = 1)
head(regions)

intersect(dmls$name, rownames(regions))
# "Itgb8"(1)






# Targeted DMRplots ----------------------------------------------------------------

## 5mCpG ----

rm(list=ls())

annoTrack <- dmrseq::getAnnot("mm10")
#save(annoTrack, file = "annoTrack_mm10.RData")

BSobj <- readRDS("manuscript/BSobj.m_uncalled.rds")
pData(BSobj) <- data.frame(group = sampleNames(BSobj))

targets <- read.csv("manuscript/tables/DMRs_5mC.uncalled_regions.139_Th17.KO.vs.WT.csv", row.names = 1)
head(targets)

# uncalled-targeted regions
regions <- read.csv("targeted.regions.csv", row.names = 1)
head(regions)

# Cas9-targeted regions
regions2 <- read.csv("targeted.regions.csv", row.names = 1)
intersect(targets$name, rownames(regions2))
# "Zfp362" "Tbx21"  "Ifng"   "Il10"   "Rorc"   "Gata3"  "Itgb8"  "Il4"
targets2 <- targets[targets$name %in% rownames(regions2), ]

sel.num <- 3 # selected region
plotDMRs(BSobj, regions=targets2[sel.num,], testCovariate="group", addRegions = targets2, 
         main = paste0(targets2$name[sel.num], " Targeted Region_5mC"),
         extend = 1000, annoTrack = annoTrack2, highlightMain = F, qval = F, stat = F, horizLegend = T)

for(i in 1:nrow(targets2)){
  jpeg(filename = paste0("manuscript/uncalled/uncalled_DMR_5mCpG_", targets2$name[i], ".nCG.", targets2$nCG[i], ".jpeg" ), height = 800, width = 1200, res = 150)
  plotDMRs(BSobj, regions=targets2[i,], testCovariate="group", addRegions = targets, 
           main = paste0(rownames(targets2)[i], ", nCG= ", targets2$nCG[i], ", Uncalled DMR 5mCpG"),
           extend = 1000, annoTrack = annoTrack, highlightMain = F, qval = F, stat = F, horizLegend = T)
  dev.off()
}


### customized DMR plots ----

rm(list=ls())

detach("package:dmrseq", unload = TRUE)

# edited dmrseq::plotDMRs functions to include chromatin status
source("plotDMRs.R")
source("helper_functions.R")

#library(annotatr)
library(bsseq)
library(RColorBrewer)
library(scales)

# load CpG and gene default annotations
# load(file = "annoTrack_mm10.RData")

# loaded chromatin annotations
# ENCODE Candidate Cis-Regulatory Elements (cCREs) combined from all cell types
#cCRE <- import.bb("/home/hh/Dropbox/HH/code/common/encode/encodeCcreCombined_mm10.bb")
#table(cCRE$ucscLabel)
#CTCF   enhD   enhP   K4m3   prom 
#24072 211185  73822  10538  24114 
#show_col(levels(factor(cCRE$itemRgb)))
"#00b0f0" "#ff0000" "#ffa700" "#ffaaaa" "#ffcd00"
#cCRE <- keepStandardChromosomes(cCRE, pruning.mode = "coarse")

# create new annoTrack with one extra layer
#annoTrack2 <- annoTrack
#annoTrack2$cCRE <- compareTrack$cCRE
#names(annoTrack2)
#save(annoTrack2, file = "annoTrack2_mm10.RData")

# optional: use fantom enhancer annotations and compareTrack parameter
#builtin_annotations()[grep("mm10", builtin_annotations())]
#fantom = build_annotations(genome = 'mm10', annotations = c("mm10_enhancers_fantom"))
#fantom <- keepStandardChromosomes(fantom, pruning.mode = "coarse")
#compareTrack = GRangesList(Fantom = fantom, cCRE = cCRE, compress = F)
#save(compareTrack, file = "compareTrack_mm10.RData")

load("annoTrack2_mm10.RData")

BSobj <- readRDS("manuscript/BSobj.m_uncalled.rds")
pData(BSobj) <- data.frame(group = sampleNames(BSobj))

targets <- read.csv("manuscript/tables/DMRs_5mC.uncalled_regions.139_Th17.KO.vs.WT.csv", row.names = 1)
head(targets)

# uncalled-targeted regions
#regions <- read.csv("targeted.regions.csv", row.names = 1)
#head(regions)

# Cas9-targeted regions
regions2 <- read.csv("targeted.regions.csv", row.names = 1)
intersect(targets$name, rownames(regions2))
# "Zfp362" "Tbx21"  "Ifng"   "Il10"   "Rorc"   "Gata3"  "Itgb8"  "Il4"
targets2 <- targets[targets$name %in% rownames(regions2), ]

table(targets2$name)
# Gata3   Ifng   Il10    Il4  Itgb8   Rorc  Tbx21 Zfp362 
#  2      2      2      1      1      2      8      3 

sel.genes <- unique(targets2$name)
i = 8
sel.targets <- targets2[targets2$name == sel.genes[i], ]
region.cols <- c(rep(alpha("orange", 0.3), nrow(sel.targets)))
region.extended <- GRanges(paste0(sel.targets$chr[1], ":", min(sel.targets$start)-1000 , "-",  max(sel.targets$end)+5000))
#region.extended <- GRanges(paste0(sel.targets$chr[1], ":", min(sel.targets$start)+6000 , "-",  max(sel.targets$end)+2000))
plotDMRs(BSobj, regions= region.extended, regionCol = region.cols, addRegions = sel.targets, testCovariate="group",
         extend = 0, annoTrack = annoTrack2, compareTrack = NULL, highlightMain = F, qval = F, stat = F, horizLegend = T, addLines = T)





## 5hmCpG ----

rm(list=ls())

detach("package:dmrseq", unload = TRUE)

# edited dmrseq::plotDMRs functions to include chromatin status
source("plotDMRs.R")
source("helper_functions.R")

#library(annotatr)
library(bsseq)
library(RColorBrewer)
library(scales)

load("annoTrack2_mm10.RData")

BSobj.m <- readRDS("manuscript/BSobj.m_uncalled.rds")
BSobj.h <- readRDS("manuscript/BSobj.h_uncalled.rds")
pData(BSobj.m) <- data.frame(group = sampleNames(BSobj.m))
pData(BSobj.h) <- data.frame(group = sampleNames(BSobj.h))

targets.m <- read.csv("manuscript/tables/DMRs_5mC.uncalled_regions.139_Th17.KO.vs.WT.csv", row.names = 1)
head(targets.m)
targets.h <- read.csv("manuscript/tables/DMRs_5hmC.uncalled_regions.139_Th17.KO.vs.WT.csv", row.names = 1)
head(targets.h)

intersect(targets.h$name, targets.m$name)
targets.m <- targets.m[targets.m$name %in% c("Zeb2","Pparg"), ]

region.cols <- alpha("orange", 0.3)

region.extended <- GRanges(targets.h[1, ])
#start(region.extended) <- start(region.extended)-20000

source("helper_functions_5hmC.R")
plotDMRs(BSobj.h, regions= region.extended, regionCol = region.cols, addRegions = targets.h, testCovariate="group",
         extend = 3000, annoTrack = annoTrack2, compareTrack = NULL, highlightMain = F, qval = F, stat = F, horizLegend = T, addLines = T)

source("helper_functions.R")
plotDMRs(BSobj.m, regions= region.extended, regionCol = region.cols, addRegions = targets.m, testCovariate="group",
         extend = 3000, annoTrack = annoTrack2, compareTrack = NULL, highlightMain = F, qval = F, stat = F, horizLegend = T, addLines = T)


sel.genes <- unique(targets2$name)
i = 8
sel.targets <- targets2[targets2$name == sel.genes[i], ]
region.cols <- c(rep(alpha("orange", 0.3), nrow(sel.targets)))
region.extended <- GRanges(paste0(sel.targets$chr[1], ":", min(sel.targets$start)-1000 , "-",  max(sel.targets$end)+5000))
#region.extended <- GRanges(paste0(sel.targets$chr[1], ":", min(sel.targets$start)+6000 , "-",  max(sel.targets$end)+2000))
plotDMRs(BSobj, regions= region.extended, regionCol = region.cols, addRegions = sel.targets, testCovariate="group",
         extend = 0, annoTrack = annoTrack2, compareTrack = NULL, highlightMain = F, qval = F, stat = F, horizLegend = T, addLines = T)




#

annoTrack <- getAnnot("mm10")

BSobj <- readRDS("manuscript/BSobj.h_uncalled.rds")
pData(BSobj) <- data.frame(group = sampleNames(BSobj))

targets <- read.csv("manuscript/tables/DMRs_5hmC.uncalled_regions.139_Th17.KO.vs.WT.csv", row.names = 1)
head(targets)

sel.num <- 1 # selected region
plotDMRs(BSobj, regions=targets[sel.num,], testCovariate="group", addRegions = targets, 
         main = paste0(targets$name[sel.num], ", nCG= ", targets$nCG[sel.num], " Uncalled DMR 5hmCpG"),
         extend = 1000, annoTrack = annoTrack, highlightMain = F, qval = F, stat = F, horizLegend = T)





# end ---------------------------------------------------------------------

sessionInfo()
