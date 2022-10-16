

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
  library("TxDb.Mmusculus.UCSC.mm10.knownGene")
  
})

genome <- BSgenome.Mmusculus.UCSC.mm10


setwd("~/Dropbox/BioInfo/Lab/Tcell_ID")



# load bed file list ------------------------------------------------------

load("bed.list.5mC.RData")
pdata <- read.csv("pdata.5mC.csv", row.names = 1)



# prepare BSseq object ----------------------------------------------------

length(bed.list)
pre.BSseq <- lapply(bed.list, function(x) as.data.frame(x))
head(pre.BSseq[[1]])
pre.BSseq <- lapply(pre.BSseq, function(x) x[,c("seqnames","start","blockCount","blockSizes")])
pre.BSseq <- lapply(pre.BSseq, function(x) {colnames(x)<-c("chr", "pos", "N", "X");x})
pre.BSseq <- lapply(pre.BSseq, function(x)  transform(x, X = X*N/100))
#pre.BSseq <- lapply(pre.BSseq, function(x) {x <- x[x$N > 0, ]})
#pre.BSseq <- lapply(pre.BSseq, function(x) {x$N <- as.numeric(x$N);x})
#pre.BSseq <- lapply(pre.BSseq, function(x) {x$N <- as.numeric(x$N);x})
#pre.BSseq <- lapply(pre.BSseq, function(x) {x$X <- round(x$X, digits = 0);x})
lapply(pre.BSseq, dim)

head(pre.BSseq[[1]])

## make BSseq objects
BSobj <- makeBSseqData(pre.BSseq, sampleNames = rownames(pdata))
sampleNames(BSobj)
pData(BSobj) <- pdata
save(BSobj, file= "BSobj_Th_5mCpG.RData")



# fit DSS multi-factor --------------------------------------------------------

load("BSobj_Th_5mCpG.RData")

genes <- annotateTranscripts(TxDb.Mmusculus.UCSC.mm10.knownGene)

pheno <- as.data.frame(pData(BSobj))
formula=~0+pheno$group+pheno$assay
model.matrix(formula)

DMLfit = DMLfit.multiFactor(BSobj, design= pheno, formula=formula, smoothing=F)

save(DMLfit, file="DMLfit.RData")



# Th1 vs Th0 --------------------------------------------------------------

load("DMLfit.RData")

head(DMLfit$X)

Contrast = matrix(c(-1,1,0,0,0,0), ncol=1)
DMLtest.Th1.vs.Th0 = DMLtest.multiFactor(DMLfit, Contrast=Contrast)
head(DMLtest.Th1.vs.Th0)
head(DMLtest.Th1.vs.Th0[order(DMLtest.Th1.vs.Th0$fdrs, decreasing = F), ])
table(DMLtest.Th1.vs.Th0$fdrs < 0.05) # 3

# Treg.vs.nTreg.sites <- callDML(DMLtest.Treg.vs.nTreg, delta=0, p.threshold=1e-5) # delta cannot be specified in multifactor test
DMLtest.Th1.vs.Th0.sites <- DMLtest.Th1.vs.Th0[DMLtest.Th1.vs.Th0$fdrs<0.05, ]
DMLtest.Th1.vs.Th0.sites <- na.omit(DMLtest.Th1.vs.Th0.sites)
gr1 <- DMLtest.Th1.vs.Th0.sites
colnames(gr1)[2] <- "start"
gr1$end <- gr1$start
gr1 <- GRanges(gr1)
match1 <- matchGenes(gr1, genes, type = "fiveprime", promoterDist = 2500, skipExons = FALSE, verbose = TRUE)
DMLtest.Th1.vs.Th0.sites <- cbind(DMLtest.Th1.vs.Th0.sites, match1)
head(DMLtest.Th1.vs.Th0.sites)

DMLtest.Th1.vs.Th0.regions <- callDMR(DMLtest.Th1.vs.Th0, p.threshold = 0.05, minCG = 3, dis.merge = 50, minlen = 10, pct.sig = 0.5)
match2 <- matchGenes(DMLtest.Th1.vs.Th0.regions, genes, type = "fiveprime", promoterDist = 2500, skipExons = FALSE, verbose = TRUE)
DMLtest.Th1.vs.Th0.regions <- cbind(DMLtest.Th1.vs.Th0.regions, match2)
head(DMLtest.Th1.vs.Th0.regions) # 2

save(DMLtest.Th1.vs.Th0, file="DMLtest.Th1.vs.Th0.rData")
write.csv(DMLtest.Th1.vs.Th0.sites, file="DMLtest.Th1.vs.Th0.sites.csv")
write.csv(DMLtest.Th1.vs.Th0.regions, file="DMLtest.Th1.vs.Th0.regions.csv")

## look at distributions of test statistics and p-values
par(mfrow=c(2,1), mar = c(5,5,5,5))
hist(DMLtest.Th1.vs.Th0$stat, 100, main="test statistics")
hist(DMLtest.Th1.vs.Th0$pvals, 100, main="P values") # not great



# Th2 vs Th0 --------------------------------------------------------------

load("DMLfit.RData")

head(DMLfit$X)

Contrast = matrix(c(-1,0,0,1,0,0), ncol=1)
DMLtest.Th1.vs.Th0 = DMLtest.multiFactor(DMLfit, Contrast=Contrast)
head(DMLtest.Th1.vs.Th0)
head(DMLtest.Th1.vs.Th0[order(DMLtest.Th1.vs.Th0$fdrs, decreasing = F), ])
table(DMLtest.Th1.vs.Th0$fdrs < 0.05) # 103

# Treg.vs.nTreg.sites <- callDML(DMLtest.Treg.vs.nTreg, delta=0, p.threshold=1e-5) # delta cannot be specified in multifactor test
DMLtest.Th1.vs.Th0.sites <- DMLtest.Th1.vs.Th0[DMLtest.Th1.vs.Th0$fdrs<0.05, ]
DMLtest.Th1.vs.Th0.sites <- na.omit(DMLtest.Th1.vs.Th0.sites)
gr1 <- DMLtest.Th1.vs.Th0.sites
colnames(gr1)[2] <- "start"
gr1$end <- gr1$start
gr1 <- GRanges(gr1)
match1 <- matchGenes(gr1, genes, type = "fiveprime", promoterDist = 2500, skipExons = FALSE, verbose = TRUE)
DMLtest.Th1.vs.Th0.sites <- cbind(DMLtest.Th1.vs.Th0.sites, match1)
head(DMLtest.Th1.vs.Th0.sites)

DMLtest.Th1.vs.Th0.regions <- callDMR(DMLtest.Th1.vs.Th0, p.threshold = 0.05, minCG = 3, dis.merge = 50, minlen = 10, pct.sig = 0.5)
match2 <- matchGenes(DMLtest.Th1.vs.Th0.regions, genes, type = "fiveprime", promoterDist = 2500, skipExons = FALSE, verbose = TRUE)
DMLtest.Th1.vs.Th0.regions <- cbind(DMLtest.Th1.vs.Th0.regions, match2)
head(DMLtest.Th1.vs.Th0.regions) # 20

save(DMLtest.Th1.vs.Th0, file="DMLtest.Th2.vs.Th0.rData")
write.csv(DMLtest.Th1.vs.Th0.sites, file="DMLtest.Th2.vs.Th0.sites.csv")
write.csv(DMLtest.Th1.vs.Th0.regions, file="DMLtest.Th2.vs.Th0.regions.csv")

## look at distributions of test statistics and p-values
par(mfrow=c(2,1), mar = c(5,5,5,5))
hist(DMLtest.Th1.vs.Th0$stat, 100, main="test statistics")
hist(DMLtest.Th1.vs.Th0$pvals, 100, main="P values") # not great



# Th17 vs Th0 --------------------------------------------------------------

load("DMLfit.RData")

head(DMLfit$X)

Contrast = matrix(c(-1,0,1,0,0,0), ncol=1)
DMLtest.Th1.vs.Th0 = DMLtest.multiFactor(DMLfit, Contrast=Contrast)
head(DMLtest.Th1.vs.Th0)
head(DMLtest.Th1.vs.Th0[order(DMLtest.Th1.vs.Th0$fdrs, decreasing = F), ])
table(DMLtest.Th1.vs.Th0$fdrs < 0.05) # 28

# Treg.vs.nTreg.sites <- callDML(DMLtest.Treg.vs.nTreg, delta=0, p.threshold=1e-5) # delta cannot be specified in multifactor test
DMLtest.Th1.vs.Th0.sites <- DMLtest.Th1.vs.Th0[DMLtest.Th1.vs.Th0$fdrs<0.05, ]
DMLtest.Th1.vs.Th0.sites <- na.omit(DMLtest.Th1.vs.Th0.sites)
gr1 <- DMLtest.Th1.vs.Th0.sites
colnames(gr1)[2] <- "start"
gr1$end <- gr1$start
gr1 <- GRanges(gr1)
match1 <- matchGenes(gr1, genes, type = "fiveprime", promoterDist = 2500, skipExons = FALSE, verbose = TRUE)
DMLtest.Th1.vs.Th0.sites <- cbind(DMLtest.Th1.vs.Th0.sites, match1)
head(DMLtest.Th1.vs.Th0.sites)

DMLtest.Th1.vs.Th0.regions <- callDMR(DMLtest.Th1.vs.Th0, p.threshold = 0.05, minCG = 3, dis.merge = 50, minlen = 10, pct.sig = 0.5)
match2 <- matchGenes(DMLtest.Th1.vs.Th0.regions, genes, type = "fiveprime", promoterDist = 2500, skipExons = FALSE, verbose = TRUE)
DMLtest.Th1.vs.Th0.regions <- cbind(DMLtest.Th1.vs.Th0.regions, match2)
head(DMLtest.Th1.vs.Th0.regions) # 12

save(DMLtest.Th1.vs.Th0, file="DMLtest.Th17.vs.Th0.rData")
write.csv(DMLtest.Th1.vs.Th0.sites, file="DMLtest.Th17.vs.Th0.sites.csv")
write.csv(DMLtest.Th1.vs.Th0.regions, file="DMLtest.Th17.vs.Th0.regions.csv")

## look at distributions of test statistics and p-values
par(mfrow=c(2,1), mar = c(5,5,5,5))
hist(DMLtest.Th1.vs.Th0$stat, 100, main="test statistics")
hist(DMLtest.Th1.vs.Th0$pvals, 100, main="P values") # not great



# Th17 vs Th1 --------------------------------------------------------------

load("DMLfit.RData")

head(DMLfit$X)

Contrast = matrix(c(0,-1,1,0,0,0), ncol=1)
DMLtest.Th1.vs.Th0 = DMLtest.multiFactor(DMLfit, Contrast=Contrast)
head(DMLtest.Th1.vs.Th0)
head(DMLtest.Th1.vs.Th0[order(DMLtest.Th1.vs.Th0$fdrs, decreasing = F), ])
table(DMLtest.Th1.vs.Th0$fdrs < 0.05) # 42

# Treg.vs.nTreg.sites <- callDML(DMLtest.Treg.vs.nTreg, delta=0, p.threshold=1e-5) # delta cannot be specified in multifactor test
DMLtest.Th1.vs.Th0.sites <- DMLtest.Th1.vs.Th0[DMLtest.Th1.vs.Th0$fdrs<0.05, ]
DMLtest.Th1.vs.Th0.sites <- na.omit(DMLtest.Th1.vs.Th0.sites)
gr1 <- DMLtest.Th1.vs.Th0.sites
colnames(gr1)[2] <- "start"
gr1$end <- gr1$start
gr1 <- GRanges(gr1)
match1 <- matchGenes(gr1, genes, type = "fiveprime", promoterDist = 2500, skipExons = FALSE, verbose = TRUE)
DMLtest.Th1.vs.Th0.sites <- cbind(DMLtest.Th1.vs.Th0.sites, match1)
head(DMLtest.Th1.vs.Th0.sites)

DMLtest.Th1.vs.Th0.regions <- callDMR(DMLtest.Th1.vs.Th0, p.threshold = 0.05, minCG = 3, dis.merge = 50, minlen = 10, pct.sig = 0.5)
match2 <- matchGenes(DMLtest.Th1.vs.Th0.regions, genes, type = "fiveprime", promoterDist = 2500, skipExons = FALSE, verbose = TRUE)
DMLtest.Th1.vs.Th0.regions <- cbind(DMLtest.Th1.vs.Th0.regions, match2)
head(DMLtest.Th1.vs.Th0.regions) # 13

save(DMLtest.Th1.vs.Th0, file="DMLtest.Th17.vs.Th1.rData")
write.csv(DMLtest.Th1.vs.Th0.sites, file="DMLtest.Th17.vs.Th1.sites.csv")
write.csv(DMLtest.Th1.vs.Th0.regions, file="DMLtest.Th17.vs.Th1.regions.csv")

## look at distributions of test statistics and p-values
par(mfrow=c(2,1), mar = c(5,5,5,5))
hist(DMLtest.Th1.vs.Th0$stat, 100, main="test statistics")
hist(DMLtest.Th1.vs.Th0$pvals, 100, main="P values") # not great



# Th17 vs Th2 --------------------------------------------------------------

load("DMLfit.RData")

head(DMLfit$X)

Contrast = matrix(c(0,0,1,-1,0,0), ncol=1)
DMLtest.Th1.vs.Th0 = DMLtest.multiFactor(DMLfit, Contrast=Contrast)
head(DMLtest.Th1.vs.Th0)
head(DMLtest.Th1.vs.Th0[order(DMLtest.Th1.vs.Th0$fdrs, decreasing = F), ])
table(DMLtest.Th1.vs.Th0$fdrs < 0.05) # 146

# Treg.vs.nTreg.sites <- callDML(DMLtest.Treg.vs.nTreg, delta=0, p.threshold=1e-5) # delta cannot be specified in multifactor test
DMLtest.Th1.vs.Th0.sites <- DMLtest.Th1.vs.Th0[DMLtest.Th1.vs.Th0$fdrs<0.05, ]
DMLtest.Th1.vs.Th0.sites <- na.omit(DMLtest.Th1.vs.Th0.sites)
gr1 <- DMLtest.Th1.vs.Th0.sites
colnames(gr1)[2] <- "start"
gr1$end <- gr1$start
gr1 <- GRanges(gr1)
match1 <- matchGenes(gr1, genes, type = "fiveprime", promoterDist = 2500, skipExons = FALSE, verbose = TRUE)
DMLtest.Th1.vs.Th0.sites <- cbind(DMLtest.Th1.vs.Th0.sites, match1)
head(DMLtest.Th1.vs.Th0.sites)

DMLtest.Th1.vs.Th0.regions <- callDMR(DMLtest.Th1.vs.Th0, p.threshold = 0.05, minCG = 3, dis.merge = 50, minlen = 10, pct.sig = 0.5)
match2 <- matchGenes(DMLtest.Th1.vs.Th0.regions, genes, type = "fiveprime", promoterDist = 2500, skipExons = FALSE, verbose = TRUE)
DMLtest.Th1.vs.Th0.regions <- cbind(DMLtest.Th1.vs.Th0.regions, match2)
head(DMLtest.Th1.vs.Th0.regions) # 32

save(DMLtest.Th1.vs.Th0, file="DMLtest.Th17.vs.Th2.rData")
write.csv(DMLtest.Th1.vs.Th0.sites, file="DMLtest.Th17.vs.Th2.sites.csv")
write.csv(DMLtest.Th1.vs.Th0.regions, file="DMLtest.Th17.vs.Th2.regions.csv")

## look at distributions of test statistics and p-values
par(mfrow=c(2,1), mar = c(5,5,5,5))
hist(DMLtest.Th1.vs.Th0$stat, 100, main="test statistics")
hist(DMLtest.Th1.vs.Th0$pvals, 100, main="P values") # better



# Th1 vs Th2 --------------------------------------------------------------

load("DMLfit.RData")

head(DMLfit$X)

Contrast = matrix(c(0,1,0,-1,0,0), ncol=1)
DMLtest.Th1.vs.Th0 = DMLtest.multiFactor(DMLfit, Contrast=Contrast)
head(DMLtest.Th1.vs.Th0)
head(DMLtest.Th1.vs.Th0[order(DMLtest.Th1.vs.Th0$fdrs, decreasing = F), ])
table(DMLtest.Th1.vs.Th0$fdrs < 0.05) # 114

# Treg.vs.nTreg.sites <- callDML(DMLtest.Treg.vs.nTreg, delta=0, p.threshold=1e-5) # delta cannot be specified in multifactor test
DMLtest.Th1.vs.Th0.sites <- DMLtest.Th1.vs.Th0[DMLtest.Th1.vs.Th0$fdrs<0.05, ]
DMLtest.Th1.vs.Th0.sites <- na.omit(DMLtest.Th1.vs.Th0.sites)
gr1 <- DMLtest.Th1.vs.Th0.sites
colnames(gr1)[2] <- "start"
gr1$end <- gr1$start
gr1 <- GRanges(gr1)
match1 <- matchGenes(gr1, genes, type = "fiveprime", promoterDist = 2500, skipExons = FALSE, verbose = TRUE)
DMLtest.Th1.vs.Th0.sites <- cbind(DMLtest.Th1.vs.Th0.sites, match1)
head(DMLtest.Th1.vs.Th0.sites)

DMLtest.Th1.vs.Th0.regions <- callDMR(DMLtest.Th1.vs.Th0, p.threshold = 0.05, minCG = 3, dis.merge = 50, minlen = 10, pct.sig = 0.5)
match2 <- matchGenes(DMLtest.Th1.vs.Th0.regions, genes, type = "fiveprime", promoterDist = 2500, skipExons = FALSE, verbose = TRUE)
DMLtest.Th1.vs.Th0.regions <- cbind(DMLtest.Th1.vs.Th0.regions, match2)
head(DMLtest.Th1.vs.Th0.regions) # 16

save(DMLtest.Th1.vs.Th0, file="DMLtest.Th1.vs.Th2.rData")
write.csv(DMLtest.Th1.vs.Th0.sites, file="DMLtest.Th1.vs.Th2.sites.csv")
write.csv(DMLtest.Th1.vs.Th0.regions, file="DMLtest.Th1.vs.Th2.regions.csv")

## look at distributions of test statistics and p-values
par(mfrow=c(2,1), mar = c(5,5,5,5))
hist(DMLtest.Th1.vs.Th0$stat, 100, main="test statistics")
hist(DMLtest.Th1.vs.Th0$pvals, 100, main="P values") # better





# end ---------------------------------------------------------------------
sessionInfo()

