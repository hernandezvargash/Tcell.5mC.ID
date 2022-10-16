
#######################################################

# 190721
# Code for:
# Regulatory T cell differentiation is controlled by αKG-induced alterations in mitochondrial metabolism and lipid homeostasis
# https://pubmed.ncbi.nlm.nih.gov/34731632/

#######################################################

# selected runs corresponding to two experiments, CAR157 and CAR158

# CAR157 (submitted to GEO as Assay 2)
# Car157_merged is made of 20201215 (IDs 29, 32, and 40) and the combined 20201207/08/09 (IDs 28 and 31).
########

# 20201207/08/09 "nightmare run"
# Crispr targeted LIGATION sequencing SQK-LSK109
# EXP-NBD104
# RB1=ID28, RB2=ID31, RB3=ID34, RB4=ID37
# Imune identity panel + UNCALLED samples from Valeries Mice (second set of samples sent) 
# 3 runs using this flow cell and the same samples (20201207-09). Tried uncalled and it wasn't workign well, so switche to not using for 08. 

# 20201215
# Crispr targeted LIGATION sequencing SQK-LSK109
# EXP-NBD104
# RB8=ID29, Treg
# RB9=ID32, DMK
# RB10=ID40, tTreg
# Imune identity panel samples from Valeries Mice (second set of samples sent) 


# CAR158 (submitted to GEO as Assay 1)
# Imune identity panel samples from Valeries Mice (second set of samples sent) 
########

# 20201210
# Crispr targeted LIGATION sequencing SQK-LSK109
# EXP-NBD104
# RB5=ID43, Treg
# RB6=ID46, DMK
# RB7=ID49, tTreg


#######################################################

# these two runs were not used ("third replicate") because they were from a mix of male and female mice

#20200907
#Crispr targeted LIGATION sequencing SQK-LSK109
#Cas9 enrichment for immune identity panel Sample ID CAR152 DMK (pooled samples 10 and 11) Valerie

#20200910
#Crispr targeted LIGATION sequencing SQK-LSK109
#Cas9 enrichment for immune identity panel Sample ID CAR152 TR (pooled samples 7 and 8) Valerie

#######################################################

# Rerio model:

# res_dna_r941_min_modbases_5mC_5hmC_v001.cfg

#######################################################

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
})

genome <- BSgenome.Mmusculus.UCSC.mm10



#######################################################

setwd("~/Dropbox/BioInfo/Tcell_ID/")
list.files()


# prepare pData

assay <- c(rep("CAR158", 3),rep("CAR157", 3))
group <- c("Treg","DMK","nTreg","DMK","nTreg","Treg")
sample_ID <- paste0(assay, "_", group)

pheno <- data.frame(sample_ID=sample_ID,assay=assay,group =group)
rownames(pheno) <- pheno$sample_ID
head(pheno)


# get bedmethyl files

setwd("~/Dropbox/BioInfo/Lab/Tcell_ID/")
bed.files <- list.files(pattern = "5mC.bed", recursive = T)
bed.list <- list()

#region = GRanges("chr1:1-248956422")
#region =  GRanges("chrM:1-16569")

for(i in 1:length(bed.files)){
  bed.list[[i]] <- import.bed(bed.files[i], which = NULL)
  bed.list[[i]] <- keepStandardChromosomes(bed.list[[i]], pruning.mode = "coarse")
  bed.list[[i]] <- sort(bed.list[[i]])
  bed.list[[i]] <- resize(bed.list[[i]], width = 3, fix = "start")
  bed.list[[i]]$context <- getSeq(genome, bed.list[[i]])
  bed.list[[i]]$blockSizes <- as.numeric(as.character(bed.list[[i]]$blockSizes))
}

names(bed.list) <- pheno$sample_ID
save(bed.list, file="bed.list.RData")

load("bed.list.RData")

lapply(bed.list, length)
# $CAR158_Treg
# [1] 6145824

# $CAR158_DMK
# [1] 2261327

# $CAR158_nTreg
# [1] 2177426

# $CAR157_DMK
# [1] 1260564

# $CAR157_nTreg
# [1] 1552107

# $CAR157_Treg
# [1] 1926103

par(mar = c(10,10,8,10))
barplot(as.numeric(lapply(bed.list, length)), las = 2,
        names.arg = pheno$sample_ID,
        main = "Total reads per sample",
        col = viridis(6))

blocks <- bed.list
blocks <- lapply(blocks, function(x) as.numeric(x$blockSizes))

# check by context
CpG <- c("CGA","CGC","CGG","CGT")
bed.list.cg <- lapply(bed.list, function(x) x[x$context%in%CpG, ])
lapply(bed.list.cg, function(x) mean(as.numeric(x$blockSizes)))
#head(bed.list.cg)
blocks.cg <- lapply(bed.list.cg, function(x) as.numeric(x$blockSizes))

bed.list.noncg <- lapply(bed.list, function(x) x[!x$context%in%CpG, ])
lapply(bed.list.noncg, function(x) mean(as.numeric(x$blockSizes)))
#head(bed.list.noncg)
blocks.noncg <- lapply(bed.list.noncg, function(x) as.numeric(x$blockSizes))

par(mar=c(10,5,5,5), mfrow=c(1,1))
#colors <- rainbow(length(names(bed.list)))
colors <- viridis(length(names(bed.list)))
colors <- c("chartreuse","blue","red")

vioplot(blocks, 
        col = colors,
        ylab = "5mC [BlockSizes]", las= 2,
        main = "5mC distribution\nall CpG sites")


#lapply(bed.list, function(x) summary(as.numeric(x$blockSizes)))
#lapply(bed.list.cg, function(x) summary(as.numeric(x$blockSizes)))
#lapply(bed.list.noncg, function(x) summary(as.numeric(x$blockSizes)))

means.nTreg <- unlist(lapply(blocks.cg, function(x) mean(x)))[c(3,5)]
#CAR158_nTreg  CAR157_nTreg 
#75.20696     75.77404 
means.Treg <- unlist(lapply(blocks.cg, function(x) mean(x)))[c(1,6)]
#CAR158_Treg CAR157_Treg 
#75.15997    76.50365 
means.DMK <- unlist(lapply(blocks.cg, function(x) mean(x)))[c(2,4)]
#CAR158_DMK CAR157_DMK 
#74.60982   74.71242

boxplot(means.nTreg, means.Treg, means.DMK,
        col = colors,
        names = c("nTreg","Treg","DMK"),
        ylim = c(72,78),
        main = "mean methylation per condition",
        ylab = "5mCpG methylation [%]")
#boxplot(means.nTreg, means.Treg, means.DMK, ylim = c(0,100))


t.test(means.nTreg, means.Treg, paired = F)
# p-value = 0.7043
t.test(unlist(blocks.cg[c(3,5)]), unlist(blocks.cg[c(1,6)]), paired = F)
# p-value = 0.1575

t.test(means.nTreg, means.DMK, paired = F)
# p-value = 0.2004
t.test(unlist(blocks.cg[c(3,5)]), unlist(blocks.cg[c(2,4)]), paired = F)
# p-value = < 2.2e-16

t.test(means.Treg, means.DMK, paired = F)
# p-value = 0.3303
t.test(unlist(blocks.cg[c(1,6)]), unlist(blocks.cg[c(2,4)]), paired = F)
# p-value  < 2.2e-16




## Differential methylation analysis using DSS - CpG level -multi-factor model
#############################################################

# DSS (Dispersion shrinkage for sequencing data) is a variance shrinkage method
# adapted to counts data and assuming a beta-binomial distribution for methylation data
# it calculates dispersion parameters, and this is followed by a Wald statistical test

## Preparing the methylation_frequency tables to apear the same as the example tabes from bismak
length(bed.list)
pre.BSseq <- lapply(bed.list, function(x) as.data.frame(x))
head(pre.BSseq[[1]])
#pre.BSseq <- lapply(pre.BSseq, function(x) as.data.frame(x))
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
BSobj <- makeBSseqData(pre.BSseq, sampleNames = pheno$sample_ID)
sampleNames(BSobj)
pData(BSobj) <- pheno
save(BSobj, file= "BSobj_5mCpG.RData")



#####

load("BSobj_5mCpG.RData")

library(bumphunter)
library("TxDb.Mmusculus.UCSC.mm10.knownGene")
genes <- annotateTranscripts(TxDb.Mmusculus.UCSC.mm10.knownGene)

pheno <- as.data.frame(pData(BSobj))
formula=~0+pheno$group+pheno$assay
model.matrix(formula)

DMLfit = DMLfit.multiFactor(BSobj, design= pheno, formula=formula, smoothing=F)
# save(DMLfit, file="DMLfit.RData")

load("DMLfit.RData")

# Treg vs nTreg

Contrast = matrix(c(0,-1,1,0), ncol=1)
DMLtest.Treg.vs.nTreg = DMLtest.multiFactor(DMLfit, Contrast=Contrast)
head(DMLtest.Treg.vs.nTreg)
head(DMLtest.Treg.vs.nTreg[order(DMLtest.Treg.vs.nTreg$fdrs, decreasing = F), ])
table(DMLtest.Treg.vs.nTreg$fdrs < 0.05) # 142

# Treg.vs.nTreg.sites <- callDML(DMLtest.Treg.vs.nTreg, delta=0, p.threshold=1e-5) # delta cannot be specified in multifactor test
Treg.vs.nTreg.sites <- DMLtest.Treg.vs.nTreg[DMLtest.Treg.vs.nTreg$fdrs<0.05, ]
Treg.vs.nTreg.sites <- na.omit(Treg.vs.nTreg.sites)
gr1 <- Treg.vs.nTreg.sites
colnames(gr1)[2] <- "start"
gr1$end <- gr1$start
gr1 <- GRanges(gr1)
match1 <- matchGenes(gr1, genes, type = "fiveprime", promoterDist = 2500, skipExons = FALSE, verbose = TRUE)
Treg.vs.nTreg.sites <- cbind(Treg.vs.nTreg.sites, match1)
head(Treg.vs.nTreg.sites)

Treg.vs.nTreg.regions <- callDMR(DMLtest.Treg.vs.nTreg, p.threshold = 0.05, minCG = 3, dis.merge = 50, minlen = 10, pct.sig = 0.5)
match2 <- matchGenes(Treg.vs.nTreg.regions, genes, type = "fiveprime", promoterDist = 2500, skipExons = FALSE, verbose = TRUE)
Treg.vs.nTreg.regions <- cbind(Treg.vs.nTreg.regions, match2)
head(Treg.vs.nTreg.regions)

save(DMLtest.Treg.vs.nTreg, file="DMLtest.Treg.vs.nTreg.rData")
write.csv(Treg.vs.nTreg.sites, file="Treg.vs.nTreg.sites.csv")
write.xlsx(Treg.vs.nTreg.sites, file="Treg.vs.nTreg.sites.xlsx")
write.csv(Treg.vs.nTreg.regions, file="Treg.vs.nTreg.regions.csv")
write.xlsx(Treg.vs.nTreg.regions, file="Treg.vs.nTreg.regions.xlsx")


# Treg vs DMK

Contrast = matrix(c(-1,0,1,0), ncol=1)
DMLtest.Treg.vs.DMK = DMLtest.multiFactor(DMLfit, Contrast=Contrast)
head(DMLtest.Treg.vs.DMK)
head(DMLtest.Treg.vs.DMK[order(DMLtest.Treg.vs.DMK$fdrs, decreasing = F), ])
table(DMLtest.Treg.vs.DMK$fdrs < 0.05) # 8

Treg.vs.DMK.sites <- DMLtest.Treg.vs.DMK[DMLtest.Treg.vs.DMK$fdrs<0.05, ]
Treg.vs.DMK.sites <- na.omit(Treg.vs.DMK.sites)
gr1b <- Treg.vs.DMK.sites
colnames(gr1b)[2] <- "start"
gr1b$end <- gr1b$start
gr1b <- GRanges(gr1b)
match1b <- matchGenes(gr1b, genes, type = "fiveprime", promoterDist = 2500, skipExons = FALSE, verbose = TRUE)
Treg.vs.DMK.sites <- cbind(Treg.vs.DMK.sites, match1b)
head(Treg.vs.DMK.sites)

Treg.vs.DMK.regions <- callDMR(DMLtest.Treg.vs.DMK, p.threshold = 0.05, minCG = 3, dis.merge = 50, minlen = 10, pct.sig = 0.5)
match2b <- matchGenes(Treg.vs.DMK.regions, genes, type = "fiveprime", promoterDist = 2500, skipExons = FALSE, verbose = TRUE)
Treg.vs.DMK.regions <- cbind(Treg.vs.DMK.regions, match2b)
head(Treg.vs.DMK.regions)

save(DMLtest.Treg.vs.DMK, file="DMLtest.Treg.vs.DMK.rData")
write.csv(Treg.vs.DMK.sites, file="Treg.vs.DMK.sites.csv")
write.xlsx(Treg.vs.DMK.sites, file="Treg.vs.DMK.sites.xlsx")
write.csv(Treg.vs.DMK.regions, file="Treg.vs.DMK.regions.csv")
write.xlsx(Treg.vs.DMK.regions, file="Treg.vs.DMK.regions.xlsx")


## look at distributions of test statistics and p-values
par(mfrow=c(1,2))
hist(DMLtest.Treg.vs.nTreg$stat, 100, main="test statistics")
hist(DMLtest.Treg.vs.nTreg$pvals, 100, main="P values")
par(mfrow=c(1,2))
hist(DMLtest.Treg.vs.DMK$stat, 100, main="test statistics")
hist(DMLtest.Treg.vs.DMK$pvals, 100, main="P values")



## Visualizations
#############################################################

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

Treg.vs.DMK.regions <- read.csv(file="Treg.vs.DMK.regions.csv", row.names = 1)
head(Treg.vs.DMK.regions)
Treg.vs.nTreg.regions <- read.csv(file="Treg.vs.nTreg.regions.csv", row.names = 1)
head(Treg.vs.nTreg.regions)



# using dmrseq
plotDMRs(BSobj, regions=Treg.vs.nTreg.regions[3,], testCovariate="group", addRegions = Treg.vs.nTreg.regions,
         extend = 5000, annoTrack = annoTrack, highlightMain = F, qval = F, stat = F, horizLegend = T)

plotDMRs(BSobj, regions=Treg.vs.DMK.regions[1,], testCovariate="group", addRegions = Treg.vs.DMK.regions,
         extend = 5000, annoTrack = annoTrack, highlightMain = F, qval = F, stat = F, horizLegend = T)


#pheno <- pData(BSobj)
pheno$col <- pheno$group
pheno$col <- gsub("nTreg", "green", pheno$col)
pheno$col <- gsub("Treg", "blue", pheno$col)
pheno$col <- gsub("DMK", "red", pheno$col)

pData(BSobj) <- pheno


for(i in 1:nrow(Treg.vs.nTreg.regions)){
  jpeg(filename = paste0("Treg.vs.nTreg_DMR_", i, ".jpeg") , width = 960, height = 960, quality = 100)
  plotDMRs(BSobj, regions=Treg.vs.nTreg.regions[i,], testCovariate="group", addRegions = Treg.vs.nTreg.regions,
           extend = 5000, annoTrack = annoTrack, highlightMain = F, qval = F, stat = F, horizLegend = T)
  dev.off()
}

for(i in 1:nrow(Treg.vs.DMK.regions)){
  jpeg(filename = paste0("Treg.vs.DMK_DMR_", i, ".jpeg") , width = 960, height = 960, quality = 100)
  plotDMRs(BSobj, regions=Treg.vs.DMK.regions[i,], testCovariate="group", addRegions = Treg.vs.DMK.regions,
           extend = 5000, annoTrack = annoTrack, highlightMain = F, qval = F, stat = F, horizLegend = T)
  dev.off()
}


#getMeth(BSobj)
meth <- getBSseq(BSobj, "M" )
head(meth)
colMeans(meth)


# Foxp3

# full gene: 7,579,567-7,593,799
foxp3 <- data.frame(chr="chrX",start=7579567, end=7593799)

# CNS0: chrX:7,565,986-7,567,542  chrX:7,571,383-7,576,032
cns0a <- data.frame(chr="chrX",start=7565986, end=7567542)
cns0b <- data.frame(chr="chrX",start=7571383, end=7576032)

# Promoter chrX:7,579,204-7,579,623
promoter <- data.frame(chr="chrX",start=7577567, end=7580567)

# CNS1 chrX:7,581,695-7,582,035
cns1 <- data.frame(chr="chrX",start=7581695, end=7582035)

# CNS2 (mm10) chrX:7,583,960-7,584,385
cns2 <- data.frame(chr="chrX",start=7583960, end=7584385)

# CNS3 chrX:7,586,562-7,586,795
cns3 <- data.frame(chr="chrX",start=7586562, end=7586795)

foxp3.all <- rbind(cns0a, cns0b, promoter, cns1, cns2, cns3)

all.regions <- data.frame(chr="chrX",start=7565986, end=7593799)

plotDMRs(BSobj, regions=all.regions, testCovariate="group", addRegions = foxp3.all, regionCol = alpha("orange", 0.2),
         extend = 1000, annoTrack = annoTrack, highlightMain = F, qval = F, stat = F, horizLegend = T, addLines = T)

plotDMRs(BSobj, regions=cns2, testCovariate="group", addRegions = foxp3.all, regionCol = alpha("orange", 0.2),
         extend = 5000, annoTrack = annoTrack, highlightMain = F, qval = F, stat = F, horizLegend = T, addLines = T)







#######################################################
sessionInfo()
#######################################################
















