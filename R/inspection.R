

# main --------------------------------------------------------------------

# in vitro differentiation experiments:

# 20210113: second set of samples, replicates from the AK28 set
# 20210312: second set of samples sent, replicates from the AKXX set
# 20210317: second set of samples sent, replicates from the ARJP28 set

# Crispr targeted LIGATION sequencing SQK-LSK109
# EXP-NBD104 barcoding
# RB1=Th0, RB2=Th1, RB3=Th2, RB4=TH17

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
  
})

genome <- BSgenome.Mmusculus.UCSC.mm10


setwd("~/Dropbox/BioInfo/Lab/Tcell_ID")



# prepare phenotype data -----------------------------------------------------------

bed.files <- list.files(path = "./remora/", pattern = "5mC.bed", recursive = T, full.names = F)

pdata <- as.data.frame(str_split(bed.files, "/", simplify = T))
colnames(pdata) <- c("run","group","basename")
pdata$bed.file <- bed.files
pdata$assay <- rep(paste0("Assay.", 1:3), each = 4)
rownames(pdata) <- paste0(pdata$group, "_", pdata$assay)

head(pdata)

write.csv(pdata, file = "pdata.5mC.csv")


# load bed files ----------------------------------------------------------

bed.files.full <- list.files(path = "./remora/", pattern = "5mC.bed", recursive = T, full.names = T)


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
save(bed.list, file="bed.list.5mC.RData")



# inspection --------------------------------------------------------------

load("bed.list.5mC.RData")
pdata <- read.csv("pdata.5mC.csv", row.names = 1)

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
        ylim = c(67,72),
        main = "mean methylation per condition",
        ylab = "5mCpG methylation [%]")


t.test(means.Th0, means.Th1, paired = F)
# p-value = 0.6949


# methylation distribution

block.sizes.df <- data.frame(Reduce(cbind, block.sizes))
colnames(block.sizes.df) <- names(block.sizes)
block.sizes.df <- as.matrix(block.sizes.df)
block.sizes.df <- block.sizes.df/100 # to get beta values
head(block.sizes.df)

pdata <- pdata[order(pdata$group), ]
block.sizes.df <- block.sizes.df[, rownames(pdata)]

par(mar = c(10,10,8,10), mfrow = c(1,1))
densityBeanPlot(block.sizes.df, 
                sampGroups = pdata$group,
                sampNames = rownames(pdata),
                numPositions = 10000)


# check distribution for high coverage loci (>5)

test <- lapply(bed.list, function(x)  {x <- x[x$blockCount > 5, ]})

barplot(as.numeric(lapply(test, length)), las = 2,
        names.arg = names(bed.list),
        main = "Total reads per sample (coverage > 5)",
        col = viridis(4))

test2 <- lapply(test, function(x) as.numeric(x$blockSizes))
test2 <- data.frame(Reduce(cbind, test2))
colnames(test2) <- names(block.sizes)
test2 <- as.matrix(test2)
test2 <- test2/100 # to get beta values
test2 <- test2[, rownames(pdata)]
head(test2)

densityBeanPlot(test2, 
                sampGroups = pdata$group,
                main = "5mC distribution (coverage > 5)",
                numPositions = 1000)




# end ---------------------------------------------------------------------
sessionInfo()

