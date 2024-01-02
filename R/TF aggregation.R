
# main --------------------------------------------------------------------

# TF aggregation

# The goal is to select regions (from off-target data) based on TF binding motifs and aggregate 5mC and 5hmC values for each TF
# Cas9-targeted loci: Foxp3, Il17a, Rorc, Ifng, Il10, Maf, Il4, Gata3, Zfp362, Tbx21, Itgb8
# TFs of interest: Foxp3, Rorc, Stat3, Maf, Gata3, Zfp362, Tbx21, Bcl6, Runx3, Hif1alpha
# including known lineage-specific master TFs (T-BET, GATA3, and ROR-γt) and STATs (STAT1 and STAT4, STAT6, and STAT3)

# cistromes from: https://bmcresnotes.biomedcentral.com/articles/10.1186/s13104-018-3856-x
# https://doi.org/10.6084/m9.figshare.7087697
# cistromes also used in this paper about pioneer TFs: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8184831/


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
  library(pals)
  library(scales)
  library(ggpubr)
  library(pheatmap)
  
})

# The data set contains the cistrome genomic map (human, mm10 genome assembly) as standard genomic regions files (.bed). 
# Inside the zip archive, there are two folders (mm10_cistrome , mm10_cismotifs) corresponding to the base cistrome and the motif-annotated subset. 
# Inside each folder there are *.bed files names as <TF>.<reliability>.bed 
# where <TF> is the UniProt ID of a particular transcription factor and <reliability> is the reliability category from A (the highest) to D (limited). 
# The bed-files within mm10_cismotifs contain an additional column which lists the maximal HOCOMOCO motif score reachable for a particular cistrome segment 
# (-log10 scale, the filtering was performed using the score of 4 which corresponds to the motif P-value of 0.0001).



# loop 5mC --------------------------------------------------------------------

rm(list = ls())

bed.list <- readRDS(file="manuscript/bedlist.m.rds")

pdata <- read.csv(file = "manuscript/pdata.csv", row.names = 1)



# TF data

#TFs.of.interest <- c("GATA3","TBX21","FOXP3","STAT3","BCL6","RUNX3","STAT1")

# list all available mm10 TFs
all.TFs <- list.files("cistrome_mm10/mm10_cistrome/", full.names = T)
# select motifs with highest confidence
all.TFs <- all.TFs[grep("MOUSE.A", all.TFs)]

nrows <- sapply(all.TFs, function(f) nrow(read.csv(f)) )
#for high-confidence analyses use >10k sites and >100 coverage
#all.TFs <- all.TFs[nrows > 10000] # 87 TFs
all.TFs <- all.TFs[nrows > 500] # 163 TFs

sort(all.TFs)
all.TFs <- str_sub(all.TFs, 30, -13)

# check for available TFs of interest
#sel.TFs <- intersect(TFs.of.interest, all.TFs)

sel.TFs <- all.TFs

group.colors <- pals::alphabet(n=20)
group.colors <- group.colors[c(5,18,15,19,3,7)]
names(group.colors) <- levels(factor(pdata$group))

#t=41
#i=1
for(t in 1:length(sel.TFs)){
  
  regions <- import.bed(paste0("cistrome_mm10/mm10_cistrome/", sel.TFs[t], "_MOUSE.A.bed"))
  expanded.regions <- resize(regions, 4000, fix="center")
  
  bed.list.temp <- lapply(bed.list, function(x) subsetByOverlaps(x, expanded.regions))
  #lapply(bed.list.temp, length)
  
  pre.DT <- lapply(bed.list.temp, function(x) as.data.frame(x))
  pre.DT <- lapply(pre.DT, function(x) x[,c("seqnames","start","freq","coverage")])
  pre.DT <- lapply(pre.DT, function(x) {colnames(x)<-c("chr", "start", "methylProp", "coverage");x})
  pre.DT <- lapply(pre.DT, function(x)  transform(x, methylProp = methylProp/100))

  BSDTList <- lapply(pre.DT, data.table)
  
#  bigBin <- lapply(X=BSDTList, FUN=aggregateMethyl, GRList=expanded.regions, binNum=21, minBaseCovPerBin = 100) # used 100 coverage for high-confidence analysis
  bigBin <- lapply(X=BSDTList, FUN=aggregateMethyl, GRList=expanded.regions, binNum=21, minBaseCovPerBin = 1) # used 100 coverage for high-confidence analysis
  bigBin <- lapply(bigBin, function(x) transform(x, featureID = sel.TFs[t]) )
  
  if(all(lapply(bigBin, length)>1) == TRUE){
    
    bigBinDT <- rbindNamedList(bigBin)
    setkey(bigBinDT, sampleName)
    annotDT <- data.table(sampleName = names(BSDTList), pdata)
    setkey(annotDT, sampleName)
    
    bigBinDT2 <- merge(bigBinDT, annotDT, all.x=TRUE)
    #plotMIRAProfiles(binnedRegDT=bigBinDT2)
    
    # normalize
    bigBinDT2[, methylProp := methylProp - min(methylProp) + .05, by=.(featureID, sampleName)]
    
    # prepare plots
    
    binnedRegDT = bigBinDT2
    featID = unique(binnedRegDT[, featureID])
    sampleTypeColName="group"
    binNum <- max(binnedRegDT[, bin])
    setkey(binnedRegDT, featureID)
    binPlot <- ggplot(data = binnedRegDT[featID], 
                      mapping = aes(x = factor(bin), 
                                    y = methylProp * 100,
                                    col = group)) +
      theme_classic() + ylim(c(0, 100)) +
      geom_hline(yintercept=c(0), alpha=.2) +
      ylab("Normalized DNA Methylation (5mCpG %)") + 
      xlab("Genome Regions Surrounding Sites") +
      #  scale_x_discrete(labels=xAxisForRegionPlots(binNum)) +
      geom_jitter(alpha = .8, size = 3) + 
      scale_color_manual(values= group.colors) +
      ggtitle(paste0(featID, " [5mCpG]")) +
      theme(
        plot.title = element_text(color="black", size=24, face="bold", hjust = 0.5),
        axis.title.x = element_text(color="black", size=14, face="bold"),
        axis.title.y = element_text(color="black", size=14, face="bold")
      ) +
      theme(legend.position="bottom", legend.text = element_text(size=16))

    sampleScores <- calcMIRAScore(bigBinDT2,
                                  shoulderShift="auto",
                                  regionSetIDColName="featureID",
                                  sampleIDColName="sampleName")
    sampleScores
    sampleScores$group <- annotDT$group
    
    scoreDT <- sampleScores
    sampleTypeNum <- length(unique(scoreDT[, group]))
    setkey(scoreDT, featureID)
    scorePlot <- ggplot(data = scoreDT[featID], 
                        mapping = aes(x = group, 
                                      y = score)) + 
      theme_classic() +
      ylab("MIRA Score") + xlab("Sample Type") +
      geom_boxplot(aes(fill = group), alpha = 0.8) + 
      geom_jitter(data = scoreDT[featID], mapping = aes(x = group, y = score)) + 
      scale_fill_manual(values= group.colors) +
      ggtitle("") +
      theme(
        plot.title = element_text(color="black", size=24, face="bold", hjust = 0.5),
        axis.title.x = element_text(color="black", size=14, face="bold"),
        axis.title.y = element_text(color="black", size=14, face="bold")
      ) +
      theme(legend.position="bottom", legend.text = element_text(size=16)) +
      stat_compare_means(method = "anova", label.y = min(sampleScores$score) - 0.02, size = 7) + # Add global p-value
      stat_compare_means(aes(label = sprintf("p = %5.2f", as.numeric(..p.format..))), ref.group = "Th0", method = "t.test", size = 5)
    
    
#    tiff(filename = paste0("manuscript/cistromes/5mC/10k.sites_100.coverage_4k.regions/MIRA_5mCpG_Profiles_and_Scores_", sel.TFs[t], ".tiff"), height = 500, width = 1000)
    tiff(filename = paste0("manuscript/cistromes/5mC/500.sites_1.coverage_4k.regions/MIRA_5mCpG_Profiles_and_Scores_", sel.TFs[t], ".tiff"), height = 500, width = 1000)
    p1 <- grid.arrange(binPlot, scorePlot, ncol = 2)
    plot(p1)
    dev.off()
    
#    write.csv(sampleScores, file = paste0("manuscript/cistromes/5mC/10k.sites_100.coverage_4k.regions/MIRA_5mCpG_Scores_", sel.TFs[t], ".csv") )
    write.csv(sampleScores, file = paste0("manuscript/cistromes/5mC/500.sites_1.coverage_4k.regions/MIRA_5mCpG_Scores_", sel.TFs[t], ".csv") )
    
  }
  
  rm(regions, expanded.regions, bigBinDT2, sampleScores, p1)
    
}



# loop 5hmC --------------------------------------------------------------------


rm(list = ls())

bed.list <- readRDS(file="manuscript/bedlist.h.rds")

pdata <- read.csv(file = "manuscript/pdata.csv", row.names = 1)



# TF data

#TFs.of.interest <- c("GATA3","TBX21","FOXP3","STAT3","BCL6","RUNX3","STAT1")

# list all available mm10 TFs
all.TFs <- list.files("cistrome_mm10/mm10_cistrome/", full.names = T)
# select motifs with highest confidence
all.TFs <- all.TFs[grep("MOUSE.A", all.TFs)]

nrows <- sapply(all.TFs, function(f) nrow(read.csv(f)) )
all.TFs <- all.TFs[nrows > 500]

sort(all.TFs)
all.TFs <- str_sub(all.TFs, 30, -13)

# check for available TFs of interest
#sel.TFs <- intersect(TFs.of.interest, all.TFs)

sel.TFs <- all.TFs

group.colors <- pals::alphabet(20)
group.colors <- group.colors[c(5,18,15,19,3,7)]
names(group.colors) <- levels(factor(pdata$group))

#t=78
#i=1
for(t in 1:length(sel.TFs)){
  
  regions <- import.bed(paste0("cistrome_mm10/mm10_cistrome/", sel.TFs[t], "_MOUSE.A.bed"))
  expanded.regions <- resize(regions, 2000, fix="center")
  
  bed.list.temp <- lapply(bed.list, function(x) subsetByOverlaps(x, expanded.regions))
  #lapply(bed.list.temp, length)
  
  pre.DT <- lapply(bed.list.temp, function(x) as.data.frame(x))
  pre.DT <- lapply(pre.DT, function(x) x[,c("seqnames","start","freq","coverage")])
  pre.DT <- lapply(pre.DT, function(x) {colnames(x)<-c("chr", "start", "methylProp", "coverage");x})
  pre.DT <- lapply(pre.DT, function(x)  transform(x, methylProp = methylProp/100))
  
  BSDTList <- lapply(pre.DT, data.table)
  
  bigBin <- lapply(X=BSDTList, FUN=aggregateMethyl, GRList=expanded.regions, binNum=21, minBaseCovPerBin = 1)
  bigBin <- lapply(bigBin, function(x) transform(x, featureID = sel.TFs[t]) )
  
  if(all(lapply(bigBin, length)>1) == TRUE){
    
    bigBinDT <- rbindNamedList(bigBin)
    setkey(bigBinDT, sampleName)
    annotDT <- data.table(sampleName = names(BSDTList), pdata)
    setkey(annotDT, sampleName)
    
    bigBinDT2 <- merge(bigBinDT, annotDT, all.x=TRUE)
    
    # normalize
    bigBinDT2[, methylProp := methylProp - min(methylProp) + .05, by=.(featureID, sampleName)]
    
    # prepare plots
    
    binnedRegDT = bigBinDT2
    featID = unique(binnedRegDT[, featureID])
    sampleTypeColName="group"
    binNum <- max(binnedRegDT[, bin])
    setkey(binnedRegDT, featureID)
    binPlot <- ggplot(data = binnedRegDT[featID], 
                      mapping = aes(x = factor(bin), 
                                    y = methylProp * 100,
                                    col = group)) +
      theme_classic() + ylim(c(0, 100)) +
      geom_hline(yintercept=c(0), alpha=.2) +
      ylab("Normalized DNA Methylation (5hmCpG %)") + 
      xlab("Genome Regions Surrounding Sites") +
      #  scale_x_discrete(labels=xAxisForRegionPlots(binNum)) +
      geom_jitter(alpha = .8, size = 3) + 
      scale_color_manual(values= group.colors) +
      ggtitle(paste0(featID, " [5hmCpG]")) +
      theme(
        plot.title = element_text(color="black", size=18, face="bold", hjust = 0.5),
        axis.title.x = element_text(color="black", size=14, face="bold"),
        axis.title.y = element_text(color="black", size=14, face="bold")
      )
    
    sampleScores <- calcMIRAScore(bigBinDT2,
                                  shoulderShift="auto",
                                  regionSetIDColName="featureID",
                                  sampleIDColName="sampleName")
    sampleScores
    sampleScores$group <- annotDT$group
    
    scoreDT <- sampleScores
    sampleTypeNum <- length(unique(scoreDT[, group]))
    setkey(scoreDT, featureID)
    scorePlot <- ggplot(data = scoreDT[featID], 
                        mapping = aes(x = group, 
                                      y = score)) + 
      theme_classic() +
      ylab("MIRA Score") + xlab("Sample Type") +
      geom_boxplot(aes(fill = group), alpha = 0.8) + 
      geom_jitter(data = scoreDT[featID], mapping = aes(x = group, y = score)) + 
      scale_fill_manual(values= group.colors) +
      ggtitle(paste0(featID, " [5hmCpG]")) +
      theme(
        plot.title = element_text(color="black", size=18, face="bold", hjust = 0.5),
        axis.title.x = element_text(color="black", size=14, face="bold"),
        axis.title.y = element_text(color="black", size=14, face="bold")
      ) +
      stat_compare_means(method = "anova", label.y = 0.8)+ # Add global p-value
      stat_compare_means(aes(label = sprintf("p = %5.2f", as.numeric(..p.format..))), ref.group = "Th0", method = "wilcox")
    
    
    jpeg(filename = paste0("MIRA_5hmCpG_Profiles_and_Scores_", sel.TFs[t], ".jpeg"), height = 600, width = 1200, quality = 100)
    p1 <- grid.arrange(binPlot, scorePlot, ncol = 2)
    plot(p1)
    dev.off()
    
    write.csv(sampleScores, file = paste0("MIRA_5hmCpG_Scores_", sel.TFs[t], ".csv") )
    
  }
  
  rm(regions, expanded.regions, bigBinDT2, sampleScores, p1)
  
}




# stats -------------------------------------------------------------------

rm(list=ls())

#score.files <- list.files(path = "manuscript/cistromes/5mC/500.sites_1.coverage_4k.regions/", pattern = ".csv", full.names = T)
score.files <- list.files(path = "manuscript/cistromes/5mC/10k.sites_100.coverage_4k.regions/", pattern = ".csv", full.names = T)

TF.names <- str_sub(str_split_i(score.files, "[_]", 6), 1, -5)

test <- read.csv(score.files[1], row.names = 1)
head(test)

scores.table <- data.frame(matrix(ncol = length(score.files), nrow = length(test$group)))
colnames(scores.table) <- TF.names
rownames(scores.table) <- test$sampleName
head(scores.table)

for(i in 1:length(score.files)){
  
  temp1 <- read.csv(score.files[i], row.names = 1)
  
  if (all(temp1$sampleName == rownames(scores.table))) {
    scores.table[ , i] <- temp1$score
  } else {
    print("samples do not match")
  }
  
  rm(temp1)
  
}

head(scores.table)
dim(scores.table)
scores.table$group <- as.factor(test$group)


# add stats

#test <- rstatix::kruskal_test(PAX5 ~ group, data = scores.table)
#test$p
#overall.p <- sapply(scores.table[, 1:(ncol(scores.table)-1)], function(x) rstatix::kruskal_test(x ~ group, data = scores.table)$p)

overall.p <- sapply(scores.table[, 1:(ncol(scores.table)-1)], function(x) summary(aov(x ~ scores.table$group))[[1]][["Pr(>F)"]][1])
colnames(scores.table)[overall.p < 0.05]

fdr <- p.adjust(overall.p, method = "fdr", n = length(overall.p))
colnames(scores.table)[fdr < 0.05]
# "IKZF1" "PAX5" 

scores <- scores.table[, overall.p < 0.05]
scores$sampleName <- rownames(scores)
scores$sampleType <- scores.table$group



#write.csv(scores, "manuscript/cistromes/MIRA_Scores_5mC_500.sites_1.coverage_4k.regions.csv")
write.csv(scores, "manuscript/cistromes/MIRA_Scores_5mC_10k.sites_100.coverage_4k.regions.csv")



# same results using a linear model

overall_p <- function(my_model) {
  f <- summary(my_model)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

#extract overall p-value of model
overall.p <- sapply(scores.table[, 1:(ncol(scores.table)-1)], function(x) overall_p(lm(x ~ group, data = scores.table)))
colnames(scores.table)[overall.p < 0.05]

fdr <- p.adjust(overall.p, method = "fdr", n = length(overall.p))
colnames(scores.table)[fdr < 0.05]


#




# heatmaps --------------------------------------------------------------------

rm(list = ls())

#scores <- read.csv("manuscript/cistromes/MIRA_Scores_5mC_500.sites_1.coverage_4k.regions.csv", row.names = 1)
scores <- read.csv("manuscript/cistromes/MIRA_Scores_5mC_10k.sites_100.coverage_4k.regions.csv", row.names = 1)

head(scores)

#heatmap.data <- t(scores[, 1:27])
heatmap.data <- t(scores[, 1:17])

annotation_col = scores[ , "sampleType" , drop = F]
colnames(annotation_col) <- "group"
annotation_col$group <- as.factor(annotation_col$group)

group.colors <- alphabet(20)
group.colors <- group.colors[c(5,7,19,3,18,15)]

names(group.colors) <- levels(annotation_col$group)
group.colors =  list(group = group.colors)

pheatmap::pheatmap(heatmap.data,
                   scale = "row",
                   cluster_cols = T,
                   annotation_col = annotation_col,
                   annotation_colors = group.colors,
                   show_rownames = T, fontsize_row = 12,
                   show_colnames = T,
                   annotation_legend = T,
                   angle_col = 45,
                   border_color = "grey",
                   main = "MIRA Scores 5mC",
                   fontsize = 12,
                   cellwidth = 25,
                   cellheight = 20)


# 5hmC TF target mean methylation ----

# as MIRA is not well adapted to induced-5hmC activity
# this analyses simply takes TF binding sites and produces a global average
# without infering any score

rm(list = ls())

bed.list <- readRDS(file="manuscript/bedlist.h.rds")

pdata <- read.csv(file = "manuscript/pdata.csv", row.names = 1)


# list all available mm10 TFs
all.TFs <- list.files("cistrome_mm10/mm10_cistrome/", full.names = T) # 877 TFs

# select motifs with highest confidence
all.TFs <- all.TFs[grep("MOUSE.A", all.TFs)] # 181 TFs

nrows <- sapply(all.TFs, function(f) nrow(read.csv(f)) )
all.TFs <- all.TFs[nrows > 500] # 163 TFs

sort(all.TFs)
all.TFs <- str_sub(all.TFs, 30, -13) 

sel.TFs <- all.TFs

group.colors <- pals::alphabet(20)
group.colors <- group.colors[c(5,18,15,19,3,7)]
names(group.colors) <- levels(factor(pdata$group))

#t=1
for(t in 1:length(sel.TFs)){
  
  regions <- import.bed(paste0("cistrome_mm10/mm10_cistrome/", sel.TFs[t], "_MOUSE.A.bed"))
  #summary(width(regions))
  #expanded.regions <- resize(regions, 2000, fix="center") # not extending regions for 5hmC
  expanded.regions <- regions
  
  bed.list.temp <- lapply(bed.list, function(x) subsetByOverlaps(x, expanded.regions))
  #lapply(bed.list.temp, length)
  tab1 <- as.data.frame(unlist(lapply(bed.list.temp, function(x) mean(x$freq))))
  colnames(tab1) <- "mean.meth"
  tab1$group <- pdata$group
  
  methPlot <- ggplot(data = tab1, 
                      mapping = aes(x = group, 
                                    y = mean.meth)) + 
    theme_classic() +
    ylab("Average methylation %") + xlab("") +
    geom_boxplot(aes(fill = group), alpha = 0.8) + 
    scale_fill_manual(values= group.colors) +
    ggtitle(paste0(sel.TFs[t], " [5hmCpG]")) +
    theme(
      plot.title = element_text(color="black", size=24, face="bold", hjust = 0.5),
      axis.title.x = element_text(color="black", size=14, face="bold"),
      axis.title.y = element_text(color="black", size=14, face="bold")
    ) +
    theme(legend.position="bottom", legend.text = element_text(size=16)) +
    stat_compare_means(method = "anova", label.y = 0.08, size = 7) + # Add global p-value
    stat_compare_means(aes(label = sprintf("p = %5.2f", as.numeric(..p.format..))), ref.group = "Th0", method = "t.test", size = 5)
  
#  tiff(filename = paste0("manuscript/cistromes/5hmC/MIRA_5hmCpG_Mean.Methylation_", sel.TFs[t], ".tiff"), height = 500, width = 500)
  jpeg(filename = paste0("manuscript/cistromes/5hmC/MIRA_5hmCpG_Mean.Methylation_", sel.TFs[t], ".jpeg"), height = 500, width = 500, quality = 100)
  plot(methPlot)
  dev.off()

  rm(regions, expanded.regions, tab1, methPlot)
  
}




# 5mC TF target mean methylation ----

# this analyses simply takes TF binding sites and produces a global average
# without inferring any score

rm(list = ls())

bed.list <- readRDS(file="manuscript/bedlist.m.rds")

pdata <- read.csv(file = "manuscript/pdata.csv", row.names = 1)


# list all available mm10 TFs
all.TFs <- list.files("cistrome_mm10/mm10_cistrome/", full.names = T) # 877 TFs

# select motifs with highest confidence
all.TFs <- all.TFs[grep("MOUSE.A", all.TFs)] # 181 TFs

nrows <- sapply(all.TFs, function(f) nrow(read.csv(f)) )
#all.TFs <- all.TFs[nrows > 500] # 163 TFs
all.TFs <- all.TFs[nrows > 10000]

sort(all.TFs)
all.TFs <- str_sub(all.TFs, 30, -13) 

sel.TFs <- all.TFs

group.colors <- pals::alphabet(20)
group.colors <- group.colors[c(5,18,15,19,3,7)]
names(group.colors) <- levels(factor(pdata$group))

#t=1
for(t in 1:length(sel.TFs)){
  
  regions <- import.bed(paste0("cistrome_mm10/mm10_cistrome/", sel.TFs[t], "_MOUSE.A.bed"))
  #summary(width(regions))
  expanded.regions <- resize(regions, 2000, fix="center") # not extending regions for 5hmC
  #expanded.regions <- regions
  
  bed.list.temp <- lapply(bed.list, function(x) subsetByOverlaps(x, expanded.regions))
  #lapply(bed.list.temp, length)
  tab1 <- as.data.frame(unlist(lapply(bed.list.temp, function(x) mean(x$freq))))
  colnames(tab1) <- "mean.meth"
  tab1$group <- pdata$group
  
  methPlot <- ggplot(data = tab1, 
                     mapping = aes(x = group, 
                                   y = mean.meth)) + 
    theme_classic() +
    ylab("Average methylation %") + xlab("") +
    geom_boxplot(aes(fill = group), alpha = 0.8) + 
    scale_fill_manual(values= group.colors) +
    ggtitle(paste0(sel.TFs[t], " [5mCpG]")) +
    theme(
      plot.title = element_text(color="black", size=24, face="bold", hjust = 0.5),
      axis.title.x = element_text(color="black", size=14, face="bold"),
      axis.title.y = element_text(color="black", size=14, face="bold")
    ) +
    theme(legend.position="bottom", legend.text = element_text(size=16)) +
    stat_compare_means(method = "anova", label.y = 0.08, size = 7) + # Add global p-value
    stat_compare_means(aes(label = sprintf("p = %5.2f", as.numeric(..p.format..))), ref.group = "Th0", method = "t.test", size = 5)
  
  #  tiff(filename = paste0("manuscript/cistromes/5hmC/MIRA_5hmCpG_Mean.Methylation_", sel.TFs[t], ".tiff"), height = 500, width = 500)
  jpeg(filename = paste0("manuscript/cistromes/5mC/mean_methylation/MIRA_5mCpG_Mean.Methylation_", sel.TFs[t], ".jpeg"), height = 500, width = 500, quality = 100)
  plot(methPlot)
  dev.off()
  
  rm(regions, expanded.regions, tab1, methPlot)
  
}




# end ---------------------------------------------------------------------
sessionInfo()
