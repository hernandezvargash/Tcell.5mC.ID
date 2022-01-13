

######################################
# 180521
######################################

# code from Chloe:

#uncalled data 
dat1.1 <- read.csv("/home/tgfb/Downloads/regions_uncalled_mouse3.csv", sep="", header=FALSE)
colnames(dat1.1) <- c("chr", "start", "end", "Gene.names")
# uncalled performed on the laptop with 400ng DNA and 10KB flanking region
dat2.1 <-readGAlignments("/media/tgfb/DATA/Data_store/20210205/Analysed2/merged.sorted.bam")
genes <- GRanges((dat1.1), names = (dat1.1$Gene.names))
polish <- subsetByOverlaps(GRanges(dat2.1), genes)
ol <- genes[subjectHits(findOverlaps(polish, genes))]


# basecalling with megalodon
######################################

conda activate megalodon

DIR=/mnt/sdb/results/Th17/WT
FAST5PATH=/mnt/sdb/projects/Th17/uncalled/20210419/fast5
REFERENCE=/mnt/sdb/refs/GRCm38.p6.genome.fa
RERIO_PATH=/mnt/sdb/scripts/rerio/basecall_models/
  CONFIG="res_dna_r941_min_modbases_5mC_v001.cfg"

megalodon $FAST5PATH \
--guppy-params "-d $RERIO_PATH" \
--guppy-config $CONFIG \
--guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server  \
--outputs basecalls mods per_read_mods mod_basecalls mappings mod_mappings  \
--write-mod-log-probs \
--write-mods-text \
--mappings-format bam \
--mod-motif m CN 0 \
--sort-mappings \
--mod-map-emulate-bisulfite \
--mappings-format bam \
--mod-map-base-conv C T \
--mod-map-base-conv m C \
--mod-binary-threshold 0.8 \
--reference ${REFERENCE} \
--mod-output-formats bedmethyl modvcf wiggle \
--output-directory ${DIR} \
--devices "cuda:all:100%" \
--processes 8

conda deactivate


######################################


conda activate megalodon

DIR=/mnt/sdb/results/Th17/KO
FAST5PATH=/mnt/sdb/projects/Th17/uncalled/20210402/fast5
REFERENCE=/mnt/sdb/refs/GRCm38.p6.genome.fa
RERIO_PATH=/mnt/sdb/scripts/rerio/basecall_models/
  CONFIG="res_dna_r941_min_modbases_5mC_v001.cfg"

megalodon $FAST5PATH \
--guppy-params "-d $RERIO_PATH" \
--guppy-config $CONFIG \
--guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server  \
--outputs basecalls mods per_read_mods mod_basecalls mappings mod_mappings  \
--write-mod-log-probs \
--write-mods-text \
--mappings-format bam \
--mod-motif m CN 0 \
--sort-mappings \
--mod-map-emulate-bisulfite \
--mappings-format bam \
--mod-map-base-conv C T \
--mod-map-base-conv m C \
--mod-binary-threshold 0.8 \
--reference ${REFERENCE} \
--mod-output-formats bedmethyl modvcf wiggle \
--output-directory ${DIR} \
--devices "cuda:all:100%" \
--processes 8

conda deactivate


######################################

rm(list=ls())

suppressPackageStartupMessages({
  library(rtracklayer)
  library(DSS)
  library(bsseq)
  library(DMRcate)
  library(BSgenome.Mmusculus.UCSC.mm10)
  library(rtracklayer)
  library(gridExtra)
  library(reshape2)
  library(minfi)
  library(wateRmelon)
  library(NMF)
  library(mgcv)
  library(ggplot2)
  library(vioplot)
  library(corrplot)
  library(viridis)
  library(EnrichedHeatmap)
  library(circlize)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(scales)
  library(annotatr)
  library(xlsx)
  
})

genome <- BSgenome.Mmusculus.UCSC.mm10

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene



# prepare regions and BedMethyl files
######################################

setwd("/mnt/sdb/results/Th17")
list.files()

regions <- read.csv("regions.csv", row.names = 1, header = F)
head(regions)
colnames(regions) <- c("chr","start","end")
table(regions$chr)
# chr1 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19  chr2  chr3  chr4  chr5  chr6  chr7  chr8  chr9  chrX 
# 9     9    21     7     7     8     5     3     5     3     2     7     8    10     7     9     6     4     8     1 
regions <- GRanges(regions)
# 139 ranges
seqlevels(regions)
summary(width(regions))
hist(width(regions), breaks = 100)

WT.bed <- import.bed("/mnt/sdb/results/Th17/WT/modified_bases.5mC.bed", which = regions)
save(WT.bed, file="WT.bed.regions.RData")
rm(WT.bed)
KO.bed <- import.bed("/mnt/sdb/results/Th17/KO/modified_bases.5mC.bed", which = regions)
save(KO.bed, file="KO.bed.regions.RData")
rm(KO.bed)


# Inspection
######################################

setwd("/mnt/sdb/results/Th17")

load("WT.bed.regions.RData")
WT.bed
# 2502022 ranges
wt <- WT.bed
load("KO.bed.regions.RData")
KO.bed
# 3087885 ranges
ko <- KO.bed

#par(mfrow=c(1,2))
vioplot(as.numeric(wt$blockSizes), as.numeric(ko$blockSizes),
        col = c("blue","orange"), ylab = "% methylation",
        names = c("WT","KO"),
        main = "Global 5mC distribution")
t.test(as.numeric(wt$blockSizes), as.numeric(ko$blockSizes))
# p-value < 2.2e-16
# mean of x mean of y 
# 3.450839  3.317017 
wilcox.test(as.numeric(wt$blockSizes), as.numeric(ko$blockSizes))
# p-value < 2.2e-16


# distribution and global differences
######################################

# get genomic context
bed.list <- list(wt,ko)
names(bed.list) <- c("wt","ko")
bed.list <- lapply(bed.list, function(x) sort(x))
bed.list <- lapply(bed.list, function(x) resize(x, width = 3, fix = "start"))
bed.list <- lapply(bed.list, function(x) {x$context <- getSeq(genome, x); x})
bed.list <- lapply(bed.list, function(x) {x$blockSizes <- as.numeric(as.character(x$blockSizes)); x})
head(bed.list)
lapply(bed.list, length)
par(mar=c(5,5,5,5), mfrow=c(1,2))
lapply(bed.list, function(x) hist(x$blockCount))
#save(bed.list, file="bed.list.RData")

# create a data.frame to extract beta-like values
bed.df <- lapply(bed.list, as.data.frame)
head(bed.df[[1]])
for(d in 1:length(bed.df)){
  rownames(bed.df[[d]]) <- paste0("pos.",bed.df[[d]]$seqnames, "_", bed.df[[d]]$start, "_", bed.df[[d]]$strand)
  #  bed.df[[d]] <- bed.df[[d]][order(rownames(bed.df[[d]])), ]
}

bed.df <- lapply(bed.df, function(x) {x<-x[x$blockCount >= 1, ]})
lapply(bed.df, dim)

save(bed.df, file = "bed.df.RData")

# as some sites are missing in some of the samples...
common.pos <- Reduce(intersect, list(rownames(bed.df[[1]]), rownames(bed.df[[2]])))
length(common.pos) # 1979101
lapply(bed.df, dim)
bed.df.short <- bed.df
for(d in 1:length(bed.df.short)){
  bed.df.short[[d]] <- bed.df.short[[d]][common.pos, ]
  bed.df.short[[d]]$blockSizes <- as.numeric(bed.df.short[[d]]$blockSizes)
}
lapply(bed.df.short, dim)
head(bed.df.short[[1]])


# prepare a matrix of beta values (blockSizes)
betas.all = do.call(cbind, bed.df.short)
betas.all <- betas.all[, c(2, 5, seq(10, length(colnames(betas.all)), by=14), length(colnames(betas.all)))]
colnames(betas.all) <- c("start", "strand", names(bed.list), "context")
betas.all <- betas.all[order(betas.all$start), ]
head(betas.all)

table(betas.all$context)
# CAA    CAC    CAG    CAT    CCA    CCC    CCG    CCT    CGA    CGC    CGG    CGT    CTA    CTC    CTG    CTT 
# 148624 159063 212454 146057 176765 151570  31170 183485  20812  23358  29305  24949 107385 174756 211359 177989 
CpG <- c("CGA","CGC","CGG","CGT")

rows.cg <- rownames(betas.all[betas.all$context%in%CpG, ])
rows.noncg <- rownames(betas.all[!betas.all$context%in%CpG, ])
rows.L.all <- rownames(betas.all[betas.all$strand=="+", ])
rows.L.cg <- rownames(betas.all[betas.all$strand=="+" & betas.all$context%in%CpG, ])
rows.L.noncg <- rownames(betas.all[betas.all$strand=="+" & !betas.all$context%in%CpG, ])
rows.H.all <- rownames(betas.all[betas.all$strand=="-", ])
rows.H.cg <- rownames(betas.all[betas.all$strand=="-" & betas.all$context%in%CpG, ])
rows.H.noncg <- rownames(betas.all[betas.all$strand=="-" & !betas.all$context%in%CpG, ])

betas <- betas.all[,3:4]
#colnames(betas) <- pheno.filt$sample_ID
head(betas)

#save(betas, file = "betas.RData")


# distribution and global differences
#############################################################

#load("betas.RData")

par(mar=c(5,5,5,5), mfrow=c(3,3))
#par(mar=c(5,5,5,5), mfrow=c(1,3))
colors <- c("blue","orange")
vioplot(betas, col = colors,
        ylab = "5mC [BlockSizes]", las= 1, 
        cex.axis = 1.5, cex.main = 2, ylim = c(0,100),
        main = "5mC distribution\nall cytosines")
vioplot(betas[rows.cg, ], col = colors,
        ylab = "5mC [BlockSizes]", las= 1,
        cex.axis = 1.5, cex.main = 2, ylim = c(0,100),
        main = "5mC distribution\nonly CpG sites")
vioplot(betas[rows.noncg, ], col = colors,
        ylab = "5mC [BlockSizes]", las= 1,
        cex.axis = 1.5, cex.main = 2, ylim = c(0,100),
        main = "5mC distribution\nnon-CpG sites")

# by strand : H

vioplot(betas[rows.H.all, ], col = colors,
        ylab = "5mC [BlockSizes]", las= 1,
        cex.axis = 1.5, cex.main = 2, ylim = c(0,100),
        main = "5mC distribution\n- strand: all cytosines")
vioplot(betas[rows.H.cg, ], col = colors,
        ylab = "5mC [BlockSizes]", las= 1,
        cex.axis = 1.5, cex.main = 2, ylim = c(0,100),
        main = "5mC distribution\n- strand: only CpG sites")
vioplot(betas[rows.H.noncg, ], col = colors,
        ylab = "5mC [BlockSizes]", las= 1,
        cex.axis = 1.5, cex.main = 2, ylim = c(0,100),
        main = "5mC distribution\n- strand: non-CpG sites")

# by strand : L

vioplot(betas[rows.L.all, ], col = colors,
        ylab = "5mC [BlockSizes]", las= 1,
        cex.axis = 1.5, cex.main = 2, ylim = c(0,100),
        main = "5mC distribution\n+ strand: all cytosines")
vioplot(betas[rows.L.cg, ], col = colors,
        ylab = "5mC [BlockSizes]", las= 1,
        cex.axis = 1.5, cex.main = 2, ylim = c(0,100),
        main = "5mC distribution\n+ strand: only CpG sites")
vioplot(betas[rows.L.noncg, ], col = colors,
        ylab = "5mC [BlockSizes]", las= 1,
        cex.axis = 1.5, cex.main = 2, ylim = c(0,100),
        main = "5mC distribution\n+ strand: non-CpG sites")



# stats (after selecting common sites)
#############

betas.cg <- betas[rows.cg, ]
# 98424
head(betas.cg)
dim(betas.cg)
summary(betas.cg)
# mean wt: 67.26, mean ko: 62.8
save(betas.cg, file="betas.cg.RData")

cor = round(cor(betas.cg$wt, betas.cg$ko, method = "spearman"), digits = 3)
# 0.845 for all commmon cytosines
# 0.555 for only CpG sites

ggplot(betas.cg, aes(x= wt, y= ko) ) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  scale_fill_viridis() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(title = 'beta value correlation'
       ,subtitle = paste0('(spearman = ', cor, ")")) +
  theme(text = element_text(family = 'Gill Sans', color = "#444444")
        ,panel.background = element_rect(fill = "white")
        ,panel.grid.minor = element_line(color = '#4d5566')
        ,panel.grid.major = element_line(color = '#586174')
        ,plot.title = element_text(size = 22)
        ,plot.subtitle = element_text(size = 16)
        ,axis.title = element_text(size = 18, color = '#555555')
        ,axis.title.y = element_text(vjust = .5, angle = 0)
        ,axis.title.x = element_text(hjust = .5)
  )


t.test(betas.cg$wt, betas.cg$ko)
# p-value < 2.2e-16
# mean of x mean of y 
#  67.25902  62.80205 
wilcox.test(betas.cg$wt, betas.cg$ko)
# p-value < 2.2e-16



# Differential methylation at selected TGFb regions
#############################################################

head(betas)

## Preparing the methylation_frequency tables to apear the same as the example tabes from bismak
length(bed.df.short)
lapply(bed.df.short, dim)
# $wt
# [1] 1979101      14
# $ko
# [1] 1979101      14

pre.BSseq <- bed.df.short

head(pre.BSseq[[1]])
#pre.BSseq <- lapply(pre.BSseq, function(x) as.data.frame(x))
pre.BSseq <- lapply(pre.BSseq, function(x) x[,c("seqnames","start","blockCount","blockSizes")])
pre.BSseq <- lapply(pre.BSseq, function(x) {colnames(x)<-c("chr", "pos", "N", "X");x})
pre.BSseq <- lapply(pre.BSseq, function(x)  transform(x, X = X*N/100))
pre.BSseq <- lapply(pre.BSseq, function(x) x <- x[rows.cg, ])
lapply(pre.BSseq, dim)
head(pre.BSseq[[1]])

## make BSseq objects
BSobj <- makeBSseqData(pre.BSseq, sampleNames = c("wt","ko"))
sampleNames(BSobj)
# save(BSobj, file = "BSobj.RData")

dmlTest <- DMLtest(BSobj, group1="ko", group2= "wt", smoothing=TRUE, smoothing.span=500, BPPARAM = bpparam()) # default smoothing span = 500 bp
# save(dmlTest, file= "dmlTest.RData")

## call DMLs with specific thresholds and P-values
dmls <- callDML(dmlTest, delta=0.1, p.threshold=0.05)
head(dmls)
#write.csv(dmls, "DMLs.KO.vs.WT.csv")

## call DMRs
dmrs = callDMR(dmlTest, delta=0.1, p.threshold=0.05)
head(dmrs)
#write.csv(dmrs, "DMRs.KO.vs.WT.csv")


# visualizations
#############################################################

rm(list=ls())

setwd("/mnt/sdb/results/Th17")
load("BSobj.RData")
BSobj
load("betas.cg.RData")
head(betas.cg)
dmls <- read.csv("DMLs.KO.vs.WT.csv", row.names = 1)
head(dmls)
dmrs <- read.csv("DMRs.KO.vs.WT.csv", row.names = 1)
head(dmrs)

dim(dmls)
hist(dmls$diff, breaks = 100)
dmls.sel <- dmls[dmls$fdr<0.01, ]
dmls.sel <- dmls.sel[dmls.sel$diff < -0.7 | dmls.sel$diff > 0.7, ]

test <- pre.BSseq[[1]]
head(test)
test <- test[as.numeric(rownames(dmls.sel)), ]
head(test)

betas.cg.sel <- betas.cg[rownames(test), ]
head(betas.cg.sel)

par(mfrow=c(1,1), mar=c(5,5,5,5))
#vioplot(betas.cg.sel)
#plot(betas.cg.sel$wt, betas.cg.sel$ko)

aheatmap(betas.cg.sel, scale = "row", distfun = "euclidean",)

betas.cg.sel.2 <- betas.cg.sel + 0.001
head(betas.cg.sel.2)
aheatmap(log2(betas.cg.sel.2), scale = "r1")


# annotate dmrs

genes <- annotateTranscripts(TxDb.Mmusculus.UCSC.mm10.knownGene)
tab<- matchGenes(dmrs,genes)
head(tab)
dmrs2 <- cbind(dmrs,tab)
head(dmrs2)
rownames(dmrs2) <- paste0("dmr.", 1:nrow(dmrs2))

#write.xlsx(dmrs2, "DMRs.KO.vs.WT.xlsx")

unique(dmrs2$name)
# 98 names
#[1] "Zfp362"        "Gata3"         "Abi3"          "Tbx21"         "Lef1"          "Runx3"         "E2f7"          "Il22"          "Gadd45b"       "Zeb2"          "Rorc"          "Anxa2"        
#[13] "Bhlhe40"       "Espnl"         "Itga2"         "Smad3"         "Ppp1r12a"      "Lmnb2"         "Ifng"          "Itgae"         "Espn"          "Cdkn2b"        "Ldhd"          "Chn2"         
#[25] "Pmepa1"        "Gls"           "S100a9"        "Rgs10"         "Furin"         "Serpine1"      "Gzmg"          "Il12b"         "1700016G22Rik" "Cpd"           "Stat1"         "Tbkbp1"       
#[37] "Ccr7"          "Tnfrsf25"      "Alox5ap"       "Pparg"         "0610040F04Rik" "Emid1"         "Gm15270"       "Id2"           "Ccl12"         "Agap1"         "Samhd1"        "Smad6"        
#[49] "Tnfsf11"       "Gzma"          "Areg"          "Grasp"         "Ltf"           "Gzmd"          "Ccr8"          "Itga1"         "Ppp2r2c"       "Mcpt2"         "Skil"          "1700084D21Rik"
#[61] "Cdh1"          "Tnfaip3"       "Cd83"          "Gpnmb"         "Cd7"           "Tnfrsf1b"      "Timp2"         "Gzmb"          "Ccl5"          "Klhl30"        "Myc"           "Casp4"        
#[73] "Maf"           "Atf3"          "Csnk1d"        "Il4"           "Il10"          "Nr4a1"         "Klf6"          "Etfb"          "Park7"         "Ereg"          "Ablim2"        "Lasp1"        
#[85] "Id1"           "Tnc"           "Tgfbr3"        "Lag3"          "Ccl8"          "Ccl1"          "Itgb8"         "E2f2"          "Il2ra"         "Odc1"          "Nfkbia"        "Mdm1"         
#[97] "Znrf1"         "Mapk1"

table(dmrs2$name)
#0610040F04Rik 1700016G22Rik 1700084D21Rik          Abi3        Ablim2         Agap1       Alox5ap         Anxa2          Areg          Atf3       Bhlhe40         Casp4          Ccl1         Ccl12 
#1             1             1             4             2             7             2             4             2             1             5             1             1             1 
#Ccl5          Ccl8          Ccr7          Ccr8           Cd7          Cd83          Cdh1        Cdkn2b          Chn2           Cpd        Csnk1d          E2f2          E2f7         Emid1 
#1             1             4             1             3             1             2             2             6             4             1             1             1             2 
#Ereg          Espn         Espnl          Etfb         Furin       Gadd45b         Gata3           Gls       Gm15270         Gpnmb         Grasp          Gzma          Gzmb          Gzmd 
#1             2             1             1             2             3             1             1             1             1             1             1             2             1 
#Gzmg           Id1           Id2          Ifng          Il10         Il12b          Il22         Il2ra           Il4         Itga1         Itga2         Itgae         Itgb8          Klf6 
#1             2             1             1             1             2             4             1             1             4             2             3             1             1 
#Klhl30          Lag3         Lasp1          Ldhd          Lef1         Lmnb2           Ltf           Maf         Mapk1         Mcpt2          Mdm1           Myc        Nfkbia         Nr4a1 
#1             1             2             2            10             1             1             1             1             2             1             1             1             1 
#Odc1         Park7        Pmepa1         Pparg      Ppp1r12a       Ppp2r2c         Rgs10          Rorc         Runx3        S100a9        Samhd1      Serpine1          Skil         Smad3 
#1             1             2             3             3             1             2             4             8             1             1             2             2             8 
#Smad6         Stat1        Tbkbp1         Tbx21        Tgfbr3         Timp2           Tnc       Tnfaip3      Tnfrsf1b      Tnfrsf25       Tnfsf11          Zeb2        Zfp362         Znrf1 
#2             1             3             4             2             1             3             2             1             1             2             8             3             1 

# DMR plots

head(dmrs)
table(dmrs$chr)

showOneDMR(dmrs[4,], BSobj)

sampleNames(BSobj)
groups <- c(wt="blue",ko="orange")
cols <- groups[as.character(c("wt","ko"))]
options(ucscChromosomeNames=FALSE, localHub = FALSE)
DMR.plot(GRanges(dmrs), dmr = 4, CpGs = BSobj, genome = "mm10", phen.col = cols)

dmrs2[dmrs2$name=="Tbx21", ]
DMR.plot(GRanges(dmrs2), dmr = 5, CpGs = BSobj, genome = "mm10", phen.col = cols)
dmrs2[dmrs2$name=="Ifng", ]
DMR.plot(GRanges(dmrs2), dmr = 24, CpGs = BSobj, genome = "mm10", phen.col = cols)
dmrs2[dmrs2$name=="Smad3", ]
DMR.plot(GRanges(dmrs2), dmr = 95, CpGs = BSobj, genome = "mm10", phen.col = cols)
dmrs2[dmrs2$name=="Runx3", ]
DMR.plot(GRanges(dmrs2), dmr = 7, CpGs = BSobj, genome = "mm10", phen.col = cols)
dmrs2[dmrs2$name=="Lef1", ]
DMR.plot(GRanges(dmrs2), dmr = 6, CpGs = BSobj, genome = "mm10", phen.col = cols)
DMR.plot(GRanges(dmrs2), dmr = 11, CpGs = BSobj, genome = "mm10", phen.col = cols)


pData(BSobj) <- data.frame(groups = c("wt","ko"))

plotDMRs(BSobj, regions= GRanges(dmrs2[5, ]), regionCol = alpha("orange", 0.2), addRegions = GRanges(dmrs2), testCovariate="groups", 
         extend = 2000, annoTrack = annoTrack, highlightMain = F, qval = F, stat = F, horizLegend = T, addLines = T)





# Enriched Heatmaps
#############################################################

setwd("/mnt/sdb/results/Th17")
load("bed.list.RData")
bed.list


# 5mC near TSS

genes <- genes(txdb)
sel.genes.1 <- subsetByOverlaps(genes, bed.list[[1]]) # 236

mat.wt.1 = normalizeToMatrix(bed.list[[1]], sel.genes.1, extend = c(500,1500), mean_mode = "absolute", w = 50, smooth = F, target_ratio = 0.6, value_column = "blockSizes")
mat.ko.1 = normalizeToMatrix(bed.list[[2]], sel.genes.1, extend = c(500,1500), mean_mode = "absolute", w = 50, smooth = F, target_ratio = 0.6, value_column = "blockSizes")

#show_col(plasma(5))

ht_list.1 <- EnrichedHeatmap(mat.wt.1, name = "5mC WT", col = c("whitesmoke", plasma(5)[1]),
                             column_title = "5mC in WT cells\ngene bodies",
                             axis_name_rot = 90,
                             top_annotation = HeatmapAnnotation(height = unit(5,"cm"),
                                                                enrich = anno_enriched(ylim = c(0, 10), 
                                                                                       gp = gpar(col = plasma(5)[1], lty = 1, lwd = 4),
                                                                                       axis_param = list(side = "left")))) +
  EnrichedHeatmap(mat.ko.1, name = "5mC KO", col = c("whitesmoke", plasma(5)[3]),
                  column_title = "5mC in KO cells\ngene bodies",
                  axis_name_rot = 90,
                  top_annotation = HeatmapAnnotation(height = unit(5,"cm"),
                                                     enrich = anno_enriched(ylim = c(0, 10), 
                                                                            gp = gpar(col = plasma(5)[3], lty = 1, lwd = 4),
                                                                            axis_param = list(side = "left"))))

jpeg("bodies.jpeg", height = 600, width = 350, quality = 100)
draw(ht_list.1, ht_gap = unit(c(10), "mm"))
dev.off()


# 5mC near CGI

session<-browserSession()
genome(session) <- "mm10"
#trackNames(session)

cgi <- session[["CpG Islands"]] # 16023
cgi <- keepStandardChromosomes(cgi, pruning.mode="coarse") # 27949

sel.genes.2 <- subsetByOverlaps(cgi, bed.list[[1]]) # 146

mat.wt.2 = normalizeToMatrix(bed.list[[1]], sel.genes.2, extend = 500, mean_mode = "absolute", w = 50, smooth = F, target_ratio = 0.3, value_column = "blockSizes")
mat.ko.2 = normalizeToMatrix(bed.list[[2]], sel.genes.2, extend = 500, mean_mode = "absolute", w = 50, smooth = F, target_ratio = 0.3, value_column = "blockSizes")

ht_list.2 <- EnrichedHeatmap(mat.wt.2, name = "5mC WT", col = c("whitesmoke", plasma(5)[1]),
                             column_title = "5mC in WT cells\nCpG Islands",
                             axis_name_rot = 90,
                             top_annotation = HeatmapAnnotation(height = unit(5,"cm"),
                                                                enrich = anno_enriched(ylim = c(0, 10), 
                                                                                       gp = gpar(col = plasma(5)[1], lty = 1, lwd = 4),
                                                                                       axis_param = list(side = "left")))) +
  EnrichedHeatmap(mat.ko.2, name = "5mC KO", col = c("whitesmoke", plasma(5)[3]),
                               column_title = "5mC in KO cells\nCpG Islands",
                               axis_name_rot = 90,
                               top_annotation = HeatmapAnnotation(height = unit(5,"cm"),
                                                                  enrich = anno_enriched(ylim = c(0, 10), 
                                                                                         gp = gpar(col = plasma(5)[3], lty = 1, lwd = 4),
                                                                                         axis_param = list(side = "left"))))

jpeg("islands.jpeg", height = 960, width = 800, quality = 100)
draw(ht_list.2, ht_gap = unit(c(10), "mm"))
dev.off()


# enhancers

builtin_annotations()
annots = c('mm10_enhancers_fantom')
annotations = build_annotations(genome = 'mm10', annotations = annots) # 44465

sel.genes.3 <- subsetByOverlaps(annotations, bed.list[[1]]) # 527

mat.wt.3 = normalizeToMatrix(bed.list[[1]], sel.genes.3, extend = 5000, mean_mode = "absolute", w = 50, smooth = F, target_ratio = 0.3, value_column = "blockSizes")
mat.ko.3 = normalizeToMatrix(bed.list[[2]], sel.genes.3, extend = 5000, mean_mode = "absolute", w = 50, smooth = F, target_ratio = 0.3, value_column = "blockSizes")

ht_list.3 <- EnrichedHeatmap(mat.wt.3, name = "5mC WT", col = c("whitesmoke",  plasma(5)[1]),
                             column_title = "5mC in WT cells\nEnhancers",
                             axis_name_rot = 90,
                             top_annotation = HeatmapAnnotation(height = unit(5,"cm"),
                                                                enrich = anno_enriched(ylim = c(0, 10), 
                                                                                       gp = gpar(col =  plasma(5)[1], lty = 1, lwd = 4),
                                                                                       axis_param = list(side = "left")))) +
  EnrichedHeatmap(mat.ko.3, name = "5mC KO", col = c("whitesmoke",  plasma(5)[3]),
                  column_title = "5mC in KO cells\nEnhancers",
                  axis_name_rot = 90,
                  top_annotation = HeatmapAnnotation(height = unit(5,"cm"),
                                                     enrich = anno_enriched(ylim = c(0, 10), 
                                                                            gp = gpar(col =  plasma(5)[3], lty = 1, lwd = 4),
                                                                            axis_param = list(side = "left"))))

jpeg("enhancers.jpeg", height = 960, width = 800, quality = 100)
draw(ht_list.3, ht_gap = unit(c(10), "mm"))
dev.off()








######################################
sessionInfo()
######################################


























