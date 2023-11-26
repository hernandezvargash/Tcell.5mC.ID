

# based on targeted 5mCpG ------------------------------------------------------------


library(fmsb)

rm(list=ls())

setwd("~/Dropbox/BioInfo/Lab/Tcell_ID")


data <- read.csv(file = "manuscript/heatmaps/5mC.targeted.heatmap.data.csv")
#data <- read.csv(file = "manuscript/heatmaps/5mC.heatmap.data.csv")

# select one DMR per targeted gene name
sel.DMRs <- c(1,5,7,10,11,17,22)

data <- data.frame(t(data[sel.DMRs, ]))
colnames(data) <- data[1, ]
data <- data[-1, ]
rownames(data) 
colnames(data)[1] <- "Foxp3"
head(data)

max_min <- as.data.frame(matrix(nrow = 2, ncol = ncol(data)))
rownames(max_min) <- c("Max", "Min")
colnames(max_min) <- colnames(data)
max_min[1, ] <- 1
max_min[2, ] <- 0

library(pals)
group.colors <- pals::alphabet(20)
group.colors <- group.colors[c(5,18,15,19,3,7)]
names(group.colors) <- c("Th0","Th1","Th17","Th17_noTgfb","Th2","Treg") 
#group.colors =  list(group = group.colors)

par(mfrow=c(3,3), mar=c(2,2,2,2))

plot.data <- rbind(max_min, as.numeric(data[3, ]))
radarchart(plot.data, pfcol = group.colors[1], axistype = 3, vlcex = 2, title = "Th0", cex.main = 2)
plot.data <- rbind(max_min, as.numeric(data[7, ]))
radarchart(plot.data, pfcol = group.colors[1], axistype = 3, vlcex = 2, title = "Th0", cex.main = 2)
plot.data <- rbind(max_min, as.numeric(data[11, ]))
radarchart(plot.data, pfcol = group.colors[1], axistype = 3, vlcex = 2, title = "Th0", cex.main = 2)

plot.data <- rbind(max_min, as.numeric(data[5, ]))
radarchart(plot.data, pfcol = group.colors[5], axistype = 3, vlcex = 2, title = "Th2", cex.main = 2)
plot.data <- rbind(max_min, as.numeric(data[9, ]))
radarchart(plot.data, pfcol = group.colors[5], axistype = 3, vlcex = 2, title = "Th2", cex.main = 2)
plot.data <- rbind(max_min, as.numeric(data[13, ]))
radarchart(plot.data, pfcol = group.colors[5], axistype = 3, vlcex = 2, title = "Th2", cex.main = 2)

plot.data <- rbind(max_min, as.numeric(data[1, ]))
radarchart(plot.data, pfcol = group.colors[6], axistype = 3, vlcex = 2, title = "Treg", cex.main = 2)
plot.data <- rbind(max_min, as.numeric(data[2, ]))
radarchart(plot.data, pfcol = group.colors[6], axistype = 3, vlcex = 2, title = "Treg", cex.main = 2)

par(mfrow=c(3,3))

plot.data <- rbind(max_min, as.numeric(data[4, ]))
radarchart(plot.data, pfcol = group.colors[2], axistype = 3, vlcex = 2, title = "Th1", cex.main = 2)
plot.data <- rbind(max_min, as.numeric(data[8, ]))
radarchart(plot.data, pfcol = group.colors[2], axistype = 3, vlcex = 2, title = "Th1", cex.main = 2)
plot.data <- rbind(max_min, as.numeric(data[12, ]))
radarchart(plot.data, pfcol = group.colors[2], axistype = 3, vlcex = 2, title = "Th1", cex.main = 2)

plot.data <- rbind(max_min, as.numeric(data[6, ]))
radarchart(plot.data, pfcol = group.colors[3], axistype = 3, vlcex = 2, title = "Th17", cex.main = 2)
plot.data <- rbind(max_min, as.numeric(data[10, ]))
radarchart(plot.data, pfcol = group.colors[3], axistype = 3, vlcex = 2, title = "Th17", cex.main = 2)
plot.data <- rbind(max_min, as.numeric(data[14, ]))
radarchart(plot.data, pfcol = group.colors[3], axistype = 3, vlcex = 2, title = "Th17", cex.main = 2)

plot.data <- rbind(max_min, as.numeric(data[15, ]))
radarchart(plot.data, pfcol = group.colors[4], axistype = 3, vlcex = 2, title = "Th17/Th1", cex.main = 2)
plot.data <- rbind(max_min, as.numeric(data[16, ]))
radarchart(plot.data, pfcol = group.colors[4], axistype = 3, vlcex = 2, title = "Th17/Th1", cex.main = 2)


## grouped spidernet plots ----

par(mfrow=c(1,2), mar=c(2,2,5,2))

plot.data <- rbind(max_min, as.numeric(data[3, ]))
plot.data <- rbind(plot.data, as.numeric(data[13, ]))
plot.data <- rbind(plot.data, as.numeric(data[2, ]))
radarchart(plot.data, pcol = group.colors[c(1,5,6)], axistype = 3, 
           vlcex = 2, plty = 3, plwd = 5, cglcol = "navy",
           title = "Th0 vs Th2 vs Treg", cex.main = 2)

plot.data <- rbind(max_min, as.numeric(data[4, ]))
plot.data <- rbind(plot.data, as.numeric(data[10, ]))
plot.data <- rbind(plot.data, as.numeric(data[15, ]))
radarchart(plot.data, pcol = group.colors[c(2,3,4)], axistype = 3, 
           vlcex = 2, plty = 3, plwd = 5, cglcol = "navy",
           title = "Th1 vs Th17 vs Th17/Th1", cex.main = 2)



par(mfrow=c(1,1), mar=c(2,2,5,2))

plot.data <- rbind(max_min, as.numeric(data[3, ]))
plot.data <- rbind(plot.data, as.numeric(data[13, ]))
plot.data <- rbind(plot.data, as.numeric(data[2, ]))
plot.data <- rbind(plot.data, as.numeric(data[4, ]))
plot.data <- rbind(plot.data, as.numeric(data[10, ]))
plot.data <- rbind(plot.data, as.numeric(data[15, ]))
radarchart(plot.data, pcol = group.colors[c(1,5,6,2,3,4)], axistype = 3, 
           vlcex = 2, plty = 3, plwd = 5, cglcol = "navy",
           title = "Target activity scores", cex.main = 2)

legend(1.2, 1.2, legend = names(group.colors), 
       fill = group.colors, cex = 2,
       border = group.colors, bty = "n", horiz = F)




# based on selected TF targets 5mCpG ------------------------------------------------------------


library(fmsb)

rm(list=ls())

setwd("~/Dropbox/BioInfo/Lab/Tcell_ID")

#sel.TFs <- c("FOXP3","ETS1","RUNX3","SIX2","JUND","TBX21","SP7","TET1")
sel.TFs <- c("ETS1","IRF8","JUND","RUNX3","SIX2","SP7","TET1","ZFP57", # high confidence TF
             "FOXP3") # "BMAL1","CTCFL","LYL1","NKX21","REST","STAT4","UTF1") # plus significant not high-confidence TFs

#i = 1
for(i in 1:length(sel.TFs)){
  if(i == 1){
    
    data.temp <- read.csv(file = paste0("manuscript/cistromes/5mC/MIRA_5mCpG_Scores_", sel.TFs[i], ".csv"), row.names = 1)
    colnames(data.temp)[3] <- sel.TFs[i]
    
  } else {
    
    data.temp2 <- read.csv(file = paste0("manuscript/cistromes/5mC/MIRA_5mCpG_Scores_", sel.TFs[i], ".csv"), row.names = 1)
    data.temp[, ncol(data.temp)+1] <- data.temp2$score
    colnames(data.temp)[ncol(data.temp)] <- sel.TFs[i]
    
    rm(data.temp2)
  }
  
}

rownames(data.temp) <- data.temp$sampleName

data <- data.temp[, -c(1,2,4)]

max_min <- as.data.frame(matrix(nrow = 2, ncol = ncol(data)))
rownames(max_min) <- c("Max", "Min")
colnames(max_min) <- colnames(data)
max_min[1, ] <- apply(data, 2, max)
max_min[2, ] <- apply(data, 2, min)

library(pals)
group.colors <- pals::alphabet(20)
group.colors <- group.colors[c(5,18,15,19,3,7)]
names(group.colors) <- c("Th0","Th1","Th17","Th17_noTgfb","Th2","Treg") 
#group.colors =  list(group = group.colors)

par(mfrow=c(3,3), mar=c(2,2,2,2))

plot.data <- rbind(max_min, as.numeric(data[1, ]))
radarchart(plot.data, pfcol = group.colors[1], axistype = 0, vlcex = 2, title = "Th0", cex.main = 3)
plot.data <- rbind(max_min, as.numeric(data[2, ]))
radarchart(plot.data, pfcol = group.colors[1], axistype = 0, vlcex = 2, title = "Th0", cex.main = 3)
plot.data <- rbind(max_min, as.numeric(data[3, ]))
radarchart(plot.data, pfcol = group.colors[1], axistype = 0, vlcex = 2, title = "Th0", cex.main = 3)

plot.data <- rbind(max_min, as.numeric(data[12, ]))
radarchart(plot.data, pfcol = group.colors[5], axistype = 0, vlcex = 2, title = "Th2", cex.main = 3)
plot.data <- rbind(max_min, as.numeric(data[13, ]))
radarchart(plot.data, pfcol = group.colors[5], axistype = 0, vlcex = 2, title = "Th2", cex.main = 3)
plot.data <- rbind(max_min, as.numeric(data[14, ]))
radarchart(plot.data, pfcol = group.colors[5], axistype = 0, vlcex = 2, title = "Th2", cex.main = 3)

plot.data <- rbind(max_min, as.numeric(data[15, ]))
radarchart(plot.data, pfcol = group.colors[6], axistype = 0, vlcex = 2, title = "Treg", cex.main = 3)
plot.data <- rbind(max_min, as.numeric(data[16, ]))
radarchart(plot.data, pfcol = group.colors[6], axistype = 0, vlcex = 2, title = "Treg", cex.main = 3)

par(mfrow=c(3,3))

plot.data <- rbind(max_min, as.numeric(data[9, ]))
radarchart(plot.data, pfcol = group.colors[2], axistype = 0, vlcex = 2, title = "Th1", cex.main = 3)
plot.data <- rbind(max_min, as.numeric(data[10, ]))
radarchart(plot.data, pfcol = group.colors[2], axistype = 0, vlcex = 2, title = "Th1", cex.main = 3)
plot.data <- rbind(max_min, as.numeric(data[11, ]))
radarchart(plot.data, pfcol = group.colors[2], axistype = 0, vlcex = 2, title = "Th1", cex.main = 3)

plot.data <- rbind(max_min, as.numeric(data[6, ]))
radarchart(plot.data, pfcol = group.colors[3], axistype = 0, vlcex = 2, title = "Th17", cex.main = 3)
plot.data <- rbind(max_min, as.numeric(data[7, ]))
radarchart(plot.data, pfcol = group.colors[3], axistype = 0, vlcex = 2, title = "Th17", cex.main = 3)
plot.data <- rbind(max_min, as.numeric(data[8, ]))
radarchart(plot.data, pfcol = group.colors[3], axistype = 0, vlcex = 2, title = "Th17", cex.main = 3)

plot.data <- rbind(max_min, as.numeric(data[4, ]))
radarchart(plot.data, pfcol = group.colors[4], axistype = 0, vlcex = 2, title = "Th17/Th1", cex.main = 3)
plot.data <- rbind(max_min, as.numeric(data[5, ]))
radarchart(plot.data, pfcol = group.colors[4], axistype = 0, vlcex = 2, title = "Th17/Th1", cex.main = 3)



## grouped spidernet plots ----

par(mfrow=c(1,2), mar=c(2,2,5,2))

plot.data <- rbind(max_min, as.numeric(data[2, ]))
plot.data <- rbind(plot.data, as.numeric(data[14, ]))
plot.data <- rbind(plot.data, as.numeric(data[15, ]))
radarchart(plot.data, pcol = group.colors[c(1,5,6)], axistype = 0, 
           vlcex = 2, plty = 3, plwd = 5, cglcol = "navy",
           title = "Th0 vs Th2 vs Treg", cex.main = 2)

plot.data <- rbind(max_min, as.numeric(data[10, ]))
plot.data <- rbind(plot.data, as.numeric(data[7, ]))
plot.data <- rbind(plot.data, as.numeric(data[4, ]))
radarchart(plot.data, pcol = group.colors[c(2,3,4)], axistype = 0, 
           vlcex = 2, plty = 3, plwd = 5, cglcol = "navy",
           title = "Th1 vs Th17 vs Th17/Th1", cex.main = 2)


par(mfrow=c(1,1), mar=c(2,2,5,2))

plot.data <- rbind(max_min, as.numeric(data[2, ]))
plot.data <- rbind(plot.data, as.numeric(data[14, ]))
plot.data <- rbind(plot.data, as.numeric(data[15, ]))
plot.data <- rbind(plot.data, as.numeric(data[10, ]))
plot.data <- rbind(plot.data, as.numeric(data[7, ]))
plot.data <- rbind(plot.data, as.numeric(data[4, ]))
radarchart(plot.data, pcol = group.colors[c(1,5,6,2,3,4)], axistype = 0, 
           vlcex = 2, plty = 3, plwd = 5, cglcol = "navy",
           title = "TF activity scores", cex.main = 2)

legend(1.5, 0, legend = names(group.colors), 
       fill = group.colors, cex = 2,
       border = group.colors, bty = "n", horiz = F)


# end ---------------------------------------------------------------------

sessionInfo()
