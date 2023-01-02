

# based on 5mC ------------------------------------------------------------


library(fmsb)

rm(list=ls())

setwd("~/Dropbox/BioInfo/Lab/Tcell_ID")


data <- read.csv(file = "5mC.targeted.heatmap.data.csv")
data <- data.frame(t(data[c(1,3,4,8,15,16,17), ]))
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

group.colors <- alphabet(20)
group.colors <- group.colors[c(5,18,15,19,3,7)]
names(group.colors) <- c("Th0","Th1","Th17","Th17_noTgfb","Th2","Treg") 
#group.colors =  list(group = group.colors)

par(mfrow=c(3,3))

plot.data <- rbind(max_min, as.numeric(data[1, ]))
radarchart(plot.data, pfcol = group.colors[1], axistype = 3, vlcex = 2, title = "Th0", cex.main = 2)
plot.data <- rbind(max_min, as.numeric(data[2, ]))
radarchart(plot.data, pfcol = group.colors[1], axistype = 3, vlcex = 2, title = "Th0", cex.main = 2)
plot.data <- rbind(max_min, as.numeric(data[3, ]))
radarchart(plot.data, pfcol = group.colors[1], axistype = 3, vlcex = 2, title = "Th0", cex.main = 2)

plot.data <- rbind(max_min, as.numeric(data[12, ]))
radarchart(plot.data, pfcol = group.colors[5], axistype = 3, vlcex = 2, title = "Th2", cex.main = 2)
plot.data <- rbind(max_min, as.numeric(data[13, ]))
radarchart(plot.data, pfcol = group.colors[5], axistype = 3, vlcex = 2, title = "Th2", cex.main = 2)
plot.data <- rbind(max_min, as.numeric(data[14, ]))
radarchart(plot.data, pfcol = group.colors[5], axistype = 3, vlcex = 2, title = "Th2", cex.main = 2)

plot.data <- rbind(max_min, as.numeric(data[15, ]))
radarchart(plot.data, pfcol = group.colors[6], axistype = 3, vlcex = 2, title = "Treg", cex.main = 2)
plot.data <- rbind(max_min, as.numeric(data[16, ]))
radarchart(plot.data, pfcol = group.colors[6], axistype = 3, vlcex = 2, title = "Treg", cex.main = 2)

par(mfrow=c(3,3))

plot.data <- rbind(max_min, as.numeric(data[4, ]))
radarchart(plot.data, pfcol = group.colors[2], axistype = 3, vlcex = 2, title = "Th1", cex.main = 2)
plot.data <- rbind(max_min, as.numeric(data[5, ]))
radarchart(plot.data, pfcol = group.colors[2], axistype = 3, vlcex = 2, title = "Th1", cex.main = 2)
plot.data <- rbind(max_min, as.numeric(data[6, ]))
radarchart(plot.data, pfcol = group.colors[2], axistype = 3, vlcex = 2, title = "Th1", cex.main = 2)

plot.data <- rbind(max_min, as.numeric(data[7, ]))
radarchart(plot.data, pfcol = group.colors[3], axistype = 3, vlcex = 2, title = "Th17", cex.main = 2)
plot.data <- rbind(max_min, as.numeric(data[8, ]))
radarchart(plot.data, pfcol = group.colors[3], axistype = 3, vlcex = 2, title = "Th17", cex.main = 2)
plot.data <- rbind(max_min, as.numeric(data[9, ]))
radarchart(plot.data, pfcol = group.colors[3], axistype = 3, vlcex = 2, title = "Th17", cex.main = 2)

plot.data <- rbind(max_min, as.numeric(data[10, ]))
radarchart(plot.data, pfcol = group.colors[4], axistype = 3, vlcex = 2, title = "Th17/Th1", cex.main = 2)
plot.data <- rbind(max_min, as.numeric(data[11, ]))
radarchart(plot.data, pfcol = group.colors[4], axistype = 3, vlcex = 2, title = "Th17/Th1", cex.main = 2)


# grouped spidernet plots

par(mfrow=c(1,2))

plot.data <- rbind(max_min, as.numeric(data[2, ]))
plot.data <- rbind(plot.data, as.numeric(data[12, ]))
plot.data <- rbind(plot.data, as.numeric(data[15, ]))
radarchart(plot.data, pcol = group.colors[c(1,5,6)], axistype = 3, 
           vlcex = 2, plty = 3, plwd = 5, cglcol = "navy",
           title = "Th0 vs Th2 vs Treg", cex.main = 2)

plot.data <- rbind(max_min, as.numeric(data[5, ]))
plot.data <- rbind(plot.data, as.numeric(data[8, ]))
plot.data <- rbind(plot.data, as.numeric(data[10, ]))
radarchart(plot.data, pcol = group.colors[c(2,3,4)], axistype = 3, 
           vlcex = 2, plty = 3, plwd = 5, cglcol = "navy",
           title = "Th1 vs Th17 vs Th17/Th1", cex.main = 2)



# end ---------------------------------------------------------------------

sessionInfo()
