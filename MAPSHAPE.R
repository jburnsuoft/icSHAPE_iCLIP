library(pheatmap)

#Define plFold and icSHAPE dataset folders
#plFolder <- 
#vitroShapeFolder <- 
#vivoShapeFolder <- 
#pname <- 


loadplfold <- function(shapePath) {
  shape <- read.table(shapePath, 
                      skip = 2,
                      header = FALSE)
  return(shape)
}

loadShapeFile <- function(shapePath) {
  shape <- read.table(shapePath, 
                      sep = "\t",
                      header = FALSE)
  return(shape)
}

loadshuffleFile <- function(shapePath) {
  shape <- read.table(shapePath, 
                      skip = 2,
                      header = FALSE)
  return(shape)
}

plfiles <- list.files(path = plFolder, 
                      all.files = FALSE, 
                      full.names = TRUE,
                      recursive = FALSE)
vitroFiles <- list.files(path = vitroShapeFolder, 
                         all.files = FALSE, 
                         full.names = TRUE,
                         recursive = FALSE)
vivoFiles <- list.files(path = vivoShapeFolder, 
                        all.files = FALSE, 
                        full.names = TRUE,
                        recursive = FALSE)
shufflefiles <- list.files(path = plshuffle, 
                           all.files = FALSE, 
                           full.names = TRUE,
                           recursive = FALSE)


plTabs <- lapply(plfiles, loadplfold)
vitroShapeTabs <- lapply(vitroFiles, loadShapeFile)
vivoShapeTabs <- lapply(vivoFiles, loadShapeFile)
shuffletabs <- lapply(shufflefiles, loadplfold)

combineCols <- function(shapeDfs) {
  shape <- data.frame()
  for (item in shapeDfs) {
    shape <- rbind(shape, item$"V2")
  }
  nums <- 1:length(shapeDfs)
  rownames(shape) <- 
    return(shape)
}

plCombine <- combineCols(plTabs)
vitroShape <- combineCols(vitroShapeTabs)
vivoShape <- combineCols(vivoShapeTabs)
shufflecombine <- combineCols(shuffletabs)
#rownames(shape) <- vitroSiteNames
#colnames(shape) <- 1:41

plMatrix <- data.matrix(plCombine)
vitroShapeMatrix <- data.matrix(vitroShape)
vivoShapeMatrix <- data.matrix(vivoShape)
shuffleMatrix <- data.matrix(shufflecombine)

#Plot and Save PDF
plColAvg <- apply(plMatrix, 2, median)
vitroColAvg <- apply(vitroShapeMatrix, 2, median)
vivoColAvg <- apply(vivoShapeMatrix, 2, median)
shuffleColAvg <- apply(shuffleMatrix, 2, median)
set.seed(1)
background <- sample(plColAvg)


pdf(file=paste(pname, "new_plot_shuffle.pdf", sep = ""), width = 8, height = 8)
plot(c(0, length(vitroColAvg)), c(0, 1), type="n", xlab="Position",
     ylab="Unpaired Probability", cex.lab = 1.2)



lines(1:length(plColAvg), plColAvg, col = "gray", lwd = "3")
lines(1:length(vitroColAvg), vitroColAvg, col = "blue", lwd = "3")
lines(1:length(vivoColAvg), vivoColAvg, col = "orange", lwd = "3")
lines(1:length(shuffleColAvg), shuffleColAvg, col = "green", lwd = "3")
#lines(1:length(background), background, col = "green", lwd = "3")
title(pname)
#legend("topleft", c("RNAplFold", "vitro", "vivo"), col = c("gray", "blue", "orange"), lty=1:2, cex=0.5)
dev.off()

#MEDIAN
plColAvg <- apply(plMatrix, 2, median)
vitroColAvg <- apply(vitroShapeMatrix, 2, median)
vivoColAvg <- apply(vivoShapeMatrix, 2, median)
set.seed(1)
background <- sample(vitroColAvg)

jpeg(file=paste(pname, "new_plot_median.jpeg", sep = ""), width = 463, height = 463 )
plot(c(0, length(vitroColAvg)), c(0, 1), type="n", xlab="Position",
     ylab="Unpaired Probability" )



plMatrix <- data.matrix(plCombine)
vitroShapeMatrix <- data.matrix(vitroShape)
vivoShapeMatrix <- data.matrix(vivoShape)
shuffleMatrix <- data.matrix(shufflecombine)
heatmap(vitroShapeMatrix, Rowv = NA, Colv = NA)
MAP <- pheatmap(vitroShapeMatrix, cluster_cols = FALSE, cluster_rows = TRUE, clustering_method = "complete", main = pname, kmeans_k = 10, annotation_names_row = TRUE)
CLUSTOUT <- as.data.frame(MAP$kmeans$cluster)
rownames(CLUSTOUT) <- vitroFiles
VFNAMES <- as.data.frame(vitroFiles)
rownames(VFNAMES) <- vitroFiles
TRY2 <- pheatmap(vitroShapeMatrix, annotation_names_row = TRUE)


TRY <- MAP$kmeans$cluster
MAPS <- pheatmap(vitroShapeMatrix, annotation_row = annotation_row)
as.dendrogram(MAPS)
hc <- MAP$tree_row
lbl <- cutree(hc, 10)
which(lbl==9)
MAP$kmeans$centers

hc <- as.hclust(MAP$tree_row)
MAPCLUSTS <- (MAP$kmeans$cluster)

MAPCLUST <- cbind(vitroShapeMatrix, cluster = cutree(MAPS$tree_col, k = 10))


plcrosslink <- plMatrix[,36]
write.table(plcrosslink, sep = "\t", file = paste(pname,"_plcrosslink.txt", sep=""), quote = FALSE)




#MEDIAN
plColAvg <- apply(plMatrix, 2, median)
vitroColAvg <- apply(vitroShapeMatrix, 2, median)
vivoColAvg <- apply(vivoShapeMatrix, 2, median)
shuffleColAvg <- apply(shuffleMatrix, 2, median)
set.seed(1)
background <- sample(vitroColAvg)

plot(c(0, 101), c(0, 1), type="n", xlab="Position",
     ylab="Unpaired Probability")

lines(1:101, plColAvg, col = "purple")
lines(1:101, vitroColAvg, col = "blue")
lines(1:101, vivoColAvg, col = "orange")
lines(1:101, shuffleColAvg, col = "green")
title(paste(pname,"Structure Profiles"))
legend("topleft", c("RNAplFold", "vitro", "vivo", "plFold Shuffle"), col = c("purple", "blue", "orange", "green"), lty=1, cex=0.8)


jpeg(file=paste(pname, "plot_median.jpeg", sep = ""), width = 463, height = 463 )
plot(c(0, length(vitroColAvg)), c(0, 1), type="n", xlab="Position",
     ylab="Unpaired Probability" )

lines(1:length(plColAvg), plColAvg, col = "purple")
lines(1:length(vitroColAvg), vitroColAvg, col = "blue")
lines(1:length(vivoColAvg), vivoColAvg, col = "orange")
lines(1:length(shuffleColAvg), shuffleColAvg, col = "green")
title(paste(pname,"RNAplfold_Shuffle"))
legend("topleft", c("RNAplFold", "vitro", "vivo", "plFold Shuffle"), col = c("purple", "blue", "orange", "green"), lty=1:2, cex=0.8)
dev.off()

#MEAN
plColAvg <- apply(plMatrix, 2, mean)
vitroColAvg <- apply(vitroShapeMatrix, 2, median)
vivoColAvg <- apply(vivoShapeMatrix, 2, median)
shuffleColAvg <- apply(shuffleMatrix, 2, mean)
set.seed(1)
background <- sample(vitroColAvg)

plot(c(0, 101), c(0, 1), type="n", xlab="Position",
     ylab="Unpaired Probability")

lines(1:101, plColAvg, col = "purple")
lines(1:101, vitroColAvg, col = "blue")
lines(1:101, vivoColAvg, col = "orange")
lines(1:101, shuffleColAvg, col = "green")
title(paste(pname,"Structure Profiles"))
legend("topleft", c("RNAplFold", "vitro", "vivo", "plFold Shuffle"), col = c("purple", "blue", "orange", "green"), lty=1:2, cex=0.5)


jpeg(file=paste(pname, "plot_mean.jpeg", sep = ""), width = 463, height = 463 )
plot(c(0, length(vitroColAvg)), c(0, 1), type="n", xlab="Position",
     ylab="Unpaired Probability" )

lines(1:length(plColAvg), plColAvg, col = "purple")
lines(1:length(vitroColAvg), vitroColAvg, col = "blue")
lines(1:length(vivoColAvg), vivoColAvg, col = "orange")
lines(1:length(shuffleColAvg), shuffleColAvg, col = "green")
title(paste(pname,"RNAplfold"))
legend("topleft", c("RNAplFold", "vitro", "vivo", "plFold Shuffle"), col = c("purple", "blue", "orange", "green"), lty=1:2, cex=0.5)
dev.off()




plColAvg <- apply(plMatrix, 2, mean)
vitroColAvg <- apply(vitroShapeMatrix, 2, mean)
vivoColAvg <- apply(vivoShapeMatrix, 2, mean)
set.seed(1)
background <- sample(vitroColAvg)

plot(c(0, 71), c(0, 1), type="n", xlab="Position",
     ylab="Unpaired Probability" )

lines(1:71, plColAvg, col = "gray")
lines(1:71, vitroColAvg, col = "blue")
lines(1:71, vivoColAvg, col = "orange")
lines(1:71, shuffleColAvg, col = "green")
title(paste(pname,"RNAplfold"))
legend("topleft", c("RNAplFold", "vitro", "vivo", "plFold Shuffle"), col = c("purple", "blue", "orange", "green"), lty=1:2, cex=0.5)

