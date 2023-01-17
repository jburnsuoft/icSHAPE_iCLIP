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

#Generate Data Matrix
plMatrix <- data.matrix(plCombine)
vitroShapeMatrix <- data.matrix(vitroShape)
vivoShapeMatrix <- data.matrix(vivoShape)
shuffleMatrix <- data.matrix(shufflecombine)

#Plot and Save PDF of Median icSHAPE Profile
plColAvg <- apply(plMatrix, 2, median)
vitroColAvg <- apply(vitroShapeMatrix, 2, median)
vivoColAvg <- apply(vivoShapeMatrix, 2, median)
shuffleColAvg <- apply(shuffleMatrix, 2, median)
set.seed(1)
background <- sample(plColAvg)
pdf(file=paste(pname, "new_plot_shuffle.pdf", sep = ""), width = 8, height = 8)

plot(c(0, length(vitroColAvg)), c(0, 1), type="n", xlab="Position", ylab="Unpaired Probability", cex.lab = 1.2)
lines(1:length(plColAvg), plColAvg, col = "gray", lwd = "3")
lines(1:length(vitroColAvg), vitroColAvg, col = "blue", lwd = "3")
lines(1:length(vivoColAvg), vivoColAvg, col = "orange", lwd = "3")
lines(1:length(shuffleColAvg), shuffleColAvg, col = "green", lwd = "3")
title(pname)
legend("topleft", c("RNAplFold", "vitro", "vivo"), col = c("gray", "blue", "orange"), lty=1:2, cex=0.5)
dev.off()



#Plot and Save PDF of Median icSHAPE Profile with Shuffleplfold
plColAvg <- apply(plMatrix, 2, mean)
vitroColAvg <- apply(vitroShapeMatrix, 2, median)
vivoColAvg <- apply(vivoShapeMatrix, 2, median)
shuffleColAvg <- apply(shuffleMatrix, 2, mean)
set.seed(1)
background <- sample(vitroColAvg)

jpeg(file=paste(pname, "plot_mean.jpeg", sep = ""), width = 463, height = 463 )
plot(c(0, length(vitroColAvg)), c(0, 1), type="n", xlab="Position", ylab="Unpaired Probability")
lines(1:length(plColAvg), plColAvg, col = "purple")
lines(1:length(vitroColAvg), vitroColAvg, col = "blue")
lines(1:length(vivoColAvg), vivoColAvg, col = "orange")
lines(1:length(shuffleColAvg), shuffleColAvg, col = "green")
title(paste(pname,"RNAplfold"))
legend("topleft", c("RNAplFold", "vitro", "vivo", "plFold Shuffle"), col = c("purple", "blue", "orange", "green"), lty=1:2, cex=0.5)
dev.off()
