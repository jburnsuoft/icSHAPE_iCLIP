library(pheatmap)
#### R Script Used to Average and Display Resulting Overlap of icSHAPE and iCLIP Coordinates Determined by CLIPSHAPE.R

### Define plFold and icSHAPE dataset folders

## Protein Name For Identification
#pname <- 
## Folder Containing Individual in vitro icSHAPE .SHAPE Files for Protein of Interest
#vitroShapeFolder <- 
## Folder Containing Individual in vivo icSHAPE .SHAPE Files for Protein of Interest
#vivoShapeFolder <- 
## Folder Containing Individual in silico Unpaired Probability lunp Files for Protein of Interest
#plFolder <- 


### Read Data From Defined Folders
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

### Format and Combine Data
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

### Generate Data Matrix for Each Dataset
plMatrix <- data.matrix(plCombine)
vitroShapeMatrix <- data.matrix(vitroShape)
vivoShapeMatrix <- data.matrix(vivoShape)
shuffleMatrix <- data.matrix(shufflecombine)

### Plot and Save PDF of Median icSHAPE Profile
plColAvg <- apply(plMatrix, 2, median)
vitroColAvg <- apply(vitroShapeMatrix, 2, median)
vivoColAvg <- apply(vivoShapeMatrix, 2, median)
shuffleColAvg <- apply(shuffleMatrix, 2, median)
set.seed(1)
background <- sample(plColAvg)

pdf(file=paste(pname, "median_plot.pdf", sep = ""), width = 8, height = 8)
  plot(c(0, length(vitroColAvg)), c(0, 1), type="n", xlab="Position", ylab="Unpaired Probability", cex.lab = 1.2)
  lines(1:length(plColAvg), plColAvg, col = "gray", lwd = "3")
  lines(1:length(vitroColAvg), vitroColAvg, col = "blue", lwd = "3")
  lines(1:length(vivoColAvg), vivoColAvg, col = "orange", lwd = "3")
  lines(1:length(shuffleColAvg), shuffleColAvg, col = "green", lwd = "3")
  title(pname)
  legend("topleft", c("RNAplFold", "vitro", "vivo"), col = c("gray", "blue", "orange"), lty=1:2, cex=0.5)
dev.off()
