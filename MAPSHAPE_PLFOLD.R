setwd("~/icSHAPE_Project")
library(pheatmap)

plFolder <- "ADAR_RNAplfold//"
vitroShapeFolder <- "ADAR_22_SHAPE_vitro_20_0//"
vivoShapeFolder <- "ADAR_22_SHAPE_VIVO_20_0//"

plFolder <- "ZRANB2_RNAplfold/"
vitroShapeFolder <- "ZRANB2_SHAPE_VITRO_20_0//"
vivoShapeFolder <- "ZRANB2_SHAPE_VIVO_20_0//"

plFolder <- "ZRANB2_RNAplfold/"
vitroShapeFolder <- "ZRANB2_SHAPE_VITRO_20_0//"
vivoShapeFolder <- "ZRANB2_SHAPE_VIVO_20_0//"

plFolder <- "SRSF1_plfold//"
vitroShapeFolder <- "SRSF1_A_SHAPE_VITRO_20_0//"
vivoShapeFolder <- "SRSF1_A_SHAPE_VIVO_20_0//"

plFolder <- "SP1_A_vitro_plfold/"
vitroShapeFolder <- "SP1_A_SHAPE_VITRO_20_0//"
vivoShapeFolder <- "SP1_A_SHAPE_VIVO_20_0//"

plFolder <- "WT1_plfold//"
vitroShapeFolder <- "WT1_A_SHAPE_VITRO_20_0//"
vivoShapeFolder <- "WT1_A_SHAPE_VIVO_20_0//"

plFolder <- "MAZ_plfold//"
vitroShapeFolder <- "MAZ_A_SHAPE_VITRO_20_0//"
vivoShapeFolder <- "MAZ_A_SHAPE_VIVO_20_0//"

plFolder <- "WT1_plfold_100//"
vitroShapeFolder <- "WT1_B_100_SHAPE_VITRO_50_0//"
vivoShapeFolder <- "WT1_B_100_SHAPE_VIVO_50_0//"

plFolder <- "ADAR_plfold_100//"
vitroShapeFolder <- "ADAR_50_SHAPE_vitro_20_0//"
vivoShapeFolder <- "ADAR_50_SHAPE_VIVO_20_0//"

plFolder <- "ZRANB2_plfold_100/"
vitroShapeFolder <- "ZRANB2_100_A_SHAPE_VITRO_20_0//"
vivoShapeFolder <- "ZRANB2_100_A_SHAPE_VIVO_20_0//"

plFolder <- "SRSF1_plfold_100_VITRO//"
vitroShapeFolder <- "SRSF1_A_100_SHAPE_VITRO_20_0//"
vivoShapeFolder <- "SRSF1_A_100_SHAPE_VIVO_20_0//"

plFolder <- "ZNF121_PLFOLD_100//"
vitroShapeFolder <- "ZNF121_A_100_SHAPE_VITRO_20_0//"
vivoShapeFolder <- "ZNF121_A_100_SHAPE_VIVO_20_0//"

plFolder <- "MAZ_A_PLFOLD_100///"
vitroShapeFolder <- "MAZ_A_100_SHAPE_VITRO_20_0//"
vivoShapeFolder <- "MAZ_A_100_SHAPE_VIVO_20_0//"

plFolder <- "METTL3_B_PLFOLD_100/"
vitroShapeFolder <- "METTL3_B_100_SHAPE_VITRO_20_0//"
vivoShapeFolder <- "METTL3_B_100_SHAPE_VIVO_20_0//"
pname <- "METTL3_100"

plFolder <- "YY1_B_PLFOLD_100/"
vitroShapeFolder <- "YY1_B_100_SHAPE_VITRO_20_0//"
vivoShapeFolder <- "YY1_B_100_SHAPE_VIVO_20_0//"
pname <- "YY1_100"

plFolder <- "SP1_A_PLFOLD_100/"
vitroShapeFolder <- "SP1_A_100_SHAPE_VITRO_20_0//"
vivoShapeFolder <- "SP1_A_100_SHAPE_VIVO_20_0//"
pname <- "SP1_100"

plFolder <- "ZBTB48_A_PLFOLD_100//"
vitroShapeFolder <- "ZBTB48_A_100_SHAPE_VITRO_20_0//"
vivoShapeFolder <- "ZBTB48_A_100_SHAPE_VIVO_20_0//"
pname <- "ZBTB48_100"

plFolder <- "CTCF_A_PLFOLD_100//"
vitroShapeFolder <- "CTCF_A_100_SHAPE_VITRO_20_0//"
vivoShapeFolder <- "CTCF_A_100_SHAPE_VIVO_20_0//"
pname <- "CTCF_100"


#35 (70)

plFolder <- "ZRANB2_35_PLFOLD///"
vitroShapeFolder <- "ZRANB2_35_SHAPE_VITRO_20_0///"
vivoShapeFolder <- "ZRANB2_35_SHAPE_VIVO_20_0///"
pname <- "ZRANB2_35"

plFolder <- "ADAR_35_PLFOLD////"
vitroShapeFolder <- "ADAR_35_SHAPE_VITRO_20_0////"
vivoShapeFolder <- "ADAR_35_SHAPE_VIVO_20_0////"
pname <- "ADAR_35"

plFolder <- "WT1_B_35_PLFOLD////"
vitroShapeFolder <- "WT1_B_35_SHAPE_VITRO_20_0////"
vivoShapeFolder <- "WT1_B_35_SHAPE_VIVO_20_0////"
pname <- "WT1_35"

plFolder <- "YY1_B_35_PLFOLD////"
vitroShapeFolder <- "YY1_B_35_SHAPE_VITRO_20_0////"
vivoShapeFolder <- "YY1_B_35_SHAPE_VIVO_20_0////"
pname <- "YY1_35"

plFolder <- "SRSF1_A_35_PLFOLD//////"
vitroShapeFolder <- "SRSF1_35_SHAPE_VITRO_20_0///////"
vivoShapeFolder <- "SRSF1_35_SHAPE_VIVO_20_0//////"
pname <- "SRSF1_35"

plFolder <- "CTCF_A_35_PLFOLD/////"
vitroShapeFolder <- "CTCF_A_35_SHAPE_VITRO_20_0/////"
vivoShapeFolder <- "CTCF_A_35_SHAPE_VIVO_20_0/////"
pname <- "CTCF_35"

plFolder <- "SP1_A_35_PLFOLD//////"
vitroShapeFolder <- "SP1_A_35_SHAPE_VITRO_20_0//////"
vivoShapeFolder <- "SP1_A_35_SHAPE_VIVO_20_0//////"
pname <- "SP1_35"


############### 150
plFolder <- "SP1_150_PLFOLD///////"
vitroShapeFolder <- "SP1_A_75_SHAPE_VITRO_20_0///////"
vivoShapeFolder <- "SP1_A_75_SHAPE_VIVO_20_0///////"
pname <- "SP1_75"

############### 10
plFolder <- "SP1_10_PLFOLD////////"
vitroShapeFolder <- "SP1_A_10_SHAPE_VITRO_20_0////////"
vivoShapeFolder <- "SP1_A_10_SHAPE_VIVO_20_0////////"
pname <- "SP1_10"

############### 10
plFolder <- "SP1_E_50_PLFOLD/////////"
vitroShapeFolder <- "SP1_E_SEP15_SHAPE_VITRO_20_0/////////"
vivoShapeFolder <- "SP1_E_SEP15_SHAPE_VIVO_20_0/////////"
pname <- "SP1_E_50"

############### 10
plFolder <- "SP1_A_35_PLFOLD/"
vitroShapeFolder <- "SP1_A_35_SHAPE_VITRO_20_0//////////"
vivoShapeFolder <- "SP1_A_35_SHAPE_VIVO_20_0//////////"
plshuffle <- "SP1_A_UNI_SHUFFLE_PLFOLD///////////"
pname <- "SP1 "

############### ZRANB2 MINUS
plFolder <- "ZRANB2_MERGED_MINUS_35_PLFOLD//"
vitroShapeFolder <- "ZRANB2_MERGED_MINUS_35_SEP16_SHAPE_VITRO_20_0///////////"
vivoShapeFolder <- "ZRANB2_MERGED_MINUS_35_SEP16_SHAPE_VIVO_20_0///////////"
plshuffle <- "SP1_A_UNI_SHUFFLE_PLFOLD///////////"
pname <- "ZRANB2_MERGED_MINUS_35"

############### ZRANB2 PLUS
plFolder <- "ZRANB2_MERGED_PLUS_35_PLFOLD//"
vitroShapeFolder <- "ZRANB2_MERGED_PLUS_35_SEP16_SHAPE_VITRO_20_0///////////"
vivoShapeFolder <- "ZRANB2_MERGED_PLUS_35_SEP16_SHAPE_VIVO_20_0///////////"
plshuffle <- "SP1_A_UNI_SHUFFLE_PLFOLD///////////"
pname <- "ZRANB2_MERGED_PLUS_35"

############### SP1 PLUS
plFolder <- "SP1_MERGED_PLUS_50_PLFOLD///"
vitroShapeFolder <- "SP1_MERGED_PLUS_50_SEP16_SHAPE_VITRO_20_0////////////"
vivoShapeFolder <- "SP1_MERGED_PLUS_50_SEP16_SHAPE_VIVO_20_0////////////"
#plshuffle <- "SP1_A_UNI_SHUFFLE_PLFOLD///////////"
plshuffle <- "200906_2_thUni_merged.bam.crosslink_sites_NOV25_FULL_VITRO_50_SHUFFLE_PLFOLD"
pname <- "SP1_MERGED_PLUS_50"

############### SP1 MINUS
plFolder <- "SP1_MERGED_MINUS_50_PLFOLD///"
vitroShapeFolder <- "SP1_MERGED_MINUS_50_SEP16_SHAPE_VITRO_20_0////////////"
vivoShapeFolder <- "SP1_MERGED_MINUS_50_SEP16_SHAPE_VIVO_20_0////////////"
plshuffle <- "SP1_MERGE///////////"
pname <- "SP1_MERGED_MINUS_50"

#200906_1
plFolder <- "S200906_1_thUni_merged.bam.crosslink_sites_NOV25_50_NEW_plfold"
vitroShapeFolder <- "S200906_1_thUni_merged.bam.crosslink_sites_NOV25_SHAPE_VITRO_50_NEW_0/"
vivoShapeFolder <- "S200906_1_thUni_merged.bam.crosslink_sites_NOV25_SHAPE_VIVO_50_NEW_0/"
plshuffle <- "ADAR_NEW_1_50_SHUFFLED_PLFOLD/"
pname <- "ADAR_100_SHUFFLE"

#200906_2
plFolder <- "200906_2_thUni_merged.bam.crosslink_sites_NOV25_FULL_VITRO_50_PLFOLD/"
vitroShapeFolder <- "200906_2_thUni_merged.bam.crosslink_sites_NOV25_FULL_NEWL_SHAPE_VITRO_50_20_0//"
vivoShapeFolder <- "200906_2_thUni_merged.bam.crosslink_sites_NOV25_FULL_NEWL_SHAPE_VIVO_20_0_50/"
plshuffle <- "200906_2_thUni_merged.bam.crosslink_sites_NOV25_FULL_VITRO_50_SHUFFLE_PLFOLD//"
pname <- "ADAR_2_100_PLFOLD"

plshuffle<- "ADAR_NEW_PLFOLD_SHUFFLE"


##### ADAR_A1
vitroShapeFolder <- "ADAR1_A1.sh.bam.crosslink_sites_JAN14_FULL_SHAPE_VITRO_20_0/"
vivoShapeFolder <- "ADAR1_A1.sh.bam.crosslink_sites_JAN14_FULL_SHAPE_VIVO_20_0/"
pname <- "ADAR_A1"

##### ADAR_A1
vitroShapeFolder <- "ADAR1_A1.sh.bam.crosslink_sites_JAN14_FULL_SHAPE_VITRO_20_0/"
vivoShapeFolder <- "ADAR1_A1.sh.bam.crosslink_sites_JAN14_FULL_SHAPE_VIVO_20_0/"
pname <- "ADAR_A1"

##### ADAR_REP1
vitroShapeFolder <- "ADAR_REP1_MERGED.thUni.bam.crosslink_sites_JAN14_FULL_SHAPE_VITRO_20_0/"
vivoShapeFolder <- "ADAR_REP1_MERGED.thUni.bam.crosslink_sites_JAN14_FULL_SHAPE_VIVO_20_0/"
pname <- "ADAR_REP1"

plFolder <- "SP1_B_WINPUT_PLFOLD////"
vitroShapeFolder <- "SP1_B_WINPUT_SEP16_SHAPE_VITRO_20_0/////"
vivoShapeFolder <- "SP1_B_WINPUT_SEP16_SHAPE_VIVO_20_0/////"
pname <- "SP1_B_WINPUT"

plFolder <- "YY1_A_WINPUT_35_SEPT16_PLFOLD/////"
vitroShapeFolder <- "YY1_A_WINPUT_35_SEP16_SHAPE_VITRO_20_0///////"
vivoShapeFolder <- "YY1_A_WINPUT_35_SEP16_SHAPE_VIVO_20_0//////"
pname <- "YY1_A_35_WITHINPUT"

plFolder <- "YY1_A_WINPUT_35_SEPT16_PLFOLD/////"
vitroShapeFolder <- "ZRANB2_MERGED_MINUS_20_SEP16_SHAPE_VITRO_20_0////////"
vivoShapeFolder <- "ZRANB2_MERGED_MINUS_20_SEP16_SHAPE_VIVO_20_0///////"
pname <- "ZRANB2_MERGED_MINUS_20"

plFolder <- "YY1_A_WINPUT_35_SEPT16_PLFOLD/////"
vitroShapeFolder <- "ZRANB2_MERGED_MINUS_35_SEP16_SHAPE_VITRO_20_0/////////"
vivoShapeFolder <- "ZRANB2_MERGED_MINUS_35_SEP16_SHAPE_VIVO_20_0///////"
pname <- "ZRANB2_MERGED_MINUS_35"

plFolder <- "YY1_A_WINPUT_35_SEPT16_PLFOLD/////"
vitroShapeFolder <- "ZRANB2_MERGED_PLUS_35_SEP16_SHAPE_VITRO_20_0/////////"
vivoShapeFolder <- "ZRANB2_MERGED_PLUS_35_SEP16_SHAPE_VIVO_20_0///////"
pname <- "ZRANB2_MERGED_PLUS_35"

vitroShapeFolder <- "ZRANB2_MERGED_PLUS_35_SEP16_SHAPE_VITRO_20_0/////////"
vivoShapeFolder <- "ZRANB2_MERGED_PLUS_35_SEP16_SHAPE_VIVO_20_0///////"
pname <- "ZRANB2_MERGED_PLUS_35"

############### NOSHAPE
plFolder <- "SRSF1_FULL_PLFOLD"
pname <- "SRSF1_FULL_PLFOLD"
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

#PDF
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

