library(factoextra)
library(NbClust)
library(fpc)
library(cluster)
library(vegan)
library(mclust)
library(pheatmap)
library(apcluster)

setwd("~/icSHAPE_Project")
vitroShapeFolder <- "SRSF1_A_100_SHAPE_VITRO_20_0//"
vivoShapeFolder <- "SRSF1_A_100_SHAPE_VIVO_20_0//"

loadShapeFile <- function(shapePath) {
  shape <- read.table(shapePath, 
                      sep = "\t",
                      header = FALSE)
  return(shape)
}

vitroFiles <- list.files(path = vitroShapeFolder, 
                         all.files = FALSE, 
                         full.names = TRUE,
                         recursive = FALSE,
                         pattern = "\\-.SHAPE")
vivoFiles <- list.files(path = vivoShapeFolder, 
                        all.files = FALSE, 
                        full.names = TRUE,
                        recursive = FALSE,
                        pattern = "\\-.SHAPE")

getClLocvitro <- function(file) {
  siteLoc <- substr(file, nchar(vitroShapeFolder) + 1, nchar(file) - 6)
  return(siteLoc)
}

getClLocvivo <- function(file) {
  siteLoc <- substr(file, nchar(vivoShapeFolder) + 1, nchar(file) - 6)
  return(siteLoc)
}

vitrositeNames <- lapply(vitroFiles, getClLocvitro)
vivositeNames <- lapply(vivoFiles, getClLocvivo)


vitroShapeTabs <- lapply(vitroFiles, loadShapeFile)
vivoShapeTabs <- lapply(vivoFiles, loadShapeFile)

vitroshape <- combineCols(vitroShapeTabs)
rownames(vitroshape) <- vitrositeNames
ncol(vitroshape)
colnames(vitroshape) <- 1:ncol(vitroshape)

vivoshape <- combineCols(vivoShapeTabs)
rownames(vivoshape) <- vivositeNames
ncol(vivoshape)
colnames(vivoshape) <- 1:ncol(vivoshape)

invitrovivo<-which(rownames(vitroshape) %in% rownames(vivoshape))    
vitroshapeinvivo<-vitroshape[invitrovivo, ]

shapediffs <- vitroshapeinvivo - vivoshape
shortshape <- shapediffs[,((ncol(shapediffs)/2-10.5):(ncol(shapediffs)/2+10.5))]

combineCols <- function(shapeDfs) {
  shape <- data.frame()
  for (item in shapeDfs) {
    shape <- rbind(shape, item$"V2")
  }
  nums <- 1:length(shapeDfs)
  rownames(shape) <- 
    return(shape)
}

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}




HEATMAP <- pheatmap(data.matrix(shapediffs), clustering_distance_rows = "euclidean", clustering_method = "complete", cluster_cols = FALSE, fontsize_row = 10, kmeans_k = 20, cutree_rows = 2, main = vitroShapeFolder)
SHORTMAP <- pheatmap(data.matrix(shortshape), clustering_distance_rows = "euclidean", clustering_method = "complete", cluster_cols = FALSE, fontsize_row = 20, kmeans_k = 3, cutree_rows = 2, main = vitroShapeFolder)
CLUSTERDF <- as.data.frame(HEATMAP$kmeans$cluster, header = F )
SHORTCLUSTERDF <- as.data.frame(SHORTMAP$kmeans$cluster, header = F )
HEATMAP <- pheatmap(data.matrix(shapediffs), cluster_cols = FALSE, fontsize_row = 10)
save_pheatmap_pdf(HEATMAP, paste(vitroShapeFolder, "CLUSTMAP_NCLUST_", length(HEATMAP$tree_row$labels), ".pdf", sep = ""))
save_pheatmap_pdf(SHORTMAP, paste(vitroShapeFolder, "SHORTCLUSTMAP_DIFF", ".pdf", sep = ""))
write.table(HEATMAP$kmeans$cluster, file = paste(vitroShapeFolder, "KMEANSCLUST.txt", sep = ""), quote = FALSE, sep = '\t' )
write.table(SHORTMAP$kmeans$cluster, file = paste(vitroShapeFolder, "SHORT_KMEANSCLUST.txt", sep = ""), quote = FALSE, sep = '\t' )
heatmap(data.matrix(shapediffs), Rowv = FALSE, Colv = TRUE, )

#ELBOW CLUSTERING
set.seed(13)
wss <- (nrow(shapediffs))*sum(apply(shapediffs,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(shapediffs,
                                     centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

#SILHOUETTE ANALYSIS
pamk.best2 <- pamk(shapediffs)
cat("number of clusters estimated by optimum average silhouette width:", pamk.best2$nc, "\n")
plot(pam(shapediffs, pamk.best2$nc))

#GAP STATISTICS
clusGap(shapediffs, kmeans, 10, B = 100, verbose = interactive())

#HIERARCHICAL CLUSTERING
shape_dist <- dist(as.matrix(shapediffs))   # find distance matrix
plot(hclust(shape_dist))

#CALINSKY CRITERION
cal_fit2 <- cascadeKM(shapediffs, 1, 10, iter = 1000)
plot(cal_fit2, sortg = TRUE, grpmts.plot = TRUE)
calinski.best2 <- as.numeric(which.max(cal_fit2$results[2,]))
cat("Calinski criterion optimal number of clusters:", calinski.best2, "\n")

#BAYESIAN - MAXIMIZATION
d_clust2 <- Mclust(as.matrix(shapediffs), G=1:20)
m.best2 <- dim(d_clust2$z)[2]
cat("model-based optimal number of clusters:", m.best2, "\n")
plot(d_clust2)


#HEIRARCHICAL CLUSTERING
set.seed(13)
d.apclus2 <- apcluster(negDistMat(r=2), shapediffs)
cat("affinity propogation optimal number of clusters:", length(d.apclus2@clusters), "\n")
heatmap(d.apclus2)




#ELBOW CLUSTERING
set.seed(13)
wss <- (nrow(shortshape))*sum(apply(shortshape,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(shortshape,
                                     centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

#SILHOUETTE ANALYSIS
pamk.best2 <- pamk(shortshape)
cat("number of clusters estimated by optimum average silhouette width:", pamk.best2$nc, "\n")
plot(pam(shortshape, pamk.best2$nc))

#GAP STATISTICS
clusGap(shortshape, kmeans, 10, B = 100, verbose = interactive())

#HIERARCHICAL CLUSTERING
shortshape_dist <- dist(as.matrix(shortshape))   # find distance matrix
plot(hclust(shortshape_dist))

#CALINSKY CRITERION
cal_fit2 <- cascadeKM(shortshape, 1, 10, iter = 1000)
plot(cal_fit2, sortg = TRUE, grpmts.plot = TRUE)
calinski.best2 <- as.numeric(which.max(cal_fit2$results[2,]))
cat("Calinski criterion optimal number of clusters:", calinski.best2, "\n")

#BAYESIAN - MAXIMIZATION
d_clust2 <- Mclust(as.matrix(shortshape), G=1:20)
m.best2 <- dim(d_clust2$z)[2]
cat("model-based optimal number of clusters:", m.best2, "\n")
plot(d_clust2)


#HEIRARCHICAL CLUSTERING
set.seed(13)
d.apclus2 <- apcluster(negDistMat(r=2), shortshape)
cat("affinity propogation optimal number of clusters:", length(d.apclus2@clusters), "\n")
heatmap(d.apclus2)

combineCols <- function(shapeDfs) {
  shape <- data.frame()
  for (item in shapeDfs) {
    shape <- rbind(shape, item$"V2")
  }
  nums <- 1:length(shapeDfs)
  rownames(shape) <- 
    return(shape)
}

shapeav <- combineCols(shapediffs)

write.table(shapediffs, file = "SHAPEDIFFS.txt", quote = FALSE, sep = '\t')