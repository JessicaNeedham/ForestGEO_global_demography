####################################################################
### THIS SCRIPT PROVIDES STATISTICAL JUSTIFICATION FOR THE EIGHT 
### DFTs 
###################################################################
rm(list = ls())
### LIBRARIES ###
library(RColorBrewer)
library(kernlab)
library(fpc)
library(cluster)
library(foreach)
library(doParallel)
library(xtable)
library(plyr)
library(factoextra)

#################################################################
### FUNCTIONS ### 
source('Analysis/Function_files/Cluster_functions.R')

#################################################################
### DATA ###
geog <- 'global'
load(file = sprintf('Analysis/Outputs/%s_pca_mat.RData', geog))
load(file = sprintf('Analysis/Outputs/PCA_%s.RData', geog))

sp.coords <- res$ind$coord
sp.coords <- cbind(pca.matrix[ ,c('Site', 'Mnemonic', 'Species')],
                   sp.coords, stringsAsFactors = FALSE)
pca.coords <- res$ind$coord

#######################################################################
### GO ### 
# these are the clustering algorithms to compare
method.names <- c("kmeans", "pam", "hclust", "kkmeans")
# The number of clusters to use with each algorithm
ind <- seq(5, 10)

dist.coords <- dist(pca.coords)
cluster.kmeans.1 <- kmeans(pca.coords, centers = 1, iter.max = 20, 
                           nstart = 200)
tot.ss <- cluster.kmeans.1$totss

#nproc <- length(ind)
#registerDoParallel(nproc)
ind.coord <- pca.coords

nproc <- length(ind)
registerDoParallel(nproc)
foreach(i = 1:(length(ind)))%do%{
  
  set.seed(100)
  cluster.ls <- c()
  
  k.no <- ind[i]
  cluster.kmeans  <- kmeans(pca.coords, centers = k.no, iter.max = 20, 
                            nstart = 200)
  cluster.pam     <- eclust(pca.coords, "pam", k = k.no, graph = FALSE)
  
  cluster.hclust  <- eclust(pca.coords, "hclust", k = k.no,
                            method = "complete", graph = FALSE)
  cluster.kkmeans <- kkmeans(pca.coords, centers = k.no)
  # Other stats
  cluster.ls$kkmeans$stats  <- cluster.stats(dist.coords,
                                             cluster.kkmeans@.Data)
  cluster.ls$kmeans$stats  <- cluster.stats(dist.coords,
                                            cluster.kmeans$cluster)
  cluster.ls$pam$stats  <- cluster.stats(dist.coords,
                                         cluster.pam$cluster)
  cluster.ls$hclust$stats  <- cluster.stats(dist.coords,
                                            cluster.hclust$cluster)
  
  
  cluster.ls$kmeans$withinss <- get.withinss(ind.coord,
                                             clust.id = cluster.kmeans$cluster, 
                                             centers = cluster.kmeans$centers)
  cluster.ls$pam$withinss <- get.withinss(ind.coord,
                                          clust.id = cluster.pam$cluster, 
                                          centers = cluster.pam$medoids)
  cluster.ls$hclust$withinss <- get.withinss(ind.coord,
                                             clust.id = cluster.hclust$cluster, 
                                             centers = get.centers(cluster.hclust))
  cluster.ls$kkmeans$withinss <- get.withinss(ind.coord,
                                              clust.id = cluster.kkmeans@.Data, 
                                              centers = cluster.kkmeans@centers)
  cluster.ls$kmeans$r2 <-
    get.explained.var(tot.ss, cluster.ls$kmeans$withinss)
  cluster.ls$pam$r2 <-
    get.explained.var(tot.ss, cluster.ls$pam$withinss)
  cluster.ls$hclust$r2 <-
    get.explained.var(tot.ss, cluster.ls$hclus$withinss)
  cluster.ls$kkmeans$r2 <-
    get.explained.var(tot.ss, cluster.ls$kkmeans$withinss)
  
  cluster.ls$kmeans$min.sil   <-
    get.min.sil(silhouette(cluster.kmeans$cluster, dist.coords,
                           Ordered = TRUE))
  cluster.ls$pam$min.sil   <-
    get.min.sil(silhouette(cluster.pam$cluster, dist.coords,
                           Ordered = TRUE))
  cluster.ls$hclust$min.sil   <-
    get.min.sil(silhouette(cluster.hclust$cluster, dist.coords,
                           Ordered = TRUE))
  cluster.ls$kkmeans$min.sil  <-
    get.min.sil(silhouette(cluster.kkmeans@.Data, dist.coords,
                           Ordered = TRUE))
  save(cluster.ls, cluster.kmeans, cluster.pam,
       cluster.hclust, cluster.kkmeans, 
       file = sprintf('Analysis/Outputs/clusterls_%s.RData', k.no))
  print(i)
}

#### RESTART HERE ### 
jpeg(file = "Paper_Figures_v3/silhouette_clusters.jpeg", bg = 'white')
par(mfrow = c(6,4), mar = c(1,1,1,1), oma = c(2,3,3,0))
#for(i in 1:length(ind)){
for(i in 1:length(ind)){ 
  
  k.no <- ind[i]
  load(file = sprintf('Analysis/Outputs/clusterls_%s.RData', k.no))
  
  make.silplot(silhouette(cluster.kmeans$cluster, dist.coords,
                          Ordered = TRUE), method = "kmeans")
  mtext(sprintf('No. clust: %i', i+4), side = 2, line = 1, 
        outer = FALSE, cex = 0.7)
  if(i == 1){
    mtext(method.names[1], side = 3, line = 1, outer = FALSE, 
          cex = 0.7)
  }
  make.silplot(silhouette(cluster.pam$cluster, dist.coords,
                          Ordered = TRUE), method = "pam")
  if(i == 1){
    mtext(method.names[2], side = 3, line = 1, outer = FALSE, 
          cex = 0.7)
  }
  make.silplot(silhouette(cluster.hclust$cluster, dist.coords,
                          Ordered = TRUE), method = "hclust")
  if(i == 1){
    mtext(method.names[3], side = 3, line = 1, outer = FALSE, 
          cex = 0.7)
  }
  make.silplot(silhouette(cluster.kkmeans@.Data, dist.coords,
                          Ordered = TRUE), method = "kkmeans")
  if(i == 1){
    mtext(method.names[4], side = 3, line = 1, outer = FALSE, 
          cex = 0.7)
  }
  if(i == 1){
    mtext(method.names[5], side = 3, line = 1, outer = FALSE, 
          cex = 0.7)
  }
  
}
dev.off()

jpeg(file = 'Paper_Figures_v3/cluster_method_barplots.jpeg', bg = 'white')

par(mfrow = c(length(ind), 4), oma = c(3.5,3,3,0), 
    mar = c(1,2,1,1))

for(i in 1:length(ind)){
  
  k.no <- ind[i]
  load(file = sprintf('Analysis/Outputs/clusterls_%s.RData', k.no))
  
  barplot(c(cluster.ls$kmeans$stats$dunn2, cluster.ls$pam$stats$dunn2,
            cluster.ls$hclust$stats$dunn2, cluster.ls$kkmeans$stats$dunn2),
          main = "", col = "tomato2",
          names.arg = NULL, cex.names = 0.7)
  if(i == 1){
    mtext('Dunn2 index', side = 3, line = 1, cex = 0.8)
  }
  if(i == length(ind)){
    axis(side = 1, labels = method.names, cex = 0.5, 
         at = seq(0.5, 4.5, length.out = 4), las = 2)
  }
  mtext(sprintf('No.clust: %i', i+4), side = 2, line = 2.5, 
        cex = 0.8)
  
  barplot(c(cluster.ls$kmeans$min.sil,
            cluster.ls$pam$min.sil,
            cluster.ls$hclust$min.sil,
            cluster.ls$kkmeans$min.sil),
          main = "", col = "tomato2",
          names.arg = NULL, cex.names = 0.7)
  if(i == 1){
    mtext('Number negative sil', side = 3, line = 1, cex = 0.8)
  }
  if(i == length(ind)){
    axis(side = 1, labels = method.names, cex = 0.5, 
         at = seq(0.5, 4.5, length.out = 4), las = 2)
  }
  
  barplot(c(cluster.ls$kmeans$stats$min.cluster.size,
            cluster.ls$pam$stats$min.cluster.size,
            cluster.ls$hclust$stats$min.cluster.size,
            cluster.ls$kkmeans$stats$min.cluster.size),
          main = "", col = "tomato2",
          names.arg = NULL, cex.names = 0.7)
  if(i == 1){
    mtext('Smallest cluster', side = 3, line = 1, cex = 0.8)
  }
  if(i == length(ind)){
    axis(side = 1, labels = method.names, cex = 0.5, 
         at = seq(0.5, 4.5, length.out = 4), las = 2)
  }
  
  barplot(c(cluster.ls$kmeans$r2, cluster.ls$pam$r2,
            cluster.ls$hclust$r2, cluster.ls$kkmeans$r2),
          main = "", col = "steelblue",
          names.arg = NULL, cex.names = 0.7, ylim = c(0, 100))
  if(i == 1){
    mtext('Explained variance', side = 3, line = 1, cex = 0.8)
  }
  if(i == length(ind)){
    axis(side = 1, labels = method.names, cex = 0.5, 
         at = seq(0.5, 4.5, length.out = 4), las = 2)
  }
}

dev.off()

### Save results ###
# GATHER RESULTS INTO TABLE
ls <- vector('list', length(ind))
for(c in 1:length(ind)) {
  k.no <- ind[c]
  load(file = sprintf('Analysis/Outputs/clusterls_%s.RData', k.no))
  stats <- lapply(cluster.ls, function(x) c(x$stats$dunn2, 
                                            x$r2, 
                                            x$stats$avg.silwidth, 
                                            x$min.sil))
  ls[[c]] <- do.call(rbind, stats)
}
results.df <- do.call(rbind, ls)
results.df <- cbind(rep(ind, each = length(method.names)), results.df)
methods <- rownames(results.df)
rownames(results.df) <- NULL
results.df <- cbind(methods, as.data.frame(results.df))

colnames(results.df) <- c('Method', "Cluster no.", "Dunn2", "Var. explained",
                          "Mean sil width", "Neg. sils (%)")
save(results.df, file = 'Analysis/Outputs/cluster_stats.RData')


###################
results.df[ ,'Method'] <- as.character(results.df[ ,'Method'])
