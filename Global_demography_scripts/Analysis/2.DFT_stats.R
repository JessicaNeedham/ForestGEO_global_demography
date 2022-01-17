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
source('Analysis/Function_files/DFT_clustering_functions.R')

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
       file = sprintf('PCA_Zone/Output/clusterls_%s.RData', k.no))
  print(i)
}

#### RESTART HERE ### 
jpeg(file = "Paper_Figures_v2/silhouette_clusters.jpeg", bg = 'white')
par(mfrow = c(6,4), mar = c(1,1,1,1), oma = c(2,3,3,0))
#for(i in 1:length(ind)){
for(i in 1:length(ind)){ 
  
  k.no <- ind[i]
  load(file = sprintf('PCA_Zone/Output/clusterls_%s.RData', k.no))
  
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

jpeg(file = 'Paper_Figures_v2/cluster_method_barplots.jpeg', bg = 'white')

par(mfrow = c(length(ind), 4), oma = c(3.5,3,3,0), 
    mar = c(1,2,1,1))

for(i in 1:length(ind)){
  
  k.no <- ind[i]
  load(file = sprintf('PCA_Zone/Output/clusterls_%s.RData', k.no))
  
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
  load(file = sprintf('PCA_Zone/Output/clusterls_%s.RData', k.no))
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
save(results.df, file = 'PCA_Zone/Output/cluster_stats.RData')


###################
results.df[ ,'Method'] <- as.character(results.df[ ,'Method'])

###########################################################################
# plot vital rates with each cluster method 
source('PCA_Zone/Source/PCA_functions.R')
n.clust <- 8
cols <-  c(brewer.pal(8, 'Dark2')[c(6,1,2,3,4,5,7,8)], 
           brewer.pal((n.clust - 8), 'Paired'))
cols.clst.light <- make.transp(cols)
cols.clst.vlight <- make.transp(cols, 25)
cols.light <- make.transp(cols) 

for(i in ind){
  
  load(file = sprintf('PCA_Zone/Output/clusterls_%s.RData', i))
  n.clust <- i
  
  for(j in 1:length(method.names)){
    method <- method.names[j]
    
    tmp <- eval(as.name(paste0('cluster.', method)))
    pc.tmp <- pca.matrix
    if(method == 'kkmeans'){
      pc.tmp[ ,'clust'] <- tmp@.Data
    }else{
      pc.tmp[ ,'clust'] <- tmp$cluster
    }
    ord <- as.numeric(names(sort(tapply(pc.tmp[ ,'ex.fast'], pc.tmp[ ,'clust'], median))))
    pc.tmp[ ,'clust'] <- mapvalues(pc.tmp[ ,'clust'], 
                                   from = ord, to = seq(n.clust))
    cl.lst <- split(pc.tmp, pc.tmp[ ,'clust'])
    
    # for each cluster plot the survival curve
    cl.surv.quants <- vector('list', n.clust)
    z <- seq(2000)
    # growth distributions
    xx <- seq(0.1, 50, 0.05)
    cl.gr.quants.slow <- cl.gr.quants.fast <- vector('list', n.clust)
    
    for(cl in 1:n.clust){
      #  cl <- 1 
      sp.params <- cl.lst[[cl]]
      s.params <- sp.params[ ,c('K', 'p1', 'r1', 'p2', 'r2', 'thresh')]
      g.params <- sp.params[ ,c('alpha1', 'alpha2', 'beta1', 'beta2', 'incr.thresh')]
      pars.samp <- cbind(s.params, g.params)
      
      surv.par.samp <- pars.samp[ ,colnames(s.params)]
      grow.par.samp <- pars.samp[ ,colnames(g.params)]
      
      # now with each set of these parameters get survival and growth - then 
      # take quantiles. 
      s.curves <- t(apply(surv.par.samp, 1, s_z_vec, z = z))
      cl.surv.quants[[cl]] <- apply(s.curves, 2, quantile, 
                                    prob = c(0.025, 0.25, 0.5,
                                             0.75, 0.975))
      
      # GROWTH 
      slow.ps <- grow.par.samp[ ,c('alpha1', 'beta1')]
      fast.ps <- grow.par.samp[ ,c('alpha2', 'beta2')]
      
      slow.dists <- t(apply(slow.ps, 1, function(ps) 
        dgamma(xx, shape = ps[1], rate = ps[2])))
      fast.dists <- t(apply(fast.ps, 1, function(ps)
        dgamma(xx, shape = ps[1], rate = ps[2])))
      
      cl.gr.quants.slow[[cl]] <- apply(slow.dists, 2, quantile, 
                                       prob = c(0.025, 0.25, 0.5, 0.75, 0.975))
      cl.gr.quants.fast[[cl]] <- apply(fast.dists, 2, quantile, 
                                       prob = c(0.025, 0.25, 0.5, 0.75, 0.975))
    }
    
    samp.clst <- as.numeric(tapply(pc.tmp[ ,'Species'], pc.tmp[ ,'clust'], 
                                   function(x) length(unique(x))))
    
    pdf(file = sprintf('PCA_Zone/Figures/cluster_vitals_%s_%d.pdf', 
                       method, i),  bg = 'white')
    par(mfrow = c(n.clust, 2), mar = c(1,4,1,1), 
        oma = c(4,1,1.2, 0))
    
    for(cl in 1:n.clust){
      sp.params <- cl.lst[[cl]]
      p2mu <- median(sp.params[ ,'p2'])
      if(cl == n.clust){
        plot(z, cl.surv.quants[[cl]][3, ], 
             #meds[[cl]], 
             type = 'l', bty = 'l', 
             col = cols[cl], lwd = 2, ylim = c(0,1),
             ylab = '', xlab = '', las = 1, cex.axis = 1)
        mtext('DBH (mm)', side = 1, line = 2, cex = 1.2, outer = TRUE, 
              adj = 0.25)
      }else{
        plot(z, cl.surv.quants[[cl]][3, ],
             #meds[[cl]],
             type = 'l', bty = 'l', 
             col = cols[cl], lwd = 2, xaxt = 'n',
             ylab = '', xlab = '', las = 1, cex.axis = 1,
             ylim = c(0, 1))
      }
      
      polygon(c(z, rev(z)), c(cl.surv.quants[[cl]][2, ], 
                              rev(cl.surv.quants[[cl]][4, ])), 
              col = cols.clst.light[cl], border = NA)
      polygon(c(z, rev(z)), c(cl.surv.quants[[cl]][1, ], 
                              rev(cl.surv.quants[[cl]][5, ])), 
              col = cols.clst.vlight[cl], border = NA)
      #text(1400, 0.85, sprintf('N: %d', samp.clst[cl]), cex = 1.2)
      if(cl == round(i/2)){
        mtext('Survival probability', side = 2, line = 3, 
              cex = 1.2, adj = 0.6)
      }
      
      if(cl == n.clust){
        plot(xx, cl.gr.quants.slow[[cl]][3, ], type = 'l', bty = 'l', 
             col = cols[[cl]], lwd = 2, las = 1, 
             ylab = '', xlab = '', yaxt = 'n')
        mtext('DBH increment (mm)', side = 1, line = 2, cex = 1.2, 
              outer = TRUE, adj = 0.9)
      }else{
        plot(xx, cl.gr.quants.slow[[cl]][3, ], type = 'l', bty = 'l', 
             col = cols[[cl]], lwd = 2, las = 1, 
             ylab = '', xlab = '', xaxt = 'n', yaxt = 'n')
      }
      points(xx, cl.gr.quants.fast[[cl]][3, ], type = 'l', 
             lwd = 2, col = cols[[cl]])
      text(20, max(cl.gr.quants.slow[[cl]][3,]*0.88),
           sprintf('N: %d', samp.clst[cl]), cex = 1.2)
      
      polygon(c(xx, rev(xx)), 
              c(cl.gr.quants.slow[[cl]][2, ], 
                rev(cl.gr.quants.slow[[cl]][4, ])), 
              col = cols.clst.light[cl], border = NA)
      polygon(c(xx, rev(xx)), 
              c(cl.gr.quants.fast[[cl]][2, ], 
                rev(cl.gr.quants.fast[[cl]][4, ])), 
              col = cols.clst.light[cl], border = NA)
      polygon(c(xx, rev(xx)), 
              c(cl.gr.quants.slow[[cl]][1, ], 
                rev(cl.gr.quants.slow[[cl]][5, ])), 
              col = cols.clst.vlight[cl], border = NA)
      polygon(c(xx, rev(xx)), 
              c(cl.gr.quants.fast[[cl]][1, ], 
                rev(cl.gr.quants.fast[[cl]][5, ])), 
              col = cols.clst.vlight[cl], border = NA)
      if(cl == round(i/2)){
        mtext('Growth density', side = 2, line = 2, 
              cex = 1.2, adj = 0.6)
      }
    }
    mtext('Survival', side = 3, line = -1, outer = TRUE, 
          cex = 1.2, adj = 0.25)
    mtext('Growth', side = 3, line = -1, outer = TRUE, 
          cex = 1.2, adj = 0.8)
    mtext(paste(method, ' ', i), side = 3, line = -1, outer = TRUE, cex = 1.4)
    dev.off()
  }
}


############################################################################
rm(list = ls())
### LIBRARIES ###
library(FactoMineR)
library(factoextra)
library(RColorBrewer)
library(kernlab)
library(fpc)
library(cluster)
library(foreach)
library(doParallel)
library(xtable)
library(plyr)
library(stringi)
library(plotrix)
library(tmvtnorm)
library(corpcor)
library(kernlab)
library(fpc)
library(viridis)
library(plot3D)
library(caTools)
library(ade4)
library(grid)
library(gridBase)
library(ggmap)

source('PCA_Zone/Source/PCA_functions.R')
source('PCA_Zone/Source/Ade4_functions.R')

geog <- 'global'
load(file = sprintf('PCA_Zone/Output/%s_pca_mat.RData', geog))
load(file = sprintf('PCA_Zone/Output/PCA_%s.RData', geog))

site.names <- as.character(unique(pca.matrix[ ,'Site']))
site.names <- ifelse(nchar(site.names) <= 4, toupper(site.names), 
                     stri_trans_totitle(site.names))
site.names[which(site.names == 'Khaochong')] <- 'Khao Chong'
site.names[which(site.names == 'Laplanada')] <- 'La Planada'
site.names[which(site.names == 'Windriver')] <- 'Wind River'
site.names[which(site.names == 'XTBG')] <- 'XSBN'
var.names <- c('MS', 'JS', 'ST', 
               '95G', '5G')
var.full.names <- c('Maximum survival', 'Juvenile survival', 
                    'Stature', '95% slowest growth', '5% fastest growth')
arrow.cols <- c('red', 'red', 'darkorange', 'darkblue', 'darkblue')

n.clust <- 8
cols <- c(brewer.pal(7, 'Dark2')[c(6,1,2,3,4,5,7)], brewer.pal(3, 'Set1')[2])
cols.light <- make.transp(cols)
cols.v.light <- make.transp(cols, rep(30, n.clust))

cols.points <- cols.v.light
cols.lines <- make.transp(cols, c(20,20,20,20,20,20,20,20))

ind.coords <- res$ind$coord
var.coords <- res$var$coord


load(file = sprintf('PCA_Zone/Output/clusterls_%s.RData', n.clust))
method.names <- c("kmeans", "pam", "hclust", "kkmeans")

pdf(file="Paper_Figures/Cluster_methods_urchin_plots.pdf", 
    bg = 'white', width = 12, height = 10)
par(mfrow = c(2,2), oma = c(0,0,0,0), mar = c(5,4,4,2))
  for(j in 1:length(method.names)){
  #j <- 1 
   method <- method.names[j]
    
    tmp <- eval(as.name(paste0('cluster.', method)))
    pc.tmp <- pca.matrix
    if(method == 'kkmeans'){
      pc.tmp[ ,'clust'] <- tmp@.Data
    }else{
      pc.tmp[ ,'clust'] <- tmp$cluster
    }
    ord <- as.numeric(names(sort(tapply(pc.tmp[ ,'ex.fast'], pc.tmp[ ,'clust'], median))))
    pc.tmp[ ,'clust'] <- mapvalues(pc.tmp[ ,'clust'], 
                                   from = ord, to = seq(n.clust))
    clusters <- pc.tmp[ ,'clust']
    Ns <- pca.matrix[ ,'N']
    
    custom_biplot(res, axes = c(1,2), cols = cols.points, clusters = clusters, 
                  var.names = var.names, pch = 16, lwd = 1.8, cex.text = 1.8, 
                  plot.axes = TRUE, plot.axes.labs = TRUE, exclude.outliers = FALSE, 
                  cex = sqrt(Ns)*0.02, lab.pos = c(3,3,3,4,4), arrow.colours = arrow.cols, 
                  xlim = c(-3, 10), ylim = c(-4, 4))
    scatterutil.grid(1)
    abline(h = 0, v = 0, lty = 1)
    mtext(method, side = 3, line = 0.5, outer = FALSE, cex = 2)
    
    dfxy <- data.frame(ind.coords[ ,1:2])
    dfdistri <- fac2disj(factor(clusters)) * rep(1, length(clusters))
    coul <- col
    w1 <- unlist(lapply(dfdistri, sum))
    dfdistri <- t(t(dfdistri)/w1)
    coox <- as.matrix(t(dfdistri)) %*% dfxy[,1]
    cooy <- as.matrix(t(dfdistri)) %*% dfxy[,2]
    
    for (i in 1:ncol(dfdistri)) {
      scatterutil.star(dfxy[ ,1], dfxy[ ,2], dfdistri[, i], cstar = 1, 
                       cols.lines[i])
    }
    #scatterutil.eti(coox, cooy, label = seq(length(cols)), clabel = 0.6, 
    #               coul = cols)
    text(coox, cooy, label = seq(length(cols)), cex = 1)
    if(j == 1){
    legend('bottomright', cex = 1, col = arrow.cols, lwd = 1.6,
           legend = paste(var.names, var.full.names, sep = ': '))
    }
  }
    
dev.off()
    
    
        
