get.anova.letters <- function(trait, cluster) {
  out <- posthoc.kruskal.dunn.test(trait ~ cluster, p.adjust = "bonf")
  out.p <- get.pvalues(out)
  out.mcV <- multcompLetters(out.p, threshold = 0.05)
  Rij <- rank(trait)
  Rj.mean <- tapply(Rij, cluster, mean)
  return(out.mcV$Letters)
}

get.row.wt <- function(census) {
  row.wt <- rep((1 / max(census)), length = length(census))
  return(as.numeric(row.wt))
}

ss <- function(x) sum((x - mean(x, na.rm = TRUE)) ^ 2, na.rm = TRUE)

get.par.ss <- function(time.vars) {
  scale.pars <- apply(time.vars[, -c(1:2)], MAR = 2, scale)
  ss.out <- matrix(NA, ncol = 3, nrow = length(unique(time.vars$Species)))
  for(i in 1:dim(scale.pars)[2]) {
    ss.out[, i] <- tapply(scale.pars[, i], time.vars$Species, FUN = ss)
  }
  return(rep(apply(ss.out, 1, sum), each = max(time.vars$Census_int)))
}

get.dr1 <- function(params, dbh = 10) {
  return((pred.s(dbh + 0.01, params) -
            pred.s(dbh, params)) / 0.01)
}

get.r1.int <- function(params, dbh = 10) {
  return(pred.s(dbh, params))
}

get.min.sil <- function(sil.data, percent = TRUE) {
  neg.no <- length(which(sil.data[ ,3] < 0))
  if (percent == TRUE) {
    return (neg.no / dim(sil.data)[1])
  } else {
    return(neg.no)
  }
}

get.centers <- function(cluster.hclust) {
  clus.list <- split(as.data.frame(ind.coord), cluster.hclust$cluster)
  return(t(sapply(clus.list, colMeans)))
}

get.ss <- function(cls.coord, center) {
  dist.mat <- dist(rbind(center, cls.coord))
  return(sum(as.matrix(dist.mat)[, 1] ^ 2))
}

get.withinss <- function(ind.coord, clust.id, centers) {
  # uses coordinates of data, cluster id, and centers of clusters
  # to estimate within-cluster sum of squares.
  
  # # Input includes:
  # ind.coord -- the matrix of PCA coordinates for every
  # # species / census.
  # clust.id -- the cluster assignment
  # centers -- the coordinates of the cluster centers
  # Output is vector of silhouette values.
  
  clus.list  <- split(as.data.frame(ind.coord), clust.id)
  within.ss <- c()
  for(i in 1:length(clus.list)) {
    within.ss[i] <- get.ss(clus.list[[i]], centers[i, ])
  }
  return(sum(within.ss))
}

get.explained.var <- function(tot.ss, within.ss) {
  return(explained.var <- 100 * (tot.ss - within.ss) / tot.ss)
}

make.silplot <- function(silinfo, method ="", make.legend = FALSE) {
  sil.sort <- sortSilhouette(silinfo)
  cols <- brewer.pal(max(silinfo[, 2]), "Set3")
  
  barplot(rev(sil.sort[, "sil_width"]), col = cols[rev(sil.sort[, "cluster"])],
          border =  NA, horiz = TRUE, names.arg = NA, space = 0,
          main = '')
  abline(v = mean(sil.sort[, 3]), lty = 2)
  if(make.legend == TRUE) {
    legend("left", legend = paste("Cluster", seq(length(cols)), sep = " "),
           col = cols, lwd = 6, cex = 0.4)
  }
}