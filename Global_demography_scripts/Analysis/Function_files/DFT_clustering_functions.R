
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

pred.s <- function(dbh, pars) {
  pars <- as.numeric(pars)
  p1  <- pars[1]
  r1   <- pars[2]
  K    <- pars[3]
  p2  <- pars[4]
  r2   <- pars[5]
  thresh <- pars[6]
  prob.surv <- ifelse(dbh < thresh,
    K / (1 + exp(-r1 * (dbh - p1))),
    K / (1 + exp(-r2 * (dbh - p2))))
  return(prob.surv)
}

make.poly.vecs <- function(surv.par.samp, dbh = seq(10, 1500),
  CI = c(0.05, 0.5, 0.95)) {
  surv.mat <- apply(surv.par.samp, MAR = 1, pred.s, dbh = dbh)
  CI.surv <- apply(surv.mat, MARGIN = 1, FUN =  quantile, probs = CI)
  return(CI.surv)
}

make.grow.vecs <- function(grow.par.samp, incr.vec = seq(0, 15, length = 100),
  CI = c(0.05, 0.5, 0.95)) {
  CI.grow <- vector("list", 2)
  grow.mat <- apply(grow.par.samp, MAR = 1, pred.g, incr.vec = incr.vec)
  incr.mat.1 <- matrix(c(lapply(grow.mat, "[", c(1)), recursive = TRUE),
    ncol = length(incr.vec), byrow = TRUE)
  incr.mat.2 <- matrix(c(lapply(grow.mat, "[", c(2)), recursive = TRUE),
    ncol = length(incr.vec), byrow = TRUE)
  CI.grow[[1]] <- apply(incr.mat.1, MARGIN = 2, FUN =  quantile, probs = CI,
    na.rm = TRUE)
  CI.grow[[2]] <- apply(incr.mat.2, MARGIN = 2, FUN =  quantile, probs = CI,
    na.rm = TRUE)
  return(CI.grow)
}

pred.g <- function(incr.vec, grow.pars) {
  a1 <- grow.pars[1]
  a2 <- grow.pars[2]
  b1 <- grow.pars[3]
  b2 <- grow.pars[4]
  y1 <- dgamma(incr.vec, shape = a1, rate = b1)
  y2 <- dgamma(incr.vec, shape = a2, rate = b2)
  y1 <- ifelse(y1 == Inf, 0, y1)
  y2 <- ifelse(y2 == Inf, 0, y2)
  return(list(y1, y2))
}


mix.gamma <- function(x, gam.pars) {
  slow.prob <- dgamma(x, gam.pars[1], gam.pars[3])
  fast.prob <- dgamma(x, gam.pars[2], gam.pars[4])
  return(tot.prob <- 0.95 * slow.prob + 0.05 * fast.prob)
}

spline_function <- function(x, y){
   # Returns a spline function (made up of cubic polynomials) that interpolates the
   # points given by the x and y vectors. The function has range [min(x), max(x)].
   fn <- splinefun(x, y, method="natural")
   min_x <- min(x)
   max_x <- max(x)
   f <- function(x){
      y <- fn(x)
      y[x < min_x | max_x < x] <- NA
      return(y)
   }
   return(f)
}

cdf <- function(fn, min_x, max_x){
   # Returns a cumulative distribution function for a non-negative function over a given range.
   f <- function(x){
      y <- rep(NA, length(x))
      total_area_under_curve <- integrate(fn, min_x, max_x)$value
      apply_fn <- function(xval){integrate(fn, min_x, min(xval, max_x))$value}
      y[min_x <= x & x <= max_x] <-
      sapply(x[min_x <= x & x <= max_x], FUN=apply_fn) /
      total_area_under_curve
      y[x < min_x] <- 0
      y[x > max_x] <- 1
      return(y)
   }
   return(f)
}

inverse <- function(fn, min_x, max_x){
   # Returns the inverse of a function for a given range.
   # E.g. inverse(sin, 0, pi/2)(sin(pi/4)) equals pi/4 because 0 <= pi/4 <= pi/2
   fn_inv <- function(y){
      uniroot((function(x){fn(x) - y}), lower=min_x, upper=max_x)[1]$root
   }
   return(Vectorize(fn_inv))
}
