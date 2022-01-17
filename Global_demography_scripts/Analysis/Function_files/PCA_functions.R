########################################
### Functions used in the clustering script 
########################################

get.clst.mid <- function(ind.coords, cluster.mids){
  ssqs <- apply(ind.coords, 1, function(x)
    sum((x-cluster.mids)^2))
  return(which.min(ssqs))
}

# logistic function
logistic <- function(x, a = 0, b = 1, alpha = 0, beta = 1) {
  return(y <- alpha + beta * log((x - a) / (b - x)))
}


# make nice colour palettes
make.ramp <- function(tmpcol, levels = 10, transp = NA) {
  col.out <- colorRampPalette(tmpcol)(levels)
  if(!is.na(transp)) {
    col.out <- paste(col.out, transp, sep = "")
  }
  return(col.out)
}

make.transp <- function(tmpcol, transp = 50) {
  col.out <- paste(tmpcol, transp, sep = "")
  return(col.out)
}

# predict survival for size z
s_z <- function(z, params.s, time = 1) {
  K <- params.s[1]
  p <- params.s[2]
  r <- params.s[3]
  # Linear predictor parameters
  pred.y <-  (K / 
                (1 + exp(-(r * ((z - p) )))))  ^ (time)
  return(pred.y)
}

# make a survival vector
s_z_vec <- function(z, m.par, time = 1){
  
  thresh <- as.numeric(m.par['thresh'])
  small.ps <- as.numeric(m.par[1:3])
  big.ps <- as.numeric(m.par[c(1,4,5)])
  
  small <- z[which(z < thresh)]
  big <- z[which(z >= thresh)]
  
  p.small <- s_z(small, params.s = small.ps)
  p.big <- s_z(big, big.ps)
  
  surv <- c(p.small, p.big)
  
  return(surv)
}

# get sum of squares - each species/census combo to cluster center
get.ss <- function(cls.coord, center) {
  dist.mat <- dist(rbind(center, cls.coord))
  return(sum(as.matrix(dist.mat)[, 1] ^ 2))
}

### Make a biplot showing individual coordinates in PCA space 
### colour colded by cluster
### and variable loadings 
custom_biplot <- function(X, axes = c(1,2), cols, clusters, var.names, 
                          cex = 0.6, pch = 1, lwd = 1.2, 
                          cex.text = 1, 
                          plot.axes = TRUE,
                          plot.axes.labs = TRUE, 
                          xlims = NULL, ylims = NULL, 
                          exclude.outliers = TRUE,
                          lab.pos = 3, 
                          cex.axis.lab = 1,
                          arrow.colours = 'black'){
  
  
  # 1. Get the scaling factor for variable loadings
  var <- facto_summarize(X, element = "var", 
                         result = c("coord", 
                                    "contrib", "cos2"), axes = axes)
  colnames(var)[2:3] <- c("x", "y")
  pca.ind <- get_pca_ind(X)
  ind <- data.frame(pca.ind$coord[, axes, drop = FALSE])
  colnames(ind) <- c("x", "y")
  r <- min((max(ind[, "x"]) - 
              min(ind[, "x"])/(max(var[, "x"]) - 
                                 min(var[, "x"]))), (max(ind[, "y"]) - 
                                                       min(ind[, "y"])/(max(var[, "y"])
                                                                        - min(var[, "y"]))))
  # 2. 
  ### The individuals ###
  element <- 'ind'
  summary.res <- c("coord", "contrib", "cos2")
  df <- facto_summarize(X, element = element, axes = axes, 
                        result = summary.res)
  colnames(df)[2:3] <- c("x", "y")
  
  #####################
  # 3. The variables ###
  # (scaled so they fit)
  element <- 'var'
  geom <- c("arrow", "text")
  summary.res <- c("coord", "contrib", "cos2")
  df.var <- facto_summarize(X, element = element, axes = axes, 
                            result = summary.res)
  colnames(df.var)[2:3] <- c("x", "y")
  scale. = r * 0.3
  df.var[ ,c("x", "y")] <- df.var[ ,c("x", "y")] * scale.
  
  percents <- X$eig[,2][axes]
  
  # if remove outliers = take out coordinates of points outside 95th middle
  if(exclude.outliers == TRUE){
    out <- which(df[ ,'x'] < quantile(df[ ,'x'], prob = 0.005) | 
                   df[ ,'x'] > quantile(df[ ,'x'], prob = 0.995) |
                   df[ ,'y'] < quantile(df[ ,'y'], prob = 0.005) |
                   df[ ,'y'] > quantile(df[ ,'y'], prob = 0.995))
    if(length(out) > 0){
      df <- df[-out, ]
    clusters <- clusters[-out]
    pch <- pch[-out]
  }
  }
  xmn <- min(c(df.var[ ,'x'], df[ ,'x']))
  xmx <- max(c(df.var[ ,'x'], df[ ,'x']))
  ymn <- min(c(df.var[ ,'y'], df[ ,'y']))
  ymx <- max(c(df.var[ ,'y'], df[ ,'y']))
  
  xmin <- ifelse(xmn < 0, xmn*1.25, xmn*0.75)
  ymin <- ifelse(ymn < 0, ymn*1.25, ymn*0.75)
  xmax <- ifelse(xmx < 0, xmx*0.75, xmx*1.25)
  ymax <- ifelse(ymx < 0, ymx*0.75, ymx*1.25)
  
  if(is.null(xlims)){
    xlims <- c(xmin, xmax)
  }
  if(is.null(ylims)){
    ylims <- c(ymin, ymax)
  }
  
  if(plot.axes == TRUE){
  plot(df[ ,'x'], df[ ,'y'], col = cols[clusters], 
       xlab = '', ylab = '', cex.axis = 1.5, 
       las = 1, cex = cex, pch = pch, 
       xaxt = 'n', yaxt = 'n', #asp = 1, 
       xlim = xlims, ylim = ylims, bty = 'o')
    axis(side=1, at = NULL, lwd.ticks = 1.5,# padj=-2,
         cex.axis=1, tck=0)
    axis(side=2, at = NULL, lwd.ticks = 1.5, padj=1.5,
         cex.axis=1, tck = 0)
  }else{
    plot(df[ ,'x'], df[ ,'y'], col = cols[clusters], 
         xlab = '', ylab = '', cex.axis = 1.5, 
         las = 1, cex = cex, pch = pch, 
         xlim = xlims, ylim = ylims, 
         xaxt = 'n', yaxt = 'n', asp = 1, bty = 'o')
  }
  
  abline(v = 0, lty = 1)
  abline(h = 0, lty = 1)
  
  if(plot.axes.labs == TRUE){
    mtext(paste0('PC ', axes[1], ': ', round(percents[1], 1), ' %'), 
          side = 1, line = 2.5, cex = cex.axis.lab)
    mtext(paste0('PC ', axes[2], ': ', round(percents[2], 1), ' %'),
          side = 2, line = 2.5, cex = cex.axis.lab)
  }
  
  ## Arrows
  #apply(df.var, 1, function(x)
    #arrows(0,0, 
       #    as.numeric(x['x']), as.numeric(x['y']), 
      #     length = 0.1, angle = 30, lwd = lwd))
  
  arrow.mat <- as.matrix(df.var[ ,c('x', 'y')])
  mapply(function(x,acol){
    arrows(0,0, 
           as.numeric(x[1]), as.numeric(x[2]), 
           length = 0.1, angle = 30, lwd = 1.5, col = acol)
  }, 
  x = split(arrow.mat, row(arrow.mat)), 
  acol = arrow.colours)
  
  ## Labels
  text(as.numeric(df.var[,'x']), as.numeric(df.var[,'y']), 
       labels = var.names, pos = lab.pos,
       cex = cex.text)
}

##########################################################
mix.gamma <- function(x, gam.pars, q) {
  slow.prob <- dgamma(x, gam.pars[1], gam.pars[3])
  fast.prob <- dgamma(x, gam.pars[2], gam.pars[4])
  return(tot.prob <- (q*slow.prob) + ((1-q)*fast.prob))
}

cdf <- function(fn, min_x, max_x){
  # Returns a cumulative distribution function for a non-negative function over a given range.
  f <- function(x){
    y <- rep(NA, length(x))
    total_area_under_curve <- integrate(fn, min_x, max_x)$value
    apply_fn <- function(xval){integrate(fn, min_x, min(xval, max_x))$value}
    y <- sapply(x, FUN=apply_fn) / total_area_under_curve
    y <- cummax(y)
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


scatterutil.grid <- function (cgrid, col = 'lightgrey') 
{

  lty <- 1
  xaxp <- par("xaxp")
  ax <- (xaxp[2] - xaxp[1])/xaxp[3]
  yaxp <- par("yaxp")
  ay <- (yaxp[2] - yaxp[1])/yaxp[3]
  a <- min(ax, ay)
  v0 <- seq(xaxp[1], xaxp[2], by = a)
  h0 <- seq(yaxp[1], yaxp[2], by = a)
  abline(v = v0, col = col, lty = lty)
  abline(h = h0, col = col, lty = lty)
  if (cgrid <= 0) 
    return(invisible())
  cha <- paste(" d = ", a, " ", sep = "")
  cex0 <- par("cex") * cgrid
  xh <- strwidth(cha, cex = cex0)
  yh <- strheight(cha, cex = cex0) * 5/3
  x1 <- par("usr")[2]
  y1 <- par("usr")[4]
  #rect(x1 - xh, y1 - yh, x1 + xh, y1 + yh, col = "white", border = 0)
  #text(x1 - xh/2, y1 - yh/2, cha, cex = cex0)
}

scatterutil.star.hack <- function (x, y, z, cstar, coul = rep(1, length(x))) 
{
  z <- z/sum(z)
  x1 <- sum(x * z)
  y1 <- sum(y * z)
  for (i in which(z > 0)) {
    hx <- cstar * (x[i] - x1)
    hy <- cstar * (y[i] - y1)
    segments(x1, y1, x1 + hx, y1 + hy, col = coul[i])
  }
}

s.class.hack <- function (dfxy, fac, wt = rep(1, length(fac)), xax = 1, yax = 2, 
                          cstar = 1, cellipse = 1.5, axesell = TRUE, label = levels(fac), 
                          clabel = 1, cpoint = 1, pch = 20, col.points, col.lines, 
                          xlim = NULL, ylim = NULL, grid = TRUE, addaxes = TRUE,
                          origin = c(0,0), include.origin = TRUE, sub = "", csub = 1,
                          possub = "bottomleft", cex.point, 
                          cgrid = 0, add.plot = FALSE) 
{
  opar <- par(mar = par("mar"))
  par(mar = c(0.1, 0.1, 0.1, 0.1))
  on.exit(par(opar))
  dfxy <- data.frame(dfxy)
  dfdistri <- fac2disj(fac) * wt
  coul <- col
  w1 <- unlist(lapply(dfdistri, sum))
  dfdistri <- t(t(dfdistri)/w1)
  coox <- as.matrix(t(dfdistri)) %*% dfxy[, xax]
  cooy <- as.matrix(t(dfdistri)) %*% dfxy[, yax]
  
  coo <- scatterutil.base.hack(dfxy = dfxy, xax = xax, yax = yax, 
                          xlim = xlim, ylim = ylim, grid = grid, addaxes = addaxes, 
                          cgrid = cgrid, include.origin = include.origin, origin = origin, 
                          sub = sub, csub = csub, possub = possub, add.plot = add.plot)
  if (cpoint > 0) 
    for (i in 1:ncol(dfdistri)) {
      pch <- rep(pch, length = nrow(dfxy))
      points(coo$x[dfdistri[, i] > 0], coo$y[dfdistri[,i] > 0], 
             pch = pch[dfdistri[, i] > 0], cex = cex.point,
             col = col.points)
    }
  if (cstar > 0) 
    for (i in 1:ncol(dfdistri)) {
      scatterutil.star.hack(coo$x, coo$y, dfdistri[, i], cstar = cstar, 
                            col.lines)
    }
  if (cellipse > 0) 
    for (i in 1:ncol(dfdistri)) {
      scatterutil.ellipse(coo$x, coo$y, dfdistri[, i], 
                          cellipse = cellipse, axesell = axesell, coul = 'red')
    }
  if (clabel > 0) 
    scatterutil.eti(coox, cooy, label, clabel, coul = col)
  box()
  invisible(match.call())
}

scatterutil.base.hack <- function(dfxy, xax, yax, xlim, 
         ylim, grid, addaxes, cgrid, include.origin, 
         origin, sub, csub, possub, add.plot){
  df <- data.frame(dfxy)
  x <- df[, xax]
  y <- df[, yax]
  if (is.null(xlim)) {
    x1 <- x
    if (include.origin) 
      x1 <- c(x1, origin[1])
    x1 <- c(x1 - diff(range(x1)/10), x1 + diff(range(x1))/10)
    xlim <- range(x1)
  }
  if (is.null(ylim)) {
    y1 <- y
    if (include.origin) 
      y1 <- c(y1, origin[2])
    y1 <- c(y1 - diff(range(y1)/10), y1 + diff(range(y1))/10)
    ylim <- range(y1)
  }
  
  if (!add.plot) 
    plot.default(0, 0, type = "n", asp = 1,
                 xlab = "", ylab = "", 
                 xaxt = "n", yaxt = "n", xlim = xlim, ylim = ylim, 
                 xaxs = "i", yaxs = "i", frame.plot = FALSE)
  if (grid & !add.plot) 
    scatterutil.grid(cgrid)
  if (addaxes & !add.plot) 
    abline(h = 0, v = 0, lty = 1)
  if (csub > 0) 
    scatterutil.sub(sub, csub, possub)
  return(list(x = x, y = y))
}