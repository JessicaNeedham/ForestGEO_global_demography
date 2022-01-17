fac2disj<- function(fac, drop = FALSE) {
  ## Returns the disjunctive table corrseponding to a factor
  n <- length(fac)
  fac <- as.factor(fac)
  if(drop)
    fac <- factor(fac)
  x <- matrix(0, n, nlevels(fac))
  x[(1:n) + n * (unclass(fac) - 1)] <- 1
  dimnames(x) <- list(names(fac), as.character(levels(fac)))
  return(data.frame(x, check.names = FALSE))
}


s.class <- function (dfxy, fac, wt = rep(1, length(fac)), xax = 1, yax = 2, 
          cstar = 1, cellipse = 1.5, axesell = TRUE, label = levels(fac), 
          clabel = 1, cpoint = 1, pch = 20, col = rep(1, length(levels(fac))), 
          xlim = NULL, ylim = NULL, grid = TRUE, addaxes = TRUE,
          origin = c(0,0), include.origin = TRUE, sub = "", csub = 1, possub = "bottomleft", 
          cgrid = 1, add.plot = FALSE) 
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
  
  coo <- scatterutil.base(dfxy = dfxy, xax = xax, yax = yax, 
                          xlim = xlim, ylim = ylim, grid = grid, addaxes = addaxes, 
                          cgrid = cgrid, include.origin = include.origin, origin = origin, 
                          sub = sub, csub = csub, possub = possub, add.plot = add.plot)
  if (cpoint > 0) 
    for (i in 1:ncol(dfdistri)) {
      pch <- rep(pch, length = nrow(dfxy))
      points(coo$x[dfdistri[, i] > 0], coo$y[dfdistri[,i] > 0], 
             pch = pch[dfdistri[, i] > 0], cex = par("cex") * 
               cpoint, col = coul[i])
    }
  if (cstar > 0) 
    for (i in 1:ncol(dfdistri)) {
      scatterutil.star(coo$x, coo$y, dfdistri[, i], cstar = cstar, 
                       coul[i])
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

scatterutil.base <- function(dfxy, xax, yax, xlim, 
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


scatterutil.eti <- function (x, y, label, clabel, boxes = TRUE, coul = rep(1, length(x)), 
          horizontal = TRUE, bg = "white") 
{
  if (length(label) == 0) 
    return(invisible())
  if (is.null(label)) 
    return(invisible())
  if (any(label == "")) 
    return(invisible())
  cex0 <- par("cex") * clabel
  for (i in 1:(length(x))) {
    cha <- as.character(label[i])
    cha <- paste(" ", cha, " ", sep = "")
    x1 <- x[i]
    y1 <- y[i]
    xh <- strwidth(cha, cex = cex0)
    yh <- strheight(cha, cex = cex0) * 5/3
    if (!horizontal) {
      tmp <- scatterutil.convrot90(xh, yh)
      xh <- tmp[1]
      yh <- tmp[2]
    }
    if (boxes) {
      rect(x1 - xh/2, y1 - yh/2, x1 + xh/2, y1 + yh/2, 
           col = bg, border = coul[i])
    }
    if (horizontal) {
      text(x1, y1, cha, cex = cex0, col = coul[i])
    }
    else {
      text(x1, y1, cha, cex = cex0, col = coul[i], srt = 90)
    }
  }
}
