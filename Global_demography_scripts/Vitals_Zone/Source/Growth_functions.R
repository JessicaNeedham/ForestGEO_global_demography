################################################################
### GROWTH FUNCTIONS ###
################################################################
mix.gamma <- function(x, gam.pars, q) {
  slow.prob <- dgamma(x, gam.pars[1], gam.pars[3])
  fast.prob <- dgamma(x, gam.pars[2], gam.pars[4])
  return(tot.prob <- q * slow.prob + (1-q) * fast.prob)
}

cdf <- function(fn, min_x, max_x){
  # Returns a cumulative distribution function for a non-negative function over a given range.
  f <- function(x){
    y <- rep(NA, length(x))
    total_area_under_curve <- integrate(fn, min_x, max_x)$value
    apply_fn <- function(xval){integrate(fn, min_x, min(xval, max_x))$value}
    y <- sapply(x, 
                FUN=apply_fn) / total_area_under_curve
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


# make nice colour palettes
make.ramp <- function(tmpcol, levels = 10, transp = NA) {
  col.out <- colorRampPalette(tmpcol)(levels)
  if(!is.na(transp)) {
    col.out <- paste(col.out, transp, sep = "")
  }
  return(col.out)
}
