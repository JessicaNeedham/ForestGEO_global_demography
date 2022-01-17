###########################################################
### IBM FUNCTIONS ###
############################################################
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

transition.trees <- function(z1, f.m, f.c){
  p.fast <- z1*f.m + f.c
  fast <- rbinom(length(z1), 1, prob = p.fast)
  return(fast+1)
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
  
  surv.thresh <- m.par$surv.thresh
  small.params <- m.par$surv.small
  big.params <- m.par$surv.big
  
  small <- z[which(z < surv.thresh)]
  big <- z[which(z >= surv.thresh)]
  
  p.small <- s_z(small, small.params)
  p.big <- s_z(big, big.params)
  
  surv <- c(p.small, p.big)
  
  return(surv)
}

t_col <- function(color, percent = 50, name = NULL) {
  #	  color = color name
  #	percent = % transparency
  #	   name = an optional name for the color
  ## Get RGB values for named color
  rgb.val <- col2rgb(color)
  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100-percent)*255/100,
               names = name)
  ## Save the color
  invisible(t.col)
}

# get passage times and longevities
get.passage.time.real <- function(recruits, threshold, times){
  
  # find which year if any the recruits made it to the threshold
  pts <- apply(recruits, 1, function(x) which(x >= threshold)[1])
  ind <- which(!is.na(pts))  # just the ones who made it
  col1 <- recruits[ind, 1] # starting size
  
  # get time to the threshold
  col2 <-  mapply(function(x,y) x[y],
                  x = split(times[ind, ], row(times[ind, ])),  
                  y = pts[ind])
  # size at size above threshold
  col3 <- apply(recruits, 1, function(x) x[x>=threshold][1])
  col3 <- col3[ind]
  out <- cbind(col1, col2, col3)
  return(out)
}

get.le.real <- function(recruits, recruit.times){
  
  # find recruits that died
  died <- recruits[which(is.na(recruits[ ,ncol(recruits)])), ]
  died.times <- recruit.times[which(is.na(recruits[ ,ncol(recruits)])), ]
  # find age at death
  les <- as.numeric(mapply(function(x,y){z <- which(is.na(x))[1];
  return(y[z])}, 
  x = split(died, row(died)), 
  y = split(died.times, row(died.times))))
  return(les)
}


get.passage.time.ibm <- function(recruits, threshold){
  
  # find which year if any the recruits made it to the threshold
  pts <- apply(recruits, 1, function(x) which(x >= threshold)[1])
  pts <- pts[!is.na(pts)]
  return(pts)
}

get.lifespan.ibm <- function(recruits){
  
  # find which year if any the recruits died
  yr.died <- apply(recruits, 1, function(x) which(is.na(x))[1])
  yr.died <- yr.died[!is.na(yr.died)]
  return(yr.died)
}

# make nice colour palettes
make.ramp <- function(tmpcol, levels = 10, transp = NA) {
  col.out <- colorRampPalette(tmpcol)(levels)
  if(!is.na(transp)) {
    col.out <- paste(col.out, transp, sep = "")
  }
  return(col.out)
}


# make transparent colours
make.transp <- function(tmpcol, transp = 50) {
  col.out <- paste(tmpcol, transp, sep = "")
  return(col.out)
}

