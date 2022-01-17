####################################################################
### Fits growth models ###
### using STAN MCMC ###
####################################################################
rm(list = ls())
### LIBRARIES ###
library(foreach)
library(doParallel)
library(rstan)
library(RColorBrewer)
library(xtable)
library(MASS)

### FUNCTIONS ###
source('Vitals_Zone/Source/Growth_functions.R')
rstan_options(auto_write = TRUE)
### DATA ###
## load the data
sites <- c('amacayacu','bci', 'changbaishan', 'fushan', 'hkk', 
           'ituri', 'khaochong', 'korup', 'lambir', 'laplanada', 
           'luquillo', 'pasoh', 'palanan',
           'scbi', 'serc', 'windriver', 
           'wytham', 'yasuni', 'xtbg')

n.sites <- length(sites)

for(ii in 1:n.sites){
 #ii <- 3
  site <- sites[ii]
  data.ls <- readRDS(file = sprintf('Data_Zone/Output/%s.RData', site))
  
  # use try since not all sites have all components 
  try(data.mat <- data.ls$data.mat, silent = TRUE)
  try(sp.names <- data.ls$sp.names, silent = TRUE)
  try(sp.latins <- data.ls$sp.latins, silent = TRUE)
  try(years <- data.ls$years, silent = TRUE)
  sp.latins <- as.character(sp.latins[ ,1])
  n.census <- length(years)
  n.intervals <- n.census - 1  
  n.sp <- length(sp.names)  # number of species
  
  sp.list <- lapply(split(data.mat, data.mat[,'sp.id']),
                    matrix, ncol = ncol(data.mat), 
                    dimnames = list(NULL, colnames(data.mat)))
  
  
  # dummy data so we can compile model outside of the 
  # foreach loop
  N <- 100
  incr <- rnorm(N, 1, 0.2)
  incr.groups <- rbinom(N, 1, 0.05) 
  
  incr.data <- list(N = length(incr), 
                    incr = incr,
                    groups = incr.groups)
  q <- 0.95
  
  # compile the model - choose appropriate stan script for the no. distributions
  # this may take about 40-60 secs
  growth.model <- stan(file = 'Vitals_Zone/Source/gamma_growth.stan',
                       data = incr.data, chains = 0, verbose = FALSE)
  
  #n.cores <- detectCores()-1  # for quantlab
  n.cores <- 5
  registerDoParallel(n.cores)
  n.iter <- 2500
  
  foreach(sp = 1:n.sp, 
          .errorhandling = 'pass', 
          .packages = 'rstan')%dopar%{
            
            data.mat <- sp.list[[sp]]
            sp.name <- sp.names[sp]
            sizes <- data.mat[ ,grep('dbh', colnames(data.mat))]
            
            g.params <- vector('list', 4)
            g.params[[1]] <- rep(NA, 4)  # parameters
            g.params[[2]] <- NULL  # HMC chains
            g.params[[3]] <- NULL 
            g.params[[4]] <- NULL
            
            # for each census get growth mode
            size <- data.mat[ ,grep('dbh.1', colnames(data.mat))]  # size start of census
            incr <- data.mat[ ,grep('incr', colnames(data.mat))]  # incr
            time <- data.mat[ ,grep('time', colnames(data.mat))]  # time
            treeID <- data.mat[ ,grep('treeID', colnames(data.mat))]
            # for each tree that lived get the largest stem that survived
            
            # remove stems dead/not recruited in either of the focal census years
            alive <- which(!is.na(incr))
            size <- size[alive]
            incr <- incr[alive]
            time <- time[alive]
            treeID <- treeID[alive]
            
            # for full not stem tables
            if(length(treeID) > 0){
              
              # index of largest stem per tree
              max.size <- tapply(size, treeID, which.max)
              # corresponding to rows... 
              size <- tapply(size, treeID, max)  # max size of each tree
              incr <- mapply(function(x,y) x[y], 
                             x = split(incr, treeID), 
                             y = max.size)
              time <- mapply(function(x,y) x[y], 
                             x = split(time, treeID), 
                             y = max.size)
            }
            
            time[is.na(time)] <- sample(time[!is.na(time)], 
                                        length(which(is.na(time))))
            
            # if there are less than 100 individuals in a census interval then 
            # put NAs and move on
            # this is now redundent since we only use species with 200 individuals
            if(length(size) < 100){
              stop()
            }
            
            # MAKE INCREMENT ANNUAL #
            incr <- incr/time
            q <- 0.95
            g.params[[3]] <- quantile(incr, prob = q)
            ############################################################################
            ### MCMC ###
            G <- 2
            dists <- rep('gamma', G)
            
            # INCREMENT GROUP
            quant.incr <- g.params[[3]]
            # divide the data into increment groups based on these increment thresholds
            incr.groups <- cut(incr, 
                               breaks = c(0, quant.incr, Inf),
                               labels = FALSE)
            
            # put the data into a list object
            incr.data <- list(N = length(incr), 
                              incr = incr,
                              groups = incr.groups)
            
            sp.1 <- as.numeric(fitdistr(incr[incr.groups == 1],
                                        'gamma')$estimate)
            sp.2 <- as.numeric(fitdistr(incr[incr.groups == 2],
                                        'gamma')$estimate)
            
            start.params <- list(alpha1 = sp.1[1], alpha2 = sp.2[1], 
                                 beta1 = sp.1[2], beta2 = sp.2[2])
            
            # parameter sampling
            growth.fit <- stan(fit = growth.model, data = incr.data, verbose = FALSE,
                               iter = n.iter, chains = 3, 
                               init = list(start.params, start.params, start.params),
                               control = list(adapt_delta = 0.99, 
                                              max_treedepth = 13))
            
            # check the MCMC results
            summary(growth.fit)
            
            # extract the parameters
            growth.mcmc <- unlist(lapply(extract(growth.fit, permuted = TRUE), median))
            g.params[[1]] <- growth.mcmc[1:4]
            g.params[[2]] <- do.call(cbind, extract(growth.fit))
            g.params[[4]] <- length(incr)
            
            names(g.params[[1]]) <- c('alpha1', 'alpha2', 'beta1', 'beta2')
            
            save(g.params, file = sprintf('Vitals_Zone/Output/%s/%s_growthps_%s.RData', 
                                          site, site, sp.name))
          }
  
  #############################
  # Make parameter tables #
  ##############################
  # Latex version for appendices and a version to 
  # work with in PCA zone #
  growth.tab <- vector('list', n.sp)
  intervals <- na.omit(expand.grid(years, years)[(seq(5)*(n.census+1)), ])
  intervals <- matrix(paste(intervals[,1], intervals[,2], sep = '-'), 
                      nrow = n.intervals)
  intervals <- intervals[n.intervals, 1]
  xx <- seq(0.1, 50, 0.1)
  xlim <- c(min(xx), max(xx))
  q <- 0.95
  
  for(sp in 1:n.sp){
    
    sp.name <- sp.names[sp]
    sp.latin <- sp.latins[sp]
    
    trycheck <- try(load(file = sprintf('Vitals_Zone/Output/%s/%s_growthps_%s.RData', 
                                        site,   site,  sp.name)), silent = TRUE)
    
    if(class(trycheck) == 'try-error'){
      next
    }
    if(all(is.na(g.params[[1]]))){
      next
    }
    
    g.params[[1]] <- c(g.params[[1]], g.params[[3]])
    names(g.params[[1]]) <- c('alpha1', 'alpha2', 'beta1', 'beta2', 'incr.thresh')
    gr.ps <- g.params[[1]]
    names(gr.ps) <- names(g.params[[1]])
    # for each set of parameters make the growth distributions and get the 
    # expectation of growth - then take the median of these
    chains <- g.params[[2]]
    mix.dist <- mix.gamma(xx, apply(chains, 2, median), q)
    thresh <- g.params[[3]]
    x1 <- xx[xx < thresh]
    x2 <- xx[xx >= thresh]
    fn <- splinefun(xx, mix.dist)
    cdf_fn.slow <- cdf(fn, min(xx), thresh)
    cdf_fn.fast <- cdf(fn, thresh, max(xx))
    
    cdf_inv.slow <- inverse(cdf_fn.slow, min(xx), thresh)
    cdf_inv.fast <- inverse(cdf_fn.fast, thresh, max(xx))
    ex.slow <- cdf_inv.slow(0.5)
    ex.fast <- try(cdf_inv.fast(0.5), silent = TRUE)
    if(class(ex.fast) == 'try-error'){
      next()
    }
    
    gr.ps <- c(gr.ps, ex.slow, ex.fast)
    gr.ps <- matrix(gr.ps, nrow = 1)
    gr.ps <- as.data.frame(gr.ps)
    gr.ps <- cbind(intervals, gr.ps, stringsAsFactors = FALSE)
    gr.ps <- cbind(sp.latin, gr.ps, stringsAsFactors = FALSE)
    
    colnames(gr.ps) <- c('Species', 'Census', 'Alpha1', 'Alpha2', 'Beta1', 'Beta2', 
                         'Distribution threshold', 'DBH increment slow', 'DBH increment fast')
    growth.tab[[sp]] <- gr.ps
    print(sp)
  }
  
  growth.tab <- do.call(rbind, growth.tab)
  
  save(growth.tab, 
       file = sprintf('Vitals_Zone/Output/latex_%s_growthp_table.RData',
                      site))
  
  ## repeat but make a version that is better for working in PCA zone
  growth.tab <- vector('list', n.sp)
  xx <- seq(0.1, 50, 0.1)
  xlim <- c(min(xx), max(xx))
  
  for(sp in 1:n.sp){
    
    sp.name <- sp.names[sp]
    sp.latin <- sp.latins[sp]
    
    trycheck <- try(load(file = sprintf('Vitals_Zone/Output/%s/%s_growthps_%s.RData', 
                                        site,   site,  sp.name)), silent = TRUE)
    
    if(class(trycheck) == 'try-error'){
      next
    } 
    if(all(is.na(g.params[[1]]))){
      next
    }
    
    g.params[[1]] <- c(g.params[[1]], g.params[[3]])
    names(g.params[[1]]) <- c('alpha1', 'alpha2', 'beta1', 'beta2', 'incr.thresh')
    gr.ps <- matrix(g.params[[1]], 
                    nrow = 1)
    colnames(gr.ps) <- names(g.params[[1]])
    # for each set of parameters make the growth distributions and get the 
    # expectation of growth - then take the median of these
    ex.slow <- ex.fast <- rep(NA, n.intervals)
    
    chains <- g.params[[2]]
    mix.dist <- mix.gamma(xx, apply(chains, 2, median), q)
    thresh <- g.params[[3]]
    x1 <- xx[xx < thresh]
    x2 <- xx[xx >= thresh]
    fn <- splinefun(xx, mix.dist)
    cdf_fn.slow <- cdf(fn, min(xx), thresh)
    cdf_fn.fast <- cdf(fn, thresh, max(xx))
    cdf_inv.slow <- inverse(cdf_fn.slow, min(xx), thresh)
    cdf_inv.fast <- inverse(cdf_fn.fast, thresh, max(xx))
    ex.slow <- cdf_inv.slow(0.5)
    ex.fast <- try(cdf_inv.fast(0.5), silent = TRUE)
    if(class(ex.fast) == 'try-error'){
      next()
    }
    
    gr.ps <- as.data.frame(gr.ps)
    gr.ps <- cbind(gr.ps, ex.slow, ex.fast)
    # gr.ps <- cbind(intervals, gr.ps, stringsAsFactors = FALSE)
    gr.ps <- cbind(sp.latin,  gr.ps,
                   stringsAsFactors = FALSE)
    gr.ps <- cbind(sp.name, gr.ps, stringsAsFactors = FALSE)
    colnames(gr.ps) <- c('sp', 'Latin', 'Alpha1', 'Alpha2', 'Beta1', 'Beta2', 
                         'incr.thresh', 'ex.slow', 'ex.fast')
    rownames(gr.ps) <- NULL
    gr.ps <- gr.ps[!all(is.na(gr.ps[ ,4:ncol(gr.ps)]))]
    growth.tab[[sp]] <- gr.ps
  }
  
  growth.tab <- do.call(rbind, growth.tab)
  save(growth.tab, 
       file = sprintf('Vitals_Zone/Output/%s_growthp_table.RData',
                      site))
  
  ###############################
  ### PLOTS ###
  ###############################
  cols <- brewer.pal(9, 'YlGn')[4]
  cols.u <- brewer.pal(9, 'YlOrRd')[4]
  cols.l <- brewer.pal(9, 'YlGnBu')[4]
  
  cols.u.light <- make.ramp(cols.u, transp = 40,
                            levels = length(cols.u))
  cols.l.light <- make.ramp(cols.l, transp = 40,
                            levels = length(cols.l))
  
  G <- 2
  dists <- rep('gamma', G)
  keep.sp <- seq(n.sp)
  
  ##########################################################################
  # plot diagnostics 
  pdf(file = sprintf('Vitals_Zone/Figures/%s_Growth_diagnostics.pdf', site), 
      bg = 'white')
  
  par(mfrow = c(3,3), mar = c(5,4,4,2), oma = c(2,3,0,0))
  
  
  for(sp in 1:n.sp){
    
    sp.name <- sp.names[sp]
    sp.latin <- sp.latins[sp]
    
    trycheck <- try(load(file = sprintf('Vitals_Zone/Output/%s/%s_growthps_%s.RData', 
                                        site,   site,  sp.name)), silent = TRUE)
    
    if(class(trycheck) == 'try-error' | 
       all(is.na(g.params[[1]]))){
      next
    }
    
    chains <- g.params[[2]]
    plot(chains[ ,ncol(chains)], main = bquote(''~italic(.(sp.latin))), 
         ylab = '')
    print(sp)
  }
  
  dev.off()
  
  print(site)
}
