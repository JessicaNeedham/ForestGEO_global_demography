####################################################################
### Fits two logistic regressions to each species that meet at 
### a size threshold - captures size-dependent survival probability
### using STAN MCMC ###
####################################################################
rm(list = ls())
### LIBRARIES ###
library(foreach)
library(doParallel)
library(rstan)
library(RColorBrewer)
library(xtable)

rstan_options(auto_write = TRUE)

### FUNCTIONS ###
## four parameter survival curve
predict.surv <- function(x, params, time = 1) {
  K <- params[1]
  ip <- params[2]
  r <- params[3]
  # Linear predictor parameters
  pred.y <-  (K / 
                (1 + exp(-(r * ((x - ip) )))))  ^ (time)
  return(pred.y)
}

# make nice colour palettes
make.ramp <- function(tmpcol, levels = 10, transp = NA) {
  col.out <- colorRampPalette(tmpcol)(levels)
  if(!is.na(transp)) {
    col.out <- paste(col.out, transp, sep = "")
  }
  return(col.out)
}

### DATA ###
## load the data
sites <- c('amacayacu', 'bci', 'changbaishan', 'fushan', 'hkk', 
           'ituri', 'khaochong', 'korup', 'lambir', 'laplanada', 
           'luquillo', 'pasoh', 'palanan',
           'scbi', 'serc', 'windriver', 
           'wytham', 'yasuni', 'xtbg')

n.sites <- length(sites)

for(ii in 1:n.sites){
  #ii <- 1 
  site <- sites[ii]
  
  data.ls <- readRDS(file = sprintf('Data_Zone/Output/%s.RData', site))
  # use try since not all sites have all components 
  try(data.mat <- data.ls$data.mat, silent = TRUE)
  try(sp.names <- data.ls$sp.names, silent = TRUE)
  try(sp.latins <- data.ls$sp.latins, silent = TRUE)
  try(years <- data.ls$years, silent = TRUE)
  
  sp.latins <- as.character(sp.latins[ ,1])
  mat.colnames <- colnames(data.mat)
  
  n.census <- length(grep('dbh', colnames(data.mat)))
  n.intervals <- n.census - 1  
  n.sp <- length(sp.names)  # number of species
  sp.list <- lapply(split(data.mat, data.mat[,grep('sp.id', 
                                                   colnames(data.mat))]),
                    matrix, ncol = ncol(data.mat), 
                    dimnames = list(NULL, mat.colnames))
  
  # here we just generate some fake data in order to pass it to stan()
  # function so that we can compile the stan model outside of the 
  # foreach loop
  N <- 100
  size <- rnorm(N, 50, 5)
  surv <- rbinom(N, 1, 0.1)
  time <- rep(5, N)
  thresh <- 30
  lowr2 <- 20
  upr2 <- 200
  
  # make the data list to run stan model
  surv.data <- list(N = length(size), 
                    dbh = size, 
                    surv = surv, 
                    time = time,
                    thresh = thresh, 
                    lowr2 = lowr2,
                    upr2 = upr2)    
  
  # compile the model (it may take 40-60 secs)
  surv.model <- stan(file = 'Vitals_Zone/Source/surv.stan', 
                     data = surv.data, chains = 0)
  
  #n.cores <- detectCores()-1  # for quantlab
  n.cores <- 16 #
  registerDoParallel(n.cores)
  
  n.iter <- 2500
  
  foreach(sp = 1:n.sp, 
          .errorhandling = 'pass', 
          .packages = 'rstan')%dopar%{
            data.mat <- sp.list[[sp]]
            sp.name <- sp.names[sp]
            
            # matrix of just sizes
            sizes <- data.mat[ ,grep('dbh', colnames(data.mat))]
            
            # to hold outputs
            s.params <- vector('list', 3)
            s.params[[1]] <- rep(NA, 7)
            s.params[[2]] <- NULL
            s.params[[3]] <- NULL
            
            # get time in years between censuses 
            time.vec <- data.mat[ ,grep('time', colnames(data.mat))]
            
            # get vectors of size survival and time
            # because we are working with stem data need to make sure that
            # all stems of a tree die before we call it dead
            size <- data.mat[ ,grep('dbh.1', colnames(data.mat))]
            surv <- data.mat[ ,grep('surv', colnames(data.mat),
                                    ignore.case = TRUE)]
            treeID <- data.mat[ ,grep('treeID', colnames(data.mat))]
            
            # only keep stems that are alive in the first census
            no.sizes <- which(is.na(size))
            
            if(length(no.sizes)>0){
              size <- size[-no.sizes]
              surv <- surv[-no.sizes]
              time.vec <- time.vec[-no.sizes]
              treeID <- treeID[-no.sizes]
            }
            
            # for those sites with treeIDs 
            if(length(treeID) > 0){
              
              # index of largest stem per tree
              max.size <- tapply(size, treeID, which.max)
              # corresponding to rows... 
              size <- tapply(size, treeID, max)  # max size of each tree
              
              time <- mapply(function(x,y) x[y], 
                             x = split(time.vec, treeID), 
                             y = max.size)
              
              # whole trees alive or dead?
              # vector length of no. trees with alive dead for whole tree
              tree.survs <- tapply(surv, treeID, 
                                   mean, na.rm = TRUE)
              # only those which are all dead (i.e. a 0) get called dead
              tree.survs <- ifelse(tree.survs == 0, 0, 1)
              surv <- tree.survs
            }
            time[is.na(time)] <- sample(time[!is.na(time)], 
                                        length(which(is.na(time))))
            
            ### if there are less than 100 individuals in a census interval then 
            # put NAs and move on
            if(length(size) < 20){
              next('Too few trees')
            }
            
            max.size <- max(size, na.rm = TRUE)
            lowr2 <- -0.05
            upr2 <- -0.00001
            thresh <- max.size*0.2
            
            # make the data list to run stan model
            surv.data <- list(N = length(size), 
                              dbh = size, 
                              surv = as.numeric(surv), 
                              time = time,
                              thresh = thresh, 
                              lowr2 = lowr2,
                              upr2 = upr2)    
            
            # run the MCMC
            surv.fit <- stan(fit = surv.model, data = surv.data, iter = n.iter, 
                             chains = 3, control = list(adapt_delta = 0.98))
            
            # extract the sampling parameters
            surv.mcmc <- lapply(extract(surv.fit), median)
            surv.chains <- do.call(cbind, extract(surv.fit))
            s.params[[2]] <- surv.chains
            s.params[[1]][1:3] <- unlist(surv.mcmc)[c(1,2,3)]
            s.params[[1]][4:6] <- unlist(surv.mcmc)[c(1,4,5)]
            s.params[[1]][7] <- thresh
            s.params[[3]] <- length(size)
            
            
            save(max.size, 
                 file = sprintf('Vitals_Zone/Output/%s/%s_maxsizes_%s.RData',
                                site, site, sp.name))
            save(s.params, 
                 file = sprintf('Vitals_Zone/Output/%s/%s_survps_%s.RData',
                                site, site, sp.name))
          }
  
  
  surv.tab <- vector('list', n.sp)
  intervals <- na.omit(expand.grid(years, years)[(seq(5)*(n.census+1)), ])
  intervals <- matrix(paste(intervals[,1], intervals[,2], sep = '-'), 
                      nrow = n.intervals)
  intervals <- intervals[n.intervals, 1]
  
  for(sp in 1:n.sp){
    
    sp.name <- sp.names[sp]
    sp.latin <- sp.latins[sp]
    data.mat <- sp.list[[sp]]
    dbhs <- data.mat[ ,grep('dbh', colnames(data.mat))]
    
    trycheck <- try(load(file = sprintf('Vitals_Zone/Output/%s/%s_survps_%s.RData', 
                                        site, site, sp.name)), silent = TRUE)
    # skip species where model didn't run in 
    # any of the census intervals
    if(class(trycheck) == 'try-error'){
      next
    }
    if(all(is.na(s.params[[1]]))){
      next
    }
    
    load(file = sprintf('Vitals_Zone/Output/%s/%s_maxsizes_%s.RData',
                        site, site, sp.name))
    
    # median parameters in each census interval
    surv.ps <- matrix(c( s.params[[3]],
                         s.params[[1]][c(1,2,3,5,6,7)]), nrow = 1)
    surv.ps <- as.data.frame(surv.ps)
    surv.ps <- cbind(intervals, surv.ps, stringsAsFactors = FALSE)
    surv.ps <- cbind(c(sp.latin, rep(' ', nrow(surv.ps)-1)), surv.ps, 
                     stringsAsFactors = FALSE)
    colnames(surv.ps) <- c('Species', 'Census', 'N', 'K', 'p1', 
                           'r1', 'p2', 'r2', 'thresh')
    rownames(surv.ps) <- NULL
    # remove any row with no parameters - i.e.
    # a census interal where pop size was too small to fit model
    surv.ps <- surv.ps[!all(is.na(surv.ps[ ,3:8])), ]
    
    U <- max.size
    xx <- seq(U)
    
    spost <- s.params[[2]]
    extras <- matrix(NA, nrow = nrow(spost), ncol = 5)
    for(tt in 1:nrow(spost)){
      yy <- predict.surv(xx,
                         params = spost[tt, 1:3])
      extras[tt, 1] <- max(yy, na.rm = TRUE)
      derivs <- diff(yy)/diff(xx)
      # rate at 10mm drawn as tangent to curve
      extras[tt, 2] <- yy[10]
      extras[tt, 3] <- derivs[10]
      yy <- predict.surv(xx, params = spost[tt, c(1,4,5)])
      derivs <- diff(yy)/diff(xx)
      extras[tt, 4] <- yy[length(yy)]
      extras[tt, 5] <- derivs[length(derivs)]
    }
    
    minmax <- apply(extras, 2, 
                    quantile, prob = 0.5, na.rm = TRUE)
    minmax <- matrix(minmax, nrow =1)
    colnames(minmax) <- c('Max surv', 'Surv 10', 'Rate 10',
                          'Surv max', 'Rate max')
    surv.ps <- cbind(surv.ps, minmax, stringsAsFactors = FALSE)
    surv.tab[[sp]] <- surv.ps
    print(sp)
  }
  surv.tab <- do.call(rbind, surv.tab)
  
  save(surv.tab, 
       file = sprintf('Vitals_Zone/Output/latex_%s_survp_table.RData',
                      site))
  
  
  
  surv.tab <- vector('list', n.sp)
  intervals <- na.omit(expand.grid(years, years)[(seq(5)*(n.census+1)), ])
  intervals <- matrix(paste(intervals[,1], intervals[,2], sep = '-'), 
                      nrow = n.intervals)
  intervals <- intervals[n.intervals, 1]
  
  ####
  ## repeat but make the table suitable for PCA zone
  surv.tab <- vector('list', n.sp)
  
  for(sp in 1:n.sp){
    
    sp.name <- sp.names[sp]
    sp.latin <- sp.latins[sp]
    data.mat <- sp.list[[sp]]
    dbhs <- data.mat[ ,grep('dbh', colnames(data.mat))]
    
    trycheck <- try(load(file = sprintf('Vitals_Zone/Output/%s/%s_survps_%s.RData', 
                                        site, site, sp.name)), silent = TRUE)
    
    # skip species where model didn't run in 
    # any of the census intervals
    if(class(trycheck) == 'try-error'){
      next
    }
    if(all(is.na(s.params[[1]]))){
      next
    }
    
    load(file = sprintf('Vitals_Zone/Output/%s/%s_maxsizes_%s.RData',
                        site, site, sp.name))
    
    # median parameters in each census interval
    surv.ps <- matrix(c(s.params[[3]], 
                        s.params[[1]][c(1,2,3,5,6,7)]), nrow = 1)
    surv.ps <- as.data.frame(surv.ps)
    surv.ps <- cbind(sp.latin,  surv.ps, stringsAsFactors = FALSE)
    surv.ps <- cbind(sp.name, surv.ps, stringsAsFactors = FALSE)
    colnames(surv.ps) <- c('sp', 'Latin', 'N', 'K', 'p1', 
                           'r1', 'p2', 'r2', 'thresh')
    rownames(surv.ps) <- NULL
    # remove any row with no parameters - i.e.
    # a census interal where pop size was too small to fit model
    surv.ps <- surv.ps[!all(is.na(surv.ps[ ,3:8])), ]
    
    U <- max.size
    xx <- seq(U)
    
    spost <- s.params[[2]]
    
    extras <- matrix(NA, nrow = nrow(spost), ncol = 5)
    for(tt in 1:nrow(spost)){
      yy <- predict.surv(xx,
                         params = spost[tt, 1:3])
      extras[tt, 1] <- max(yy, na.rm = TRUE)
      derivs <- diff(yy)/diff(xx)
      # rate at 10mm drawn as tangent to curve
      extras[tt, 2] <- yy[10]
      extras[tt, 3] <- derivs[10]
      yy <- predict.surv(xx, params = spost[tt, c(1,4,5)])
      derivs <- diff(yy)/diff(xx)
      extras[tt, 4] <- yy[length(yy)]
      extras[tt, 5] <- derivs[length(derivs)]
    }
    
    minmax <- apply(extras, 2, 
                    quantile, prob = 0.5, na.rm = TRUE)
    minmax <- matrix(minmax, nrow = 1)
    colnames(minmax) <- c('max.surv', 'surv.10', 'rate.10',
                          'surv.max', 'rate.max')
    surv.ps <- cbind(surv.ps, minmax, stringsAsFactors = FALSE)
    surv.tab[[sp]] <- surv.ps
    print(sp)
  }
  
  surv.tab <- do.call(rbind, surv.tab)
  
  save(surv.tab, 
       file = sprintf('Vitals_Zone/Output/%s_survp_table.RData',
                      site))
  
  
  ################
  ##############################################################
  ### plot diagnostics
  pdf(file = sprintf('Vitals_Zone/Figures/%s_Survival_diagnostics.pdf', site), 
      bg = 'white')
  
  par(mfrow = c(3,3), mar = c(5,4,4,2), oma = c(2,3,0,0))
  
  
  for(sp in 1:n.sp){
    sp.name <- sp.names[sp]
    sp.latin <- sp.latins[sp]
    
    trycheck <- try(load(file = sprintf('Vitals_Zone/Output/%s/%s_survps_%s.RData',
                                        site, site, sp.name)), silent = TRUE)
    
    if(class(trycheck) == 'try-error' | 
       all(is.na(s.params[[1]]))){
      next
    }
    
    chains <- s.params[[2]]
    plot(chains[ ,ncol(chains)], main = bquote(''~italic(.(sp.latin))), 
         ylab = '')
    
    print(sp)
  }
  
  dev.off()
  
  print(site)
}
