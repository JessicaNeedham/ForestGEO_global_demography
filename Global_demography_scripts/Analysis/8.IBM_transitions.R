####################################################################################
### Simulation code WITH TRANSITION PROBABILITIES BETWEEN GROWTH DISTRIBUTIONS ###
# 1. project each cluster forward starting with 10 mm stems
# and get passage times and life expectancies
# This version does not include parameter uncertainty. Use 
# the parameter medians for each cluster = variation comes 
# from demographic stochasticity

## AGB dynamics use site specific info - but are based on COHORTS from 
## simulations - i.e. not an estimate of AGB at each site. That is 
# estimated in Biomass.R
########################################################################
rm(list = ls())
####################
### LIBRARIES ###
####################
library(Matrix)
library(RColorBrewer) 
library(plyr)
library(doParallel)
library(foreach)
set.seed(1)

#####################
### FUNCTIONS ###
#####################
source('Analysis/Function_Files/IBM_functions.R')

######################
### DATA ###
######################
# load the table with parameter values and clusters
geog <- 'global'
load(file = sprintf('Analysis/Outputs/%s_pca_mat.RData', geog))
load(file = sprintf('Analysis/Outputs/PCA_%s.RData', geog))
param.mat <- pca.matrix

# split it up into clusters
cl.list <- split(param.mat, param.mat[ ,'clust'])
n.clust <- length(cl.list)

q <- 0.95

nproc <- 1
registerDoParallel(nproc)

############################
### SIMULATIONS ###
############################
foreach(cl = 1:n.clust)%do%{
  #cl <- 5
  sp.params <- cl.list[[cl]]
  
  # if there are size dependent transition probabilities between 
  # growth distributions then we need an indicator of the max size
  # we can use p2 * 1.5 - doesn't need to be accurate - just needs
  # to be big enough so that trees of all sizes can transition
  L <- 10
  U <- max(sp.params[ ,'p2'])*1.5
  S <- 1000
  h <- (U-L)/S 
  G <- 2
  meshpts <-  L + (1:S)*h -h/2 # midpoints
  boundary <-  L + c(0:S)*h  # boundary points of mesh cells
  
  # 10000 new recruits for each cluster
  start.size <- rep(10, 25000)
  
  # transition probabilities between growth distributions
  p.F.lowest <- 0.1
  p.F.highest <- 0.8
  f.m <- (p.F.lowest - p.F.highest)/(L-U)
  f.c <- p.F.highest - (f.m*U)
  # 2. EACH ITERATION MULTIPLY BY GROWTH, AND SURVIVAL - THEN MOVE
  # INDIVIDUALS BETWEEN GROWTH DISTRIBUTIONS
  n.yrs <- 1500
  
  # uncertianty from covariance matrices
  # and you need to weight them by the number of census intervals
  # so that one species from long time series plot doesn't have 
  # disproportionate influence
  s.params <- sp.params[ ,c('K', 'p1', 'r1', 'p2', 'r2', 'thresh')]
  g.params <- sp.params[ ,c('alpha1', 'alpha2', 'beta1', 'beta2')]
  
  # use weighted means - weighted by number of censuses
  s.params <- apply(s.params, 2, weighted.mean, 
                    w = sp.params[ ,'n.interval'])
  g.params <- apply(g.params, 2, weighted.mean, 
                    w = sp.params[ ,'n.interval'])
  incr.thresh <- weighted.mean(sp.params[ ,'incr.thresh'], 
                               w = sp.params[ ,'n.interval'], 
                               na.rm = TRUE)
  
  g.slow <- as.numeric(g.params[c('alpha1', 'beta1')])
  g.fast <- as.numeric(g.params[c('alpha2', 'beta2')])
  small.surv <- as.numeric(s.params[c('K', 'p1', 'r1')])
  big.surv <- as.numeric(s.params[c('K', 'p2', 'r2')])
  surv.thresh <- as.numeric(s.params['thresh'])  
  
  # make the inverse CDFs of growth parameters so that
  # we can sample
  gam.pars <- as.numeric(g.params)
  x.1 <- seq(0.001, incr.thresh, length = 100)
  x.2 <- seq(incr.thresh, 40, length = 100)
  
  y.1 <- mix.gamma(x.1, gam.pars, q)
  y.2 <- mix.gamma(x.2, gam.pars, q)
  xlim <- c(min(x.1), max(x.2))
  
  fn.1 <- splinefun(x.1, y.1)
  fn.2 <- splinefun(x.2, y.2)
  # Turn the spline function into a cumulative distribution 
  # function 
  cdf_fn.slow <- cdf(fn.1, x.1[1], incr.thresh)
  cdf_fn.fast <- cdf(fn.2, incr.thresh, x.2[length(x.2)])
  cdf_inv.lo <- inverse(cdf_fn.slow, min(x.1), incr.thresh)
  cdf_inv.hi <- inverse(cdf_fn.fast, incr.thresh, max(x.2))
  
  ### sample here!!! 
  slow.incrs <- cdf_inv.lo(runif(10000))
  fast.incrs <- cdf_inv.hi(runif(10000))
  
  # matrix to hold the population each year
  # Gmat holds which distribution each tree was in each year
  rmat <- Gmat <-  matrix(NA, ncol = n.yrs,
                          nrow = length(start.size))
  rmat[ ,1] <- start.size
  # 50 % start fast  -
  
  ### NB For actual populations this should be 5 and 95% but here we are just tryign 
  # to demonstrate what the passage times and life expectancies are for fast and slow
  # in each cluster - so sample sizes should be the same
  start.fast <- sample(length(start.size), 
                       round(length(start.size)*0.05), replace = FALSE)
  Gmat[ ,1] <- 1
  Gmat[start.fast, 1] <- 2
  yr <- 2
  
  while(yr <= n.yrs){
    
    # size vector from previous year
    z <- rmat[ ,yr-1]
    # growth distribution vector from previous year
    g <- Gmat[ ,yr-1]  
    
    # GROWTH
    slow <- which(g == 1)
    fast <- which(g == 2)
    
    ### 
    z1 <- rep(NA, length(z))
    
    if(length(slow) > 0){
      z1[slow] <- z[slow] + sample(slow.incrs, length(slow), replace = TRUE)
    }
    if(length(fast) > 0){
      z1[fast] <- z[fast] + sample(fast.incrs, length(fast), replace = TRUE)
    }
    
    # SURVIVAL
    s1 <- which(z < surv.thresh)
    s2 <- which(z >= surv.thresh)
    
    s1.surv <- rbinom(length(s1), 1, 
                      s_z(z[s1], small.surv))
    s2.surv <- rbinom(length(s2), 1, 
                      s_z(z[s2], big.surv))
    z1[s1[which(s1.surv == 0)]] <- NA
    z1[s2[which(s2.surv == 0)]] <- NA
    
    # GROWTH CLASS TRANSITIONS
    gs <- transition.trees(z1, f.m, f.c)
    Gmat[ ,yr] <- gs
    
    # add everything to the master matrices
    rmat[ ,yr] <- z1
    Gmat[ ,yr] <- gs
    
    # print(yr)
    yr <- yr + 1
  }
  
  # now separate the population into those that started fast and those 
  # that started slow
  slow.mat <- rmat[which(Gmat[ ,1] == 1), ]
  fast.mat <- rmat[which(Gmat[ ,1] == 2), ]
  
  # size summary each year
  slow.size.summary <- apply(slow.mat, 2, 
                             quantile, prob = c(0, 0.025, 0.25, 0.5, 
                                                0.75, 0.975, 1), na.rm = TRUE)
  fast.size.summary <- apply(fast.mat, 2, 
                             quantile, prob = c(0, 0.025, 0.25, 0.5, 
                                                0.75, 0.975, 1), na.rm = TRUE)
  
  # size at death 
  dbh.at.death.slow <- apply(slow.mat, 1, max, na.rm = TRUE)
  dbh.at.death.fast <- apply(fast.mat, 1, max, na.rm = TRUE)
  
  # percent of time in slow
  g.slow.mat <- Gmat[which(Gmat[ ,1] == 1), ]
  g.fast.mat <- Gmat[which(Gmat[ ,1] == 2), ]
  
  percent.slows <- apply(g.slow.mat, 1, function(x){
    x <- x[!is.na(x)];
    s <- length(which(x == 1));
    percent.s <- s/length(x)*100;
    return(percent.s)
  })
  pc.slow <- mean(unlist(percent.slows))
  pc.fast <- mean(100 - unlist(percent.slows))
  
  # put them together to get pts and les
  pts <- get.passage.time.ibm(rmat, threshold = 100)
  les <- get.lifespan.ibm(rmat)
  
  # population size through time
  pop.N <- apply(rmat, 2, function(x) length(which(!is.na(x))))
  plot(pop.N, type = 'l')
  start <- length(start.size)
  # year that 99% of pop is extinct
  pop.99.ex.yr <- which.min(abs(pop.N - (start/100)))
  # size 99th quantile of size at death
  dbh.at.death <- apply(rmat, 1, max, na.rm = TRUE)
  quantile(dbh.at.death, 0.99)
  
  
  save(les, pts, slow.size.summary, fast.size.summary, 
       dbh.at.death.slow, dbh.at.death.fast, 
       pc.slow, pc.fast, percent.slows, Gmat, rmat, 
       file = sprintf('Analysis/Outputs/%d_pts_les_transitions.RData', cl))
  
}  # end of clusters


#######################################################
### Now for each mode at each site get the biomass

source('Analysis/Function_Files/Biomass_functions.R')

##################
### DATA ###
##################
# The pca matrix and the site dataframe with meta data for each site
geog <- 'global'
load(file = sprintf('Analysis/Outputs/%s_pca_mat.RData', geog))
load('Analysis/Outputs/site_df.RData')
# remove sites with no data for now! - but add back later 
pca.matrix <- pca.matrix[-which(pca.matrix[ ,'Site'] %in% c('mudumalai')), ]
pca.matrix[ ,'Site'] <- as.character(pca.matrix[ ,'Site'])
site.df <- site.df[-which(site.df[ ,'Site'] == 'Mudumalai'), ]
sites <- as.character(unique(pca.matrix[ ,'Site']))
n.sites <- length(sites)
n.clust <- length(unique(pca.matrix[, 'clust']))

# load the dft mean wsg
load(file = 'Analysis/Outputs/clust_wsgs.RData')

# load Es - the climate variables needed for allometries
Es <- read.table('Data_Zone/Data/EforBiomass.txt')
full.site.names <- c('Amacayacu', 'Barro Colorado Island', 'Changbaishan', 'Fushan', 'Huai Kha Khaeng', 
                     'Ituri', 'Khao Chong', 'Korup', 'Lambir', 'La Planada', 'Luquillo', 'Palanan', 'Pasoh', 'SCBI', 'SERC', 
                     'Wytham Woods', 'Wind River', 'Xishuangbanna ', 'Yasuni')
E_biomass <- Es[match(full.site.names, Es[ ,'Site']), 'E']


pca.mats <- split(pca.matrix, pca.matrix[ ,'Site'])

n.yrs <- 1500
#nproc <- length(sites)
#registerDoParallel(nproc)

# use Amacayacu as an example 
ii <- 1 
site <- sites[ii]
E <- E_biomass[ii]

agb.dfts <- vector('list', n.clust)

for(cl in 1:n.clust){
  #cl <- 2
  
  load(sprintf("Analysis/Outputs/%s_pts_les_transitions.RData", cl))
  
  # ok - now do the AGB for cohorts of 1000 individuals
  slows <- rmat[which(Gmat[ ,1] == 1), ]
  fasts <- rmat[which(Gmat[ ,1] == 2), ]
  slow.samp.size <- 950
  fast.samp.size <- 50
  
  slow.samps <- replicate(250, sample(seq(nrow(slows)), 
                                      size = slow.samp.size, replace = TRUE), 
                          simplify = FALSE)
  fast.samps <- replicate(250, sample(seq(nrow(fasts)), 
                                      size = fast.samp.size, replace = TRUE), 
                          simplify = FALSE)
  
  agbs <- matrix(NA, length(slow.samps), ncol = n.yrs)
  
  for(ss in 1:length(slow.samps)){
    slowrsamp <- matrix(slows[slow.samps[[ss]], ], ncol = n.yrs)
    fastrsamp <- matrix(fasts[fast.samps[[ss]], ], ncol = n.yrs)
    rsamp <- rbind(slowrsamp, fastrsamp)
    
    # get sizes at each time step and calculate agb for each year
    # AGB each year of each stem 
    stems.AGB <- apply(rsamp, 2, function(x)  # Chave's tropical equations
      agb.est(x/10, E, clust.wsgs[cl]))
    
    stems.AGB <- stems.AGB/1000  # make it in Mg 
    dft.AGB <- colSums(stems.AGB, na.rm = TRUE)  # total AGB each year
    agbs[ss, ] <- dft.AGB
  }
  
  agb.dfts[[cl]] <- agbs
  
}

save(agb.dfts, file = 'Analysis/Outputs/agb_dfts.RData')
