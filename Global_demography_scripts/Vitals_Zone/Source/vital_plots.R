rm(list = ls())		
### LIBRARIES ###		
library(foreach)		
library(doParallel)		
library(rstan)		
library(RColorBrewer)		
library(xtable)		


### FUNCTIONS ###		
# make nice colour palettes		
make.ramp <- function(tmpcol, levels = 10, transp = NA) {		
  col.out <- colorRampPalette(tmpcol)(levels)		
  if(!is.na(transp)) {		
    col.out <- paste(col.out, transp, sep = "")		
  }		
  return(col.out)		
}		

## four parameter survival curve		
predict.surv <- function(x, params, time = 1){		
  K <- params[1]		
  ip <- params[2]		
  r <- params[3]		
  # Linear predictor parameters		
  pred.y <-  (K /(1 + exp(-(r * ((x - ip) ))))) ^ (time)		
  return(pred.y)		
}		


### DATA ###		
## load the data		
sites <- c('amacayacu', 'bci', 	
           'changbaishan', 'fushan', 'hkk', 'ituri', 'khaochong',
           'korup', 'lambir', 		
           'laplanada', 'luquillo', 
           'mudumalai', 
           'palanan', 'pasoh',
           'scbi', 'serc', 
           'windriver', 'wytham', 'xtbg', 'yasuni')		

n.sites <- length(sites)		

### GROWTH ###
nproc <- 15
registerDoParallel(nproc)

foreach(jj = 1:n.sites)%dopar%{		
  site <- sites[jj]		
  try(data.ls <- readRDS(file = sprintf('Data_Zone/Output/%s.RData', site)), 
      silent = TRUE)
  if(site == 'mudumalai'){
    data.ls <-readRDS(file = sprintf('Data_Zone/Output/%s_spnames_years.RData', site)) 
  }
  
  # use try since not all sites have all components 		
  try(data.mat <- data.ls$data.mat, silent = TRUE)		
  try(sp.names <- data.ls$sp.names, silent = TRUE)		
  try(sp.latins <- data.ls$sp.latins, silent = TRUE)		
  try(years <- data.ls$years, silent = TRUE)		
  sp.latins <- as.character(sp.latins[ ,1])
  
  n.census <- length(years)	
  n.intervals <- n.census - 1  		
  n.sp <- length(sp.names)  # number of species		
  n.censuses <- n.census		
  colgs <- brewer.pal(9, 'YlGn')[6]		
  colgf <- brewer.pal(9, 'YlOrRd')[6]		
  colgs.light <- make.ramp(colgs, transp = 40,		
                           levels = length(colgs))		
  colgf.light <- make.ramp(colgf, transp = 40,		
                           levels = length(colgs))		
  G <- 2		
  dists <- rep('gamma', G)		
  
  pdf(file = sprintf('Vitals_Zone/Figures/%s_growth_sp.pdf', site), 		
      width = 10, bg = 'white')		
  par(mfrow = c(3,5), mar = c(3,2,2,0.5), oma = c(3,3,1,1))		
  
  load(file = sprintf('Vitals_Zone/Output/%s_growthp_table.RData',
                      site))
  
  for(sp in 1:nrow(growth.tab)){		
    
    sp.name <- growth.tab[sp, 'sp']	
    sp.latin <- growth.tab[sp, 'Latin']	
    load(file = sprintf('Vitals_Zone/Output/%s/%s_growthps_%s.RData', 
                        site, site, sp.name))
    N <- try(g.params[[4]])
   
    if(all(par()$mfg == c(3,5,3,5))){		
      par(mar = c(0,0,0,2))		
      plot.new()		
      legend('left', lwd = 2,  cex = 1,		
             legend = c('Slow', 'Fast'),		
             col = c(colgs, colgf), 		
             bty = 'n')		
      par(mar = c(3,2,2,0.5))	
      mtext('Increment (mm)', side = 1, line = 1, 		
            cex = 1.4, outer = TRUE)	
      mtext('Density', side = 2, line = 1, 		
            cex = 1.4, outer = TRUE)		
    }		
    
    xx <- seq(0.01, 35, 0.05)		
    census <- n.intervals	
    slow.cens <- dgamma(xx, shape = growth.tab[sp, 'Alpha1'],		
                        rate = growth.tab[sp, 'Beta1'])		
    fast.cens <- dgamma(xx, shape = growth.tab[sp, 'Alpha2'],		
                        rate = growth.tab[sp, 'Beta2'])		
    
    ymax <- max(c(fast.cens, slow.cens), na.rm = TRUE)		
    tmp <- strsplit(sp.latin, ' ')		
    tmp1 <- tmp[[1]][1]		
    tmp2 <- tmp[[1]][2]		
    
    plot(xx, slow.cens, type = 'l', col = colgs,		
         lwd = 2, ylim = c(0, ymax),		
         las = 1, cex.axis = 1.2,		
         xlab = '', ylab = '', bty = 'l',		
         main = '', yaxt = 'n')		
    
    points(xx, fast.cens, type = 'l', col = colgf, lwd = 2)	
    mtext(bquote(''~italic(.(tmp1))), side = 3, line = 0.8, 		
          cex = 0.8)		
    mtext(bquote(''~italic(.(tmp2))), side = 3, line = -0.05, 		
          cex = 0.8)		
    if(class(N) != 'try-error'){
    mtext(paste0('N:', N), side = 3, 
          line = -1, cex = 0.8)
    }
    # uncertainty
    slow.ps <- g.params[[2]][ ,c('alpha1', 'beta1')]
    fast.ps <- g.params[[2]][ ,c('alpha2', 'beta2')]
    slow.dists <- t(apply(slow.ps, 1, function(ps) 
      dgamma(xx, shape = ps[1], rate = ps[2])))
    fast.dists <- t(apply(fast.ps, 1, function(ps)
      dgamma(xx, shape = ps[1], rate = ps[2])))
    
    gr.quants.slow <- apply(slow.dists, 2, quantile, 
                            prob = c(0.025, 0.25, 0.5, 0.75, 0.975))
    gr.quants.fast <- apply(fast.dists, 2, quantile, 
                            prob = c(0.025, 0.25, 0.5, 0.75, 0.975))
    polygon(c(xx, rev(xx)), c(gr.quants.slow[2, ], 
                              rev(gr.quants.slow[4, ])), 
            col = colgs.light, border = NA)
    polygon(c(xx, rev(xx)), c(gr.quants.fast[2, ], 
                              rev(gr.quants.fast[4, ])), 
            col = colgf.light, border = NA)
    print(sp)
  }		
  dev.off()		
  print(site)
}


### SURVIVAL ###
foreach(jj = 1:n.sites)%dopar%{		
  
  site <- sites[jj]		
  try(data.ls <- readRDS(file = sprintf('Data_Zone/Output/%s.RData', site)), 
      silent = TRUE)
  
  # use try since not all sites have all components 		
  try(data.mat <- data.ls$data.mat, silent = TRUE)		
  try(sp.names <- data.ls$sp.names, silent = TRUE)		
  try(sp.latins <- data.ls$sp.latins, silent = TRUE)		
  try(years <- data.ls$years, silent = TRUE)		
  
  sp.latins <- as.character(sp.latins[ ,1])		
  n.census <- length(years)		
  n.intervals <- n.census - 1  		
  
  load(file = sprintf('Vitals_Zone/Output/%s_survp_table.RData',
                      site))
  
  n.sp <- nrow(surv.tab)  # number of species		
  cols <- brewer.pal(9, 'PuBu')[6]		
  cols.light <- make.ramp(cols, transp = 40,		
                          levels = length(cols))	
  
  pdf(file = sprintf('Vitals_Zone/Figures/%s_survival.pdf', site), 		
      width = 10, bg = 'white')		
  par(mfrow = c(3,5), mar = c(3,2,2,0.5), oma = c(3,3,1,1))		
  
  for(sp in 1:nrow(surv.tab)){		
    
    sp.name <- surv.tab[sp, 'sp']	
    sp.latin <- surv.tab[sp, 'Latin']	
    load(file = sprintf('Vitals_Zone/Output/%s/%s_maxsizes_%s.RData',
                        site, site, sp.name))
    load(file = sprintf('Vitals_Zone/Output/%s/%s_survps_%s.RData',
                        site, site, sp.name)) 
    N <- s.params[[3]]
    
    # plots of SURVIVAL in each census		
    thresh <- surv.tab[sp, 'thresh']
    small.surv <- surv.tab[sp, c('K', 'p1', 'r1')]
    big.surv <- surv.tab[sp, c('K', 'p2', 'r2')]		
    U <- 1900
    
    size.1a <- seq(thresh) # to the threshold		
    size.2a <- seq(max(size.1a)+1, U) # ensures no overlap 		
    s.survs <- predict.surv(size.1a, 
                            params = as.numeric(small.surv)) 		
    b.survs <- predict.surv(size.2a, 
                            params = as.numeric(big.surv))		
    
    tmp <- strsplit(sp.latin, ' ')		
    tmp1 <- tmp[[1]][1]		
    tmp2 <- tmp[[1]][2]		
    
    if(par()$mfg[2] == 1 & par()$mfg[1] == 1 | par()$mfg[2] == 5){		
      plot(seq(U), c(s.survs, b.survs), 		
           main = '',		
           las = 1, bty = 'l',		
           lwd = 2, col = cols, ylim = c(0, 1), type = "l",		
           xlab = '', ylab = '')		
    }else{		
      plot(seq(U), c(s.survs, b.survs), 		
           main = '', 		
           yaxt = 'n',		
           las = 1, bty = 'l',		
           lwd = 2, col = cols, ylim = c(0, 1), type = "l",		
           xlab = '', ylab = '')		
    }		
    
    mtext(bquote(''~italic(.(tmp1))), side = 3, line = 0.8, 		
          cex = 0.8)		
    mtext(bquote(''~italic(.(tmp2))), side = 3, line = -0.05, 		
          cex = 0.8)		
    mtext(paste0('N:', N), side = 3, line = 0.5, cex = 0.8, 
          adj = -0.25)
    if(all(par()$mfg == c(3,5,3,5))){		
      mtext('DBH (mm)', side = 1, line = 1, 		
            cex = 1.4, outer = TRUE)		
      mtext('Survival Probability', side = 2, line = 1, 		
            cex = 1.4, outer = TRUE)		
    }		
    
    thresh.tmp <- surv.tab[sp, 'thresh']
    abline(v = thresh.tmp, lty = 3)		
    
    abline(v = max.size, lty = 3, lwd = 0.8)
    chains <- s.params[[2]]
    uncert.mat <- matrix(NA, nrow(chains), U)		
    # these have already been thinned		
    for(tt in 1:dim(chains)[1]){		
      s.survs <- predict.surv(size.1a, params = chains[tt, c('K', 'p1', 'r1')]) 		
      b.survs <- predict.surv(size.2a, params = chains[tt, c('K', 'p2', 'r2')])
      uncert.mat[tt, ] <- c(s.survs, b.survs)		
    }		
    low.uncert <- apply(uncert.mat, 2, quantile, 0.25)		
    high.uncert <- apply(uncert.mat, 2, quantile, 0.75)		
    xs <- seq(U)		
    polygon(x = c(xs, rev(xs)), y = c(low.uncert, rev(high.uncert)), 		
            col = cols.light, border = cols.light)		
    print(sp) 		
  }		
  
  dev.off()		
  print(site)
}  

#############################
### SURVIVAL SCALED ###
foreach(jj = 1:n.sites)%dopar%{		
 # jj <- 1
  site <- sites[jj]		
  try(data.ls <- readRDS(file = sprintf('Data_Zone/Output/%s.RData', site)), 
      silent = TRUE)
  
  # use try since not all sites have all components 		
  try(data.mat <- data.ls$data.mat, silent = TRUE)		
  try(sp.names <- data.ls$sp.names, silent = TRUE)		
  try(sp.latins <- data.ls$sp.latins, silent = TRUE)		
  try(years <- data.ls$years, silent = TRUE)		
  
  sp.latins <- as.character(sp.latins[ ,1])		
  n.census <- length(years)		
  n.intervals <- n.census - 1  
  load(file = sprintf('Vitals_Zone/Output/%s_survp_table.RData',
                      site))
  keep <- which(sp.names %in% surv.tab[ ,'sp'])
  
  if(site != 'mudumalai'){
  sp.list <- lapply(split(data.mat, data.mat[,grep('sp.id', 
                                                   colnames(data.mat))]),
                    matrix, ncol = ncol(data.mat), 
                    dimnames = list(NULL, colnames(data.mat)))
  sp.list <- sp.list[keep]
  }
  
  n.sp <- nrow(surv.tab)  # number of species		
  cols <- brewer.pal(9, 'PuBu')[6]		
  cols.light <- make.ramp(cols, transp = 40,		
                          levels = length(cols))	
  
  pdf(file = sprintf('Vitals_Zone/Figures/%s_survival_scaled.pdf', site), 		
      width = 10, bg = 'white')		
  par(mfrow = c(3,5), mar = c(3,2,2,0.5), oma = c(3,3,1,1))		
  
  for(sp in 1:nrow(surv.tab)){		
   #for(sp in 1:24){ 
    sp.name <- surv.tab[sp, 'sp']	
    sp.latin <- surv.tab[sp, 'Latin']	
    load(file = sprintf('Vitals_Zone/Output/%s/%s_maxsizes_%s.RData',
                        site, site, sp.name))
    load(file = sprintf('Vitals_Zone/Output/%s/%s_survps_%s.RData',
                        site, site, sp.name))
    N <- try(s.params[[3]])
    
    # plots of SURVIVAL in each census		
    thresh <- surv.tab[sp, 'thresh']
    small.surv <- surv.tab[sp, c('K', 'p1', 'r1')]
    big.surv <- surv.tab[sp, c('K', 'p2', 'r2')]		
    U <- as.numeric(big.surv['p2']*1.4)
    
    size.1a <- seq(thresh) # to the threshold		
    size.2a <- seq(max(size.1a)+1, U) # ensures no overlap 		
    s.survs <- predict.surv(size.1a, 
                            params = as.numeric(small.surv)) 		
    b.survs <- predict.surv(size.2a, 
                            params = as.numeric(big.surv))		
    
    tmp <- strsplit(sp.latin, ' ')		
    tmp1 <- tmp[[1]][1]		
    tmp2 <- tmp[[1]][2]		
    
    if(par()$mfg[2] == 1 & par()$mfg[1] == 1 | par()$mfg[2] == 5){		
      plot(seq(U), c(s.survs, b.survs), 		
           main = '',		
           las = 1, bty = 'l',		
           lwd = 2, col = cols, ylim = c(0, 1), type = "l",		
           xlab = '', ylab = '')		
    }else{		
      plot(seq(U), c(s.survs, b.survs), 		
           main = '', 		
           yaxt = 'n',		
           las = 1, bty = 'l',		
           lwd = 2, col = cols, ylim = c(0, 1), type = "l",		
           xlab = '', ylab = '')		
    }		
    if(all(par()$mfg == c(3,5,3,5))){		
      mtext('DBH (mm)', side = 1, line = 1, 		
            cex = 1.4, outer = TRUE)		
      mtext('Survival Probability', side = 2, line = 1, 		
            cex = 1.4, outer = TRUE)		
    }		
    mtext(bquote(''~italic(.(tmp1))), side = 3, line = 0.8, 		
          cex = 0.8)		
    mtext(bquote(''~italic(.(tmp2))), side = 3, line = -0.05, 		
          cex = 0.8)	
    if(class(N) != 'try-error'){
    mtext(paste0('N:', N), side = 3, line = 0.5, cex = 0.8, 
          adj = -0.25)
    }
    thresh.tmp <- surv.tab[sp, 'thresh']
    abline(v = thresh.tmp, lty = 3)		
    
    load(file = sprintf('Vitals_Zone/Output/%s/%s_maxsizes_%s.RData',
                        site, site, sp.name))
    load(file = sprintf('Vitals_Zone/Output/%s/%s_survps_%s.RData',
                        site, site, sp.name))
    
    #abline(v = max.size, lty = 3, lwd = 0.8)
    chains <- s.params[[2]]
    uncert.mat <- matrix(NA, nrow(chains), U)		
    # these have already been thinned		
    for(tt in 1:dim(chains)[1]){		
      s.survs <- predict.surv(size.1a, params = chains[tt, c('K', 'p1', 'r1')]) 		
      b.survs <- predict.surv(size.2a, params = chains[tt, c('K', 'p2', 'r2')])
      uncert.mat[tt, ] <- c(s.survs, b.survs)		
    }		
    low.uncert <- apply(uncert.mat, 2, quantile, 0.25)		
    high.uncert <- apply(uncert.mat, 2, quantile, 0.75)		
    xs <- seq(U)		
    polygon(x = c(xs, rev(xs)), y = c(low.uncert, rev(high.uncert)), 		
            col = cols.light, border = cols.light)
    
    # now show the data
    if(site != 'mudumalai'){
    data.mat <- sp.list[[sp]]
    sizes <- data.mat[ ,grep('dbh', colnames(data.mat))]
    n.census <- length(grep('dbh', colnames(data.mat)))
    int <- n.census - 1  
    size <- data.mat[ ,grep('dbh', colnames(data.mat))[int]]
    surv <- data.mat[ ,grep('surv', colnames(data.mat),
                            ignore.case = TRUE)[(int)]]
    treeID <- data.mat[ ,grep('treeID', colnames(data.mat))]
    
    # only keep stems that are alive in the first census
    no.sizes <- which(is.na(size))
    if(length(no.sizes)>0){
      size <- size[-no.sizes]
      surv <- surv[-no.sizes]
      treeID <- treeID[-no.sizes]
    }
    
    # for those sites with treeIDs 
    if(length(treeID) > 0){
      # whole trees alive or dead?
      # vector length of no. trees with alive dead for whole tree
      tree.survs <- tapply(surv, treeID, 
                           mean, na.rm = TRUE)
      # only those which are all dead (i.e. a 0) get called dead
      tree.survs <- ifelse(tree.survs == 0, 0, 1)
      
      # stem surv variable = whole tree level 
      stem.tree.surv <- tree.survs[match(treeID, names(tree.survs))]
      
      # for trees that lived get largest stem that lived
      # for trees that died largest stem before death
      tree.size <- rep(NA, length(tree.survs))
      surv.mat <- cbind(treeID, size, stem.tree.surv, surv)
      trees.lived <- matrix(surv.mat[surv.mat[ ,3] == 1, ], 
                            ncol = 4)
      # stems that died get an NA for size ...
      trees.lived[,2] <- ifelse(trees.lived[,4] == 0, NA, trees.lived[,2])
      # ... so that the max size for the tree is the largest stem that lived
      tree.size[tree.survs == 1] <- tapply(trees.lived[,2],
                                           trees.lived[,1], max, na.rm = TRUE)
      # largest stem of trees that died
      trees.died <- matrix(surv.mat[surv.mat[,3] == 0, ], ncol = 4)
      tree.size[tree.survs == 0] <- tapply(trees.died[,2], 
                                           trees.died[,1], max, na.rm = TRUE)
    
      size <- tree.size
      surv <- tree.survs
    }
    
    rug(size[surv == 0], side = 1, col = 'red')
    rug(size[surv == 1], side = 3, col = 'black')
    }
    print(sp) 		
  }		
  
  dev.off()		
  print(site)
}  
