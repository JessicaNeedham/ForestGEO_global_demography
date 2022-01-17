#################################################
#################################################
rm(list = ls())
set.seed(1)
########################
### LIBRARIES ###
library(plotrix)
library(RColorBrewer)
library(FactoMineR)
library(plyr)
library(stringi)
library(factoextra)
library(tmvtnorm)
library(corpcor)
library(xtable)
###########################
### FUNCTIONS ###
###########################
#options(device = 'quartz')
source('Analysis/Function_files/PCA_functions.R')

#############################
### DATA ### 
#############################
global.sites <- c('amacayacu', 'bci',
                  'changbaishan', 'fushan',
                  'hkk', 'ituri',
                  'khaochong', 
                  'korup',
                  'lambir', 
                  'laplanada', 'luquillo', 
                  'mudumalai', 
                  'palanan', 'pasoh', 
                  'scbi', 'serc',
                  'windriver', 'wytham', 
                  'xtbg', 'yasuni')
###############################################
###
n.sites <- length(global.sites)
p.list <- vector('list', n.sites)

for(jj in 1:n.sites){
  #jj <- 2 
  site <- global.sites[jj]
  load(file = sprintf('Vitals_Zone/Output/%s_survp_table.RData', site))
  load(file = sprintf('Vitals_Zone/Output/%s_growthp_table.RData', site))
  data.ls <- try(readRDS(file = sprintf('Data_Zone/Output/%s.RData', site)), 
                 silent = TRUE)
  
  if(site == 'mudumalai'){
    data.ls <- readRDS(file = sprintf('Data_Zone/Output/%s_spnames_years.RData', 
                                      site))
  }
  
  # use try since not all sites have all components 
  try(data.mat <- data.ls$data.mat, silent = TRUE)
  try(sp.names <- data.ls$sp.names, silent = TRUE)
  try(sp.latins <- data.ls$sp.latins, silent = TRUE)
  try(years <- data.ls$years, silent = TRUE)
  
  n.census <- length(years)
  
  if(site == 'mudumalai'){
    surv.tab <- cbind(surv.tab[ ,1:2], 1, 
                      surv.tab[ ,3:ncol(surv.tab)])
    colnames(surv.tab)[3] <- 'N'
  }
  
  all.tmp <- merge(surv.tab, growth.tab, 
                   by = intersect(names(surv.tab)[1:2], 
                                  names(growth.tab)[1:2]), 
                   all = TRUE)
  
  # get rid of species with NAs for any parameter - cannot have them in the PCA
  na.rows <- apply(all.tmp, 1, function(x) any(is.na(x)))
  all.tmp <- all.tmp[which(na.rows == FALSE), ]
  # add a site column
  all.tmp <- cbind(site, all.tmp)
  p.list[[jj]] <- all.tmp
  print(jj)
}

# one big matrix for all species, sites and census intervals
pca.matrix <- do.call(rbind, p.list)

### PROCESSING ### 
### THE PCA ###
colnames(pca.matrix) <- c('Site', 'Mnemonic', 'Species', 'N',
                          'K', 'p1', 'r1', 'p2', 'r2', 
                          'thresh', 'maxsurv', 'surv.sm', 'rate.sm', 'surv.lg', 
                          'rate.lg', 'alpha1', 'alpha2', 'beta1',
                          'beta2', 'incr.thresh', 
                          'ex.slow', 'ex.fast')

var.cols <- as.matrix(pca.matrix[ ,c('maxsurv', 'surv.sm', 
                                     'ex.slow', 'ex.fast')])

# logistic the bounded variables
surv.probs <- var.cols[ ,c('maxsurv', 'surv.sm')]
logistic.surv.probs <- apply(surv.probs, 2, function(x) logistic(x))
colnames(logistic.surv.probs) <- c('logistic.maxsurv', 'logistic.surv.sm')
pca.matrix <- cbind(pca.matrix, logistic.surv.probs, 
                    stringsAsFactors = FALSE)


# for each site load the species table - check IDlevel
for(ii in 1:length(global.sites)){
  #for(ii in 1:10){ 
  site <- global.sites[ii]
  
  if(site %in% c('amacayacu', 'changbaishan', 'mudumalai', 'khaochong',
                 'scbi', 'serc', 'windriver', 'wytham')){
    next
  }
  
  # load the spp table
  tmp <- readRDS(sprintf('Data_Zone/Data/%s/%s.spptable_V2.rdata', site, site))
  
  if('IDLevel' %in% colnames(tmp) == FALSE){
    next
  }
  uni <- tmp[which(tmp[ ,'IDLevel'] == 'multiple'), 'sp']
  
  # remove them from the pcamatrix
  pca.matrix <- pca.matrix[(pca.matrix[ ,'Mnemonic'] %in% uni == FALSE), ]
  
  print(ii)
}

# remove species with unidentified
unidentified <- grep('unidentified', pca.matrix[ ,'Species'])
pca.matrix <- pca.matrix[-unidentified, ]

#### Weighting species by the number of sites they occur in -
# this only works for known binomials - NOT THE MORPHO SPECIES

# find number of sites of each species
n.intervals <- table(pca.matrix[ ,'Species'])
# find morpho species
morphs <- grep('sp.', names(n.intervals), fixed = TRUE)
# find duplicates
dupes <- as.numeric(which(n.intervals > 1))
# find the intersection  - i.e. morpho species in two sites
rm.sp <- intersect(morphs, dupes)
# set them to 1
n.intervals[rm.sp] <- 1
# make a vector with number of sites for each species
n.ints <- n.intervals[match(pca.matrix[ ,'Species'], names(n.intervals))]
# make a column for the number of sites a species is in (not the morphs)
pca.matrix <- cbind(as.numeric(n.ints), pca.matrix)
colnames(pca.matrix)[1] <- 'n.interval' 

# load max size for each species.
pca.matrix[ ,'maxsize'] <- rep(NA, nrow(pca.matrix))

for(sp in 1:nrow(pca.matrix)){
  site <- as.character(pca.matrix[sp, 'Site'])
  
  sp.name <- pca.matrix[sp, 'Mnemonic']
  tmp <- try(load(file = sprintf('Vitals_Zone/Output/%s/%s_maxsizes_%s.RData',
                                 site, site, sp.name)))
  if(class(tmp) == 'try-error'){next}
  
  pca.matrix[sp ,'maxsize'] <- max.size
}

na.rows <- apply(pca.matrix, 1, function(x) any(is.na(x)))
pca.matrix <- pca.matrix[which(na.rows == FALSE), ]

res <- PCA(pca.matrix[ ,c('logistic.maxsurv', 
                          'logistic.surv.sm', 
                          'maxsize',
                          'ex.slow', 'ex.fast')],
           scale.unit = TRUE, 
           row.w = 1/pca.matrix[ ,'n.interval'],  
           graph = FALSE, 
           ncp = 5)
geog <- 'global'
save(res, file = sprintf('Analysis/Outputs/PCA_%s.RData', geog))

##################################
### THE CLUSTERS ###
##################################
set.seed(100)
n.clust <- 8
kmeans.res <- kmeans(x = res$ind$coord, centers = n.clust, 
                     nstart = 200, iter.max = 20)

# add clusters to main matrix
pca.matrix <- cbind(pca.matrix, kmeans.res$cluster)
colnames(pca.matrix)[ncol(pca.matrix)] <- 'clust'
# order by ex.fast
ord <- as.numeric(names(sort(tapply(pca.matrix[ ,'ex.fast'], pca.matrix[ ,'clust'], median))))
pca.matrix[ ,'clust'] <- mapvalues(pca.matrix[ ,'clust'], 
                                   from = ord, to = seq(n.clust))

save(pca.matrix, kmeans.res, 
     file = sprintf('Analysis/Outputs/%s_pca_mat.RData', geog))

params.coords.matrix <- cbind(pca.matrix, res$ind$coord)
save(params.coords.matrix, file = 'Analysis/Outputs/params_coords_matrix.RData')


########################################################

######################################################
### EXAMPLE SPECIES FOR EACH DEMOGRAPHIC MODE ###

# Given the coordinates of each species/census combo
# which ones are closest to cluster mid-points
geog <- 'global'

load(file = sprintf('Analysis/Outputs/%s_pca_mat.RData', geog))
load(file = sprintf('Analysis/Outputs/PCA_%s.RData', geog))

# which species are closest to this?
# split the individual loadings into clusters
clst.list <- lapply(split(res$ind$coord, pca.matrix[ ,'clust']), 
                    matrix, ncol = ncol(res$ind$coord))

# pca matrix clsuters have been sorted by ex.fast of each cluster
# but kmeans.res has not be sorted = do that now 
ord <- match(table(pca.matrix[ ,'clust']), kmeans.res$size)
cluster.mids <- kmeans.res$centers[ord, ]

# index of species closest to the cluster centre
ex.sps.inds <- mapply(get.clst.mid, 
                      ind.coords = clst.list, 
                      cluster.mids = split(cluster.mids, row(cluster.mids)))

# get the coordinates of these species
sp.example.pca.coords <- mapply(function(x,y) x[y, ], 
                                x = clst.list, 
                                y = ex.sps.inds)

# What are these species? 
# Get them from the full pca matrix
pca.matrix.clsts <- split(pca.matrix, pca.matrix[ ,'clust'])
example.sp <- do.call(rbind, (mapply(function(x, y) x[y, ], 
                                     x = pca.matrix.clsts, 
                                     y = ex.sps.inds, 
                                     SIMPLIFY = FALSE)))

save(example.sp, cluster.mids, sp.example.pca.coords, 
     file = 'Analysis/Outputs/example_species.RData')



