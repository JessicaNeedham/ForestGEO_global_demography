################################################
################################################
rm(list = ls())
# LIBRARIES #
library(foreach)
library(doParallel)
library(stringr)
library(lubridate)
##############################################################################
# User defined variables      
site.name <- 'windriver'  # CTFS site name used in r table 
censuses <- c(1,2)      # which censuses to use
n.census <- length(censuses) # number of census to include in analysis
unit.dbh <-   "cm"             # input the units for dbh either "cm" or "mm"
dbh.factor <- ifelse(unit.dbh=="cm", 10, 1)   # converstion factor for dbh
nproc <- 1  # no. cores to dedicate to running this script
registerDoParallel(nproc)
source('Data_Zone/Source/CTFS_data_processing_functions.R')

##############################################################################
set.seed(1)
# load in all censuses - in wind river there are those from 2011, 
# living trees in 2016 and dead trees 2016

# original census
cens.1 <- read.csv('Data_Zone/Data/windriver/WFDP_Live2011.csv', 
                  stringsAsFactors = FALSE)
# survived
cens.2.lived <- read.csv('Data_Zone/Data/windriver/WFDP_Live2011_Live2016.csv', 
                         stringsAsFactors = FALSE)
# died
#cens.2.died <- read.csv('Data_Zone/Data/windriver/WFDP_Live2011_Dead2016.csv', 
 #                      stringsAsFactors = FALSE)
# don't need to load this as the ones not alive can get DBH NA and 0 for survival

cns.list <- list(cens.1, cens.2.lived)

### DBH ###
# convert DBH to mm (but first make it numeric)
cns.list <- lapply(cns.list, function(x){
  ind <- grep('DBH', colnames(x))
  if(length(ind) == 1){
  tmp <- matrix(x[ ,ind], nrow = nrow(x))
  }else{
    tmp <- x[ ,ind]
  };
  tmp <- apply(tmp, 2, function(y) as.numeric(y)*dbh.factor);
  x[ ,ind] <- tmp;
  return(x)
})

### SURVIVAL COLUMNS ###
# This has effectively been done by the way the data were sent
cns.list[[2]]$surv.1 <- 1

# now merge them into a single dataframe
colnames(cns.list[[2]])[2] <- 'species'  # ensure we can merge on same colnames
data.mat <- merge(cns.list[[1]], cns.list[[2]], 
                  by = c('STEM_TAG', 'species', 'PLOT_X', 'PLOT_Y'), 
                  all.x = TRUE)
# put 0 in survival for those without a 1
data.mat$surv.1 <- ifelse(is.na(data.mat$surv.1), 0, 1)

# remove unnecessary columns
data.mat <- data.mat[ ,!names(data.mat) %in% 
                        c('BIOMASS_DBH1', 'X', 'X.1', 'X.2', 'X.3')]

### STEMTAGS ###
# give each stem a column corresponding to species and coordinates
data.mat$uniqID <-  paste(as.character(data.mat$PLOT_X), 
                          as.character(data.mat$PLOT_Y), 
                          data.mat$species, 
                          sep = '.')

# now give each stem a tree ID 
data.mat$treeID <- rep(NA, nrow(data.mat))
# unique treeID for first of each unique coord/species combo
data.mat$treeID[which(duplicated(data.mat$uniqID) == FALSE)] <- 
  seq(length(which(duplicated(data.mat$uniqID) == FALSE)))
# now match coord/species combos for multistems
data.mat$treeID[which(duplicated(data.mat$uniqID) == TRUE)] <- 
  data.mat$treeID[sapply(data.mat$uniqID[which(duplicated(data.mat$uniqID) == TRUE)],
                         function(y){match(y, data.mat$uniqID)})]

# remove unnecessary columns
data.mat <- data.mat[ ,-which(colnames(data.mat) == 'uniqID')]
# split by tree ID
tree.split <- split(data.mat, data.mat$treeID)
# each stem gets a unique number
tree.split <- lapply(tree.split, function(y){
  y$stemTag <- seq(nrow(y)); y
})
data.mat <- do.call(rbind, tree.split)

# Housekeeping
data.mat <- data.mat[ ,-which(colnames(data.mat) == 'STEM_TAG')]
colnames(data.mat) <- c('sp', 'gx', 'gy', 'dbh.1', 'dbh.2', 'date.1', 'date.2', 
                        'adjusted', 'surv.1', 'treeID', 'stemTag')
data.mat <- data.mat[ ,c('treeID', 'stemTag', 'sp', 'gx', 'gy', 'dbh.1', 'dbh.2',
                         'date.1', 'date.2', 'surv.1', 'adjusted')]

### DATES etc. ###
# get the years that the censuses were conducteded
# but first get years - for graphs etc
years <- apply(data.mat[, grep('date', colnames(data.mat))], 2,
               function(x) str_sub(x, 6,11))
years <- apply(years, 2, table)
years <- lapply(years, function(x) names(which(x > 5000)))
years <- as.character(unlist(lapply(years, function(x) 
  ifelse(length(x)>1, paste(x[1], x[length(x)], sep = ':'), x))))
years
# sort the dates
# get the years over which the censuses were conducted
dates <- data.mat[ ,grep('date', colnames(data.mat))]
dates <- apply(dates, 2, mdy)

# get the years over which the censuses were conducted
# convert dates to time in years since first measurment
start.d <- min(dates, na.rm = TRUE)
dates <- apply(dates, 2, function(x){(x - start.d)/365})

data.mat[ ,grep('date', colnames(data.mat))] <- dates

# add wsg info 
sp.table <- read.csv('Data_Zone/Data/windriver/wfdp.spptable.csv', 
                     stringsAsFactors = FALSE)
tmp <- sp.table
site.name
colnames(tmp) <- c('sp', 'Genus', 'Species', 'variety', 'Family', 'Authority')
save(tmp, file = sprintf('Data_Zone/Data/%s/%s.spptable.RData', site.name, site.name))

sp.mat <- tmp
sp.mat$Latin <- paste(sp.mat$Genus, sp.mat$Species, sep = ' ')
all.sp <- data.mat
all.sp$Genus <- sp.mat[match(all.sp$sp, sp.mat$sp) ,'Genus']
all.sp$Latin <- sp.mat[match(all.sp$sp, sp.mat$sp), 'Latin']

all.sp$wsg <- NA
all.sp$wsg.level <- NA

# fill in some of the others with genus level info
wsg.mat <- read.csv('Data_Zone/Data/GlobalWoodDensityDatabase.csv')
colnames(wsg.mat)[4] <- 'WSG'
all.sp$Family <- NA
all.sp$Family <- wsg.mat[match(all.sp$Latin, wsg.mat$Binomial), 'Family']

# there are multiple entries for some species  - take the means
wsgs <- tapply(wsg.mat[ ,'WSG'], list(wsg.mat[ ,'Binomial']), mean)
# genus level
wsg.mat$Genus <- word(wsg.mat$Binomial,1,1)
genus.wsgs <- tapply(wsg.mat$WSG, wsg.mat$Genus, mean)
# family level
family.wsgs <- tapply(wsg.mat$WSG, wsg.mat$Family, mean)

g.wsgs <- genus.wsgs[sapply(all.sp[ ,'Genus'], function(x) match(x, names(genus.wsgs)))]
f.wsgs <- family.wsgs[sapply(all.sp[ ,'Family'], function(x) match(x, names(family.wsgs)))]
all.sp$wsg <- wsgs[match(all.sp$Latin, names(wsgs))]
all.sp$wsg.level <- ifelse(is.na(all.sp$wsg.level) & !is.na(all.sp[ ,'wsg']),
                           'species', all.sp[ ,'wsg.level'])
all.sp$wsg <- ifelse(is.na(all.sp$wsg), g.wsgs, all.sp$wsg)
all.sp$wsg.level <- ifelse(is.na(all.sp$wsg.level) & !is.na(all.sp[ ,'wsg']),
                           'genus', all.sp[ ,'wsg.level'])
all.sp$wsg <- ifelse(is.na(all.sp$wsg), f.wsgs, all.sp$wsg)
all.sp$wsg.level <- ifelse(is.na(all.sp$wsg.level) & !is.na(all.sp[ ,'wsg']),
                           'family', all.sp[ ,'wsg.level'])

# and global for everything else
global <- mean(wsgs, na.rm = TRUE)
all.sp$wsg <- ifelse(is.na(all.sp$wsg), global, all.sp$wsg)
all.sp$wsg.level <- ifelse(is.na(all.sp$wsg.level),
                           'global', all.sp[ ,'wsg.level'])

data.mat <- all.sp

data.mat$dbh.1 <- as.numeric(data.mat$dbh.1)

# save
site.name <- 'windriver'
saveRDS(data.mat, 
        file = sprintf('Data_Zone/Output/%s_all.sp.RData', site.name))

## CHOOSE SPECIES ##
#### load and format the species names
load(sprintf('Data_Zone/Data/%s/%s.spptable.RData',
             site.name, site.name))

# exclude tree ferns and palms - these all fall belong 
# to the families listed below
exclude.fam<-c("Arecaceae", "Thyrsopteridaceae", "Loxsomataceae",
               "Culcitatceae", 
               "Plagiohyriaceae", "Cibotiaceae", "Cyatheaceae",
               "Dicksoniaceae",
               "Meetaxyaceae", "Heliconiaceae", 
               'Poaceae')
exclude.sp <- tmp[which(tmp$Family %in% exclude.fam), 'sp']

if(length(exclude.sp)>0){
  all.sp <- all.sp[-which(all.sp$sp %in% exclude.sp), ]
}

## CODE SNIPPET TO REMOVE GRAPE VINES, 
## POISON IVY AND UNIDENTIFIED SPECIES
vines <- c('Vitis', 'Unidentified')
rm.sp <- as.character(unlist(sapply(vines, function(x)
  tmp[grep(x, tmp$Genus), 'sp'])))
poison.ivy <- tmp[which(tmp$Genus == 'Toxicodendron' & 
                          tmp$Species == 'radicans'), 'sp']
rm.sp <- c(rm.sp, poison.ivy)

if(length(rm.sp)>0){
  all.sp <- all.sp[-which(all.sp$sp %in% rm.sp), ]
}

all.sp$sp <- factor(all.sp$sp)  # refactor species codes

# split into species specific matrices
sp.list <- split(all.sp, all.sp$sp)

tmp <- lapply(sp.list, function(x) apply(x, 2, function(z) 
  length(which(!is.na(z))))[grep('dbh', colnames(x))])

keep <- which(lapply(tmp, function(x) ifelse(any(x > 200), 1, 0))
              == 1)

all.sp <- all.sp[all.sp$sp %in% names(keep) == TRUE, ]
all.sp$sp <- factor(all.sp$sp)
# split into species specific matrices
sp.list <- split(all.sp, all.sp$sp)

# vector of names
sp.names <- na.omit(as.character(unlist(lapply(sp.list, function(x) x$sp[1]))))
# save the names and list of species to model
site <- site.name

# GROWTH CORRECTION # 

# for each species get rid of those with individuals more than twice the
# 99th quantile of size. 
source('Data_Zone/Source/CTFS_data_processing_functions.R')

# bind species back together to process all at once
data.mat <- do.call(rbind, sp.list)

# make it a matrix to speed things up
sizes <- as.matrix(data.mat[ ,grep('dbh', colnames(data.mat))])    
incr <- t(t(apply(sizes, 1, diff)))  # NOT ANNUAL

# get time in years between censuses 
times <- data.mat[ ,grep('date', colnames(data.mat))]
time <- t(t(apply(times, 1, diff)))

# survival
survival <- as.matrix(data.mat[ ,grep('surv', colnames(data.mat))],
                      ncol = n.census -1)

# coordinates
gx <- as.numeric(data.mat$x)
gy <- as.numeric(data.mat$y)
# make it a matrix - much quicker to process
sp.id <- as.numeric(data.mat$sp)  # first make species a numeric
data.mat[ ,grep('tag', colnames(data.mat), ignore.case = TRUE)] <- 
 as.numeric(data.mat[ ,grep('tag', colnames(data.mat), ignore.case = TRUE)])
wsg <- data.mat[ ,'wsg']
wsg.level <- factor(data.mat$wsg.level, levels = c('species','genus','family','global'))
wsg.level <- as.numeric(wsg.level)

data.mat <- as.matrix(cbind(as.matrix(data.mat[ ,1:2]), 
                            sp.id, gx, gy, sizes, survival, incr, time, wsg, wsg.level))
colnames(data.mat) <- c('treeID','stemID','sp.id', 'gx', 'gy', 'dbh.1', 'dbh.2', 
                        'surv', 'incr', 'time', 'wsg', 'wsg.level')


# finding errors and correcting them  - 
# this is based on annual increments- stems that shrank more than 25% of original size 
# or grew more than 75 mm in a single year
pc.change <- sizes[ ,2]/sizes[ ,1]*100
data.mat <- data.mat[-which(pc.change < 75), ]
annual.incr <- data.mat[ ,'incr']/data.mat[ ,'time']
if(length(which(annual.incr > 75))){
  data.mat <- data.mat[-which(annual.incr > 75), ]
}


#  Remove stems dead in all censuses
bb <- data.mat[ ,grep('dbh', colnames(data.mat))]  
data.mat <- data.mat[which(apply(bb, 1, 
                                 function(x) all(is.na(x))) 
                           == FALSE), ]

give.up <- 1000  #  attempts at getting positive increments 

dbhs <- data.mat[ ,grep('dbh', 
                        colnames(data.mat))]  # matrix of sizes
times <- as.matrix(data.mat[ ,grep('time', 
                                   colnames(data.mat))], ncol = n.census-1)  # matrix of times
# split dbhs and times into chuncks for parallel processing
chunk.size <- dim(dbhs)[1]/nproc  
d <- seq(dim(dbhs)[1])
e <- split(d, ceiling(seq_along(d)/chunk.size))
if(length(e) > nproc){
  e[[nproc]] <- c(e[[nproc]], e[[length(e)]])
  e[[nproc +1]] <- NULL
}
dbh.list <- lapply(e, function(x, y) y[x, ], y = dbhs)
times.list <- lapply(e, function(x, y) y[x, ], y = times)

out <- foreach(proc = 1:nproc)%dopar% {
  # growth correction code
  resampled.dbhs <- correct.dbh(dbh.list[[proc]], 
                                times = times.list[[proc]], 
                                give.up = give.up)
  
  sampled <- resampled.dbhs$n.sampled
  dd <- resampled.dbhs$dbh.stars
  return(list(med.correct = dd, 
              still.neg = sampled))
}

# number of trees that the algorithm failed on - these trees 
# had increments sampled from positive increments of other trees 
sampled <- unlist(lapply(out, '[[', 2))  
print(sprintf('total sampled: %d', sampled))

# corrected sizes
med.corrections <- do.call(rbind, lapply(out, '[[', 1))  

# replace the sizes with these corrected ones
data.mat[ ,grep('dbh', colnames(data.mat))] <- med.corrections
# make sure they are all still monotonically increasing
incrs.lll <- t(apply(med.corrections, 1, diff))

# replace the increment columns with the new positive increments
data.mat[ ,grep('incr', colnames(data.mat))] <- incrs.lll
site.name <- 'windriver'
#### load and format the species names
load(sprintf('Data_Zone/Data/%s/%s.spptable.RData',
             site.name, site.name))
tmp$Latin <- paste(tmp$Genus, tmp$Species, sep = ' ')
sps <- match(sp.names, tmp$sp)
sp.latins <- cbind(tmp$Latin[sps], tmp$sp[sps])

rownames(data.mat) <- NULL

# save the output
data.ls <- list(data.mat = data.mat,
                sp.names = sp.names,
                sp.latins = sp.latins,
                years = years)
saveRDS(data.ls, 
        file = sprintf('Data_Zone/Output/%s.RData', site))

