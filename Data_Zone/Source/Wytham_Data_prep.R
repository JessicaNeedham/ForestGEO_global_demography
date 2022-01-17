##############################################
rm(list = ls())
# LIBRARIES #
library(foreach)
library(doParallel)
library(SpatialTools)
library(plyr)
library(stringi)
library(dplyr)
library(lubridate)

##############################################################################
# User defined variables      
site.name <- 'Wytham'  # CTFS site name used in r table 
censuses <- c(1,2,3)      # which censuses to use
n.census <- length(censuses) # number of census to include in analysis
unit.dbh <-   "mm"             # input the units for dbh either "cm" or "mm"
dbh.factor <- ifelse(unit.dbh=="cm", 10, 1)   # converstion factor for dbh
nproc <- 20  # no. cores to dedicate to running this script
registerDoParallel(nproc)
# End user defined portion of script
##############################################################################

### FUNCTIONS ###
source('Data_Zone/Source/CTFS_data_processing_functions.R')

### DATA ###
# 2008 - 2010 censuses
data.mat <-read.csv(file = sprintf("Data_Zone/Data/%s/WythamCompleteCensuses0810.editednumbers.csv",
                                   site.name), header=TRUE,stringsAsFactors = FALSE)
# 2016 census
data.16 <- read.csv(file = sprintf('Data_Zone/Data/%s/Wytham_2016.csv', site.name))
data.16[ ,'ly'] <- as.numeric(as.character(data.16[ ,'ly']))
# trim down to the columns worth keeping
data.mat <- data.mat[ ,1:21]

### Merge the two together based on stem tags.

# stem tag system in wytham is confusing. each unique stem is identified 
# by tag, stemtag and stemtag2. merge files based on these 
data.mat <- merge(data.mat, data.16, 
                  by = c('tag', 'stemtag', 'stemtag2'), all = TRUE)
rm(data.16)
# Since 2016 has most uptodate coordinates, use metadata from 2016 and 
# just add dates, dbhs, codes from earlier censuses
data.mat$lx <- ifelse(is.na(data.mat$lx.y), data.mat$lx.x, data.mat$lx.y)
data.mat$ly <- ifelse(is.na(data.mat$ly.y), data.mat$ly.x, data.mat$ly.y)
data.mat$spcode.y <- as.character(data.mat$spcode.y)
data.mat$spcode <- ifelse(is.na(data.mat$spcode.y), data.mat$spcode.x, data.mat$spcode.y)
data.mat$date <- as.character(data.mat$date)
data.mat$subplot <- ifelse(is.na(data.mat$subplot.y), data.mat$subplot.x, 
                           data.mat$subplot.y)

data.mat <- data.mat[ ,!names(data.mat) %in% c('lx.y', 'ly.y', 'lx..old.', 'ly..old.', 
                                               'lx.x', 'ly.x', 'plotnum.x', 'plotnum.y', 
                                               'spcode.x', 'spcode.y', 'subplotX', 'subplotY',
                                               'notes.16', 'plotX', 'plotY', "subplot.x", 'subplot.y')]

############## 
### choose species 
#sp.list <- split(data.mat, data.mat$spcode)
#tmp <- lapply(sp.list, function(x) apply(x, 2, function(z) 
#  length(which(!is.na(z))))[grep('dbh', colnames(x))])
#keep <- which(lapply(tmp, function(x) ifelse(any(x > 20), 1, 0))
#              == 1)
#data.mat <- data.mat[data.mat$spcode %in% names(keep) == TRUE, ]
data.mat$spcode <- factor(data.mat$spcode)
###################################################################
#####################################################################
# COORDINATES
qdrts <- read.table(file = sprintf("Data_Zone/Data/%s/Wytham_2010quadrat.txt", 
                                   site.name), header = TRUE)
gx <- gy <- rep(NA, times = length(data.mat$tag))

for (ii in 1:length(data.mat$tag)) {
  tmp <- match(data.mat$quadrat[ii], qdrts$quadrat)
  gx[ii] <- qdrts$startx[tmp] + data.mat$lx[ii]
  gy[ii] <- qdrts$starty[tmp] + data.mat$ly[ii]
}

data.mat$gx <- gx
data.mat$gy <- gy

# DEALING WITH FALLEN TREES AND MAKING A SURVIVAL VARIABLE
size.mat <- data.mat[ ,grep('dbh', colnames(data.mat))]
code.mat <- data.mat[ ,grep('codes', colnames(data.mat))]
for(i in 1:ncol(size.mat)){
  size.mat[ ,i] <- ifelse(code.mat[ ,i] == 'D' | code.mat[ ,i] == 'R',
                          NA, size.mat[ ,i])
  size.mat[which(size.mat[,i] == 0) ,i] <- NA
}

# survival based on NAs in sizes
surv.mat <- matrix(NA, nrow(size.mat), ncol(size.mat)-1)
for(i in 1:n.census-1){
  surv.mat[ ,i] <- ifelse(is.na(size.mat[ ,i+1]) & 
                            !is.na(size.mat[ ,i]), 0, NA)
  
  surv.mat[ ,i] <- ifelse(!is.na(size.mat[ ,i+1]) &
                            !is.na(size.mat[ ,i]), 1, 
                          surv.mat[ ,i])
}

tmp <- cbind(size.mat, surv.mat)
tmpll <- tmp
tmpll[,1] <- ifelse(is.na(tmp[,1]), NA, 5)
tmpll[ ,2] <- ifelse(is.na(tmp[,2]), NA, 5)
tmpll[,3] <- ifelse(is.na(tmp[,3]), NA, 5)
unique(tmpll)

# get dates and times
# but first get years - for graphs etc
years <- apply(data.mat[, grep('date', colnames(data.mat))], 2,
               function(x) stri_sub(x, 1, 4))
years <- apply(years, 2, table)
years <- lapply(years, function(x) names(which(x > 200)))
years <- as.character(unlist(lapply(years, function(x) 
  ifelse(length(x)>1, paste(x[1], x[2], sep = '/'), x))))

dates <- data.mat[ ,grep('date', colnames(data.mat))]
# give trees missing dates a random date from that year
dates <- apply(dates, 2, function(x){ifelse(is.na(x), sample(na.omit(x), 1), x)})
# but new recruits get an NA
dates <- apply(dates, 2, function(x) ifelse(x == '', NA, x))
dates <- apply(dates, 2, function(x) sapply(x, as.Date, format = '%Y-%m-%d'))
data.mat[ ,grep('date', colnames(data.mat))] <- dates

colnames(data.mat)[grep('date', colnames(data.mat))] <- 
  sapply(seq(n.census), function(x) paste('date', x, sep = '.'))

# convert dates to time in years since first measurment
start.d <- min(data.mat$date.1, na.rm=TRUE)
for(i in 1:n.census){
  data.mat[ ,grep(paste('date', i, sep = '.'), colnames(data.mat))] <- 
    (data.mat[ ,grep(paste('date', i, sep ='.'), colnames(data.mat))] - start.d)/365.25
}

# get time in years between censuses 
times <- data.mat[ ,grep('date', colnames(data.mat))]
times <- t(apply(times, 1, diff))

if(n.census == 2){
  times <- t(times)
}
colnames(times) <- paste('time', seq(n.census-1), sep = '.')

# make spcode numeric so that the whole thing can be a matrix - 
# much faster for most operations
sp.names <- levels(data.mat$spcode)
spcode <- as.numeric(data.mat$spcode)
data <- cbind(data.mat$tag, data.mat$stemtag,data.mat$stemtag2, spcode,
              data.mat$gx, data.mat$gy, 
              size.mat, surv.mat, times)
colnames(data) <- c('tag', 'stemtag', 'stemtag2', 'spcode', 'x', 'y', 
                    'size.1', 'size.2', 'size.3', 'surv.1', 'surv.2', 
                    'time.1', 'time.2')

# get tag system in order 
# new stems in 2016 - multistems
# get same tag as tree and stemtag is 1 more than previous highest

# new stems in 2016 - new trees
# tag and stemtag same as nearby tree but stemtag2 is A or next letter

# unique coords
# assuming that multistems have same coords and no two trees have same coords...
data$unicoord <- paste(data$x, data$y, data$spcode, sep = '.')
# split by coordinates and species
tree.split <- split(data, f = data$unicoord)
# give each stem a treeID
tree.split <- mapply(function(x,y) cbind(x,y),
                     x = tree.split, y = as.list(seq(length(tree.split))), 
                     SIMPLIFY = FALSE)
# reorder the stemtag column so that largest stem is 1 etc
tree.split <- lapply(tree.split, function(x) 
  cbind(x, rev(order(x$size.1, na.last = FALSE))))

## NB there are 36 cases where species and coordinates are the same but 
# two stems had 'main' in code08. These are going to be called multistems
# because it is impossible to assign stems to each individual based on the
# coding system
data.mat <- do.call(rbind, tree.split)
colnames(data.mat) <- c('tag', 'stemtag', 'stemtag2', 'spcode', 'x', 'y', 
                        'dbh.1', 'dbh.2', 'dbh.3', 'surv.1', 'surv.2', 
                        'time.1', 'time.2', 'unicoord', 'treeID', 'stem.tag')

data.mat <- data.mat[ ,!names(data.mat) %in% c('tag', 'stemtag',
                                               'stemtag2','unicoord')]

# SUBPLOTS
data.mat$sub <- rep(NA, nrow(data.mat))
counter <- 1
for(j in 1:30){
  for(i in 1:15){
    tmp <- which((data.mat$x >= i*20-20) & (data.mat$x < i*20) & 
                   (data.mat$y >= j*20 - 20) & (data.mat$y < j*20))
    data.mat$sub[tmp] <- counter
    counter <- counter+1
  }
}

# get rid of trees outside of the 18ha plot. 
data.mat <- data.mat[!is.na(data.mat$sub),]

# keep only the most recent census
data.mat <- data.mat[ ,-grep('dbh.1', colnames(data.mat))]
data.mat <- data.mat[ ,-grep('surv.1', colnames(data.mat))]
data.mat <- data.mat[ ,-grep('time.1', colnames(data.mat))]


############################################################
# Now add wsg info 
data.mat$spcode <- sp.names[data.mat$spcode]
# add wsg info 
# save the output
sp.latins <- c('Acer campestre', 'Acer pseudoplatanus', 'Betula sp.', 'Carpinus betulus',
               'Castanea sativa', 'Cornus sanguinea', 
               'Corylus avellana', 'Crataegus monogyna', 'Euonymus europaeus',
               'Fagus sylvatica', 'Franal', 'Fraxinus excelsior', 'Ilex aquifolium', 
               'Larix decidua', 'Picea abies', 'Pinus sylvestris', 'Prunus avium', 'Prunus sp.',
               'Pseudotsuga menziesii',
               'Quercus robur', 'Salix sp.', 'Sambucus nigra', 'Taxus baccata', 'Tilia sp.', 
               'Ulmus minor')
sp.names <- sp.names[-1]
tmp <- as.data.frame(cbind(sp.latins, sp.names))
tmp$Genus <- word(tmp$sp.latins)
colnames(tmp) <- c('Latin', 'sp', 'Genus')
sp.mat <- tmp
all.sp <- data.mat
all.sp$Genus <- sp.mat[match(all.sp$sp, sp.mat$sp) ,'Genus']
all.sp$Latin <- sp.mat[match(all.sp$sp, sp.mat$sp), 'Latin']

all.sp$wsg <- NA
all.sp$wsg.level <- NA

# fill in some of the others with genus level info
wsg.mat <- read.csv('Data_Zone/Data/GlobalWoodDensityDatabase.csv')
colnames(wsg.mat)[4] <- 'WSG'
all.sp$Family <- wsg.mat[match(all.sp$Latin, wsg.mat$Binomial) ,'Family']
sp.mat$Family <- wsg.mat[match(sp.mat$Latin, wsg.mat$Binomial), 'Family']
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
data.mat <- data.mat[ ,-grep('sub', colnames(data.mat))]
colnames(data.mat) <- c("sp", 'gx', 'gy', 'dbh.1', 'dbh.2', 'surv.1', 'time.1', 
                        'treeID', 'stem.tag', 'Genus', 'Latin', 'wsg', 'wsg.level', 'Family')
sizes <- data.mat[ ,grep('dbh', colnames(data.mat))]
all.dead <- which(is.na(rowMeans(sizes, na.rm = TRUE)))
data.mat <- data.mat[-all.dead, ]
# save
site.name <- 'wytham'
saveRDS(data.mat, 
        file = sprintf('Data_Zone/Output/%s_all.sp.RData', site.name))

save(sp.mat, file = sprintf('Data_Zone/Data/%s/%s.spptable.rdata',
             site.name, site.name))

#####################################################################
##### load and format the species names
load(sprintf('Data_Zone/Data/%s/%s.spptable.rdata',
             site.name, site.name))
tmp <-sp.mat
all.sp <- data.mat

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

keep <- which(lapply(tmp, function(x) ifelse(any(x > 100), 1, 0))
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
time <- data.mat[ ,'time.1']

# survival
survival <- as.matrix(data.mat[ ,grep('urv', colnames(data.mat))],
                      ncol = n.census -1)

# coordinates
gx <- as.numeric(data.mat$gx)
gy <- as.numeric(data.mat$gy)
# make it a matrix - much quicker to process
sp.id <- as.numeric(data.mat$sp)  # first make species a numeric
data.mat[ ,grep('tag', colnames(data.mat), ignore.case = TRUE)] <- 
  as.numeric(data.mat[ ,grep('tag', colnames(data.mat), ignore.case = TRUE)])
wsg <- data.mat[ ,'wsg']
wsg.level <- factor(data.mat$wsg.level, levels = c('species','genus','family','global'))
wsg.level <- as.numeric(wsg.level)

data.mat <- as.matrix(cbind(as.matrix(data.mat[ ,c('treeID', 'stem.tag')]), 
                            sp.id, gx, gy, sizes, survival, incr, time, wsg, wsg.level))
colnames(data.mat) <- c('treeID','stemTag', 'sp.id', 'gx', 'gy', 'dbh.1', 'dbh.2', 
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

give.up <- 10  #  attempts at getting positive increments 

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

#### load and format the species names
load(sprintf('Data_Zone/Data/%s/%s.spptable.rdata',
             site.name, site.name))
tmp <- sp.mat
# match the species codes from the chosen species iwth the codes from 
# species tables
tmp <- as.data.frame(tmp)

sps <- match(sp.names, tmp$sp)
sp.latins <- tmp[sps, ]

rownames(data.mat) <- NULL

# save the output
data.ls <- list(data.mat = data.mat,
                sp.names = sp.names,
                sp.latins = sp.latins,
                years = years)
saveRDS(data.ls, 
        file = sprintf('Data_Zone/Output/%s.RData', site))

