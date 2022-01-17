############################################################
### AMACAYACU DATA PREP ###
##########################################################
rm(list = ls())
###################
### LIBRARIES ###
###################
library(foreach)
library(doParallel)
#library(SpatialTools)
#library(stringi)
library(stringr)
library(lubridate)

#####################
### FUNCTIONS ###
#####################
source('Data_Zone/Source/CTFS_data_processing_functions.R')
nproc <- 1  # for parallel processing
registerDoParallel(nproc)
#####################
### DATA ###
#####################
site <- 'amacayacu'
data.mat <- read.csv(sprintf('Data_Zone/Data/%s/%s.csv', 
                             site, site),
                     stringsAsFactors = FALSE)

# only keep relevant columns
data.mat <- data.mat[ ,c('tag', 'StemTag', 
                         'sp', 'Species','Genus', 'Family', 'gx', 'gy', 
                         'date.1', 'date.2', 
                         'dbh.1..cm.', 'dbh.2..cm.')]



### TREE and STEM ID coloumns
data.mat$treeID <- seq(nrow(data.mat))  # unique number
# find main stem of multi stems and give them the same treeID
# first strip - from start of tag column
data.mat$tag <- str_sub(data.mat$tag, 2, 7)

multis <- which(data.mat$StemTag != '')
data.mat$treeID[multis] <- sapply(data.mat$StemTag[multis], function(x)
  data.mat$treeID[match(x, data.mat$tag)])
# now make stemTag to be numbers 1 through to number of stems
tree.split <- split(data.mat, data.mat$treeID)
tree.split <- lapply(tree.split, function(x){
  x$stemTag <- seq(nrow(x)); x})
data.mat <- do.call(rbind, tree.split)

# add wsg info and then save  - for when we do biomass calculations later
wsg.mat <- read.csv('Data_Zone/Data/GlobalWoodDensityDatabase.csv')
colnames(wsg.mat)[4] <- 'WSG'
# there are multiple entries for some species  - take the means
wsgs <- tapply(wsg.mat[ ,'WSG'], list(wsg.mat[ ,'Binomial']), mean)
# genus level
wsg.mat$Genus <- word(wsg.mat$Binomial,1,1)
genus.wsgs <- tapply(wsg.mat$WSG, wsg.mat$Genus, mean)
# family level
family.wsgs <- tapply(wsg.mat$WSG, wsg.mat$Family, mean)

data.mat <- cbind(data.mat, NA, NA)
colnames(data.mat)[(ncol(data.mat)-1):ncol(data.mat)] <- c('wsg', 'wsg.level') 

data.mat[ ,'wsg'] <- wsgs[match(data.mat[ ,'Species'], names(wsgs))]
data.mat[ ,'wsg.level'] <- ifelse(is.na(data.mat[ ,'wsg']), NA, 'species')

g.wsgs <- genus.wsgs[sapply(data.mat[ ,'Genus'], function(x) match(x, names(genus.wsgs)))]
f.wsgs <- family.wsgs[sapply(data.mat[ ,'Family'], function(x) match(x, names(family.wsgs)))]

data.mat$wsg <- wsgs[match(data.mat$Species, names(wsgs))]
data.mat$wsg.level <- ifelse(is.na(data.mat$wsg.level) & !is.na(data.mat[ ,'wsg']),
                           'species', data.mat[ ,'wsg.level'])
data.mat$wsg <- ifelse(is.na(data.mat$wsg), g.wsgs, data.mat$wsg)
data.mat$wsg.level <- ifelse(is.na(data.mat$wsg.level) & !is.na(data.mat[ ,'wsg']),
                             'genus', data.mat[ ,'wsg.level'])
data.mat$wsg <- ifelse(is.na(data.mat$wsg), f.wsgs, data.mat$wsg)
data.mat$wsg.level <- ifelse(is.na(data.mat$wsg.level) & !is.na(data.mat[ ,'wsg']),
                             'family', data.mat[ ,'wsg.level'])

# and global for everything else
global <- mean(wsgs, na.rm = TRUE)
data.mat$wsg <- ifelse(is.na(data.mat$wsg), global, data.mat$wsg)
data.mat$wsg.level <- ifelse(is.na(data.mat$wsg.level),
                             'global', data.mat[ ,'wsg.level'])


### DBH 
# convert to mm
data.mat[ ,grep('dbh', colnames(data.mat))] <- 
  data.mat[ ,grep('dbh', colnames(data.mat))]*10

### DATES
# get them in the same format
data.mat$date.2 <- ymd(data.mat$date.2)
data.mat$date.1 <- dmy(data.mat$date.1)

# format years = for figure legends etc.
dates <- data.mat[ ,grep('date', colnames(data.mat))]
tmpdates <- apply(dates, 2, as.character)
years <- apply(tmpdates, 2, 
               function(x) str_sub(x, 1,4))
years <- lapply(split(years, col(years)), table)
years <- lapply(years, function(x) names(which(x > 1000)))
years <- as.character(unlist(lapply(years, function(x) 
  ifelse(length(x)>1, paste(x[1], x[length(x)], sep = ':'), x))))

start.d <- min(dates[ ,1], na.rm=TRUE)
data.mat[ ,grep('date', colnames(data.mat))] <- apply(dates, 2, function(x) 
  difftime(x, start.d, units = 'days')/365.25)

# make a survival vector based on NAs in the data
data.mat$surv <- ifelse(is.na(data.mat$dbh.2..cm.), 0, 1)

data.mat <- data.mat[ ,3:ncol(data.mat)]

colnames(data.mat) <- c('sp', 'Latin', 'Family', 'gx', 'gy', 'date.1', 'date.2', 'dbh.1', 'dbh.2', 
                        'treeID', 'stemTag', 'wsg', 'wsg.level', 'Genus', 'surv')

data.mat <- data.mat[ ,c('treeID', 'stemTag', 'sp', 'gx', 'gy', 'dbh.1','date.1','dbh.2', 'date.2', 
                      'surv', 'wsg', 'wsg.level', 'Latin', 'Family')]

# save this complete version
site.name <- 'amacayacu'
saveRDS(data.mat, 
        file = sprintf('Data_Zone/Output/%s_all.sp.RData', site.name))

# remove stems not present in census 1 - can't be used for 
# size-dependent vital rates
data.mat <- data.mat[!is.na(data.mat$dbh.1), ]

# remove species from tree fern or palm families
exclude.fam <- c("Arecaceae", "Thyrsopteridaceae", "Loxsomataceae",
                 "Culcitatceae", 
                 "Plagiohyriaceae", "Cibotiaceae", "Cyatheaceae",
                 "Dicksoniaceae",
                 "Meetaxyaceae", "Heliconiaceae", 
                 'Poaceae')

data.mat <- data.mat[-which(data.mat$Family %in% exclude.fam), ]

### SPECIES IDENTITIES - 
data.mat$sp.id <- toupper(data.mat$sp)

# Deal with the 'unknowns' 
# first we can get rid of things with no information at all
length(which(data.mat$sp.id == 'UNIDUNID'))
data.mat <- data.mat[-which(data.mat$sp.id == 'UNIDUNID'), ]
data.mat <- data.mat[-which(is.na(data.mat$sp.id)), ]
# pull out morpho species
morphs <- data.mat[which(grepl('\\d', data.mat$sp.id) == TRUE), ]
# temporarily remove from main matrix
data.mat <- data.mat[-which(grepl('\\d', data.mat$sp.id) == TRUE), ]

# now sort the species latins for these species - 
# replace 'unidentified' with genus and number
morphs[ ,'sp.id'] <- paste(str_sub(morphs$sp.id, 1, 4), 'SP', 
                   str_sub(morphs$sp, -2, -1), 
                    sep = '.')
morphs[ ,'Latin'] <- gsub(',', '.', morphs[ ,'Latin'])

# put the cleaned up morpho species back in 
data.mat <- rbind(data.mat, morphs)

# Filter by families
prob.families <- which(is.na(data.mat[ ,'Family']) |
                         data.mat[ ,'Family'] == 'Unknown' |
                         data.mat[ ,'Family'] == '')

if(length(prob.families)>=1){
  data.mat <- data.mat[-prob.families, ]
}

# get rid of species with no species identification - not even morpho species
data.mat <- data.mat[which(data.mat[ ,'sp.id'] != ''), ]

# some checks
unique(data.mat[ ,'Family'])
unique(data.mat[ ,'sp.id'])
unique(data.mat[ ,'Latin'])

# now make a species latin look up table
sp.latins <- sort(unique(data.mat$sp.id))
latin <- sapply(sp.latins, function(x) 
  data.mat[match(x, data.mat$sp.id), 'Latin'])
sp.latins <- cbind(latin, sp.latins)

sp.names <- sp.latins[,2]
sp.latins <- as.data.frame(sp.latins)
rownames(sp.latins) <- NULL
colnames(sp.latins) <- NULL

# keep only species with at least 200 individuals in census 1
keep <- names(which(table(data.mat$sp.id) > 200))
keep.ind <- which(table(data.mat$sp.id) > 200)
data.mat <- data.mat[data.mat$sp.id %in% keep, ]
sp.names <- sp.names[keep.ind]
sp.latins <- sp.latins[keep.ind, ]

# change colname order
data.mat <- data.mat[ ,c('treeID', 'stemTag', 'sp.id', 'gx', 'gy', 
                         'date.1', 'date.2', 'dbh.1', 'dbh.2',
                         'surv', 'wsg', 'wsg.level')]

sp.latins <- cbind(sp.latins, NA, NA)
sp.latins[ ,3] <- as.numeric(data.mat[match(sp.latins[ ,2], data.mat[ ,'sp.id']) ,'wsg'])
sp.latins[ ,4] <- data.mat[match(sp.latins[ ,2], data.mat[ ,'sp.id']), 'wsg.level']

data.mat$sp.id <- as.numeric(as.factor(data.mat$sp.id))
wsg <- data.mat[ ,'wsg']
wsg.level <- factor(data.mat$wsg.level, levels = c('species','genus','family','global'))
wsg.level <- as.numeric(wsg.level)

sizes <- as.matrix(data.mat[ ,grep('dbh', colnames(data.mat))])   
incr <- apply(sizes, 1, diff)  # NOT ANNUAL
# get time in years between censuses 
times <- data.mat[ ,grep('date', colnames(data.mat))]
time <- apply(times, 1, diff)

data.mat <- cbind(data.mat[ ,c('treeID', 'stemTag', 'sp.id', 'gx', 'gy', 'dbh.1', 'dbh.2', 'surv')], 
                  incr, time, wsg, wsg.level)

# finding errors and correcting them  - 
# this is based on annual increments- stems that shrank more than 25% of original size 
# or grew more than 75 mm in a single year
pc.change <- sizes[ ,2]/sizes[ ,1]*100
data.mat <- data.mat[-which(pc.change < 75), ]

annual.incr <- data.mat[ ,'incr']/data.mat[ ,'time']
if(length(which(annual.incr > 75))){
  data.mat <- data.mat[-which(annual.incr > 75), ]
}

# - propose new sizes so that size monotonically increases
give.up <- 1000#  attempts at getting positive increments 

dbhs <- as.matrix(data.mat[ ,grep('dbh', 
                        colnames(data.mat))])  # matrix of sizes
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

out <- foreach(proc = 1:nproc)%dopar%{
  set.seed(1)
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
incrs.lll <- apply(med.corrections, 1, diff)

# replace the increment columns with the new positive increments
data.mat[ ,grep('incr', colnames(data.mat))] <- incrs.lll

# make the whole thing a matrix 
data.mat <- as.matrix(data.mat)

# save the output
data.ls <- list(data.mat = data.mat,
                sp.names = sp.names,
                sp.latins = sp.latins,
                years = years)
saveRDS(data.ls, 
        file = sprintf('Data_Zone/Output/%s.RData', site))




