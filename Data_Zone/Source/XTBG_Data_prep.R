##############################################
rm(list = ls())
# LIBRARIES #
library(foreach)
library(doParallel)
library(stringr)
library(lubridate)
##############################################################################
# User defined variables      
site.name <- 'XTBG'  # CTFS site name used in r table 
censuses <- c(1,2)      # which censuses to use
n.census <- length(censuses) # number of census to include in analysis
unit.dbh <-   "mm"             # input the units for dbh either "cm" or "mm"
dbh.factor <- ifelse(unit.dbh=="cm", 10, 1)   # converstion factor for dbh
nproc <- 1  # no. cores to dedicate to running this script
registerDoParallel(nproc)
# End user defined portion of script
##############################################################################
set.seed(1)
# load in all censuses
cns <- read.csv(file = sprintf('Data_Zone/Data/%s/%s_12.csv', 
                               site.name, site.name), 
                stringsAsFactors = FALSE)

### COORDINATES ###
# four stems have wrong coordinates - make the gx and gy equal x and y
cns$gx <- ifelse(cns$wrong.coordinate == 1, cns$x, cns$gx)
cns$gy <- ifelse(cns$wrong.coordinate == 1, cns$y, cns$gy)

# now just get the columns of interest -
# note in this data set tag is treeID (same for multistems) 
# and stem refers to multistems 
cnsdata <- cns[ ,c('tag', 'tag2', 'stem', 'sp07', 'sp12', 'gx', 'gy',  
                   'dbh07', 'dbh12', 'dbh15', 'date07', 
                   'date12', 'code07', 'code12')]
rm(cns)

### TAXONOMY ###
# how many changed species between first two censuses?
length(which(cnsdata$sp07 != cnsdata$sp12))  # arghh quite a lot - 
changed <- which(cnsdata$sp07 != cnsdata$sp12)
sp.change.mat <- cbind(cnsdata$sp07[changed], cnsdata$sp12[changed])
unique(sp.change.mat)  # not just a few species that changed taxonomy
# assume that the species in 2012 is the correct one - unless 2012 is 
# unidentified and 07 isn't
sp <- ifelse(cnsdata$sp12 == 'UN', cnsdata$sp07, cnsdata$sp12)
cnsdata <- cbind(cnsdata[ ,-grep('sp', colnames(cnsdata))], sp)

# load the species info 
sp.info <- read.csv(sprintf('Data_Zone/Data/%s/XTBG_splist_2007.csv', 
                            site.name),
                    stringsAsFactors = FALSE)

tmp <- sp.info
rm(sp.info)
colnames(tmp) <- c('order', 'Family', 'Genus', 'Species', 'sp')
site.name
save(tmp, file = sprintf('Data_Zone/Data/%s/%s.spptable.RData', site.name, site.name))

### SURVIVAL AND GROWTH ### 
# whats going on with dbh15?
length(which(!is.na(cnsdata[,'dbh15'])))
cnsdata[which(!is.na(cnsdata[,'dbh15'])), ]
# ok - looks like these trees were missing in 2012 but got remeasured in 2015
# move the dbh15 measurement to the 2012 column and make the date12 2015
cnsdata[ ,'date12'] <- ifelse(!is.na(cnsdata[ ,'dbh15']),
                              '15/06/2015', 
                              cnsdata[,'date12'])
cnsdata[ ,'dbh12'] <- ifelse(!is.na(cnsdata[ ,'dbh15']), 
                             cnsdata[ ,'dbh15'], 
                             cnsdata[ ,'dbh12'])
cnsdata <- cnsdata[ ,-grep('dbh15', colnames(cnsdata))]

# Column names
colnames(cnsdata) <- c('tag', 'tag2', 'stem', 'gx', 'gy', 'dbh.1', 'dbh.2',
                       'date.1', 'date.2', 'code.1', 'code.2', 'sp')

# try and work out what is going on with the tags
multis <- cnsdata[cnsdata$stem > 0, ]
head(multis)
# ok stem refers to which stem it is in a multistem. stem tags don't change 
# between censuses - this is equivalent of a stem table from CTFS database

# convert all dbhs to mm units
cnsdata[ ,grep('dbh', colnames(cnsdata))] <- apply(cnsdata[ ,grep('dbh', colnames(cnsdata))], 
                                                   2, function(x) x*dbh.factor)


cnsdata$sp <- factor(cnsdata$sp)
length(levels(cnsdata$sp))  # two species lost

### DATES AND YEARS ###
dates <- cnsdata[ ,grep('date', colnames(cnsdata))]
# replace '' and NAs with the most common date for that year
dates <- apply(dates, 2, function(x) 
  ifelse(x == '' | is.na(x), 
         names(which.max(table(x))), x))

# some of the dates in the first census are in a different format
dates[ ,1] <- ifelse(nchar(dates[ ,1]) !=10, 
                     as.character(as.Date(as.numeric(dates[ ,1]), origin  = '1904-01-01')),
                     as.character(dmy(dates[ ,1])))

wrong <- which(str_sub(dates[ ,1], 1, 4) == '1900')
dates[wrong, 1] <- names(which.max(table(dates[,1])))

dates[ ,2] <- as.character(dmy(dates[ ,2]))
cnsdata[ ,grep('date', colnames(cnsdata))] <- dates

# years - for graphs etc
years <- apply(cnsdata[, grep('date', colnames(cnsdata))], 2,
               function(x) str_sub(x, 1,4))
years <- apply(years, 2, table)
years <- lapply(years, function(x) names(which(x > 1000)))
years <- as.character(unlist(lapply(years, function(x) 
  ifelse(length(x)>1, paste(x[1], x[length(x)], sep = ':'), x))))


# convert dates to time in years since first measurment
cnsdata[ ,grep('date', colnames(cnsdata))] <- apply(dates, 2, 
                                                    function(x)
                                                      ymd(x))
start.d <- min(cnsdata$date.1, na.rm=TRUE)
for(i in 1:n.census){
  cnsdata[ ,grep(paste('date', i, sep = '.'), colnames(cnsdata))] <- 
    (cnsdata[ ,grep(paste('date', i, sep ='.'), colnames(cnsdata))] - start.d)/365.25
}


# survival 
dbhs <- cnsdata[ ,grep('dbh', colnames(cnsdata))]

for(i in 1:(n.census-1)){
  S <- ifelse(!is.na(dbhs[,i]) & !is.na(dbhs[,i+1]), 1, 
              ifelse(!is.na(dbhs[,i]) & is.na(dbhs[,i+1]), 0, NA))
  cnsdata[,paste0("surv.",i)] <- S
}

all.sp <- cnsdata
rm(cnsdata)


# Add wsg info 
site.name
load(sprintf('Data_Zone/Data/%s/%s.spptable.RData',
             site.name, site.name))
tmp$Latin <- paste(tmp$Genus, tmp$Species, sep = ' ')
sp.mat <- tmp
all.sp$Genus <- sp.mat[match(all.sp$sp, sp.mat$sp) ,'Genus']
all.sp$Family <- sp.mat[match(all.sp$sp, sp.mat$sp), 'Family']
all.sp$Latin <- sp.mat[match(all.sp$sp, sp.mat$sp), 'Latin']

all.sp$wsg <- NA
all.sp$wsg.level <- NA

# fill in some of the others with genus level info
wsg.mat <- read.csv('Data_Zone/Data/GlobalWoodDensityDatabase.csv')
colnames(wsg.mat)[4] <- 'WSG'

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

data.mat <- data.mat[-which(data.mat$dbh.2 > 5000), ]

# save
site.name <- 'xtbg'
saveRDS(data.mat, 
        file = sprintf('Data_Zone/Output/%s_all.sp.RData', site.name))

site.name <- 'XTBG'

### REMOVE RECRUITS ###
#  since we only have two censuses we can only use trees
# present in census 1 - so remove recruits
data.mat <- data.mat[-which(is.na(data.mat$dbh.1)), ]
###########################################################################
# these families are only palms and tree ferns
exclude.fam <- c("Arecaceae", "Thyrsopteridaceae", "Loxsomataceae",
                 "Culcitatceae", 
                 "Plagiohyriaceae", "Cibotiaceae", "Cyatheaceae",
                 "Dicksoniaceae",
                 "Meetaxyaceae", "Heliconiaceae", 'Poaceae')

# which species are in the palm/tree fern families
exclude.sp <- sp.mat[which(sp.mat$Family %in% exclude.fam), 'sp']
# remove them
if(length(exclude.sp)>0){
  data.mat <- data.mat[-which(data.mat$sp %in% exclude.sp), ]
}

## CODE SNIPPET TO REMOVE GRAPE VINES, 
## POISON IVY AND UNIDENTIFIED SPECIES
vines <- c('Vitis', 'Unidentified')
rm.sp <- as.character(unlist(sapply(vines, function(x)
  sp.mat[grep(x, sp.mat$Genus), 'sp'])))
poison.ivy <- sp.mat[which(sp.mat$Genus == 'Toxicodendron' &
                      sp.mat$Species == 'radicans'), 'sp']
rm.sp <- c(rm.sp, poison.ivy, 'UN')

if(length(rm.sp)>0){
  data.mat <- data.mat[-which(data.mat$sp %in% rm.sp), ]
}


# now get rid of any species with < 200 individuals in either census
data.mat$sp <- factor(data.mat$sp)  # refactor species codes
# split into species specific matrices
sp.list <- split(data.mat, data.mat$sp)
tmp <- lapply(sp.list, function(x) apply(x, 2, function(z) 
  length(which(!is.na(z))))[grep('dbh', colnames(x))])
keep <- which(lapply(tmp, function(x) ifelse(any(x > 200), 1, 0))
              == 1)
data.mat <- data.mat[data.mat$sp %in% names(keep) == TRUE, ]
data.mat$sp <- factor(data.mat$sp)

# GROWTH CORRECTION #
# for each species get rid of those with individuals more than twice the
# 99th quantile of size. 
source('Data_Zone/Source/CTFS_data_processing_functions.R')

# make it a matrix to speed things up
sizes <- as.matrix(data.mat[ ,grep('dbh', colnames(data.mat))])    
incr <- t(t(apply(sizes, 1, diff)))

# get time in years between censuses 
times <- data.mat[ ,grep('date', colnames(data.mat))]
time <- t(t(apply(times, 1, diff)))

# survival
survival <- as.matrix(data.mat[ ,grep('surv', colnames(data.mat))],
                      ncol = n.census -1)

# coordinates
gx <- as.numeric(data.mat$gx)
gy <- as.numeric(data.mat$gy)
# make it a matrix - much quicker to process
sp.names <- unique(data.mat$sp)
sp.id <- as.numeric(data.mat$sp)  # first make species a numeric
data.mat[ ,grep('tag', colnames(data.mat), ignore.case = TRUE)] <- 
  apply(data.mat[ ,grep('tag', colnames(data.mat), ignore.case = TRUE)], 
        2, function(x) as.numeric(x))
wsg <- data.mat[ ,'wsg']
wsg.level <- factor(data.mat$wsg.level, levels = c('species','genus','family','global'))
wsg.level <- as.numeric(wsg.level)

data.mat <- as.matrix(cbind(as.matrix(data.mat[ ,1:3]), 
                            sp.id, gx, gy, sizes, survival, incr, time, wsg, wsg.level))
colnames(data.mat) <- c('treeID','stag2', 'stem', 'sp.id', 'gx', 'gy', 'dbh.1', 'dbh.2', 
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

#### load and format the species names
site.name
load(sprintf('Data_Zone/Data/%s/%s.spptable.RData',
             site.name, site.name))

# match the species codes from the chosen species with the codes from 
# species tables
tmp <- as.data.frame(tmp)
sps <- match(sp.names, tmp$sp)
tmp$Latin <- paste(tmp$Genus, tmp$Species, sep = ' ')
sp.latins <- tmp[sps, c('Latin', 'sp')]

rownames(data.mat) <- NULL

# save the output
data.ls <- list(data.mat = data.mat,
                sp.names = sp.names,
                sp.latins = sp.latins,
                years = years)
site.name <- 'xtbg'
saveRDS(data.ls, 
        file = sprintf('Data_Zone/Output/%s.RData', site.name))


