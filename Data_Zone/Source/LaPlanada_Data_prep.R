################################################
################################################
rm(list = ls())
# LIBRARIES #
library(foreach)
library(doParallel)
library(stringr)

##############################################################################
# User defined variables      
site.name <- 'laplanada'  # CTFS site name used in r table 
censuses <- c(1,2)      # which censuses to use
n.census <- length(censuses) # number of census to include in analysis
unit.dbh <-   "mm"             # input the units for dbh either "cm" or "mm"
dbh.factor <- ifelse(unit.dbh=="cm", 10, 1)   # converstion factor for dbh
nproc <- 1  # no. cores to dedicate to running this script
registerDoParallel(nproc)

##############################################################################
set.seed(1)
# load in all censuses
for (i in censuses){ 
  load(paste('Data_Zone/Data/', site.name, '/', site.name, ".stem",i,".rdata", sep=""))
}

# put them all in a list
cns <- vector('list', n.census)

for(i in 1:n.census){
  cns[[i]] = eval(as.name(paste( site.name, ".stem",censuses[i], sep="")))
}

# rearrange columns and only keep tags, species, quadrats, coordinates, 
# height of measure, status, date and exact date
cnsdata = list() 

for (i in 1:n.census){
  cnsdata[[i]] <- cns[[i]][,c('treeID', 'stemID', 'tag', 'StemTag', 
                              "sp", "gx", "gy", "dbh", "hom", "status",
                              "date", "ExactDate")]
}

# convert all dbhs to mm units
for(i in 1:n.census){ 
  cnsdata[[i]]$dbh <-cnsdata[[i]]$dbh * dbh.factor
} 

# this appends the census number to time-dependent variables in each census
# (i.e. dbh becomes dbh.1)
for(i in 1:n.census){
  colnames(cnsdata[[i]])[8:ncol(cnsdata[[i]])] <- 
    paste(colnames(cnsdata[[i]])[8:ncol(cnsdata[[i]])],".",i, sep="")
}

# Combine list items into a data frame
# select static information from first census table
indv  <- cnsdata[[1]][,1:7]
temp2 <- lapply(cnsdata,"[", -c(1:7))
# bind dynamic data to static data
tmp <- cbind(indv, do.call(cbind ,(temp2)))
rm(temp2)
rm(indv)

# get the years that the censuses were conducteded
# but first get years - for graphs etc
years <- c('1997:1998', '2002:2003')

# convert dates to time in years since first measurment
start.d <- min(tmp$date.1, na.rm=TRUE)
for(i in 1:n.census){
  tmp[ ,grep(paste('date', i, sep = '.'), colnames(tmp))] <- 
    (tmp[ ,grep(paste('date', i, sep ='.'), colnames(tmp))] - start.d)/365.25
}

# make a survival column - for the stem data this needs to be based on 
# dbh and NOT on status - since status column refers to the whole tree

# any 0 must become an NA
tmp[ ,'dbh.1'] <- ifelse(tmp[ ,'dbh.1'] == 0, NA, tmp[ ,'dbh.1'])
tmp[ ,'dbh.2'] <- ifelse(tmp[ ,'dbh.2'] == 0, NA, tmp[ ,'dbh.2'])


# check to see if there are dbh values in dead trees
which(!is.na(tmp$dbh.1) & tmp$status.1 == 'D')
which(!is.na(tmp$dbh.2) & tmp$status.2 == 'D')

# survival based on NAs in sizes
size.mat <- tmp[ ,grep('dbh', colnames(tmp))]

S <- ifelse(is.na(size.mat[ ,2]) & 
              !is.na(size.mat[ ,1]), 0, NA)
S <-  ifelse(!is.na(size.mat[ ,2]) &
               !is.na(size.mat[ ,1]), 1, S)
tmp[ ,paste0('Surv.', i)] <- S

# get rid of trees with no measurement in any census
sizes <- tmp[ ,grep('dbh', colnames(tmp))]
all.na <- which(apply(sizes, 1, function(x) all(is.na(x))))
if(length(all.na) > 0){
  tmp <- tmp[-all.na, ]
}

all.sp <- tmp[ ,!names(tmp) %in% c(paste0("hom.", seq_along(1:n.census)),
                                   paste0("ExactDate.", seq_along(1:n.census)),
                                   paste0("status.", seq_along(1:n.census)))]

# tidy the workspace
remove(cns)
remove(cnsdata)
remove(tmp)

# add wsg info 
load(sprintf('Data_Zone/Data/%s/%s.spptable.rdata',
             site.name, site.name))
tmp <- eval(parse(text = sprintf('%s.spptable', site.name)))

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
# save
site.name <- 'laplanada'
saveRDS(data.mat, 
        file = sprintf('Data_Zone/Output/%s_all.sp.RData', site.name))

############################
## CHOOSE SPECIES ##
## CHOOSE SPECIES ##
#### load and format the species names
load(sprintf('Data_Zone/Data/%s/%s.spptable.rdata',
             site.name, site.name))
tmp <- eval(parse(text = sprintf('%s.spptable', site.name)))

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

keep <- which(lapply(tmp, function(x) ifelse(any(x > 200), 1, 0)) == 1)

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
survival <- as.matrix(data.mat[ ,grep('Surv', colnames(data.mat))],
                      ncol = n.census -1)

# coordinates
gx <- as.numeric(data.mat$gx)
gy <- as.numeric(data.mat$gy)
# make it a matrix - much quicker to process
sp.id <- as.numeric(data.mat$sp)  # first make species a numeric
data.mat[ ,grep('tag', colnames(data.mat), ignore.case = TRUE)] <- 
  apply(data.mat[ ,grep('tag', colnames(data.mat), ignore.case = TRUE)], 
        2, function(x) as.numeric(x))
wsg <- data.mat[ ,'wsg']
wsg.level <- factor(data.mat$wsg.level, levels = c('species','genus','family','global'))
wsg.level <- as.numeric(wsg.level)

data.mat <- as.matrix(cbind(as.matrix(data.mat[ ,1:4]), 
                            sp.id, gx, gy, sizes, survival, incr, time, wsg, wsg.level))
colnames(data.mat) <- c('treeID','stemID', 'tag', 'StemTag', 'sp.id', 'gx', 'gy', 'dbh.1', 'dbh.2', 
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
load(sprintf('Data_Zone/Data/%s/%s.spptable.rdata',
             site.name, site.name))
tmp <- eval(parse(text = sprintf('%s.spptable', tolower(site.name))))
# match the species codes from the chosen species iwth the codes from 
# species tables
tmp <- as.data.frame(tmp)

sapply(exclude.sp, function(x)
  tmp[match(x, tmp[ ,'sp']), 2])

sps <- match(sp.names, tmp$sp)
sp.latins <- cbind(tmp$Latin[sps], tmp$sp[sps])

na.sp <- sp.latins[is.na(sp.latins[ ,1]), 2]
genus.na.sp <- tmp[match(na.sp, tmp[ ,'sp']) ,'Genus']
sp.latins[is.na(sp.latins[ ,1]) ,1] <- paste0(genus.na.sp, ' sp.')

rownames(data.mat) <- NULL

# save the output
data.ls <- list(data.mat = data.mat,
                sp.names = sp.names,
                sp.latins = sp.latins,
                years = years)
saveRDS(data.ls, 
        file = sprintf('Data_Zone/Output/%s.RData', site))

