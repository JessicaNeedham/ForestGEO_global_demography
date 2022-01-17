#############################################
## NB for this site no taxonomic data RE families
# manually went through species with large enough sample 
# sizes for vital rate models - none are in excluded families
##############################################
rm(list = ls())
# LIBRARIES #
library(foreach)
library(doParallel)
library(stringr)

##############################################################################
# User defined variables      
site.name <- 'changbaishan'  # CTFS site name used in r table 
censuses <- c(2,3)      # which censuses to use
n.census <- length(censuses) # number of census to include in analysis
unit.dbh <-   "cm"             # input the units for dbh either "cm" or "mm"
dbh.factor <- ifelse(unit.dbh=="cm", 10, 1)   # converstion factor for dbh
nproc <- 1  # no. cores to dedicate to running this script
registerDoParallel(nproc)
# End user defined portion of script
##############################################################################
set.seed(1)
# load in all censuses
cns <- read.csv(file = sprintf('Data_Zone/Data/%s/%s_data.csv', 
                               site.name, site.name))

# rearrange columns and only keep useful columns
cnsdata <- cns[ ,c('treeID', "tag", "stemID", "Latin.name", 
                   "gx", "gy",  'dbh2', 'dbh3',
                  'status2', 'status3')]
# change status to just alive, dead or NA (for not yet recruited)
statuses <- cnsdata[ ,grep('status', colnames(cnsdata))]
statuses <- apply(statuses, 2, function(x){ x[which(x == 'lean' | x == 'alive$wd' |
                                                      x == 'alive$br' |
                                                      x == "alive$omit")] <- 'alive'; x})
cnsdata[ ,grep('status', colnames(cnsdata))] <- statuses
head(cnsdata)

# make a species code variable
cnsdata <- cnsdata[-which(cnsdata$Latin.name == '_\x90\xd5\xf1\x8d\xcb'), ]

sp.latins <- as.character(cnsdata$Latin.name)
sp.latins[which(sp.latins == '_\x90\xd5\xf1\x8d\xcb')] <- NA

# The varieties are the same for each species so just 
# use the species names
vars <- grep('var', sp.latins)
sp.latins[vars] <- sapply(sp.latins[vars], function(x) word(x, 1, 4))
sp.latins[-vars] <- sapply(sp.latins[-vars], function(x)
  ifelse(str_count(x, ' ')<1, x, word(x, 1, 2)))
cnsdata$Latin.name <- sp.latins
sp.latins <- unique(sp.latins)

sp.names <- toupper(as.character(sapply(sp.latins, function(x) 
  paste(str_sub(word(x, 1, 1), 1, 4),
        str_sub(word(x, 2,2), 1, 2),
        sep = ''))))

sp.names[which(table(sp.names) > 1)] <-
  paste0(sp.names[which(table(sp.names)>1)], 1)
sp.latins <- cbind(sp.latins, sp.names)
cnsdata$sp <- sort(sp.names)[as.numeric(as.factor(cnsdata$Latin.name))]

tmp <- cnsdata[ ,c('treeID', 'tag', 'stemID', 'sp', 'Latin.name',
                   'gx', 'gy',
                   'dbh2', 'status2', 'dbh3', 'status3')]

# no dates! so just take dates from CTFS website and assume
# there are five years between each measurement 
years <- c('2009', '2014')

# check to see if there are dbh values in dead trees
length(which(!is.na(tmp$dbh2) & tmp$status2 == 'D'))
length(which(!is.na(tmp$dbh3) & tmp$status3 == 'D'))

# change dbh to NA if status is D
tmp$dbh2 <- ifelse(tmp$status2 == 'D', NA, tmp$dbh2)
tmp$dbh3 <- ifelse(tmp$status3 == 'D', NA, tmp$dbh3)

# get rid of trees dead in both
tmp <- tmp[-which(tmp$status2=='D' & tmp$status3 == 'D'), ]

S <- ifelse(tmp$status2 == 'alive' & tmp$status3 == 'alive', 1, NA)
S <- ifelse(tmp$status2 == 'alive' & tmp$status3 == 'D', 0, S)
table(S, useNA = 'ifany')


# survival based on NAs in sizes
size.mat <- tmp[ ,grep('dbh', colnames(tmp))]
size.mat <- apply(size.mat, 2, function(x) x*dbh.factor)
tmp[ ,grep('dbh', colnames(tmp))] <- size.mat

S <- ifelse(is.na(size.mat[ ,2]) & 
                !is.na(size.mat[ ,1]), 0, NA)
S <-  ifelse(!is.na(size.mat[ ,2]) &
                 !is.na(size.mat[ ,1]), 1, S)
  
tmp[ ,'Surv'] <- S

# get rid of trees with no measurement in any census
all.na <- which(apply(size.mat, 1, function(x) all(is.na(x))))
if(length(all.na)>0){
  tmp <- tmp[-all.na, ]
}

all.sp <- tmp
# tidy the workspace
remove(cns)
remove(cnsdata)
remove(tmp)
remove(statuses)

# add family - doing this manually from the species list on ForestGEO
all.sp$Genus <- word(all.sp$Latin.name,1,1)
all.sp$Family <- NA
all.sp$Family <- ifelse(all.sp$Genus == 'Abies', 'Pinaceae', all.sp$Family)
all.sp$Family <- ifelse(all.sp$Genus == 'Acer', 'Sapindaceae', all.sp$Family)
all.sp$Family <- ifelse(all.sp$Genus == 'Actinidia', 'Actinidiaceae', all.sp$Family)
all.sp$Family <- ifelse(all.sp$Genus == 'Aralia', 'Araliaceae', all.sp$Family)
all.sp$Family <- ifelse(all.sp$Genus == 'Betula', 'Betulaceae', all.sp$Family)
all.sp$Family <- ifelse(all.sp$Genus == 'Corylus', 'Betulaceae', all.sp$Family)
all.sp$Family <- ifelse(all.sp$Genus == 'Crataegus', 'Rosaceae', all.sp$Family)
all.sp$Family <- ifelse(all.sp$Genus == 'Deutzia', 'Hydrangeaceae', all.sp$Family)
all.sp$Family <- ifelse(all.sp$Genus == 'Eleutherococcus', 'Araliaceae', all.sp$Family)
all.sp$Family <- ifelse(all.sp$Genus == 'Euonymus', 'Celastraceae', all.sp$Family)
all.sp$Family <- ifelse(all.sp$Genus == 'Fraxinus', 'Oleaceae', all.sp$Family)
all.sp$Family <- ifelse(all.sp$Genus == 'Lonicera', 'Caprifoliaceae', all.sp$Family)
all.sp$Family <- ifelse(all.sp$Genus == 'Maackia', 'Fabaceae', all.sp$Family)
all.sp$Family <- ifelse(all.sp$Genus == 'Malus', 'Rosaceae', all.sp$Family)
all.sp$Family <- ifelse(all.sp$Genus == 'Phellodendron', 'Rutaceae', all.sp$Family)
all.sp$Family <- ifelse(all.sp$Genus == 'Philadelphus', 'Hydrangeaceae', all.sp$Family)
all.sp$Family <- ifelse(all.sp$Genus == 'Pinus', 'Pinaceae', all.sp$Family)
all.sp$Family <- ifelse(all.sp$Genus == 'Populus', 'Salicaceae', all.sp$Family)
all.sp$Family <- ifelse(all.sp$Genus == 'Prunus', 'Rosaceae', all.sp$Family)
all.sp$Family <- ifelse(all.sp$Genus == 'Quercus', 'Fagaceae', all.sp$Family)
all.sp$Family <- ifelse(all.sp$Genus == 'Rhamnus', 'Rhamnaceae', all.sp$Family)
all.sp$Family <- ifelse(all.sp$Genus == 'Salix', 'Salicaceae', all.sp$Family)
all.sp$Family <- ifelse(all.sp$Genus == 'Sambucus', 'Adoxaceae', all.sp$Family)
all.sp$Family <- ifelse(all.sp$Genus == 'Schisandra', 'Schisandraceae', all.sp$Family)
all.sp$Family <- ifelse(all.sp$Genus == 'Sorbaria', 'Rosaceae', all.sp$Family)
all.sp$Family <- ifelse(all.sp$Genus == 'Sorbus', 'Rosaceae', all.sp$Family)
all.sp$Family <- ifelse(all.sp$Genus == 'Syringa', 'Oleaceae', all.sp$Family)
all.sp$Family <- ifelse(all.sp$Genus == 'Tilia', 'Malvaceae', all.sp$Family)
all.sp$Family <- ifelse(all.sp$Genus == 'Ulmus', 'Ulmaceae', all.sp$Family)
all.sp$Family <- ifelse(all.sp$Genus == 'Viburnum', 'Adoxaceae', all.sp$Family)
all.sp$Family <- ifelse(all.sp$Genus == 'Vitis', 'Vitaceae', all.sp$Family)

# WSG
wsg.mat <- read.csv('Data_Zone/Data/GlobalWoodDensityDatabase.csv')
colnames(wsg.mat)[4] <- 'WSG'
# there are multiple entries for some species  - take the means
wsgs <- tapply(wsg.mat[ ,'WSG'], list(wsg.mat[ ,'Binomial']), mean)
# genus level
wsg.mat$Genus <- word(wsg.mat$Binomial,1,1)
genus.wsgs <- tapply(wsg.mat$WSG, wsg.mat$Genus, mean)
# family level
family.wsgs <- tapply(wsg.mat$WSG, wsg.mat$Family, mean)

data.mat <- all.sp
data.mat <- cbind(data.mat, NA, NA)
colnames(data.mat)[(ncol(data.mat)-1):ncol(data.mat)] <- c('wsg', 'wsg.level') 

data.mat[ ,'wsg'] <- wsgs[match(data.mat[ ,'Latin.name'], names(wsgs))]
data.mat[ ,'wsg.level'] <- ifelse(is.na(data.mat[ ,'wsg']), NA, 'species')

g.wsgs <- genus.wsgs[sapply(data.mat[ ,'Genus'], function(x) match(x, names(genus.wsgs)))]

data.mat$wsg <- ifelse(is.na(data.mat$wsg), g.wsgs, data.mat$wsg)
data.mat$wsg.level <- ifelse(is.na(data.mat$wsg.level) & !is.na(data.mat[ ,'wsg']),
                             'genus', data.mat[ ,'wsg.level'])

f.wsgs <- family.wsgs[sapply(data.mat[ ,'Family'], function(x) match(x, names(family.wsgs)))]

data.mat$wsg <- ifelse(is.na(data.mat$wsg), f.wsgs, data.mat$wsg)
data.mat$wsg.level <- ifelse(is.na(data.mat$wsg.level) & !is.na(data.mat[ ,'wsg']),
                             'family', data.mat[ ,'wsg.level'])


# and global for everything else
global <- mean(wsgs, na.rm = TRUE)
data.mat$wsg <- ifelse(is.na(data.mat$wsg), global, data.mat$wsg)
data.mat$wsg.level <- ifelse(is.na(data.mat$wsg.level),
                             'global', data.mat[ ,'wsg.level'])

colnames(data.mat) <- c('treeID', 'tag', 'stemID', 'sp', 'Latin', 'gx', 'gy', 'dbh.1', 'status.1', 
                        'dbh.2', 'status.2', 'surv', 'Genus', 'Family', 'wsg', 'wsg.level')
data.mat <- data.mat[ ,-grep('status', colnames(data.mat))]


### 

site.name <- 'changbaishan'
saveRDS(data.mat, 
        file = sprintf('Data_Zone/Output/%s_all.sp.RData', site.name))


############################

# remove recruits
data.mat <- data.mat[-which(is.na(data.mat$dbh.1) & !is.na(data.mat$dbh.2)), ]

## CHOOSE SPECIES ##
# split into species specific matrices
sp.list <- split(data.mat, data.mat$sp)

tmp <- lapply(sp.list, function(x) apply(x, 2, function(z) 
  length(which(!is.na(z))))[grep('dbh', colnames(x))])

keep <- which(lapply(tmp, function(x) ifelse(any(x > 200), 1, 0))
              == 1)

all.sp <- data.mat[data.mat$sp %in% names(keep) == TRUE, ]
all.sp$sp <- factor(all.sp$sp)
# split into species specific matrices
sp.list <- split(all.sp, all.sp$sp)

# vector of names
sp.names <- na.omit(as.character(unlist(lapply(sp.list, function(x) x$sp[1]))))
sp.latins <- sp.latins[sp.latins[ ,2] %in% sp.names, ]
sp.latins <- sp.latins[order(sp.latins[ ,1]), ]  


# GROWTH CORRECTION # 

# for each species get rid of those with individuals more than twice the
# 99th quantile of size. 
source('Data_Zone/Source/CTFS_data_processing_functions.R')

# bind species back together to process all at once
data.mat <- do.call(rbind, sp.list)

# make it a matrix to speed things up
sizes <- as.matrix(data.mat[ ,grep('dbh', colnames(data.mat))])    
incr <- t(t(apply(sizes, 1, diff)))

# get time in years between censuses 
time <- matrix(5, nrow(data.mat), ncol = 1)

# survival
survival <- as.matrix(data.mat[ ,grep('surv', colnames(data.mat))],
                      ncol = n.census -1)
# coordinates
gx <- as.numeric(data.mat$gx)
gy <- as.numeric(data.mat$gy)
# make it a matrix - much quicker to process
sp.id <- as.numeric(data.mat$sp)  # first make species a numeric 
wsg <- data.mat[ ,'wsg']
wsg.level <- factor(data.mat$wsg.level, levels = c('species','genus','family','global'))
wsg.level <- as.numeric(wsg.level)

# check they are both alphabetical
table(sp.id)
table(data.mat[ ,'Latin'])

data.mat <- as.matrix(cbind(data.mat$treeID, data.mat$stemID,
                            sp.id, gx, gy, sizes, survival, incr, time, 
                            data.mat[ ,'wsg'], wsg.level))

colnames(data.mat)[1:2] <- c('treeID', 'stemID')


sizes <- as.matrix(data.mat[ ,grep('dbh', colnames(data.mat))])   

# finding errors and correcting them  - 
# this is based on annual increments- stems that shrank more than 25% of original size 
# or grew more than 75 mm in a single year
pc.change <- sizes[ ,2]/sizes[ ,1]*100
data.mat <- data.mat[-which(pc.change < 75), ]

colnames(data.mat) <- c('treeID', 'stemID','sp.id', 'gx', 'gy', 
                        'dbh.1', 'dbh.2', 'surv', 'incr', 'time', 'wsg', 'wsg.level')

annual.incr <- data.mat[ ,'incr']/data.mat[ ,'time']
if(length(which(annual.incr > 75))){
  data.mat <- data.mat[-which(annual.incr > 75), ]
}  
 
#  Remove stems dead in all censuses
bb <- data.mat[ ,grep('dbh', colnames(data.mat))]  
data.mat <- data.mat[which(apply(bb, 1, 
                                 function(x) all(is.na(x))) 
                           == FALSE), ]

# - propose new sizes so that size monotonically increases
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
incrs.lll <- t(apply(med.corrections, 1, diff))

# replace the increment columns with the new positive increments
data.mat[ ,grep('incr', colnames(data.mat))] <- incrs.lll

rownames(data.mat) <- NULL
colnames(data.mat) <- c('treeID', 'stemID', 'sp.id', 'gx', 'gy', 'dbh.1', 'dbh.2', 
                        'surv', 'incr', 'time', 'wsg', 'wsg.level')


# save the output
data.ls <- list(data.mat = data.mat,
                sp.names = sp.names,
                sp.latins = sp.latins,
                years = years)
saveRDS(data.ls, 
        file = sprintf('Data_Zone/Output/%s.RData', site.name))
