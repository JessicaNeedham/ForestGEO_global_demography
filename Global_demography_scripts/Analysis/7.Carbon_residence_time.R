### LIBRARIES & FUNCTIONS ###
library(RColorBrewer)
library(stringi)
########################################################
rm(list = ls())
set.seed(1)

# Temperate equation
agb.est.temp <- function(dbh, b0, b1){
  agb.est.temp <- exp(b0 + b1 * log(dbh/10))
  return(agb.est.temp)
}

# tropical equation for bioground
agb.est <- function(D, E, wsg){
  agb.est <- exp(-1.803 - 0.976*E + 0.976*log(wsg)
                 + 2.673*log(D) - 0.0299*(log(D)^2))
  return(agb.est)
}
########################################################

geog <- 'global'
load(file = sprintf('Analysis/Outputs/%s_pca_mat.RData', geog))
load(file = sprintf('Analysis/Outputs/PCA_%s.RData', geog))

site.names <- as.character(unique(pca.matrix[ ,'Site']))
site.names <- ifelse(nchar(site.names) <= 4, toupper(site.names), 
                     stri_trans_totitle(site.names))
site.names[which(site.names == 'Khaochong')] <- 'Khao Chong'
site.names[which(site.names == 'Laplanada')] <- 'La Planada'
site.names[which(site.names == 'Windriver')] <- 'Wind River'
site.names[which(site.names == 'XTBG')] <- 'XSBN'

site.names <- site.names[-which(site.names == 'Mudumalai')]
sites <- as.character(unique(pca.matrix[ ,'Site']))
sites <- sites[-which(sites == 'mudumalai')]

n.sites <- length(sites)

load('Analysis/Outputs/site_df.RData')
site.df <- site.df[-which(site.df$Site=='Mudumalai'), ]
plot.sizes <- site.df$Area

# load Es - the climate variables needed for allometries
Es <- read.table('Data_Zone/Data/EforBiomass.txt')
full.site.names <- c('Amacayacu', 'Barro Colorado Island', 'Changbaishan', 'Fushan', 'Huai Kha Khaeng', 
                     'Ituri', 'Khao Chong', 'Korup', 'Lambir', 'La Planada', 'Luquillo', 'Palanan', 'Pasoh', 'SCBI', 'SERC', 
                     'Wytham Woods', 'Wind River', 'Xishuangbanna ', 'Yasuni')
E_biomass <- Es[match(full.site.names, Es[ ,'Site']), 'E']

# load jenkins parameterd
temp.ps <- read.table('Data_Zone/Data/Jenkins_eq_params.txt', header = TRUE)

temp <- which(sites %in% c('changbaishan', 'scbi', 'serc', 'windriver', 'wytham'))
trop <- seq(n.sites)[-temp]

# for each site - get change in basal area from growth, mortality and 
# recruitment on a per hectare scale
mort <- rep(NA, n.sites)
recr <- rep(NA, n.sites)
growth <- rep(NA, n.sites)

for(j in trop){
  site.name <- sites[j]
  data.mat <- readRDS(file = sprintf('Data_Zone/Output/%s_all.sp.RData', site.name))
  # remove stems not present in either census
  data.mat <- data.mat[which(!is.na(data.mat$dbh.1) | !is.na(data.mat$dbh.2)), ]
  E <- E_biomass[j]
  wsg.global <- data.mat[which(data.mat[ ,'wsg.level'] == 'global')[1] ,'wsg']
  
  dbh1 <- data.mat$dbh.1
  dbh2 <- data.mat$dbh.2
  
  # time
  time <- data.mat$date.2 - data.mat$date.1
  time[is.na(time)] <- median(time[!which(time==0)], na.rm = TRUE)
  time[which(time==0)] <- median(time[!which(time==0)])
  
  n.years <- round(mean(time, na.rm = TRUE))
  
  # index for those that died and those that recruited
  # need to account for multi stems!
  
  died <- which(is.na(dbh2) & !is.na(dbh1))
  recruits <- which(is.na(dbh1) & !is.na(dbh2))
  survived <- which(!is.na(dbh1) & !is.na(dbh2))
  area <- plot.sizes[j]
  
  # Growth - has 2 components
  # 1. sum of growth of surviving trees
  # 2. sum of unobserved growth of trees that died
  wsg <- data.mat[ ,'wsg']
  
  agb1 <- agb.est(dbh1/10, E, wsg)
  agb2 <- agb.est(dbh2/10, E, wsg)
  
  growth1 <- sum((agb2 - agb1)/time, na.rm = TRUE)/area

  # growth 2
  N0 <- length(which(!is.na(agb1)))  # number alive in census 1
  t <- mean(time, na.rm =TRUE)  # time interval
  Nst <- length(survived)  # number surviving
  Nt <- length(!is.na(dbh2))  # number alive in census 2
  ma <- 1- (Nst/N0)^(1/t) # annual per capita mortality
  
  # how many trees died in each year of the census interval?
  n.died <- rep(0, n.years)
  for(i in 1:n.years){
    n.died[i] <- ma * (N0 - sum(n.died))
  }
  
  # mean number of years trees that died would have lived before death
  Ymean <-  sum((n.died * seq(n.years)))/sum(n.died)
  # plot level median growth for each size class
  d1 <- which(dbh1 < 200)
  d2 <- which(dbh1>=200 & dbh1 < 400)
  d3 <- which(dbh1 >= 400)
  
  g1 <- median((dbh2[d1] - dbh1[d1])/time[d1],na.rm = TRUE)
  g2 <- median((dbh2[d2] - dbh1[d2])/time[d2],na.rm = TRUE)
  g3 <- median((dbh2[d3] - dbh1[d3])/time[d3],na.rm = TRUE)
  if(site.name == 'changbaishan'){
    g1 <- median((dbh2[d1] - dbh1[d1])/time,na.rm = TRUE)
    g2 <- median((dbh2[d2] - dbh1[d2])/time,na.rm = TRUE)
    g3 <- median((dbh2[d3] - dbh1[d3])/time,na.rm = TRUE)
  }
  
  
  bins <- cut(dbh1[died], breaks = c(0, 200, 400, Inf), labels = FALSE)
  G <- c(g1,g2,g3)[bins]
  
  dbh.t0 <- dbh1[died]  # starting size of trees that died
  dbh.death <- dbh.t0 + G * Ymean  # size at death of trees that died 
  
  growth2 <- sum(agb.est(dbh.death/10, E, wsg[died]) - agb.est(dbh.t0/10, E, wsg[died]))
  growth2 <- growth2/area
  
  #growth[j] <- sum(growth1, growth2)
  
  ### Recruitment
  # 1. new recruits
  recr1 <- sum(agb.est(dbh2[recruits]/10, E, wsg[recruits]))/area/n.years
  
  # 2. unobserved recruits (that died before being recorded)
  Ma <- (N0/area)*(1-(Nst/N0)^(1/t))  # per area ANNUAL mortality
  Ra <- (Ma * (Nt-Nst)/(N0-Nst)) # per area ANNUAL recruitment rate
  
  # number of recruits that are unobserved each year
  # (1-ma)^t is survival to next census
  n.unobserved.recrs <- rep(NA, n.years)
  for(i in 1:n.years){
    n.unobserved.recrs[i] <- Ra - ((1-ma)^i * Ra)
  }
  # total unobserved recruits in the census interval
  tur <- sum(n.unobserved.recrs)
  
  # mean number of years recruits that died would have lived before death
  Ymean_r <-  (n.unobserved.recrs * seq(n.years))/sum(n.unobserved.recrs)
  
  # diameter at death
  ddeath = 10 + (g1 * Ymean_r)
  
  # divide by n.years to make it annual
  recr2 <- sum(agb.est(ddeath/10, E, wsg.global) * tur)/n.years  
  
  recr[j] <- sum(recr1, recr2)
  
  growth[j] <- sum(growth1, growth2, recr1, recr2)
  
  ### Mortality 
  # sum of trees that died (plus their unobserved growth)
  mort1 <- sum(agb.est(dbh.death/10, E, wsg[died]))/area/n.years
  
  # sum of recruits that entered and died
  mort2 <- recr2
  
  mort[j] <- sum(mort1, mort2)
  
  print(j)
}

for(j in temp){
  site.name <- sites[j]
  data.mat <- readRDS(file = sprintf('Data_Zone/Output/%s_all.sp.RData', site.name))
  # remove stems not present in either census
  data.mat <- data.mat[which(!is.na(data.mat$dbh.1) | !is.na(data.mat$dbh.2)), ]
  E <- E_biomass[j]
  wsg.global <- data.mat[which(data.mat[ ,'wsg.level'] == 'global')[1] ,'wsg']
  
  dbh1 <- data.mat$dbh.1
  dbh2 <- data.mat$dbh.2
  
  # time
  time <- data.mat$date.2 - data.mat$date.1
  time[is.na(time)] <- median(time[!which(time==0)], na.rm = TRUE)
  time[which(time==0)] <- median(time[!which(time==0)])
  if(site.name == 'changbaishan'){
    time <- 5
  }
  if(site.name == 'wytham'){
    time <- data.mat$time
  }
  
  n.years <- round(mean(time, na.rm = TRUE))
  
  # index for those that died and those that recruited
  died <- which(is.na(dbh2) & !is.na(dbh1))
  recruits <- which(is.na(dbh1) & !is.na(dbh2))
  survived <- which(!is.na(dbh1) & !is.na(dbh2))
  area <- plot.sizes[j]
  
  b0s <- b1s <- rep(NA, nrow(data.mat))
  
  # Growth - has 2 components
  # 1. sum of growth of surviving trees
  # 2. sum of unobserved growth of trees that died
  for(i in 1:nrow(data.mat)){
    # find the correct parameters
    rows <- which(temp.ps[ ,'Taxa'] %in% c(data.mat[i, 'Genus'], data.mat[i ,'Family']))
    
    b0s[i] <- ifelse(data.mat[i, 'wsg'] < temp.ps[rows[1], 'high'], 
                 temp.ps[rows[1], 'b0'], temp.ps[rows[2], 'b0'])
    b1s[i] <- ifelse(data.mat[i, 'wsg'] < temp.ps[rows[1], 'high'], 
                 temp.ps[rows[1], 'b1'], temp.ps[rows[2], 'b1'])
    
    if(is.na(b0s[i])){
      b0s[i] <- -2.48 # mixed hardwood
      b1s[i] <- 2.4835
    }
  }
  
  agb1 <- agb.est.temp(dbh1, b0s, b1s)
  agb2 <- agb.est.temp(dbh2, b0s, b1s)
  
  growth1 <- sum((agb2 - agb1)/time, na.rm = TRUE)/area
  
  # growth 2
  N0 <- length(which(!is.na(agb1)))  # number alive in census 1
  t <- mean(time, na.rm =TRUE)  # time interval
  Nst <- length(survived)  # number surviving
  Nt <- length(!is.na(dbh2))  # number alive in census 2
  ma <- 1- (Nst/N0)^(1/t) # annual per capita mortality
  
  # how many trees died in each year of the census interval?
  n.died <- rep(0, n.years)
  for(i in 1:n.years){
    n.died[i] <- ma * (N0 - sum(n.died))
  }
  
  # mean number of years trees that died would have lived before death
  Ymean <-  sum((n.died * seq(n.years)))/sum(n.died)
  # plot level median growth for each size class
  d1 <- which(dbh1 < 200)
  d2 <- which(dbh1>=200 & dbh1 < 400)
  d3 <- which(dbh1 >= 400)
  
  g1 <- median((dbh2[d1] - dbh1[d1])/time[d1],na.rm = TRUE)
  g2 <- median((dbh2[d2] - dbh1[d2])/time[d2],na.rm = TRUE)
  g3 <- median((dbh2[d3] - dbh1[d3])/time[d3],na.rm = TRUE)
  if(site.name == 'changbaishan'){
    g1 <- median((dbh2[d1] - dbh1[d1])/time,na.rm = TRUE)
    g2 <- median((dbh2[d2] - dbh1[d2])/time,na.rm = TRUE)
    g3 <- median((dbh2[d3] - dbh1[d3])/time,na.rm = TRUE)
  }
  
  
  bins <- cut(dbh1[died], breaks = c(0, 200, 400, Inf), labels = FALSE)
  G <- c(g1,g2,g3)[bins]
  
  dbh.t0 <- dbh1[died]  # starting size of trees that died
  dbh.death <- dbh.t0 + G * Ymean  # size at death of trees that died 
  
  growth2 <- sum(agb.est.temp(dbh.death, b0s[died], b1s[died]) - agb.est.temp(dbh.t0, b0s[died], b1s[died]))
  growth2 <- growth2/area
  
  #growth[j] <- sum(growth1, growth2)
  
  ### Recruitment
  # 1. new recuits
  recr1 <- sum(agb.est.temp(dbh2[recruits], b0s[recruits], b1s[recruits]))/area/n.years
  
  # 2. unobserved recruits (that died before being recorded)
  Ma <- (N0/area)*(1-(Nst/N0)^(1/t))  # per area ANNUAL mortality
  Ra <- (Ma * (Nt-Nst)/(N0-Nst)) # per area ANNUAL recruitment rate
  
  # number of recruits that are unobserved each year
  # (1-ma)^t is survival to next census
  n.unobserved.recrs <- rep(NA, n.years)
  for(i in 1:n.years){
    n.unobserved.recrs[i] <- Ra - ((1-ma)^i * Ra)
  }
  # total unobserved recruits in the census interval
  tur <- sum(n.unobserved.recrs)
  
  # mean number of years recruits that died would have lived before death
  Ymean_r <-  (n.unobserved.recrs * seq(n.years))/sum(n.unobserved.recrs)
  
  # diameter at death
  ddeath = 10 + (g1 * Ymean_r)
  
  # divide by n.years to make it annual
  recr2 <- sum(agb.est.temp(ddeath, -2.48, 2.4835) * tur)/n.years  
  
  recr[j] <- sum(recr1, recr2)
  
  growth[j] <- sum(growth1, growth2, recr1, recr2)
  
  ### Mortality 
  # sum of trees that died (plus their unobserved growth)
  mort1 <- sum(agb.est.temp(dbh.death, b0s[died], b1s[died]))/area/n.years
  
  # sum of recruits that entered and died
  mort2 <- recr2
  
  mort[j] <- sum(mort1, mort2)
  
  print(j)
}


growth <- growth/1000
mort <- mort/1000
recr <- recr/1000

save(growth, file = 'Analysis/Outputs/site_gr.RData')
save(mort, file = 'Analysis/Outputs/site_mort.RData')

names(growth) <- site.names
names(mort) <- site.names
names(recr) <- site.names

load(file = 'Analysis/Outputs/AGBs.RData')

tau <- AGBs/mort

save(tau, file = 'Analysis/Outputs/taus.RData')


