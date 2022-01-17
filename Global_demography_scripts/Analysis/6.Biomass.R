########################################
### Allometries ###
########################################
##############################
### LIBRARIES ###
library(stringr)
library(stringi)
##################################
### FUNCTIONS ###
rm(list = ls())
# Get Basal Area in m2
get.ba <- function(dbh){
  dbh <- dbh/1000  # get it in m
  ba <- pi*((dbh/2)^2)
  return(ba)
}

# Temperate equation
agb.est.temp <- function(dbh, b0, b1){
  agb.est.temp <- exp(b0 + b1 * log(dbh/10))
  return(agb.est.temp)
}

# tropical equation for bioground
#agb.est <- function(D, E, wsg){
 # agb.est <- exp(-1.803 - 0.976*E + 0.976*log(wsg)
  #               + 2.673*log(D) - 0.0299*(log(D)^2))
  #return(agb.est)
#}

agb.est <- function(D,E,wsg) {
  #look up parameters based on E 
  H <-  exp(0.893 - E + 0.760*log(D) - 0.034*log(D)^2)
  return(0.0673*(wsg*(D)^2*H)^0.976)
}


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
site.names <- site.df[ ,'Site']

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

agb.clusts <- bas.clusts <- stems.clusts <- matrix(NA, nrow = n.sites,
                                                    ncol = 9)
colnames(agb.clusts) <- seq(9)
colnames(bas.clusts) <- seq(9)
colnames(stems.clusts) <- seq(9)

AGBs <- BAs <- rep(NA, n.sites)

# Run through the sites and get BA and AGB for each species 
for(ii in trop){
  
  site <- sites[ii]
  site.name <- site
  E <- E_biomass[ii]
  
  all.sp <- readRDS(file = sprintf('Data_Zone/Output/%s_all.sp.RData', site.name))
  all.sp[all.sp$wsg.level == 'global', 'wsg'] <- 0.615
  
  # matrix of just sizes
  size <- all.sp[ ,grep('dbh.2', colnames(all.sp))]
  wsg <- all.sp[ ,'wsg']
  agbs <- agb.est(size/10, E, wsg)
  
  sp.agb <- all.sp[ ,c('sp', 'Latin')]
  sp.agb$agb <- agbs
  sp.agb$BA <- get.ba(size)
  sp.agb$dbh <- size
  
  AGBs[ii] <- sum(agb.est(size/10, E, wsg), na.rm = TRUE)
  BAs[ii] <- sum(get.ba(size), na.rm = TRUE)
  
  #save(sp.agb, file = sprintf('Data_Zone/Output/%s_all_sp_agb.RData', site.name))
  
  # which are in pca matrix
  # pca matrix for that site
  site.pc.mat <- pca.matrix[which(pca.matrix[ ,'Site']==site), ]
  
  sp.agb$clust <- site.pc.mat[match(sp.agb$Latin, site.pc.mat[ ,'Species']) ,'clust']
  
  sp.agb$clust[which(is.na(sp.agb$clust))] <- 9
  
  save(sp.agb, file = sprintf('Data_Zone/Output/%s_all_sp_agb.RData', site.name))
  
  
  agbs <- tapply(sp.agb$agb, sp.agb$clust, sum, na.rm = TRUE)
  agb.clusts[ii, match(names(agbs), colnames(agb.clusts))] <- agbs
  
  bas <- tapply(sp.agb$BA, sp.agb$clust, sum, na.rm = TRUE)
  bas.clusts[ii, match(names(bas), colnames(bas.clusts))] <- bas
  
  stems <- tapply(sp.agb$sp, sp.agb$clust, length)
  stems.clusts[ii, match(names(stems), colnames(stems.clusts))] <- stems
  
    print(ii)
}

for(jj in temp){
  
  site <- sites[jj]
  site.name <- site
  all.sp <- readRDS(file = sprintf('Data_Zone/Output/%s_all.sp.RData', site.name))
  all.sp[all.sp$wsg.level == 'global', 'wsg'] <- 0.615
  
  
  # matrix of just sizes
  size <- all.sp[ ,grep('dbh.2', colnames(all.sp))]
  wsg <- all.sp[ ,'wsg']
  BAs[jj] <- sum(get.ba(size), na.rm = TRUE)
  
  agbs <- rep(NA, nrow(all.sp))
  counter <- 0
  
  for(i in 1:nrow(all.sp)){
    # find the correct parameters
    rows <- which(temp.ps[ ,'Taxa'] %in% c(all.sp[i, 'Genus'], all.sp[i ,'Family']))
    
    b0 <- ifelse(all.sp[i, 'wsg'] < temp.ps[rows[1], 'high'], 
                 temp.ps[rows[1], 'b0'], temp.ps[rows[2], 'b0'])
    b1 <- ifelse(all.sp[i, 'wsg'] < temp.ps[rows[1], 'high'], 
                 temp.ps[rows[1], 'b1'], temp.ps[rows[2], 'b1'])
    
    if(is.na(b0)){
      counter <- counter + 1
      #b0 <- -2.48 # mixed hardwood
      #b1 <- 2.4835
      b0 <- -2.3679
      b1 <- 2.445
    }
    
    size.stem <- all.sp[i, 'dbh.2']
    agbs[i] <- agb.est.temp(size.stem, b0, b1)
  }
 
  sp.agb <- all.sp[ ,c('sp', 'Latin')]
  sp.agb$agb <- agbs
  sp.agb$BA <- get.ba(size)
  AGBs[jj] <- sum(agbs, na.rm = TRUE)
  sp.agb$dbh <- size
  
  # which are in pca matrix
  # pca matrix for that site
  site.pc.mat <- pca.matrix[which(pca.matrix[ ,'Site']==site), ]
  
  sp.agb$clust <- site.pc.mat[match(sp.agb$Latin, site.pc.mat[ ,'Species']) ,'clust']
  
  sp.agb$clust[which(is.na(sp.agb$clust))] <- 9
  
  save(sp.agb, file = sprintf('Data_Zone/Output/%s_all_sp_agb.RData', site.name))
  
  
  agbs <- tapply(sp.agb$agb, sp.agb$clust, sum, na.rm = TRUE)
  agb.clusts[jj, match(names(agbs), colnames(agb.clusts))] <- agbs
  
  bas <- tapply(sp.agb$BA, sp.agb$clust, sum, na.rm = TRUE)
  bas.clusts[jj, match(names(bas), colnames(bas.clusts))] <- bas
  
  stems <- tapply(sp.agb$sp, sp.agb$clust, length)
  stems.clusts[jj, match(names(stems), colnames(stems.clusts))] <- stems
  
  print(jj)
}

percents.agb <- apply(agb.clusts, 1, function(x) sum(x[1:8],na.rm = TRUE)/sum(x, na.rm = TRUE)*100)
names(percents.agb) <- sites
  
percents.agb

percents.ba <- apply(bas.clusts, 1, function(x) sum(x[1:8],na.rm = TRUE)/sum(x, na.rm = TRUE)*100)
names(percents.ba) <- sites

percents.ba




# Divide by plot size
load(file = 'Analysis/Outputs/site_df.RData')
site.df <- site.df[-which(site.df$Site == 'Mudumalai'), ]
areas <- site.df$Area
AGBs <- AGBs/areas
# kilograms to Mg C /ha 
AGBs <- AGBs/1000


# compare to Lutz et al. 2018
lutz <- c(268, 257, 288, 224, 258, (467+375)/2, NA, 345, 495, 270, 283, 414, 
          324, 259, 299, 532, 310, 280, 261)


discr <- AGBs - lutz

pdf(file = 'Paper_Figures/lutz_discrepancies.pdf')
par(mfrow = c(1,1), mar = c(5,5,2,1))
plot(discr, ylab = '', col = 'white')
text(sites, x = seq(19), y = discr)
mtext('Discrepancy Mg C /ha', side = 2, line = 3)
dev.off()

names(AGBs) <- sites

save('AGBs', file = 'Analysis/Outputs/AGBs.RData')

par(mfrow = c(1,1))
barplot(AGBs, names.arg = sites)



# convert to Kg C m2
AGBs.kg <- AGBs*1000
AGBs.kg.m2 <- AGBs.kg/10000

geog <- 'global'
load(file = sprintf('Analysis/Outputs/%s_pca_mat.RData', geog))

load('Analysis/Outputs/site_df.RData')
site.df <- site.df[-which(site.df$Site=='Mudumalai'), ]
plot.sizes <- site.df$Area

# Now get AGB of each species at each site
for(j in 1:n.sites){
  site <- sites[j]
  load(file = sprintf('Data_Zone/Output/%s_all_sp_agb.RData', site))
  area <- plot.sizes[j]
  if(site == 'amacayacu'){
    sp.agb.sums <- tapply(sp.agb$agb, sp.agb$Latin, sum, na.rm = TRUE)
    names(sp.agb.sums) <- gsub(',', '.', names(sp.agb.sums))
    # pca matrix for that site
    site.pc.mat <- pca.matrix[which(pca.matrix[ ,'Site']==site), ]
    sp.agb.sums.pcmat <- sp.agb.sums[which(names(sp.agb.sums)%in% site.pc.mat[ ,'Species'])]
    
    }else{
    sp.agb.sums <- tapply(sp.agb$agb, sp.agb$sp, sum, na.rm = TRUE)
    # pca matrix for that site
    site.pc.mat <- pca.matrix[which(pca.matrix[ ,'Site']==site), ]
    sp.agb.sums.pcmat <- sp.agb.sums[which(names(sp.agb.sums)%in% site.pc.mat[ ,'Mnemonic'])]
  }
  
  sp.agb.sums <- sp.agb.sums/1000
  sp.agb.sums <- sp.agb.sums/area
  
  sp.agb.sums.pcmat <- sp.agb.sums.pcmat/1000
  sp.agb.sums.pcmat <- sp.agb.sums.pcmat/area
  
  
  save(sp.agb.sums, file = sprintf('Analysis/Outputs/%s_agb_sp.RData', site))
  save(sp.agb.sums.pcmat, file = sprintf('Analysis/Outputs/%s_agb_sp_pcmat.RData', site))
  
  
}

## Check that works 

rm(list = ls())
geog <- 'global'
load(file = sprintf('Analysis/Outputs/%s_pca_mat.RData', geog))
# remove sites with no data for now! - but add back later
pca.matrix <- pca.matrix[-which(pca.matrix[ ,'Site'] %in% c('mudumalai')), ]
pca.matrix[ ,'Site'] <- as.character(pca.matrix[ ,'Site'])
sites <- as.character(unique(pca.matrix[ ,'Site']))
n.sites <- length(sites)

agbs <- rep(NA, n.sites)

for(i in 1:n.sites){
  site <- sites[i]
  load(sprintf('Analysis/Outputs/%s_agb_sp.RData', site))
  agbs[i] <- sum(sp.agb.sums)
}

load(file = 'Analysis/Outputs/AGBs.RData')
plot(AGBs, agbs)

