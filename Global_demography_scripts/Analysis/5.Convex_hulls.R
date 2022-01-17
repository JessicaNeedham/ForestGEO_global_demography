############################################
### Convex hulls ###
############################################
rm(list = ls())
set.seed(1)
##############################################
### LIBRARIES ###
##############################################
library(geometry)
library(RColorBrewer)
library(stringi)
library(foreach)
library(doParallel)
library(data.table)

get.subplot <- function(x, y, max.x){
  sub <- ((max.x * y) + (x - max.x))
  return(sub)
}
###########################################
### DATA ###
###########################################
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

sites <- as.character(unique(pca.matrix[ ,'Site']))
sites <- sites[-which(sites == 'mudumalai')]
site.names <- site.names[-which(site.names == 'Mudumalai')]
n.sites <- length(sites)

n.bootstraps <- 500
##############################################

# Site by site load data - 500 bootstrap samples of 
# 16 ha - then calculate convex hull for position of
# species in those 16ha in PCA space
load(file = 'Analysis/Outputs/site_sps_16ha.RData')
load(file = 'Analysis/Outputs/site_counts.RData')
site.median.hulls_2 <- rep(NA, n.sites)

for(i in 1:n.sites){
  
  site <- sites[i] 
  
  site.pca.index <- which(pca.matrix[ ,'Site'] == site)
  site.pca.matrix <- pca.matrix[which(pca.matrix[ ,'Site'] == site), ]
  site.pca <- res$ind$coord[site.pca.index, ]
  
  data.ls <- try(readRDS(file = sprintf('Data_Zone/Output/%s.RData', site)), 
                 silent = TRUE)
  data.mat <- data.ls$data.mat
  sp.names <- data.ls$sp.names
  
  # Keep only species in both
  keep <- which(site.pca.matrix[ ,'Mnemonic'] %in% sp.names)
  site.pca.matrix <- site.pca.matrix[keep, ]
  site.pca <- site.pca[keep, ]
  
  sp.ids <- which(sp.names %in% site.pca.matrix[ ,'Mnemonic'])
  site.pca <- cbind(site.pca, sp.ids)
  
  uni <- site.sps.16ha[[i]] 
  counts <- site.counts[[i]]
  convex.hull.vols_2 <- rep(NA, length(uni))
  coords.ls <- vector('list', length(uni))
  
  for(j in 1:length(uni)){
    sp.sample <- uni[[j]]
    coords <- site.pca[match(sp.sample, site.pca[ ,'sp.ids']), 1:5]
    coords <- coords[!is.na(coords[ ,1]), ]
    
    # weight by percent variance explained for each axis
    coords <- t(t(coords)*(res$eig[ ,2]/100))
    
    coords.ls[[j]] <- coords
    convex.hull.vols_2[j] <- convhulln(coords[ ,1:2], output.options = 'FA')$vol
    
  }
  
  saveRDS(coords.ls, file = sprintf('Analysis/Outputs/%s_coords.RData', site))
  saveRDS(convex.hull.vols_2, file = sprintf('Analysis/Outputs/%s_chulls.RData', site))
  
  chull2 <- rep(convex.hull.vols_2, times = counts)
  site.median.hulls_2[i] <- median(chull2)
 
  print(i)
}

save(site.median.hulls_2, file = 'Analysis/Outputs/Site_median_hulls_2.RData')
save(site.counts, file = 'Analysis/Outputs/Site_area_counts.RData')

par(mfrow = c(2,1))
barplot(site.median.hulls_2)
#barplot(site.median.hulls_5)

############################################################################
# save the outlines in 2D for plotting 
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

sites <- as.character(unique(pca.matrix[ ,'Site']))
n.sites <- length(sites)

convex_hull_coords <- vector('list', n.sites)
counter <- 1

for(i in 1:n.sites){
  
  site <- sites[i]
  
  site.pca.index <- which(pca.matrix[ ,'Site'] == site)
  site.pca.matrix <- pca.matrix[which(pca.matrix[ ,'Site'] == site), ]
  site.pca <- res$ind$coord[site.pca.index, ]
  
  coords <- site.pca[, 1:2]
  coords <- coords[!is.na(coords[ ,1]), ]
  
  # weight by percent variance explained for each axis
  #coords <- t(t(coords)*(res$eig[ ,2]/100))
  
  chulls <- chull(coords[ ,c(1,2)])
  chulls <- c(chulls, chulls[1])
  
  tmp <- coords[chulls,c(1,2)]
  #plot(coords[,1],coords[,2])
  #lines(tmp)
  convex_hull_coords[[i]] <- tmp
  
}

save(convex_hull_coords, file = 'Analysis/Outputs/convex_hull_coords.RData')



#############################################################################
### Pairwise intersects
############################################
### Convex hulls ###
############################################
rm(list = ls())

get.subplot <- function(x, y, max.x){
  sub <- ((max.x * y) + (x - max.x))
  return(sub)
}

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
sites <- as.character(unique(pca.matrix[ ,'Site']))
sites <- sites[-which(sites == 'mudumalai')]
site.names <- site.names[-which(site.names == 'Mudumalai')]
n.sites <- length(sites)

n.bootstraps <- 500

load('Analysis/Outputs/site_df.RData')
site.df <- site.df[-which(site.df[ ,'Site'] == 'Mudumalai'), ]
ord <- order(abs(site.df[ ,'Latitude']))
sites <- sites[ord]
site.names <- site.names[ord]
load(file = 'Analysis/Outputs/Site_area_counts.RData')
site.counts <- site.counts[ord]

for(i in 1:n.sites){
  site1 <- sites[i]
  vols.med <- rep(NA, n.sites)
  
  for(j in 1:n.sites){
    
    if(i >= j){next()}
    
    site2 <- sites[j]
    coords1 <- readRDS(sprintf('Analysis/Outputs/%s_coords.RData', site1))
    coords2 <- readRDS(sprintf('Analysis/Outputs/%s_coords.RData', site2))
    
    nsamps <- length(coords1)*length(coords2)
    print(paste0('nsamps: ', nsamps))
    
    vols <- rep(NA, nsamps)
    reps <- rep(NA, nsamps)
    
    counter <- 1
    for(k in 1:length(coords1)){
      for(l in 1:length(coords2)){
        
        vol <- try(intersectn(coords1[[k]][,1:2], coords2[[l]][,1:2], tol = 0, 
                          return.chs = TRUE, options = "Tv",
                        fp = NULL, autoscale = FALSE))
      try(vols[counter] <- vol$ch$vol)
        reps[counter] <- site.counts[[i]][k] * site.counts[[j]][l]
      counter <- counter + 1
      }
    }
    
    vols <- unlist(vols)
    vols <- rep(vols, reps)
    median.vol <- median(unlist(vols))
    
    vols.med[j] <- median.vol
    
    print(paste0('i: ', i, ' j: ', j))
  }
  save(vols.med, file = sprintf('Analysis/Outputs/vol_med_intersect_%s.RData', site1))
}

overlap <- matrix(NA, n.sites, n.sites)

for(i in 1:(n.sites)){
  site <- sites[i]
  load(sprintf('Analysis/Outputs/vol_med_intersect_%s.RData', site))
  overlap[i, ] <- vols.med
}

rownames(overlap) <- site.names
colnames(overlap) <- site.names

overlap
save(overlap, file = 'Analysis/Outputs/Overlaps.RData')



##############################################################
###############################################################
### Examples for presentations

# for five example sites save not just the volume but the outline 
examples <- c('Lambir', 'Pasoh')
n.examples <- length(examples)
ex.sites <- which(site.names %in% examples) 

convex_hull_coords <- vector('list', n.examples)
counter <- 1

for(i in ex.sites){
  
  site <- sites[i]
  
  site.pca.index <- which(pca.matrix[ ,'Site'] == site)
  site.pca.matrix <- pca.matrix[which(pca.matrix[ ,'Site'] == site), ]
  site.pca <- res$ind$coord[site.pca.index, ]
  
  data.ls <- try(readRDS(file = sprintf('Data_Zone/Output/%s.RData', site)), 
                 silent = TRUE)
  data.mat <- data.ls$data.mat
  sp.names <- data.ls$sp.names
  sp.ids <- which(sp.names %in% site.pca.matrix[ ,'Mnemonic'])
  site.pca <- cbind(site.pca, sp.ids)
  
  # make subplots
  if(site == 'ituri'){
    data.mat[ ,'quadrat'] <- as.numeric(data.mat[ ,'quadrat'])
  }else{
    max.gx <- round(max(data.mat[ ,'gx'], na.rm = TRUE))
    max.gy <- round(max(data.mat[ ,'gy'], na.rm = TRUE))
    max.x <- max.gx/20
    max.y <- max.gy/20
    sub.x <- cut(data.mat[ ,'gx'], breaks = seq(0, max.gx, 20),
                 label = FALSE, include.lowest = TRUE)
    sub.y <- cut(data.mat[ ,'gy'], breaks = seq(0, max.gy, 20), 
                 label = FALSE, include.lowest = TRUE)
    sub <- get.subplot(sub.x, sub.y, max.x)
    data.mat <- cbind(data.mat, sub)
  }
  
  samples.ls <- vector('list', n.bootstraps)
  
  if(site == 'ituri'){
    quads <- sample(unique(data.mat[ ,'quadrat']), 400, replace = TRUE)
    sp.sample <- unique(data.mat[which(data.mat[ ,'quadrat'] %in% quads), 'sp.id'])
  } else if(site == 'windriver'){
    sp.sample <- seq(nrow(site.pca.matrix))
  } else{
    subs <- sample(unique(data.mat[ ,'sub']), 400, replace = TRUE)
    sp.sample <- unique(data.mat[which(data.mat[ ,'sub'] %in% subs), 'sp.id'])
  }
  
  coords <- site.pca[match(sp.sample, site.pca[ ,'sp.ids']), 1:5]
  coords <- coords[!is.na(coords[ ,1]), ]
  
  # weight by percent variance explained for each axis
  #coords <- t(t(coords)*(res$eig[ ,2]/100))
  
  chulls <- chull(coords[ ,c(1,2)])
  chulls <- c(chulls, chulls[1])
  
  tmp <- coords[chulls,c(1,2)]
  #plot(coords[,1],coords[,2])
  #lines(tmp)
  convex_hull_coords[[counter]] <- tmp
  counter <- counter + 1
}

save(convex_hull_coords, file = 'Analysis/Outputs/convex_hull_coords_presentations.RData')
