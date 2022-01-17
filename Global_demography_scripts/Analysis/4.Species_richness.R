# Species richness at each site 
############################################
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
###########################################
### FUNCTIONS AND SETUP ###
###########################################
get.subplot <- function(x, y, max.x){
  sub <- ((max.x * y) + (x - max.x))
  return(sub)
}

###########################################
### DATA ###
###########################################
geog <- 'global'
load(file = sprintf('PCA_Zone/Output/%s_pca_mat.RData', geog))
load(file = sprintf('PCA_Zone/Output/PCA_%s.RData', geog))

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
# 16 ha - 
site.counts <- vector('list', n.sites)
site.sp.rch <- rep(NA, n.sites)
site.sp.rch.all <- rep(NA, n.sites)
site.sps.16ha <- vector('list', n.sites)

for(i in 1:n.sites){
  
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
  # now sample 400 subplots (16 ha)
  
  samples.ls <- vector('list', n.bootstraps)
  
  for(j in 1:n.bootstraps){
    
    if(site == 'ituri'){
      quads <- sample(unique(data.mat[ ,'quadrat']), 400, replace = TRUE)
      sp.sample <- unique(data.mat[which(data.mat[ ,'quadrat'] %in% quads), 'sp.id'])
    } else if(site == 'windriver'){
      sp.sample <- seq(nrow(site.pca.matrix))
    } else{
      subs <- sample(unique(data.mat[ ,'sub']), 400, replace = TRUE)
      sp.sample <- unique(data.mat[which(data.mat[ ,'sub'] %in% subs), 'sp.id'])
    }
    samples.ls[[j]] <- sort(sp.sample)
  }
  
  # What are the unique species combinations and how many times
  # does each combination occur?
  
  tmp <- do.call(rbind, samples.ls)
  for(k in 1:nrow(tmp)){
    tmp[k, ][duplicated(tmp[k,])] <- NA
  }
  
  df <- data.table(tmp)
  df <- df[, .(COUNT = .N), by = names(df)]
  site.sp.rch[i] <- mean(unlist(lapply(samples.ls, length)))
  
  counts <- df$COUNT  
  site.counts[[i]] <- counts
  
  if(length(counts) >1){
    uni <- do.call(cbind, df)
    uni <- uni[ ,1:(ncol(uni)-1)]
    uni <- asplit(uni, 1)
    uni <- lapply(uni, as.numeric)
    uni <- lapply(uni, function(x) x[!is.na(x)])
  }else{
    uni <- list(tmp[1, ])
  }
  
  site.sps.16ha[[i]] <- uni
  print(i)
}

save(site.sps.16ha, file = 'Analysis/Outputs/site_sps_16ha.RData')
save(site.counts, file = 'Analysis/Outputs/site_counts.RData')
save(site.sp.rch, file = 'Analysis/Outputs/site_sp_rch.RData')

load(file = 'Analysis/Outputs/site_counts.RData')

### Repeat but with the full data set - including rare species

for(i in 1:n.sites){
  site <- sites[i] 
  data.mat <- readRDS(file = sprintf('Data_Zone/Output/%s_all.sp.RData', site))
  
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
  # now sample 400 subplots (16 ha)
  
  samples.ls <- vector('list', n.bootstraps)
  
  for(j in 1:n.bootstraps){
    
    if(site == 'ituri'){
      quads <- sample(unique(data.mat[ ,'quadrat']), 400, replace = TRUE)
      sp.sample <- unique(data.mat[which(data.mat[ ,'quadrat'] %in% quads), 'sp'])
    } else{
      subs <- sample(unique(data.mat[ ,'sub']), 400, replace = TRUE)
      sp.sample <- unique(data.mat[which(data.mat[ ,'sub'] %in% subs), 'sp'])
    }
    samples.ls[[j]] <- sort(sp.sample)
  }
  
  # What are the unique species combinations and how many times
  # does each combination occur?
  
  tmp <- do.call(rbind, samples.ls)
  for(k in 1:nrow(tmp)){
    tmp[k, ][duplicated(tmp[k,])] <- NA
  }
  
  df <- data.table(tmp)
  df <- df[, .(COUNT = .N), by = names(df)]
  site.sp.rch.all[i] <- mean(unlist(lapply(samples.ls, length)))

  print(i)
}

save(site.sp.rch.all, file = 'Analysis/Outputs/site_sp_rch_all.RData')

