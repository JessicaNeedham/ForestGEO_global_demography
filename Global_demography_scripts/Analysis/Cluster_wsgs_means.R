
rm(list = ls())

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
n.clust <- length(unique(pca.matrix[, 'clust']))

site.list <- vector('list', n.sites)
# loop through sites and add wsg to each species in pca matrix
for(jj in 1:n.sites){
  site <- sites[jj]
  load(file = sprintf('Vitals_Zone/Output/%s_survp_table.RData', site))
  load(file = sprintf('Vitals_Zone/Output/%s_growthp_table.RData', site))
  all.sp <- readRDS(file = sprintf('Data_Zone/Output/%s_all.sp.RData', site))
  site.pca <- pca.matrix[which(pca.matrix[ ,'Site'] == site), ]
  site.pca$wsg <- all.sp$wsg[match(site.pca$Species, all.sp$Latin)]
  site.list[[jj]] <- site.pca
}

all <- do.call(rbind, site.list)

clust.wsgs <- tapply(all$wsg, all$clust, mean, na.rm = TRUE)

save(clust.wsgs, file = 'Analysis/Outputs/clust_wsgs.RData')
