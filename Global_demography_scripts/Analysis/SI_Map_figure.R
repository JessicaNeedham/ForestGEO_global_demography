########################################################
#### Map #####
#######################################################
### LIBRARIES & FUNCTIONS ###
library(RColorBrewer)
library(factoextra)
library(FactoMineR)
library(plotrix)
library(plyr)
library(stringi)
library(tmvtnorm)
library(corpcor)
library(kernlab)
library(fpc)
library(cluster)
library(viridis)
library(plot3D)
library(caTools)
library(ade4)
library(grid)
library(gridBase)
library(multcompView)
library(ggmap)
library(ggrepel) 
library(maps)
library(maptools)
library(showtext)
library(ggplot2)
########################################################
rm(list = ls())
source('PCA_Zone/Source/PCA_functions.R')
source('PCA_Zone/Source/Ade4_functions.R')
#options(device = 'quartz')
########################################################
### Fig. 1 - Global PCA with the 8 modes
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


#############################################################
######################## 
### world map with sites

# load the site data frame and get coordinates
load('PCA_Zone/Output/site_df.RData')
n.sites <- nrow(site.df)

map.points <- site.df[ ,c('Longitude', 'Latitude')]
colnames(map.points) <- c('lon', 'lat')
map.world <- map_data("world")

d <- cbind(map.points, site.names)
colnames(d) <- c('Longitude', 'Latitude', 'Site')

pdf(file = 'Paper_Figures_v3/SI_world_map.pdf', bg = 'white', 
    width = 12, height = 5)
pal <- viridis(20)
mp <- NULL
label <- site.names

#quartz()
showtext_auto()
par(mar = c(0,0,0,0))
mp <- ggplot(d, aes(Longitude, Latitude), label = site.names) +  
  geom_polygon(data = map.world, aes(x = long, y = lat, group = group)) +
  geom_point(aes(x = map.points[ ,'lon'], y = map.points[ ,'lat'], size = 0), 
             size = 1.5, col = pal[15]) + 
  coord_cartesian(ylim = c(-50, 80), xlim = c(-167, 175)) +
  geom_label_repel(aes(x = map.points[ ,'lon'], y = map.points[ ,'lat'],
                       label = label,    colour = 'black'), size = 3, color = 'black')
mp
dev.off()

##################################################################

