#####################################################################
### Make site_df ###
####################################################################
rm(list=ls())

geog <- 'global'
load(file = sprintf('Analysis/Outputs/%s_pca_mat.RData', geog))
load(file = sprintf('Analysis/Outputs/PCA_%s.RData', geog))

sites <- as.character(unique(pca.matrix[ ,'Site']))
site.names <- ifelse(nchar(sites) <= 4, toupper(sites), 
                     stri_trans_totitle(sites))
site.names[which(site.names == 'Khaochong')] <- 'Khao Chong'
site.names[which(site.names == 'Laplanada')] <- 'La Planada'
site.names[which(site.names == 'Windriver')] <- 'Wind River'
site.names[which(site.names == 'XTBG')] <- 'XSBN'
n.sites <- length(site.names)

census.years <- vector('list', n.sites)

site.df <- site.names

for(i in 1:n.sites){
  site <- sites[i]
  if(site == 'mudumalai'){
    data.ls <- readRDS('Data_Zone/Output/mudumalai_spnames_years.RData')
    census.years[[i]] <- paste(data.ls$years[1], data.ls$years[2], sep = '-')
  }else{
  readRDS(file = sprintf('Data_Zone/Output/%s.RData', site))
  data.ls <- readRDS(file = sprintf('Data_Zone/Output/%s.RData', site))
  # use try since not all sites have all components 
  try(data.mat <- data.ls$data.mat, silent = TRUE)
  try(sp.names <- data.ls$sp.names, silent = TRUE)
  try(sp.latins <- data.ls$sp.latins, silent = TRUE)
  try(years <- data.ls$years, silent = TRUE)
  
  census.years[[i]] <- years
  }
}

census.years <- unlist(lapply(census.years, function(x) paste(x[1], x[2], sep = '-')))


Region <-  c('Neotropic','Neotropic','Palearctic','Indo-Malay','Indo-Malay','Afrotropic',
             'Indo-Malay','Afrotropic','Indo-Malay','Neotropic','Neotropic','Indo-Malay',
             'Indo-Malay','Indo-Malay','Nearctic' ,'Nearctic' ,'Nearctic' ,'Palearctic',
             'Indo-Malay', 'Neotropic') 
Latitude <- c(-3.8091, 9.1543,  42.3833, 24.7614, 15.6324, 1.4368,  7.54347,
              5.07389, 4.1865,  1.1558,  18.3262, 11.5989, 17.0402, 2.982,  
              38.8935, 38.8891, 45.8197, 51.7743, 21.6117, -0.6859)
Longitude <- c(-70.26780,  -79.84610,  128.08300,  121.55500,   99.21700,   28.58260,   
               99.79800,
               8.85472,  114.01700,  -77.99350,  -65.81600,   76.53380,  122.38800,  102.31300,
               -78.14540,  -76.55940, -121.95580,   -1.33790,  101.57400,  -76.39700)
Area <- c(25.0, 50.0, 25.0, 25.0, 50.0, 40.0, 24.0, 55.0, 52.0, 25.0, 16.0, 50.0, 16.0, 50.0, 
          25.6, 16.0, 27.2, 18.0, 20.0, 25.0)
No.species <- c(1133,  299,   52,  110,  251,  445,  593,  494, 1182,  240,  138, 
                72,  335,  814,   64,   79, 26,   23,  468, 1114)
Koppen.climate.zone <-  c('Af',  'Am',  'Dwb', 'Cfa', 'Aw',  'Af',  'Am',  'Am',
                          'Af',  'Cfb', 'Am',  'Aw',  'Af', 'Af',  'Cfa', 'Cfa',
                          'Csb', 'Cfb', 'Cwa', 'Af' )
MAT <-  c(25.8, 27.1,  2.9, 18.2, 23.5, 24.3, 27.1, 26.6, 26.6, 19.0, 22.8, 22.7,
          26.1, 27.9, 12.9, 13.2, 9.2, 10.0, 21.8, 28.3)
MAP <- c(3215, 2551,  700, 4271, 1476, 1682, 2611, 5272, 2664, 4087, 3548, 1255,
         3380, 1788, 1001, 1068, 2495,  717, 1493, 3081)
Disturbance.type <- c('Fl', 'D; W', '', 'H', 'Fi; D', 'W; A', 'W;L', 'W', 'L; D',
                      'W', 'H; L','Fi; A; D','H','W','W, Ic','H; W','Fi; W; In', 
                      '', 'W; D' ,'_-')

site.df <- as.data.frame(cbind(site.names, census.years, Region, Latitude, Longitude, Area,
                 No.species, Koppen.climate.zone, MAT, MAP, Disturbance.type))
site.df$Longitude <- as.numeric(site.df$Longitude)
site.df$Latitude <- as.numeric(site.df$Latitude)
site.df$Area <- as.numeric(site.df$Area)
site.df$MAT <- as.numeric(site.df$MAT)
site.df$MAP <- as.numeric(site.df$MAP)
site.df$No.species <- as.numeric(site.df$No.species)

colnames(site.df)[1] <- 'Site'
save(site.df, file = 'Analysis/Outputs/site_df.RData')
