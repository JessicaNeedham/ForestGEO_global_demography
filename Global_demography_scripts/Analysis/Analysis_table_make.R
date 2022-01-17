### One table with all values of DD and DC for analysis
rm(list = ls())

lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

geog <- 'global'
load(file = sprintf('Analysis/Outputs/%s_pca_mat.RData', geog))
load(file = sprintf('Analysis/Outputs/PCA_%s.RData', geog))
load('Analysis/Outputs/Site_df.RData')  # MAT and MAP
load(file = 'Analysis/Outputs/Site_median_hulls_2.RData')  # DC
ch <- site.median.hulls_2
load(file = 'Analysis/Outputs/site_sp_rch.RData')  # Sp richness
load(file = 'Analysis/Outputs/AGBs.RData')  # AGB
load(file = 'Analysis/Outputs/taus.RData')

n.sites <- length(unique(pca.matrix$Site))-1
n.clust <- max(pca.matrix$clust)+1

gsms <- matrix(nrow = n.sites, ncol = n.clust)
sites <- as.character(unique(pca.matrix[ ,'Site']))
sites <- sites[-which(sites == 'mudumalai')]

# proportion of GSMs
for(i in 1:n.sites){
  site <- sites[i]
  # load sum of AGB of each species
  load(file = sprintf('Analysis/Outputs/%s_agb_sp.RData', site))
  # match to GSM
  site.pca <- pca.matrix[pca.matrix$Site == site, ]
  
  agbs <- data.frame(names(sp.agb.sums))
  agbs$agb <- as.numeric(sp.agb.sums)
  if(site == 'amacayacu'){
    agbs$clust <- site.pca[match(names(sp.agb.sums), site.pca[ ,'Species']) ,'clust']
  }else{
    agbs$clust <- site.pca[match(names(sp.agb.sums), site.pca[ ,'Mnemonic']) ,'clust']
  }
  agbs$clust[which(is.na(agbs$clust))] <- 9
  sums <- tapply(agbs$agb, agbs$clust, sum)
  vec <- rep(NA, 9)
  names(vec) <- seq(9)
  vec[match(names(sums), names(vec))] <- sums
  # vec[is.na(vec)] <- 0
  gsms[i, ] <- vec
  print(i)
}

# normalise it
rownames(gsms) <- sites
gsms_norm <- gsms[ ,1:8]
gsms_norm <- apply(gsms_norm, 1, function(x) x/sum(x, na.rm=TRUE))

test.tab <- site.df[ ,c('Site', 'MAT', 'MAP')]
test.tab <- test.tab[-which(test.tab[ ,'Site'] == 'Mudumalai'), ]

test.tab <- cbind(test.tab, ch, t(gsms_norm), site.sp.rch, AGBs, tau)

colnames(test.tab) <- c('Site', 'MAT', 'MAP', 'DD', '1','2','3','4','5','6','7','8',
                        'Sp_richness', 'AGB', 'C_res')

save(test.tab, file = 'Analysis/Outputs/test_tab.RData')

###################################################################################
rm(list = ls())
load(file = 'Analysis/Outputs/test_tab.RData')

lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

# DC and MAT 
MAT <- test.tab$MAT

mat1 <- lm(test.tab[,'1'] ~ MAT)
mat1p <- lmp(mat1)
mat2 <- lm(test.tab[,'2'] ~ MAT)
mat2p <- lmp(mat2)
mat3 <- lm(test.tab[,'3'] ~ MAT)
mat3p <- lmp(mat3)
mat4 <- lm(test.tab[,'4'] ~ MAT)
mat4p <- lmp(mat4)
mat5 <- lm(test.tab[,'5'] ~ MAT)
mat5p <- lmp(mat5)
mat6 <- lm(test.tab[,'6'] ~ MAT)
mat6p <- lmp(mat6)
mat7 <- lm(test.tab[,'7'] ~ MAT)
mat7p <- lmp(mat7)
mat8 <- lm(test.tab[,'8'] ~ MAT)
mat8p <- lmp(mat8)

matps <- c(mat1p,mat2p,mat3p,mat4p,mat5p,mat6p,mat7p,mat8p)
p.adj.mat <- p.adjust(matps, 'fdr')

summary(mat2)$r.squared

# DC and MAP
MAP <- test.tab$MAP

map1 <- lm(test.tab[,'1'] ~ MAP)
map1p <- lmp(map1)
map2 <- lm(test.tab[,'2'] ~ MAP)
map2p <- lmp(map2)
map3 <- lm(test.tab[,'3'] ~ MAP)
map3p <- lmp(map3)
map4 <- lm(test.tab[,'4'] ~ MAP)
map4p <- lmp(map4)
map5 <- lm(test.tab[,'5'] ~ MAP)
map5p <- lmp(map5)
map6 <- lm(test.tab[,'6'] ~ MAP)
map6p <- lmp(map6)
map7 <- lm(test.tab[,'7'] ~ MAP)
map7p <- lmp(map7)
map8 <- lm(test.tab[,'8'] ~ MAP)
map8p <- lmp(map8)

mapps <- c(map1p,map2p,map3p,map4p,map5p,map6p,map7p,map8p)

p.adj.map <- p.adjust(mapps, 'fdr')

# DC and AGB
AGB <- test.tab$AGB

AGB1 <- lm(test.tab[,'1'] ~ AGB)
AGB1p <- lmp(AGB1)
AGB2 <- lm(test.tab[,'2'] ~ AGB)
AGB2p <- lmp(AGB2)
AGB3 <- lm(test.tab[,'3'] ~ AGB)
AGB3p <- lmp(AGB3)
AGB4 <- lm(test.tab[,'4'] ~ AGB)
AGB4p <- lmp(AGB4)
AGB5 <- lm(test.tab[,'5'] ~ AGB)
AGB5p <- lmp(AGB5)
AGB6 <- lm(test.tab[,'6'] ~ AGB)
AGB6p <- lmp(AGB6)
AGB7 <- lm(test.tab[,'7'] ~ AGB)
AGB7p <- lmp(AGB7)
AGB8 <- lm(test.tab[,'8'] ~ AGB)
AGB8p <- lmp(AGB8)


AGBps <- c(AGB5p,AGB6p)
p.adj.AGB <- p.adjust(AGBps, 'fdr')


# DC and Tau
tau <- test.tab$C_res
why.ind <- which(test.tab$Site == 'Wytham')

tau1 <- lm(test.tab[-why.ind,'1'] ~ tau[-why.ind])
tau1p <- lmp(tau1)
tau2 <- lm(test.tab[-why.ind,'2'] ~ tau[-why.ind])
tau2p <- lmp(tau2)
tau3 <- lm(test.tab[-why.ind,'3'] ~ tau[-why.ind])
tau3p <- lmp(tau3)
tau4 <- lm(test.tab[-why.ind,'4'] ~ tau[-why.ind])
tau4p <- lmp(tau4)
tau5 <- lm(test.tab[-why.ind,'5'] ~ tau[-why.ind])
tau5p <- lmp(tau5)
tau6 <- lm(test.tab[-why.ind,'6'] ~ tau[-why.ind])
tau6p <- lmp(tau6)
tau7 <- lm(test.tab[-why.ind,'7'] ~ tau[-why.ind])
tau7p <- lmp(tau7)
tau8 <- lm(test.tab[-why.ind,'8'] ~ tau[-why.ind])
tau8p <- lmp(tau8)

tau.ps <- c(tau5p,tau6p)
p.adj.tau <- p.adjust(tau.ps, 'fdr')
p.adj.tau 

tau.ps <- c(tau1p, tau2p, tau3p, tau4p, tau5p,tau6p, tau7p, tau8p)
p.adj.tau <- p.adjust(tau.ps, 'fdr')
p.adj.tau 



spr <- test.tab$Sp_richness

spr1 <- lm(test.tab[,'1'] ~ spr)
spr1p <- lmp(spr1)
spr2 <- lm(test.tab[,'2'] ~ spr)
spr2p <- lmp(spr2)
spr3 <- lm(test.tab[,'3'] ~ spr)
spr3p <- lmp(spr3)
spr4 <- lm(test.tab[,'4'] ~ spr)
spr4p <- lmp(spr4)
spr5 <- lm(test.tab[,'5'] ~ spr)
spr5p <- lmp(spr5)
spr6 <- lm(test.tab[,'6'] ~ spr)
spr6p <- lmp(spr6)
spr7 <- lm(test.tab[,'7'] ~ spr)
spr7p <- lmp(spr7)
spr8 <- lm(test.tab[,'8'] ~ spr)
spr8p <- lmp(spr8)

spr.ps <- c(spr1p,spr2p,spr3p,spr4p,spr5p,spr6p,spr7p,spr8p)

p.adj.spr <- p.adjust(spr.ps, 'fdr')
p.adj.mat
p.adj.mat
p.adj.AGB
p.adj.tau

save(p.adj.map, file = "Analysis/p_adj_map.RData")
save(p.adj.mat, file = "Analysis/p_adj_mat.RData")
save(p.adj.AGB, file = "Analysis/p_adj_AGB.RData")
save(p.adj.tau, file = "Analysis/p_adj_tau.RData")
save(p.adj.spr, file = 'Analysis/p_adj_spr.RData')

############################################################
rm(list = ls())

lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}


geog <- 'global'
load(file = sprintf('Analysis/Outputs/%s_pca_mat.RData', geog))
load(file = 'Analysis/Outputs/site_sp_rch.RData')  # Sp richness
load(file = 'Analysis/Outputs/AGBs.RData')  # AGB
load(file = 'Analysis/Outputs/taus.RData')

n.sites <- length(unique(pca.matrix$Site))-1
n.clust <- max(pca.matrix$clust)

gsms <- matrix(nrow = n.sites, ncol = n.clust)
sites <- as.character(unique(pca.matrix[ ,'Site']))
sites <- sites[-which(sites == 'mudumalai')]

for(i in 1:n.sites){
  site <- sites[i]
  site.pca <- pca.matrix[pca.matrix$Site == site, ]
  sums <- tapply(site.pca$N, site.pca$clust, sum)
  vec <- rep(NA, 8)
  names(vec) <- seq(8)
  vec[match(names(sums), names(vec))] <- sums
  gsms[i, ] <- vec
}

rownames(gsms) <- sites

load(file = 'Analysis/Outputs/test_tab.RData')


# DC and AGB
AGB <- test.tab$AGB

AGB1 <- lm(gsms[,1] ~ AGB)
AGB1p <- lmp(AGB1)
AGB2 <- lm(gsms[,2] ~ AGB)
AGB2p <- lmp(AGB2)
AGB3 <- lm(gsms[,3] ~ AGB)
AGB3p <- lmp(AGB3)
AGB4 <- lm(gsms[,4] ~ AGB)
AGB4p <- lmp(AGB4)
AGB5 <- lm(gsms[,5] ~ AGB)
AGB5p <- lmp(AGB5)
AGB6 <- lm(gsms[,6] ~ AGB)
AGB6p <- lmp(AGB6)
AGB7 <- lm(gsms[,7] ~ AGB)
AGB7p <- lmp(AGB7)
AGB8 <- lm(gsms[,8] ~ AGB)
AGB8p <- lmp(AGB8)


AGBps <- c(AGB1p, AGB2p, AGB3p, AGB4p, AGB5p,AGB6p, AGB7p, AGB8p)
p.adj.AGB <- p.adjust(AGBps, 'fdr')
p.adj.AGB

# DC and Tau
tau <- test.tab$C_res
why.ind <- which(test.tab$Site == 'Wytham')

tau1 <- lm(gsms[-why.ind,1] ~ tau[-why.ind])
tau1p <- lmp(tau1)
tau2 <- lm(gsms[-why.ind,2] ~ tau[-why.ind])
tau2p <- lmp(tau2)
tau3 <- lm(gsms[-why.ind,3] ~ tau[-why.ind])
tau3p <- lmp(tau3)
tau4 <- lm(gsms[-why.ind,4] ~ tau[-why.ind])
tau4p <- lmp(tau4)
tau5 <- lm(gsms[-why.ind,5] ~ tau[-why.ind])
tau5p <- lmp(tau5)
tau6 <- lm(gsms[-why.ind,6] ~ tau[-why.ind])
tau6p <- lmp(tau6)
tau7 <- lm(gsms[-why.ind,7] ~ tau[-why.ind])
tau7p <- lmp(tau7)
tau8 <- lm(gsms[-why.ind,8] ~ tau[-why.ind])
tau8p <- lmp(tau8)

tau.ps <- c(tau5p,tau6p)
p.adj.tau <- p.adjust(tau.ps, 'fdr')
p.adj.tau 

tau.ps <- c(tau1p, tau2p, tau3p, tau4p, tau5p,tau6p, tau7p, tau8p)
p.adj.tau <- p.adjust(tau.ps, 'fdr')
p.adj.tau 




