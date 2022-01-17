######################################################
rm(list = ls())

library(stringi)
library(RColorBrewer)
library(factoextra)
library(FactoMineR)
library(ade4)
library(corrplot)
library(npreg)

source('Analysis/Function_files/PCA_functions.R')
source('Analysis/Function_files/Ade4_functions.R')
#source('Paper_Figures_v2/Violin_plot.R')
#options(device = 'quartz')
set.seed(1)

lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

########################################################
geog <- 'global'
load(file = sprintf('Analysis/Outputs/%s_pca_mat.RData', geog))
load(file = sprintf('Analysis/Outputs/PCA_%s.RData', geog))

# Number of species
n.species <- length(unique(pca.matrix[ ,'Species']))
n.species
n.stems <- sum(pca.matrix[ ,'N'])
# Number of sites
sites <- as.character(unique(pca.matrix[ ,'Site']))
n.sites <- length(sites)
n.sites

# latitude extremes
load('Analysis/Outputs/Site_df.RData')
range(abs(site.df[ ,'Latitude']))
site.df[ ,c('Site', 'Latitude')]


sites <- as.character(unique(pca.matrix[ ,'Site']))
site.names <- ifelse(nchar(sites) <= 4, toupper(sites), 
                     stri_trans_totitle(sites))
site.names[which(site.names == 'Khaochong')] <- 'Khao Chong'
site.names[which(site.names == 'Laplanada')] <- 'La Planada'
site.names[which(site.names == 'Windriver')] <- 'Wind River'
site.names[which(site.names == 'XTBG')] <- 'XSBN'
site.names <- site.names[-which(site.names == 'Mudumalai')]

sites <- sites[-which(sites == 'mudumalai')]
n.sites <- length(sites)
incl <- rep(NA, n.sites)

# percent of stems included
# for(i in 1:n.sites){
#   site.name <- site.names[i]
#   site <- sites[i]
#   all.sp <- readRDS(file = sprintf('Data_Zone/Output/%s_all.sp.RData', site))
#   site.pca <- pca.matrix[pca.matrix[ ,'Site'] == site, ]
#   # which were in the PCA? 
#   all.sp$pca <- match(all.sp$Latin, site.pca[ ,'Species'])
#   incl[i] <- length(which(all.sp$pca > 0))/nrow(all.sp)
#   
# }
# incl
# range(incl)
# mean(incl)


# number of stems
sum(pca.matrix[ ,'N'])

# species in more than one sites
mp.sites.sp <- names(which(table(pca.matrix[ ,'Species']) >1))
# number of species in multiple sites
length(mp.sites.sp)
# as a percent
length(mp.sites.sp)/n.species*100

sites <- as.character(unique(pca.matrix[ ,'Site']))
site.names <- ifelse(nchar(sites) <= 4, toupper(sites), 
                     stri_trans_totitle(sites))
site.names[which(site.names == 'Khaochong')] <- 'Khao Chong'
site.names[which(site.names == 'Laplanada')] <- 'La Planada'
site.names[which(site.names == 'Windriver')] <- 'Wind River'
site.names[which(site.names == 'XTBG')] <- 'XSBN'

# proportion of species in slow growing high survival 
table(pca.matrix[,'clust'])

table(pca.matrix[ ,'clust'])/nrow(pca.matrix)*100

# proportion of individuals in GSSMs 5 and 6 in lambir and windriver
lambir <- pca.matrix[which(pca.matrix[ ,'Site'] == 'lambir'), ]
tmp <- tapply(lambir[ ,'N'], lambir[ ,'clust'], sum)
(tmp[5] + tmp[6])/sum(tmp)

wr <- pca.matrix[which(pca.matrix[ ,'Site'] == 'windriver'), ]
tmp <- tapply(wr[ ,'N'], wr[ ,'clust'], sum, na.rm = TRUE)
(tmp[4] + tmp[5])/sum(tmp, na.rm = TRUE)


###########################################################
# 1. growth survival state space
###########################################################
var.names <- c('MS', 'JS', 'ST', 
               '95G', '5G')
var.full.names <- c('Maximum survival', 'Juvenile survival', 
                    'Stature', '95% slowest growth', '5% fastest growth')
arrow.cols <- c('red', 'red', 'darkorange', 'darkblue', 'darkblue')

cols <- c(brewer.pal(7, 'Dark2')[c(6,1,2,3,4,5,7)], brewer.pal(3, 'Set1')[2])
cols.light <- make.transp(cols)
cols.v.light <- make.transp(cols, c(30,20,40,30,30,50,70,80))

cols.points <- cols.v.light
cols.lines <- make.transp(cols, c(25,10,20,20,20,50,70,80))

ind.coords <- res$ind$coord
var.coords <- res$var$coord

clusters <- pca.matrix[ ,'clust']
Ns <- pca.matrix[ ,'N']

# FIg 1
pdf('Paper_Figures_v4/Figure_2.pdf', width =10, height = 7)
#jpeg('Paper_Figures_v3/PCA_plot.jpeg', width = 650, height = 400)

#par(mfrow = c(1,1))
mat <- matrix(c(1,2), nrow = 1, ncol = 2)
layout(mat, width = c(0.7, 0.3))
par(mar = c(5,4,1,0), xpd = FALSE)

custom_biplot(res, axes = c(1,2), cols = cols.points, clusters = clusters, 
              var.names = var.names, pch = 16, lwd = 2, cex.text = 1, 
              plot.axes = TRUE, plot.axes.labs = TRUE, exclude.outliers = FALSE, 
              cex = sqrt(Ns)*0.02, lab.pos = c(3,3,3,4,4), arrow.colours = arrow.cols, 
              xlim = c(-3, 11), ylim = c(-5, 4))
scatterutil.grid(1)
abline(h = 0, v = 0, lty = 1)

dfxy <- data.frame(ind.coords[ ,1:2])
dfdistri <- fac2disj(factor(clusters)) * rep(1, length(clusters))
coul <- col
w1 <- unlist(lapply(dfdistri, sum))
dfdistri <- t(t(dfdistri)/w1)
coox <- as.matrix(t(dfdistri)) %*% dfxy[,1]
cooy <- as.matrix(t(dfdistri)) %*% dfxy[,2]

for (i in 1:ncol(dfdistri)) {
  scatterutil.star(dfxy[ ,1], dfxy[ ,2], dfdistri[, i], cstar = 1, 
                   cols.lines[i])
}
#scatterutil.eti(coox, cooy, label = seq(length(cols)), clabel = 0.6, 
#               coul = cols)
text(coox, cooy, label = seq(length(cols)), cex =1)
#plot.new()
legend('bottomright', cex =0.8, col = arrow.cols, lwd = 1,
       legend = paste(var.names, var.full.names, sep = ': '))

mtext('fast growth', side = 1, line = 2.5, adj = 0.95)
mtext('slow growth', side = 1, line = 2.5, adj = 0.05)
mtext(bquote(atop('small stature', 'low survival')), 
      side = 2, line = 1.3, adj = 0.05)
mtext(bquote(atop('large stature', 'high survival')), 
      side = 2, line = 1.3, adj = 0.95)

par(xpd = TRUE)
arrows(x0 = 7, x1 = 8.6, y0= -6.4, y1 = -6.4, 
       length = 0.1)
arrows(x0 = 1, x1 = -0.5, y0= -6.4, y1 = -6.4, 
       length = 0.1)
arrows(x0 = -5, x1 = -5, y0= 1, y1 = 2.1, 
       length = 0.1)
arrows(x0 = -5, x1 = -5, y0= -2, y1 = -3.2, 
       length = 0.1)

par(mar = c(0,0,0,0))
plot.new()

GSSMs <- c('very slow-growing', 
          'slow-growing, low survival', 
          'very low survival', 
          'intermediate growth and survival', 
          'large statured, high survival',
          'very large statured, fast growth', 
          'fast growth, low survival', 
          'pioneers - very fast growth')

legend('left', legend = c(seq(8), GSSMs), cex = 0.9, 
       text.col = c(cols, rep('black', 8)), ncol = 2, 
       text.width = c(0), bty = 'n')

dev.off()

############## FIG 2 ########################
pdf(file = 'Paper_Figures_v4/Figure_3.pdf')
par(mfrow = c(4,5), oma = c(0,0,0,0), mar = c(0,0,0,0))
load(file = sprintf('Analysis/Outputs/PCA_%s.RData', geog))
load(file = 'Analysis/Outputs/example_species.RData')
site.coords <- lapply(split(res$ind$coord, pca.matrix[ ,'Site']), 
                      matrix, ncol = ncol(res$ind$coord))
site.clusters <- split(pca.matrix[ ,'clust'], pca.matrix[ ,'Site'])
site.Ns <- split(pca.matrix[ ,'N'], pca.matrix[ ,'Site'])

n.sites <- length(site.coords)

load('Analysis/Outputs/site_df.RData')
ord <- order(abs(site.df[ ,'Latitude']))

bg.pt.col <- brewer.pal(9, 'Greys')[3]
bg.pt.col <- make.transp(bg.pt.col)
bg.line.col <- brewer.pal(9, 'Greys')[5]
bg.line.col <- make.transp(bg.line.col, 50)
pt.col <- make.transp(cols, c(80,80,50,80,95,80,80,80))

load(file = 'Analysis/Outputs/convex_hull_coords.RData')
site.names <- site.df[ ,'Site']

for(ss in ord){
  
  site <- site.coords[[ss]]
  clusts <- rep(1, nrow(site))
  Ns <- rep(1000, nrow(site))
  
  cp <- cols.light
  cl <- cols.light
  
  s.class.hack(res$ind$coord[ ,1:2], fac = factor(rep(1, nrow(res$ind$coord))), 
               grid = FALSE, col.points = bg.pt.col, label = '', 
               cstar = 0, cpoint = 0.25, cellipse = 0, cex.point = 0.8, 
               xlim = c(-3, 12), ylim = c(-5, 5), axesell = FALSE)
  
  s.class.hack(site, fac = factor(rep(1, nrow(site))), 
               col.points = cp[clusts], cpoint = 0, 
               col.lines = rep(bg.line.col, length(clusts)),
               grid = FALSE, 
               label = '', cellipse = 0, 
               add.plot = TRUE, axesell = FALSE, 
               cex.point = sqrt(Ns)*0.02)
  
  s.class.hack(site, fac = factor(rep(1, nrow(site))), 
               col.points = pt.col[clusts], cstar = 0, 
               col.lines = rep(invis, length(clusts)),
               grid = FALSE, 
               label = '', cellipse = 0, 
               add.plot = TRUE, axesell = FALSE, 
               cex.point = sqrt(Ns)*0.02)
  mtext(site.names[ss], side = 3, line =-2, adj = 0.8, cex = 1)
  print(ss)
  polygon(x = convex_hull_coords[[ss]][ ,1], y = convex_hull_coords[[ss]][,2], 
          col = cols.v.light[1], border = cols.light[1], lwd =2)
  
}

dev.off()

############################################################
# 1. Does DD  and DC vary among forests around the world? 
# Does DD and DC vary with species richness
############################################################
load(file = 'Analysis/Outputs/Site_median_hulls_2.RData')  # DC
ch <- site.median.hulls_2

# 2. Are forests with more species more demographically diverse?
load(file = 'Analysis/Outputs/site_sp_rch.RData')
load(file = 'Analysis/Outputs/site_sp_rch_all.RData')

# 3. Are forests similar in demographic composition? 
n.sites <- length(unique(pca.matrix$Site))-1
n.clust <- max(pca.matrix$clust)+1

GSSMs <- matrix(nrow = n.sites, ncol = n.clust)
sites <- sites[-which(sites == 'mudumalai')]

# proportion of GSSMs
for(i in 1:n.sites){
  site <- sites[i]
  # load sum of AGB of each species
  load(file = sprintf('Analysis/Outputs/%s_agb_sp.RData', site))
  # match to GSSM
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
  GSSMs[i, ] <- vec
  print(i)
}

# normalise it
rownames(GSSMs) <- sites
GSSMs_norm <- GSSMs[ ,1:8]
GSSMs_norm <- apply(GSSMs_norm, 1, function(x) x/sum(x, na.rm=TRUE))
site.names <- site.names[-which(site.names == 'Mudumalai')]

dist <- c(4,5,11,12)

# Fig 3 
pdf('Paper_Figures_v4/Figure_4.pdf', width = 7, height = 5)
par(mfrow = c(1,1), mar = c(4,4,2,1), oma = c(0,0,0,0))

spl <- ss(log(site.sp.rch), ch, nknots = 6)
plot(spl, xlab = '', ylab = '', las=1, main = NA, xaxt = 'n', 
     ylim =c(0, max(ch)), xlim = c(1.9,6.05))
points(log(site.sp.rch), ch, col = 'skyblue3', pch = 16)
#points(log(site.sp.rch)[dist], ch[dist], col = 'goldenrod2', pch = 16)
mtext('Species Richness', side = 1, line = 2.5, cex = 1.2)
mtext('Demographic diversity', side = 2, line = 2.5, cex = 1.2)

text(log(site.sp.rch), ch, cex = 0.8, labels = site.names, 
     pos = c(2, 1, 4, 1, 2,1,4,2,2,4,2,2,1,2,4,4,3,3,1))

vals <- c(10, 20, 50, 100, 200, 500)
axis(side = 1, at = log(vals), labels = vals)

# for the full species diversity 
#text(site.sp.rch.all, ch, cex = 0.6, labels = site.names, 
#    pos = c(1, 2, 4, 4, 1,1,1,2,1,3,4,4,1,2,2,1,3,3,3))
dev.off()

#####################################################
# 4. Do DD and DC vary with climate variables? 
#####################################################
load(file = 'Analysis/Outputs/site_df.RData')

site.df <- site.df[-which(site.df$Site == 'Mudumalai'), ]

pdf('Paper_Figures_v4/Figure_5.pdf', width = 8)
#png('Paper_Figures_v3/Fig_4.png', width =10)

par(mfrow = c(2,2), mar = c(5,5,2,1), oma = c(0,0,0,1))
mat <- site.df[ ,'MAT']
plot(mat, site.median.hulls_2, xlab = '', ylab = '', las = 1, pch = 16, 
     col= 'skyblue3')
#mtext('Mean Annual Temperature (degrees C)', side = 1, line = 2.5, cex = 0.8)
mtext('Demographic diversity', side = 2, line = 2.5, cex = 0.8)
mtext('(a)', side = 3, line = 0.2, adj = 0.05)
m1 <- lm(ch~mat)
p1 <- lmp(m1)
legend(1, 8, border = NULL, bty = 'n',
       legend = paste('Rsq = ', signif(summary(m1)$r.squared, 2), 
                      'p = ', signif(p1, 2)), 
       text.col = 'black')
xs <- seq(min(mat), max(mat), 1)
ys <- m1$coefficients[1] + m1$coefficients[2]*xs
points(xs, ys, type = 'l', col = 'darkgrey', lty = 1)

tmp.ind <- which(sites %in% c('changbaishan', 'serc', 'scbi', 'windriver', 'wytham'))
trp.ind <- seq(n.sites)[-tmp.ind]

m1b <- lm(ch[trp.ind]~mat[trp.ind])
p1b <- lmp(m1b)
legend(10,1.75, border = NULL, bty = 'n',
       legend = paste('Rsq = ', signif(summary(m1b)$r.squared, 2), 
                      'p = ', signif(p1b, 2)), 
       text.col = 'darkgoldenrod4')
xs <- seq(min(mat[trp.ind]), max(mat[trp.ind]), 1)
ys <- m1b$coefficients[1] + m1b$coefficients[2]*xs
points(xs, ys, type = 'l', col = 'darkgoldenrod4', lty = 1)

plot(mat, GSSMs_norm[6, ], xlab = NA, ylab = NA, 
     pch = 17, col = cols[6], las = 1)
mtext('Mean Annual Temperature (degrees C)', side = 1,adj = -3, line = 2.5, cex = 0.8)
mtext('Relative abundance of GSSM 6', side = 2, line = 2.5, cex = 0.8)
mtext('(b)', side = 3, line = 0.2, adj = 0.05)
m2 <- lm(GSSMs_norm[6,]~mat)
p2 <- lmp(m2)
legend('topright', border = NULL, bty = 'n',
       legend = paste('Rsq = ', signif(summary(m2)$r.squared, 2), 
                      'p = ', signif(p2, 2)), 
       text.col = 'black')
xs <- seq(min(mat), max(mat), 1)
ys <- m2$coefficients[1] + m2$coefficients[2]*xs
points(xs, ys, type = 'l', col = 'darkgrey', lty = 1)

map <- site.df[ ,'MAP']
plot(map, site.median.hulls_2, xlab = '', ylab = '', las = 1, pch = 16, 
     col= 'skyblue3')
#mtext('Mean Annual Precipitation (mm)', side = 1, line = 2.5, cex = 0.8)
mtext('Demographic diversity', side = 2, line = 2.5, cex = 0.8)
mtext('(c)', side = 3, line = 0.2, adj = 0.05)
m3 <- lm(ch~map)
p3 <- lmp(m3)
legend('topright', border = NULL, bty = 'n',
       legend = paste('Rsq = ', signif(summary(m3)$r.squared, 2), 
                      'p = ', signif(p3, 2)), 
       text.col = 'black')
xs <- seq(min(map), max(map), 1)
ys <- m3$coefficients[1] + m3$coefficients[2]*xs
points(xs, ys, type = 'l', col = 'darkgrey', lty = 1)

plot(map, GSSMs_norm[6, ], xlab = NA, ylab = NA, 
     pch = 17, col = cols[6], las = 1)
mtext('Mean Annual Precipitation (mm)', side = 1, line = 2.5, 
      cex = 0.8, adj = -1.5)
mtext('Relative abundance of GSSM 6', side = 2, line = 2.5, cex = 0.8)
mtext('(d)', side = 3, line = 0.2, adj = 0.05)
m4 <- lm(GSSMs_norm[6,]~map)
p4 <- lmp(m4)
legend('topright', border = NULL, bty = 'n',
       legend = paste('Rsq = ', signif(summary(m4)$r.squared, 2), 
                      'p = ', signif(p4, 2)), 
       text.col = 'black')
xs <- seq(min(map), max(map), 100)
ys <- m4$coefficients[1] + m4$coefficients[2]*xs
points(xs, ys, type = 'l', col = 'darkgrey', lty = 1)

dev.off()

##############################################################
# 5. Are DD and DC related to forest structure and function? 
##############################################################
load(file = 'Analysis/Outputs/AGBs.RData')
load(file = 'Analysis/Outputs/taus.RData')


pdf('Paper_Figures_v4/Figure_6.pdf', width = 10)

par(mfrow = c(2,3), mar = c(4,4,3,1))
plot(ch, AGBs, xlab = NA, ylab = NA, 
     pch = 16, col = 'deepskyblue3')
mtext('Demographic Diversity', side = 1, line = 2.5, cex = 0.8)
mtext('AGB (kg C ha-1)', side = 2, line = 2.5, cex = 0.8)
mtext('(a)', side = 3, line = 0.2, adj = 0.05)
m5 <- lm(AGBs~ch)
p5 <- lmp(m5)
legend('topright', border = NULL, bty = 'n',
       legend = paste('Rsq = ', signif(summary(m5)$r.squared, 2), 
                      'p = ', signif(p5, 2)), 
       text.col = 'black',cex=1.2)
xs <- seq(min(ch), max(ch), 0.1)
ys <- m5$coefficients[1] + m5$coefficients[2]*xs
points(xs, ys, type = 'l', col = 'darkgrey', lty = 1)
signif(p5, 2)

plot(GSSMs_norm[6, ], AGBs, xlab = NA, ylab = NA, 
     pch = 17, col= cols[6])
mtext('Relative abundance of GSSM 6', side = 1, line = 2.5, cex = 0.8)
#mtext('AGB (kg C ha-1)', side = 2, line = 2.5, cex = 1)
mtext('(b)', side = 3, line = 0.2, adj = 0.05)
AGB6 <- lm(AGBs~GSSMs_norm[6,])
AGB6.p <- lmp(AGB6)
legend('topright', border = NULL, bty = 'n',
       legend = paste('Rsq = ', signif(summary(AGB6)$r.squared, 2), 
                      'p = ', signif(AGB6.p, 2)), 
       text.col = 'black',cex=1.2)
xs <- seq(min(GSSMs_norm[6,]), max(GSSMs_norm[6,]), 0.1)
ys <- AGB6$coefficients[1] + AGB6$coefficients[2]*xs
points(xs, ys, type = 'l', col = 'darkgrey', lty = 1)


plot(GSSMs_norm[5, ], AGBs, xlab = NA, ylab = NA, 
     pch = 16, col= cols[5])
mtext('Relative abundance of GSSM 5', side = 1, line = 2.5, cex = 0.8)
#mtext('AGB (kg C ha-1)', side = 2, line = 2.5, cex = 1)
mtext('(c)', side = 3, line = 0.2, adj = 0.05)
AGB5 <- lm(AGBs~GSSMs_norm[5,])
AGB5.p <- lmp(AGB5)
legend('topleft', border = NULL, bty = 'n',
       legend = paste('Rsq = ', signif(summary(AGB5)$r.squared, 2), 
                      'p = ', signif(AGB5.p, 2)), 
       text.col = 'black',cex=1.2)
xs <- seq(min(GSSMs_norm[5,]), max(GSSMs_norm[5,]), 0.1)
ys <- AGB5$coefficients[1] + AGB5$coefficients[2]*xs
points(xs, ys, type = 'l', col = 'darkgrey', lty = 1)

plot(ch, tau, xlab = NA, ylab = NA,
     pch = 16, col = 'deepskyblue3')
mtext('Demographic Diversity', side = 1, line = 2.5, cex = 0.8)
mtext('Carbon residence time (years)', side = 2, line = 2.5, cex = 0.8)
mtext('(d)', side = 3, line = 0.2, adj = 0.05)
why.ind <- which(sites == 'wytham')
m7 <- lm(tau[-why.ind]~ch[-why.ind])
p7 <- lmp(m7)
legend('topright', border = NULL, bty = 'n',
       legend = paste('Rsq = ', signif(summary(m7)$r.squared, 2), 
                      'p = ', signif(p7, 2)), 
       text.col = 'black',cex=1.2)
xs <- seq(min(ch), max(ch), 0.1)
ys <- m7$coefficients[1] + m7$coefficients[2]*xs
points(xs, ys, type = 'l', col = 'darkgrey', lty = 1)
signif(p7, 2)

plot(GSSMs_norm[6, ], tau, xlab = NA, ylab = NA, 
     pch = 17, col= cols[6])
mtext('Relative abundance of GSSM 6', side = 1, line = 2.5, cex = 0.8)
#mtext('Carbon residence time (years)', side = 2, line = 2.5, cex = 1)
mtext('(e)', side = 3, line = 0.2, adj = 0.05)
m8 <- lm(tau[-why.ind]~GSSMs_norm[6,-why.ind])
p8 <- lmp(m8)
legend('topleft', border = NULL, bty = 'n',
       legend = paste('Rsq = ', signif(summary(m8)$r.squared, 2), 
                      'p = ', signif(p8, 2)), 
       text.col = 'black',cex=1.2)
xs <- seq(min(GSSMs_norm[6,]), max(GSSMs_norm[6,]), 0.1)
ys <- m8$coefficients[1] + m8$coefficients[2]*xs
points(xs, ys, type = 'l', col = 'darkgrey', lty = 1)
signif(p8, 2)

plot(GSSMs_norm[5, ], tau, xlab = NA, ylab = NA, 
     pch = 16, col= cols[5])
mtext('Relative abundance of GSSM 5', side = 1, line = 2.5, cex = 0.8)
#mtext('Carbon residence time (years)', side = 2, line = 2.5, cex = 1)
mtext('(f)', side = 3, line = 0.2, adj = 0.05)
m5 <- lm(tau[-why.ind]~GSSMs_norm[5,-why.ind])
p5 <- lmp(m5)
legend('topright', border = NULL, bty = 'n',
       legend = paste('Rsq = ', signif(summary(m5)$r.squared, 2), 
                      'p = ', signif(p5, 2)), 
       text.col = 'black',cex=1.2)
xs <- seq(min(GSSMs_norm[5,]), max(GSSMs_norm[5,]), 0.1)
ys <- m8$coefficients[1] + m8$coefficients[2]*xs
points(xs, ys, type = 'l', col = 'darkgrey', lty = 1)
signif(p8, 2)

dev.off()

##


