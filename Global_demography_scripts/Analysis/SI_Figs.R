######################################################
rm(list = ls())

library(stringi)
library(RColorBrewer)
library(factoextra)
library(FactoMineR)
library(ade4)
library(corrplot)
library(npreg)
library(multcompView)

source('Analysis/Function_files/PCA_functions.R')
source('Analysis/Function_files/Ade4_functions.R')

#source('Paper_Figures_v3/Violin_plot.R')
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

#################################################
load(file = 'Analysis/Outputs/site_sp_rch.RData')

# 3. Are forests similar in demographic composition? 
n.sites <- length(unique(pca.matrix$Site))-1
n.clust <- max(pca.matrix$clust)+1

gsms <- matrix(nrow = n.sites, ncol = n.clust)
sites <- as.character(unique(pca.matrix[ ,'Site']))
site.names <- ifelse(nchar(sites) <= 4, toupper(sites), 
                     stri_trans_totitle(sites))
site.names[which(site.names == 'Khaochong')] <- 'Khao Chong'
site.names[which(site.names == 'Laplanada')] <- 'La Planada'
site.names[which(site.names == 'Windriver')] <- 'Wind River'
site.names[which(site.names == 'XTBG')] <- 'XSBN'
site.names <- site.names[-which(site.names == 'Mudumalai')]

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

# % AGB is GSSMs 5 and 6 in lambir and windriver
(gsms_norm[5, 'lambir'] + gsms_norm[6, 'lambir'])/sum(gsms_norm[ ,'lambir'])
(gsms_norm[5, 'windriver'] + gsms_norm[6, 'windriver'])/sum(gsms_norm[ ,'windriver'],na.rm = TRUE)


cols <- c(brewer.pal(7, 'Dark2')[c(6,1,2,3,4,5,7)], brewer.pal(3, 'Set1')[2])

# SI fig.  Species richness versus demographic composition
pdf('Paper_Figures_v3/SI_Fig_DC_v_sp_rch.pdf')
par(mfrow = c(3,3), mar = c(2,2,2,1), oma = c(3,3,0,0))

plot(site.sp.rch, gsms_norm[1,], xlab = '', ylab = '', las=1, col = cols[1], 
     pch = 16, ylim = c(0,1))
m <- lm(gsms_norm[1,]~site.sp.rch)
p <- lmp(m)
legend('topright', border = NULL, bty = 'n',
       legend = paste('Rsq = ', signif(summary(m)$r.squared, 2), 
                      'p = ', signif(p, 2)), 
       text.col = 'black', cex = 1.2)
xs <- seq(min(site.sp.rch), max(site.sp.rch), 1)
ys <- m$coefficients[1] + m$coefficients[2]*xs
points(xs, ys, type = 'l', col = cols[1], lty = 1)
mtext('a)', side = 3, line = 0.2, adj = 0.05)

plot(site.sp.rch, gsms_norm[2, ], col = cols[2], pch = 16,
     xlab = '', ylab = '', 
     las=1, ylim = c(0,1))
m <- lm(gsms_norm[2,]~site.sp.rch)
p <- lmp(m)
legend('topright', border = NULL, bty = 'n',
       legend = paste('Rsq = ', signif(summary(m)$r.squared, 2),
                      'p = ', signif(p, 2)), 
       text.col = 'black', cex = 1.2)
ys <- m$coefficients[1] + m$coefficients[2]*xs
points(xs, ys, type = 'l', col = cols[2], lty = 1)
mtext('b)', side = 3, line = 0.2, adj = 0.05)

plot(site.sp.rch, gsms_norm[3, ], col = cols[3], pch = 16, 
     xlab = '', ylab = '', 
     las=1, ylim = c(0,1))
m <- lm(gsms_norm[3,]~site.sp.rch)
p <- lmp(m)
legend('topright', border = NULL, bty = 'n',
       legend = paste('Rsq = ', signif(summary(m)$r.squared, 2),
                      'p = ', signif(p, 2)), 
       text.col = 'black', cex = 1.2)
ys <- m$coefficients[1] + m$coefficients[2]*xs
points(xs, ys, type = 'l', col = cols[3], lty =1)
mtext('c)', side = 3, line = 0.2, adj = 0.05)

plot(site.sp.rch, gsms_norm[4, ], col = cols[4], pch = 16,
     xlab = '', ylab = '',
     las=1, ylim = c(0,1))
m <- lm(gsms_norm[4,]~site.sp.rch)
p <- lmp(m)
legend('topright', border = NULL, bty = 'n',
       legend = paste('Rsq = ', signif(summary(m)$r.squared, 2),
                      'p = ', signif(p, 2)), 
       text.col = 'black', cex = 1.2)
xs <- seq(min(site.sp.rch), max(site.sp.rch), 1)
ys <- m$coefficients[1] + m$coefficients[2]*xs
points(xs, ys, type = 'l', col = cols[4], lty = 1)
mtext('d)', side = 3, line = 0.2, adj = 0.05)

plot(site.sp.rch, gsms_norm[5, ], col = cols[5], pch = 16, 
     xlab = '', ylab = '', 
     las=1, ylim = c(0,1))
m <- lm(gsms_norm[5,]~site.sp.rch)
p <- lmp(m)
legend('topright', border = NULL, bty = 'n',
       legend = paste('Rsq = ', signif(summary(m)$r.squared, 2),
                      'p = ', signif(p, 2)), 
       text.col = 'black', cex = 1.2)
ys <- m$coefficients[1] + m$coefficients[2]*xs
points(xs, ys, type = 'l', col = cols[5], 1)
mtext('e)', side = 3, line = 0.2, adj = 0.05)

plot(site.sp.rch, gsms_norm[6, ], col = cols[6], pch = 16,
     xlab = '', ylab = '', 
     las=1, ylim = c(0,1))
m <- lm(gsms_norm[6,]~site.sp.rch)
p <- lmp(m)
legend('topright', border = NULL, bty = 'n',
       legend = paste('Rsq = ', signif(summary(m)$r.squared, 2),
                      'p = ', signif(p, 2)), 
       text.col = 'black', cex = 1.2)
ys <- m$coefficients[1] + m$coefficients[2]*xs
points(xs, ys, type = 'l', col = cols[6], lty = 1)
mtext('f)', side = 3, line = 0.2, adj = 0.05)

plot(site.sp.rch, gsms_norm[7, ], col = cols[7], pch = 16, 
     xlab = '', ylab = '', 
     las=1, ylim = c(0,1))
m <- lm(gsms_norm[7,]~site.sp.rch)
p <- lmp(m)
legend('topright', border = NULL, bty = 'n',
       legend = paste('Rsq = ', signif(summary(m)$r.squared, 2),
                      'p = ', signif(p, 2)), 
       text.col = 'black', cex = 1.2)
ys <- m$coefficients[1] + m$coefficients[2]*xs
points(xs, ys, type = 'l', col = cols[7], lty = 1)
mtext('g)', side = 3, line = 0.2, adj = 0.05)

plot(site.sp.rch, gsms_norm[8, ], col = cols[8], pch = 16,
     xlab = '', ylab = '',
     las=1, ylim = c(0,1))
m <- lm(gsms_norm[8,]~site.sp.rch)
p <- lmp(m)
legend('topright', border = NULL, bty = 'n',
       legend = paste('Rsq = ', signif(summary(m)$r.squared, 2),
                      'p = ', signif(p, 2)), 
       text.col = 'black', cex = 1.2)
ys <- m$coefficients[1] + m$coefficients[2]*xs
points(xs, ys, type = 'l', col = cols[8], lty = 1)
mtext('h)', side = 3, line = 0.2, adj = 0.05)

mtext('Species Richness', side = 1, line = 1, cex = 1.2, outer=TRUE)
mtext('Demographic composition', side = 2, line = 1, cex = 1.2, outer=TRUE)

dev.off()

### Without temperate plots? 
trop.ind <- seq(19)[-which(sites %in% c('changbaishan', 'wytham', 'windriver', 'scbi', 'serc'))]
plot(site.sp.rch[trop.ind], gsms_norm[6, trop.ind],
     col = cols[6], pch = 16,
     xlab = '', ylab = '', 
     las=1, ylim = c(0,1))
m <- lm(gsms_norm[6,trop.ind]~site.sp.rch[trop.ind])
p <- lmp(m)
legend('topright', border = NULL, bty = 'n',
       legend = paste('Rsq = ', signif(summary(m)$r.squared, 2)), 
       text.col = 'black', cex = 1.2)
ys <- m$coefficients[1] + m$coefficients[2]*xs
points(xs, ys, type = 'l', col = cols[6], lty = 1)
mtext('f)', side = 3, line = 0.2, adj = 0.05)
p


## ALL PCA DIMENSIONS ##
geog <- 'global'
load(sprintf('Analysis/Outputs/pca_%s.RData', geog))
load(sprintf('Analysis/Outputs/%s_pca_mat.RData', geog))
axis.pc <- round(res$eig[ ,"percentage of variance"], 2)
pc.var.names <- mapply(function(x,y) paste('PC ', x, ': ', y, '%', sep = ''), 
                       x = seq(5), y = axis.pc)
var.names <- c('MS', 'JS', 'ST', '95G', '5G')
arrow.cols <- c('red', 'red', 'darkorange', 'darkblue', 'darkblue')


pdf(file = 'Paper_Figures_v3/SI_pc_mats.pdf')
mat <- matrix(c(1,2,3,4,
                0,5,6,7,
                0,0,8,9,
                0,0,0,10), 4,4, byrow = TRUE)
layout(mat)
par(mar = c(0,0,0,0), oma= c(5,5,2,2))
counter <- 1

n.p <- length(pc.var.names)

grey.transp.col <- brewer.pal(9, 'Set1')[9]
grey.transp.col <- make.transp(grey.transp.col, 10)

for(p in 1:n.p){
  for(q in 2:n.p){
    
    if(p >= q){
      next
    }
    
    X <- res
    if(counter %in% c(1,5,8,10)){
      custom_biplot(X, axes = c(q,p), cols = grey.transp.col, 
                    clusters = rep(1, nrow(X$ind$coord)), 
                    var.names, 
                    cex = 0.6, pch = 1, lwd = 0.6, 
                    cex.text = 1, plot.axes = TRUE, 
                    plot.axes.labs = FALSE, 
                    exclude.outliers = TRUE,
                    lab.pos = 3, 
                    xlim = c(-4,5), 
                    ylim = c(-4,6),
                    cex.axis.lab = 0.6, 
                    arrow.colours = arrow.cols)
      if(p == 1 & q == 2){
        mtext(pc.var.names[p], side = 2, line = 3)
      }
      mtext(pc.var.names[q], side = 1, line = 6)
      if(p == n.p-1 & q == n.p){
        mtext(pc.var.names[q], side = 1, line = 3)
      }
    }else{
      custom_biplot(X, axes = c(q,p), cols = grey.transp.col,
                    clusters = rep(1, nrow(X$ind$coord)), 
                    var.names, 
                    cex = 0.6, pch = 1, lwd = 0.6, 
                    cex.text = 1, plot.axes = FALSE, 
                    plot.axes.labs = FALSE, 
                    xlim = c(-4,5), 
                    ylim = c(-4,6),
                    exclude.outliers = TRUE,
                    lab.pos = 3, 
                    cex.axis.lab = 0.6, 
                    arrow.colours = arrow.cols)
    }
    
    counter <- counter + 1
  }
}
dev.off()

######################################################
load('Analysis/Outputs/Overlaps.RData')
dim(overlap)
overlap[is.na(overlap)] <- 0

pdf(file = 'Paper_Figures_v3/SI_chull_overlaps.pdf', bg = 'white')
par(mfrow = c(1,1), mar = c(5,4,2,5))
corrplot(overlap, type = "upper", 
         tl.col = "black", tl.srt = 45, is.corr = FALSE, 
         diag = FALSE, hclust.method = "ward.D2")
dev.off()
#################################################

######################################
##### Stacked barplot for each plot 
######################################
cols.light <- make.transp(cols, 90)

pdf('Paper_Figures_v3/SI_barplot.pdf')
par(mfrow = c(1,1), oma = c(1,1,0,0), mar = c(4,4,1,1))
gsms_norm[which(is.na(gsms_norm))] <- 0
bp = barplot(gsms_norm, col = cols.light, names.arg = rep('', n.sites), las = 1)
text(site.names, x= bp, y = -0.01,  srt = 45, cex = 0.75, xpd = TRUE, adj = 1)
mtext('Relative abundance of GSM', side = 2, line = 3, outer = FALSE)
dev.off()


################################################
### Climate v DC
load(file = 'Analysis/Outputs/site_df.RData')
site.df <- site.df[-which(site.df$Site == 'Mudumalai'), ]
MAT <- site.df[ ,'MAT']
MAP <- site.df[ ,'MAP']

pdf('Paper_Figures_v3/SI_GSMSvClimate.pdf')
par(mfrow = c(3,3), mar = c(2,2,2,1), oma = c(3,3,0,0))

plot(MAT, gsms_norm[1,], xlab = '', ylab = '', las=1, col = cols[1], 
     pch = 15, ylim = c(0,1))
m <- lm(gsms_norm[1,]~MAT)
p <- lmp(m)
legend('topright', border = NULL, bty = 'n',
       legend = paste('Rsq = ', signif(summary(m)$r.squared, 2), 
                      'p = ', signif(p, 2)), 
       text.col = 'black', cex = 1.2)
xs <- seq(min(MAT), max(MAT), 1)
ys <- m$coefficients[1] + m$coefficients[2]*xs
points(xs, ys, type = 'l', col = cols[1], lty = 1)
mtext('a)', side = 3, line = 0.2, adj = 0.05)

plot(MAT, gsms_norm[2, ], col = cols[2], pch = 7, xlab = '', ylab = '', 
     las=1, ylim = c(0,1))
m <- lm(gsms_norm[2,]~MAT)
p <- lmp(m)
legend('topright', border = NULL, bty = 'n',
       legend = paste('Rsq = ', signif(summary(m)$r.squared, 2), 
                      'p = ', signif(p, 2)), 
       text.col = 'black', cex = 1.2)
ys <- m$coefficients[1] + m$coefficients[2]*xs
points(xs, ys, type = 'l', col = cols[2], lty = 1)
mtext('b)', side = 3, line = 0.2, adj = 0.05)

plot(MAT, gsms_norm[3, ], col = cols[3], pch = 8, xlab = '', ylab = '', 
     las=1, ylim = c(0,1))
m <- lm(gsms_norm[3,]~MAT)
p <- lmp(m)
legend('topright', border = NULL, bty = 'n',
       legend = paste('Rsq = ', signif(summary(m)$r.squared, 2), 
                      'p = ', signif(p, 2)), 
       text.col = 'black', cex = 1.2)
ys <- m$coefficients[1] + m$coefficients[2]*xs
points(xs, ys, type = 'l', col = cols[3], lty = 1)
mtext('c)', side = 3, line = 0.2, adj = 0.05)

plot(MAT, gsms_norm[4, ], col = cols[4], pch = 4, xlab = '', ylab = '',
     las=1, ylim = c(0,1))
m <- lm(gsms_norm[4,]~MAT)
p <- lmp(m)
legend('topright', border = NULL, bty = 'n',
       legend = paste('Rsq = ', signif(summary(m)$r.squared, 2), 
                      'p = ', signif(p, 2)), 
       text.col = 'black', cex = 1.2)
ys <- m$coefficients[1] + m$coefficients[2]*xs
points(xs, ys, type = 'l', col = cols[4], lty = 1)
mtext('d)', side = 3, line = 0.2, adj = 0.05)

plot(MAT, gsms_norm[5, ], col = cols[5], pch = 5, xlab = '', ylab = '', 
     las=1, ylim = c(0,1))
m <- lm(gsms_norm[5,]~MAT)
p <- lmp(m)
legend('topright', border = NULL, bty = 'n',
       legend = paste('Rsq = ', signif(summary(m)$r.squared, 2), 
                      'p = ', signif(p, 2)), 
       text.col = 'black', cex = 1.2)
ys <- m$coefficients[1] + m$coefficients[2]*xs
points(xs, ys, type = 'l', col = cols[5], lty = 1)
mtext('e)', side = 3, line = 0.2, adj = 0.05)

plot(MAT, gsms_norm[6, ], col = cols[6], pch = 16, xlab = '', ylab = '', 
     las=1, ylim = c(0,1))
m <- lm(gsms_norm[6,]~MAT)
p <- lmp(m)
legend('topright', border = NULL, bty = 'n',
       legend = paste('Rsq = ', signif(summary(m)$r.squared, 2), 
                      'p = ', signif(p, 2)), 
       text.col = 'black', cex = 1.2)
ys <- m$coefficients[1] + m$coefficients[2]*xs
points(xs, ys, type = 'l', col = cols[6], lty = 1)
mtext('f)', side = 3, line = 0.2, adj = 0.05)

plot(MAT, gsms_norm[7, ], col = cols[7], pch = 17, xlab = '', ylab = '', 
     las=1, ylim = c(0,1))
m <- lm(gsms_norm[7,]~MAT)
p <- lmp(m)
legend('topright', border = NULL, bty = 'n',
       legend = paste('Rsq = ', signif(summary(m)$r.squared, 2), 
                      'p = ', signif(p, 2)), 
       text.col = 'black', cex = 1.2)
ys <- m$coefficients[1] + m$coefficients[2]*xs
points(xs, ys, type = 'l', col = cols[7], lty = 1)
mtext('g)', side = 3, line = 0.2, adj = 0.05)

plot(MAT, gsms_norm[8, ], col = cols[8], pch = 18, xlab = '', ylab = '',
     las=1, ylim = c(0,1))
m <- lm(gsms_norm[8,]~MAT)
p <- lmp(m)
legend('topright', border = NULL, bty = 'n',
       legend = paste('Rsq = ', signif(summary(m)$r.squared, 2), 
                      'p = ', signif(p, 2)), 
       text.col = 'black', cex = 1.2)
ys <- m$coefficients[1] + m$coefficients[2]*xs
points(xs, ys, type = 'l', col = cols[8], lty = 1)
mtext('h)', side = 3, line = 0.2, adj = 0.05)

mtext('Mean Annual Temperature', side = 1, line = 1, cex = 1.2, outer=TRUE)
mtext('Demographic composition', side = 2, line = 1, cex = 1.2, outer=TRUE)

dev.off()
#######################################################################
pdf('Paper_Figures_v3/SI_GSMSvMAP.pdf')
par(mfrow = c(3,3), mar = c(2,2,2,1), oma = c(3,3,0,0))

plot(MAP, gsms_norm[1,], xlab = '', ylab = '', las=1, col = cols[1], 
     pch = 15, ylim = c(0,1))
m <- lm(gsms_norm[1,]~MAP)
p <- lmp(m)
legend('topright', border = NULL, bty = 'n',
       legend = paste('Rsq = ', signif(summary(m)$r.squared, 2), 
                      'p = ', signif(p, 2)), 
       text.col = 'black', cex = 1.2)
xs <- seq(min(MAP), max(MAP), 1)
ys <- m$coefficients[1] + m$coefficients[2]*xs
points(xs, ys, type = 'l', col = cols[1], lty = 1)
mtext('a)', side = 3, line = 0.2, adj = 0.05)

plot(MAP, gsms_norm[2, ], col = cols[2], pch = 7, xlab = '', ylab = '', 
     las=1, ylim = c(0,1))
m <- lm(gsms_norm[2,]~MAP)
p <- lmp(m)
legend('topright', border = NULL, bty = 'n',
       legend = paste('Rsq = ', signif(summary(m)$r.squared, 2), 
                      'p = ', signif(p, 2)), 
       text.col = 'black', cex = 1.2)
ys <- m$coefficients[1] + m$coefficients[2]*xs
points(xs, ys, type = 'l', col = cols[2], lty = 1)
mtext('b)', side = 3, line = 0.2, adj = 0.05)

plot(MAP, gsms_norm[3, ], col = cols[3], pch = 8, xlab = '', ylab = '', 
     las=1, ylim = c(0,1))
m <- lm(gsms_norm[3,]~MAP)
p <- lmp(m)
legend('topright', border = NULL, bty = 'n',
       legend = paste('Rsq = ', signif(summary(m)$r.squared, 2), 
                      'p = ', signif(p, 2)), 
       text.col = 'black', cex = 1.2)
ys <- m$coefficients[1] + m$coefficients[2]*xs
points(xs, ys, type = 'l', col = cols[3], lty = 1)
mtext('c)', side = 3, line = 0.2, adj = 0.05)

plot(MAP, gsms_norm[4, ], col = cols[4], pch = 4, xlab = '', ylab = '',
     las=1, ylim = c(0,1))
m <- lm(gsms_norm[4,]~MAP)
p <- lmp(m)
legend('topright', border = NULL, bty = 'n',
       legend = paste('Rsq = ', signif(summary(m)$r.squared, 2), 
                      'p = ', signif(p, 2)), 
       text.col = 'black', cex = 1.2)
ys <- m$coefficients[1] + m$coefficients[2]*xs
points(xs, ys, type = 'l', col = cols[4], lty = 1)
mtext('d)', side = 3, line = 0.2, adj = 0.05)

plot(MAP, gsms_norm[5, ], col = cols[5], pch = 5, xlab = '', ylab = '', 
     las=1, ylim = c(0,1))
m <- lm(gsms_norm[5,]~MAP)
p <- lmp(m)
legend('topright', border = NULL, bty = 'n',
       legend = paste('Rsq = ', signif(summary(m)$r.squared, 2), 
                      'p = ', signif(p, 2)), 
       text.col = 'black', cex = 1.2)
ys <- m$coefficients[1] + m$coefficients[2]*xs
points(xs, ys, type = 'l', col = cols[5], lty = 1)
mtext('e)', side = 3, line = 0.2, adj = 0.05)

plot(MAP, gsms_norm[6, ], col = cols[6], pch = 16, xlab = '', ylab = '', 
     las=1, ylim = c(0,1))
m <- lm(gsms_norm[6,]~MAP)
p <- lmp(m)
legend('topright', border = NULL, bty = 'n',
       legend = paste('Rsq = ', signif(summary(m)$r.squared, 2), 
                      'p = ', signif(p, 2)), 
       text.col = 'black', cex = 1.2)
ys <- m$coefficients[1] + m$coefficients[2]*xs
points(xs, ys, type = 'l', col = cols[6], lty = 1)
mtext('f)', side = 3, line = 0.2, adj = 0.05)

plot(MAP, gsms_norm[7, ], col = cols[7], pch = 17, xlab = '', ylab = '', 
     las=1, ylim = c(0,1))
m <- lm(gsms_norm[7,]~MAP)
p <- lmp(m)
legend('topright', border = NULL, bty = 'n',
       legend = paste('Rsq = ', signif(summary(m)$r.squared, 2), 
                      'p = ', signif(p, 2)), 
       text.col = 'black', cex = 1.2)
ys <- m$coefficients[1] + m$coefficients[2]*xs
points(xs, ys, type = 'l', col = cols[7], lty = 1)
mtext('g)', side = 3, line = 0.2, adj = 0.05)

plot(MAP, gsms_norm[8, ], col = cols[8], pch = 18, xlab = '', ylab = '',
     las=1, ylim = c(0,1))
m <- lm(gsms_norm[8,]~MAP)
p <- lmp(m)
legend('topright', border = NULL, bty = 'n',
       legend = paste('Rsq = ', signif(summary(m)$r.squared, 2), 
                      'p = ', signif(p, 2)), 
       text.col = 'black', cex = 1.2)
ys <- m$coefficients[1] + m$coefficients[2]*xs
points(xs, ys, type = 'l', col = cols[8], lty = 1)
mtext('h)', side = 3, line = 0.2, adj = 0.05)

mtext('Mean Annual Precipitation (mm)', side = 1, line = 1, cex = 1.2, outer=TRUE)
mtext('Demographic composition', side = 2, line = 1, cex = 1.2, outer=TRUE)
dev.off()
###############################################
### DC and AGB ###
###############################################
#######################################################################
load(file = 'Analysis/Outputs/AGBs.RData')
load(file = 'Analysis/Outputs/taus.RData')

pdf('Paper_Figures_v3/SI_GSMSvAGB.pdf')
par(mfrow = c(3,3), mar = c(2,2,2,1), oma = c(3,3,0,0))

plot(gsms_norm[1,], AGBs, xlab = '', ylab = '', las=1, col = cols[1], 
     pch = 15)
m <- lm(AGBs~gsms_norm[1,])
p <- lmp(m)
legend('topright', border = NULL, bty = 'n',
       legend = paste('Rsq = ', signif(summary(m)$r.squared, 2), 
                      'p = ', signif(p, 2)), 
       text.col = 'black', cex = 1.2)
xs <- seq(min(gsms_norm[1, ]), max(gsms_norm[1, ]), length.out=50)
ys <- m$coefficients[1] + m$coefficients[2]*xs
points(xs, ys, type = 'l', col = cols[1], lty = 1)
mtext('a)', side = 3, line = 0.2, adj = 0.05)

plot(gsms_norm[2, ], AGBs, col = cols[2], pch = 7, xlab = '', ylab = '', 
     las=1)
m <- lm(AGBs~gsms_norm[2,])
p <- lmp(m)
legend('topright', border = NULL, bty = 'n',
       legend = paste('Rsq = ', signif(summary(m)$r.squared, 2), 
                      'p = ', signif(p, 2)), 
       text.col = 'black', cex = 1.2)
xs <- seq(min(gsms_norm[2, ]), max(gsms_norm[2, ]), length.out=50)
ys <- m$coefficients[1] + m$coefficients[2]*xs
points(xs, ys, type = 'l', col = cols[2], lty = 1)
mtext('b)', side = 3, line = 0.2, adj = 0.05)

plot(gsms_norm[3, ], AGBs, col = cols[3], pch = 8, xlab = '', ylab = '', 
     las=1)
m <- lm(AGBs~gsms_norm[3,])
p <- lmp(m)
legend('topright', border = NULL, bty = 'n',
       legend = paste('Rsq = ', signif(summary(m)$r.squared, 2), 
                      'p = ', signif(p, 2)), 
       text.col = 'black', cex = 1.2)
xs <- seq(min(gsms_norm[3, ]), max(gsms_norm[3, ]), length.out=50)
ys <- m$coefficients[1] + m$coefficients[2]*xs
points(xs, ys, type = 'l', col = cols[3], lty = 1)
mtext('c)', side = 3, line = 0.2, adj = 0.05)

plot(gsms_norm[4, ], AGBs, col = cols[4], pch = 4, xlab = '', ylab = '',
     las=1)
m <- lm(AGBs~ gsms_norm[4,])
p <- lmp(m)
legend('topright', border = NULL, bty = 'n',
       legend = paste('Rsq = ', signif(summary(m)$r.squared, 2), 
                      'p = ', signif(p, 2)), 
       text.col = 'black', cex = 1.2)
xs <- seq(min(gsms_norm[4, ]), max(gsms_norm[4, ]), length.out=50)
ys <- m$coefficients[1] + m$coefficients[2]*xs
points(xs, ys, type = 'l', col = cols[4], lty = 1)
mtext('d)', side = 3, line = 0.2, adj = 0.05)

plot(gsms_norm[5, ], AGBs, col = cols[5], pch = 5, xlab = '', ylab = '', 
     las=1)
m <- lm(AGBs~gsms_norm[5,])
p <- lmp(m)
legend('topright', border = NULL, bty = 'n',
       legend = paste('Rsq = ', signif(summary(m)$r.squared, 2), 
                      'p = ', signif(p, 2)), 
       text.col = 'black', cex = 1.2)
xs <- seq(min(gsms_norm[5, ]), max(gsms_norm[5, ]), length.out=50)
ys <- m$coefficients[1] + m$coefficients[2]*xs
points(xs, ys, type = 'l', col = cols[5], lty = 1)
mtext('e)', side = 3, line = 0.2, adj = 0.05)

plot(gsms_norm[6, ], AGBs, col = cols[6], pch = 16, xlab = '', ylab = '', 
     las=1)
m <- lm(AGBs~gsms_norm[6,])
p <- lmp(m)
legend('topright', border = NULL, bty = 'n',
       legend = paste('Rsq = ', signif(summary(m)$r.squared, 2), 
                      'p = ', signif(p, 2)), 
       text.col = 'black', cex = 1.2)
xs <- seq(min(gsms_norm[6, ]), max(gsms_norm[6, ]), length.out=50)
ys <- m$coefficients[1] + m$coefficients[2]*xs
points(xs, ys, type = 'l', col = cols[6],  lty = 1)
mtext('f)', side = 3, line = 0.2, adj = 0.05)

plot(gsms_norm[7, ], AGBs, col = cols[7], pch = 17, xlab = '', ylab = '', 
     las=1)
m <- lm(AGBs~gsms_norm[7,])
p <- lmp(m)
legend('topright', border = NULL, bty = 'n',
       legend = paste('Rsq = ', signif(summary(m)$r.squared, 2), 
                      'p = ', signif(p, 2)), 
       text.col = 'black', cex = 1.2)
xs <- seq(min(gsms_norm[7, ]), max(gsms_norm[7, ]), length.out=50)
ys <- m$coefficients[1] + m$coefficients[2]*xs
points(xs, ys, type = 'l', col = cols[7],  lty = 1)
mtext('g)', side = 3, line = 0.2, adj = 0.05)

plot(gsms_norm[8, ], AGBs, col = cols[8], pch = 18, xlab = '', ylab = '',
     las=1)
m <- lm(AGBs~gsms_norm[8,])
p <- lmp(m)
legend('topright', border = NULL, bty = 'n',
       legend = paste('Rsq = ', signif(summary(m)$r.squared, 2), 
                      'p = ', signif(p, 2)), 
       text.col = 'black', cex = 1.2)
xs <- seq(min(gsms_norm[8, ]), max(gsms_norm[8, ]), length.out=50)
ys <- m$coefficients[1] + m$coefficients[2]*xs
points(xs, ys, type = 'l', col = cols[8],  lty = 1)
mtext('h)', side = 3, line = 0.2, adj = 0.05)

mtext('Demographic Composition', side = 1, line = 1, cex = 1.2, outer=TRUE)
mtext('AGB (kg C m2)', side = 2, line = 1, cex = 1.2, outer=TRUE)

dev.off()
##############################################
### DC and tau ###
###############################################
#######################################################################

pdf('Paper_Figures_v3/SI_GSMSvtau.pdf')
par(mfrow = c(3,3), mar = c(2,2,2,1), oma = c(3,3,0,0), xpd = FALSE)

plot(gsms_norm[1,], tau, xlab = '', ylab = '', las=1, col = cols[1], 
     pch = 15)
why.ind <- which(sites == 'wytham')
m <- lm(tau[-why.ind]~gsms_norm[1,-why.ind])
p <- lmp(m)
legend('topright', border = NULL, bty = 'n',
       legend = paste('Rsq = ', signif(summary(m)$r.squared, 2), 
                      'p = ', signif(p, 2)), 
       text.col = 'black', cex = 1.2)
xs <- seq(min(gsms_norm[1, ]), max(gsms_norm[1, ]), length.out=50)
ys <- m$coefficients[1] + m$coefficients[2]*xs
points(xs, ys, type = 'l', col = cols[1],  lty = 1)
mtext('a)', side = 3, line = 0.2, adj = 0.05)

plot(gsms_norm[2, ], tau, col = cols[2], pch = 7, xlab = '', ylab = '', 
     las=1)
m <- lm(tau[-why.ind]~gsms_norm[2,-why.ind])
p <- lmp(m)
legend('topright', border = NULL, bty = 'n',
       legend = paste('Rsq = ', signif(summary(m)$r.squared, 2), 
                      'p = ', signif(p, 2)), 
       text.col = 'black', cex = 1.2)
xs <- seq(min(gsms_norm[2, ]), max(gsms_norm[2, ]), length.out=50)
ys <- m$coefficients[1] + m$coefficients[2]*xs
points(xs, ys, type = 'l', col = cols[2],  lty = 1)
mtext('b)', side = 3, line = 0.2, adj = 0.05)

plot(gsms_norm[3, ], tau, col = cols[3], pch = 8, xlab = '', ylab = '', 
     las=1)
m <- lm(tau[-why.ind]~gsms_norm[3,-why.ind])
p <- lmp(m)
legend('topright', border = NULL, bty = 'n',
       legend = paste('Rsq = ', signif(summary(m)$r.squared, 2), 
                      'p = ', signif(p, 2)), 
       text.col = 'black', cex = 1.2)
xs <- seq(min(gsms_norm[3, ]), max(gsms_norm[3, ]), length.out=50)
ys <- m$coefficients[1] + m$coefficients[2]*xs
points(xs, ys, type = 'l', col = cols[3],  lty = 1)
mtext('c)', side = 3, line = 0.2, adj = 0.05)

plot(gsms_norm[4, ], tau, col = cols[4], pch = 4, xlab = '', ylab = '',
     las=1)
m <- lm(tau[-why.ind]~ gsms_norm[4,-why.ind])
p <- lmp(m)
legend('topright', border = NULL, bty = 'n',
       legend = paste('Rsq = ', signif(summary(m)$r.squared, 2), 
                      'p = ', signif(p, 2)), 
       text.col = 'black', cex = 1.2)
xs <- seq(min(gsms_norm[4, ]), max(gsms_norm[4, ]), length.out=50)
ys <- m$coefficients[1] + m$coefficients[2]*xs
points(xs, ys, type = 'l', col = cols[4],  lty = 1)
mtext('d)', side = 3, line = 0.2, adj = 0.05)

plot(gsms_norm[5, ], tau, col = cols[5], pch = 5, xlab = '', ylab = '', 
     las=1)
m <- lm(tau[-why.ind]~gsms_norm[5,-why.ind])
p <- lmp(m)
legend('topright', border = NULL, bty = 'n',
       legend = paste('Rsq = ', signif(summary(m)$r.squared, 2), 
                      'p = ', signif(p, 2)), 
       text.col = 'black', cex = 1.2)
xs <- seq(min(gsms_norm[5, ]), max(gsms_norm[5, ]), length.out=50)
ys <- m$coefficients[1] + m$coefficients[2]*xs
points(xs, ys, type = 'l', col = cols[5],  lty = 1)
mtext('e)', side = 3, line = 0.2, adj = 0.05)

plot(gsms_norm[6, ], tau, col = cols[6], pch = 16, xlab = '', ylab = '', 
     las=1)
m <- lm(tau[-why.ind]~gsms_norm[6,-why.ind])
p <- lmp(m)
legend('topright', border = NULL, bty = 'n',
       legend = paste('Rsq = ', signif(summary(m)$r.squared, 2), 
                      'p = ', signif(p, 2)), 
       text.col = 'black', cex = 1.2)
xs <- seq(min(gsms_norm[6, ]), max(gsms_norm[6, ]), length.out=50)
ys <- m$coefficients[1] + m$coefficients[2]*xs
points(xs, ys, type = 'l', col = cols[6],  lty = 1)
mtext('f)', side = 3, line = 0.2, adj = 0.05)

plot(gsms_norm[7, ], tau, col = cols[7], pch = 17, xlab = '', ylab = '', 
     las=1)
m <- lm(tau[-why.ind]~gsms_norm[7,-why.ind])
p <- lmp(m)
legend('topright', border = NULL, bty = 'n',
       legend = paste('Rsq = ', signif(summary(m)$r.squared, 2), 
                      'p = ', signif(p, 2)), 
       text.col = 'black', cex = 1.2)
xs <- seq(min(gsms_norm[7, ]), max(gsms_norm[7, ]), length.out=50)
ys <- m$coefficients[1] + m$coefficients[2]*xs
points(xs, ys, type = 'l', col = cols[7],  lty = 1)
mtext('g)', side = 3, line = 0.2, adj = 0.05)

plot(gsms_norm[8, ], tau, col = cols[8], pch = 18, xlab = '', ylab = '',
     las=1)
m <- lm(tau[-why.ind]~gsms_norm[8,-why.ind])
p <- lmp(m)
legend('topright', border = NULL, bty = 'n',
       legend = paste('Rsq = ', signif(summary(m)$r.squared, 2), 
                      'p = ', signif(p, 2)), 
       text.col = 'black', cex = 1.2)
xs <- seq(min(gsms_norm[8, ]), max(gsms_norm[8, ]), length.out=50)
ys <- m$coefficients[1] + m$coefficients[2]*xs
points(xs, ys, type = 'l', col = cols[8],  lty = 1)
mtext('h)', side = 3, line = 0.2, adj = 0.05)

mtext('Demographic Composition', side = 1, line = 1, cex = 1.2, outer=TRUE)
mtext('Carbon residence time ', side = 2, line = 1, cex = 1.2, outer=TRUE)

dev.off()


# for each plot load data and plot size distribution of stems included
# and stems not included
par(mfrow = c(5,4), mar = c(5,4,1,1), oma = c(0,0,1,0))

pdf(file = 'Paper_Figures_v3/size_dists_stems_inncluded.pdf')
for(jj in 1:n.sites){
  site <- sites[jj]
  site.name <- site
  load( file = sprintf('Data_Zone/Output/%s_all_sp_agb.RData', site.name))
  
  hist(sp.agb$dbh[sp.agb$clust == 9],50, xlim = c(0,1400), 
       xlab = 'DBH (mm)', main = '')
  mtext(site.names[jj], side = 3, adj = 0.2, cex = 0.8)
  
}
dev.off()



#################################################
### Now show the AGB dynamics for these modes
site <- 'amacayacu'
load(file = 'Analysis/Outputs/agb_dfts.RData')
par(mfrow = c(1,1))

agb.dfts <- lapply(agb.dfts, function(x) apply(x, 2, function(x) x*0.5))
agb.quants <- lapply(agb.dfts, function(x) apply(x, 2, quantile, 
                                                 c(0.025, 0.25, 0.5, 0.75, 0.975)))

ymx <- max(unlist(lapply(agb.quants, max)))
xmx <- max(unlist(lapply(agb.quants, function(x){which(x[2, ] == 0)[1]})))
q <- agb.quants[[1]] 

plot(q[5, ], type = 'l', col = 'white', ylim = c(0,ymx), xlim = c(0,xmx), 
     ylab = '', xlab = '', cex.axis = 1, xaxt = 'n', yaxt = 'n')
axis(side = 1, at = NULL, padj =-2, cex.axis = 1, lwd.ticks = 0.6, tck=-0.01)
axis(side = 2, at = NULL, padj =2, cex.axis = 1, lwd.ticks = 0.6, tck=-0.01)
cols <- c(brewer.pal(7, 'Dark2')[c(6,1,2,3,4,5,7)], brewer.pal(3, 'Set1')[2])
cols.light <- make.transp(cols)
cols.v.light <- make.transp(cols, c(20,20,20,20,20,20,20,20))

for(cl in 1:length(agb.quants)){
  q <- agb.quants[[cl]] 
  polygon(x = c(seq(1500), rev(seq(1500))), 
          y = c(q[1, ], rev(q[5, ])), col = cols.v.light[cl], 
          border = cols.v.light[cl], cex.axis = 1)
  polygon(x = c(seq(1500), rev(seq(1500))), 
          y = c(q[2, ], rev(q[4, ])), col = cols.light[cl], 
          border = cols.light[cl])
  points(seq(1500), q[3, ], type = 'l', lwd = 0.6, col = cols[cl])
}
mtext('Time (years)', side = 1, line = 1.5, outer = FALSE, cex = 1)
mtext(paste('AGB of 1000 ', '\n', 'individual cohorts (Mg C)'),
      side = 2, line = 0, 
      outer = TRUE, cex = 1, adj = 0.15)    
mtext('AGB dynamics', side = 3, line = 0.1, cex = 1)
legend('topright', bty = 'n', col = cols, legend = seq(n.clust), lwd = 0.6, 
       cex = 1)


##################################################################
# Passage times and life expectancies - and AGB dynamics 
# for each mode

library(plotrix)

n.clust <- 8
# colours 

geog <- 'global'
load(file = sprintf('Analysis/Outputs/%s_pca_mat.RData', geog))
load(file = sprintf('Analysis/Outputs/PCA_%s.RData', geog))
param.mat <- pca.matrix
# split it up into clusters
cl.list <- split(param.mat, param.mat[ ,'clust'])
n.clust <- length(cl.list)

pts.cl <- les.cl <- vector('list', n.clust)
for(cl in 1:n.clust){
  load(file = sprintf('Analysis/Outputs/%d_pts_les_transitions.RData', cl))
  pts.cl[[cl]] <- pts
  les.cl[[cl]] <- les
}

mat.le <- do.call(cbind, les.cl)

save(les.cl, file = 'Analysis/Outputs/les_gsm.RData')
save(pts.cl, file = 'Analysis/Outputs/pts_gsm.RData')

# Here we do the Tukey Test
dfts <- rep(seq(n.clust), lapply(pts.cl, length))
long.format <- data.frame(cbind(dfts, unlist(pts.cl)))
colnames(long.format) <- c('clust', 'pt')
long.format[ ,'clust'] <- as.factor(long.format[ ,'clust'])

dfts <- rep(seq(n.clust), each = nrow(mat.le))
long.format.le <- data.frame(cbind(dfts, matrix(mat.le, ncol = 1)))
colnames(long.format.le) <- c('clust', 'le')
long.format.le[ ,'clust'] <- as.factor(long.format.le[ ,'clust'])

mod.pt <- aov(pt ~ clust, data = long.format)
mod.le <- aov(le ~ clust, data = long.format.le)

tukey.mod.pt <- TukeyHSD(mod.pt, 'clust')
tukey.mod.le <- TukeyHSD(mod.le, 'clust')
print(tukey.mod.pt)
print(tukey.mod.le)

tukey.pt.levels <- tukey.mod.pt[['clust']][,4]
tukey.pt.labs <- data.frame(multcompLetters(tukey.pt.levels)['Letters'])
tukey.pt.labs <- as.character(tukey.pt.labs[order(rownames(tukey.pt.labs)),])

tukey.le.levels <- tukey.mod.le[['clust']][,4]
tukey.le.labs <- data.frame(multcompLetters(tukey.le.levels)['Letters'])
tukey.le.labs <- as.character(tukey.le.labs[order(rownames(tukey.le.labs)),])

pdf(file = 'Paper_Figures_v3/LEPT.pdf', bg = 'white', 
    width = 8, height = 8)
par(mfrow = c(2,1), mar = c(2,3,0,0), oma = c(2,2,0,0.1))
#par(mfrow = c(2,1), mar = c(2,2,2,2), oma = c(2,2,0,0))
fs <- violin_plot(mat.le, 
                  col = cols.light,
                  cex.axis = 0.3, las = 1, 
                  ann = FALSE, cex = 0.6, axes = FALSE, 
                  xlab = '', ylab = '', main = '', 
                  show_mean = TRUE, mean_pch = 16, 
                  mean_pch_col = 'black', 
                  median_col = 'grey', border = NULL, 
                  ylim = c(0,725))
mtext(side = 3, line = -2.5, text = 'a)', cex = 1, adj = 0.01)
axis(side = 1, at = seq(n.clust), lwd.ticks = 0.6, labels = seq(n.clust), 
     cex.axis = 1,  tck=-0.01)
ymax <- max(mat.le)
axis(side = 2, at = NULL, lwd.ticks = 0.6, cex.axis = 1, tck=-0.01)
#text(seq(n.clust), 720, tukey.le.labs, cex=1)
mtext('Time (years)', side = 2, line = 0, outer = TRUE, cex = 1)
#mtext('DFT', side = 1, line = 0.25, outer = TRUE, cex = 1.2)
mtext('Life Expectancies', side = 3, line = -2, outer = FALSE, 
      cex = 1)

fs <- violin_plot(pts.cl[[1]], xlim = c(0,9), 
                  col = cols.light[1],
                  las = 1, 
                  ann = FALSE, axes = FALSE,
                  xlab = '', ylab = '', main = '', 
                  show_mean = TRUE, mean_pch = 16,
                  mean_pch_col = 'black', 
                  median_col = 'grey', ylim = c(0,210), 
                  border = NULL, xaxt = 'n', 
                  violin_width = 0.7, x_axis_labels = TRUE)
for(ii in 2:n.clust){
  violin_plot(pts.cl[[ii]], col = cols.light[ii], 
              at = ii, add = TRUE, ann = FALSE, show_mean = TRUE, 
              mean_pch = 16, mean_pch_col = 'black', median_col = 'grey', 
              border = NULL, violin_width = 0.7, xaxt = 'n')
}
#text(seq(n.clust), 205, tukey.pt.labs, cex = 1)
axis(side = 1, at = seq(n.clust), lwd.ticks = 0.6, labels = seq(n.clust), 
     cex.axis = 1, lwd = 0.6, tck = -0.01)
ymax2 = max(unlist(pts.cl))
axis(side = 2, at = NULL, lwd.ticks = 0.6, cex.axis = 1, tck=-0.01)
mtext('b)', cex = 1,   side = 3, line = -2, adj = 0.01)
mtext('Time to 100 mm DBH', side = 3, line = -2, outer = FALSE, 
      cex = 1)

legend('right', bty = 'n', col = cols, legend = seq(n.clust), lwd = 1.5, 
       cex = 1)
dev.off()


apply(mat.le, 2, summary)
do.call(rbind, lapply(pts.cl, summary))

