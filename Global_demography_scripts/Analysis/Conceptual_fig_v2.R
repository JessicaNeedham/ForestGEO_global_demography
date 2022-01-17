##### Niche complementary v mass ratio ######
#############################################
rm(list = ls())
library(stringi)
library(RColorBrewer)
library(factoextra)
library(FactoMineR)
library(ade4)
library(corrplot)
library(npreg)
library(geometry)
set.seed(1)

source('Analysis/Function_files/PCA_functions.R')
source('Analysis/Function_files/Ade4_functions.R')

pdf(file = 'Paper_Figures_v4/Figure_1.pdf', width = 10)

par(mfrow = c(2,3), oma = c(0,2,0,0), mar = c(3,2,3,3))
par(xpd = FALSE)

plot(c(1,2,3,4,5,6,7,8,9,  6,2,8,9),
     c(1,4,2,3,1,2,4,1,3,  8,6,9,7),
     cex = c(3,2,7,3,9,2,5,1,5,    1,3,8,1),
     xlab = '', ylab = '', xaxt = 'n', yaxt = 'n',
     xlim = c(0, 10), ylim = c(0,10))
abline(h = 5)

mtext('Resource 1', side =1, line = 1, cex = 0.9)
mtext('Resource 2', side = 2, line = 1, cex = 0.9)
mtext('(a)', side = 3, line = 0, cex = 0.9, adj = 0)
corners <- par('usr')
par(xpd = TRUE)
text(x = corners[2] + 0.5, y =corners[4]-2.4,
     'Low species richness', srt = 270)
text(x = corners[2] + 0.5, y =corners[4]-8.4,
     'High species richness', srt = 270)
par(xpd = FALSE)

cols <- c(brewer.pal(7, 'Dark2')[c(6,1,2,3,4,5,7)], brewer.pal(3, 'Set1')[2])
cols.light <- make.transp(cols)
cols.v.light <- make.transp(cols, c(30,20,40,30,30,50,70,80))

cols.points <- cols.v.light
cols.lines <- make.transp(cols, c(25,10,20,20,20,50,70,80))

xs.low <- c(-6, -1,-2,  2)
ys.low <- c(3,  5,4,  3)
low <- cbind(xs.low, ys.low)
clusts <- c(1,3,3,4)
low <- cbind(low, clusts)
conv <- chull(low[,c(1,2)])
conv <- c(conv, conv[1])
coords <- low[conv,c(1,2)]
midsx <- tapply(low[,1], low[,3], mean)
midsy <- tapply(low[,2], low[,3], mean)

plot(low[,1], low[,2], xlab = '', ylab = '', 
     xaxt = 'n', yaxt = 'n', col = cols[clusts+1], 
     xlim = c(-9,7), ylim = c(-6,6), pch = 16)
abline(h = 0)
arrows(x0 = low[2:3,1], x1 = midsx[2],
       y0 = low[2:3,2], y1 = midsy[2], length = 0, 
       col = cols[4])

polygon(x = coords[,1], y = coords[,2], 
        col = cols.v.light[1], 
        border = cols.light[1], lwd =2)
mtext('Demographic axis 1', side =1, line = 1, cex = 0.9)
mtext('Demographic axis 2', side = 2, line = 1, cex = 0.9)
mtext('(b)', side = 3, line = 0, cex = 0.9, adj = 0)
corners <- par('usr')
par(xpd = TRUE)
text(x = corners[2] + 0.5, y =corners[4]-3,
     'Low species richness', srt = 270)
text(x = corners[2] + 0.5, y =corners[4]-10,
     'High species richness', srt = 270)
par(xpd = FALSE)

xs.high <- c(-5.5,-4.7,-4,-3.8,
             -5.8,-4.1,-3.9,-3.5,
             1,2,3,4,
             3.8,4.2,4.8,5.5)
ys.high <- c(-0.4,-1,-0.8,-2, 
             -3,-4,-5,-3.4,
             -4,-2.8,-3.4,-3.5,
             -0.6,-1.2,-2.4,-0.8)
high <- cbind(xs.high, ys.high)
clusts <- c(rep(1,4), rep(2,4), rep(3,4), rep(4,4))
high <- cbind(high, clusts)
conv <- chull(high[,c(1,2)])
conv <- c(conv, conv[1])
coords <- high[conv,c(1,2)]
midsx <- tapply(high[,1], high[,3], mean)
midsy <- tapply(high[,2], high[,3], mean)

points(high[,1], high[,2],col = cols[clusts+1], pch = 16)
arrows(x0 = high[1:4,1], x1 = midsx[1],
       y0 = high[1:4,2], y1 = midsy[1], length = 0, 
       col = cols[2])
arrows(x0 = high[5:8,1], x1 = midsx[2],
       y0 = high[5:8,2], y1 = midsy[2], length = 0, 
       col = cols[3])
arrows(x0 = high[9:12,1], x1 = midsx[3],
       y0 = high[9:12,2], y1 = midsy[3], length = 0, 
       col = cols[4])
arrows(x0 = high[13:16,1], x1 = midsx[4],
       y0 = high[13:16,2], y1 = midsy[4], length = 0, 
       col = cols[5])

polygon(x = coords[,1], y = coords[,2], 
        col = cols.v.light[1], 
        border = cols.light[1], lwd =2)

plot(seq(10), seq(10), xlab = '', ylab='', 
     xaxt = 'n', yaxt = 'n', type = 'l')
mtext('Species Richness', side = 1, line =1, cex = 0.9)
mtext('DD DC AGB', side = 2, line = 1, cex = 0.9)
par(xpd = TRUE)
arrows(9.5,-0.15,10.5,-0.15, length = 0.1)
arrows(-0.2, 8.5, -0.2, 9.5, length = 0.1)
par(xpd = FALSE)
mtext('(c)', side = 3, line = 0, cex = 0.9, adj = 0)



###########################
### Mass Ratio Hypotheses
###########################

par(xpd = FALSE)
plot(c(3,5,7,8), c(7.5,8,7,9), 
     cex = runif(4, 1, 18), 
     xlab = '', ylab = '', xaxt = 'n', yaxt = 'n',
     xlim = c(0, 10), ylim = c(0,10))
mtext('Resource 1', side =1, line = 1, cex = 0.9)
mtext('Resource 2', side = 2, line = 1, cex = 0.9)
mtext('(d)', side = 3, line = 0, cex = 0.9, adj = 0)
abline(h = 5)

points(c(1,2,3,4,5,6,7,8,9,5),
       c(1.5,2,3,2.5,1.2,1.1,2.7,3,3.2,2.4), 
     cex = runif(10, 1, 18))

corners <- par('usr')
par(xpd = TRUE)
text(x = corners[2] + 0.5, y =corners[4]-3,
     'Low species richness', srt = 270)
text(x = corners[2] + 0.5, y =corners[4]-8.5,
     'High species richness', srt = 270)
par(xpd = FALSE)

xs.low <- c(-6, -1,-2,  2)
ys.low <- c(3,  5,4,  3)
low <- cbind(xs.low, ys.low)
clusts <- c(1,3,3,4)
low <- cbind(low, clusts)
conv <- chull(low[,c(1,2)])
conv <- c(conv, conv[1])
coords <- low[conv,c(1,2)]
midsx <- tapply(low[,1], low[,3], mean)
midsy <- tapply(low[,2], low[,3], mean)

plot(low[,1], low[,2], xlab = '', ylab = '', 
     xaxt = 'n', yaxt = 'n', col = cols[clusts+1], 
     xlim = c(-9,7), ylim = c(-6,6), pch = 16)
abline(h = 0)
arrows(x0 = low[2:3,1], x1 = midsx[2],
       y0 = low[2:3,2], y1 = midsy[2], length = 0, 
       col = cols[4])

polygon(x = coords[,1], y = coords[,2], 
        col = cols.v.light[1], 
        border = cols.light[1], lwd =2)
mtext('Demographic axis 1', side =1, line = 1, cex = 0.9)
mtext('Demographic axis 2', side = 2, line = 1, cex = 0.9)
mtext('(e)', side = 3, line = 0, cex = 0.9, adj = 0)
corners <- par('usr')
par(xpd = TRUE)
text(x = corners[2] + 0.5, y =corners[4]-3,
     'Low species richness', srt = 270)
text(x = corners[2] + 0.5, y =corners[4]-10,
     'High species richness', srt = 270)
par(xpd = FALSE)


xs.high <- c(-6,-4,-5.2,-3,
            2,1.5,1.8,1.4,
            -2,-1,0,-1.5)
ys.high <- c(-5, -4,-4.5, -3.8,
                 -5, -3.9, -4.2, -4.5,
                 -3,-3.6,-3.4,-3.1)
high <- cbind(xs.high, ys.high)
clusts <- c(rep(1,4), rep(4,4), rep(3,4))
high <- cbind(high, clusts)
conv <- chull(high[,c(1,2)])
conv <- c(conv, conv[1])
coords <- high[conv,c(1,2)]
midsx <- tapply(high[,1], high[,3], mean)
midsy <- tapply(high[,2], high[,3], mean)

points(high[,1], high[,2], col = cols[clusts +1], pch = 16)
arrows(x0 = high[1:4,1], x1 = midsx[1],
       y0 = high[1:4,2], y1 = midsy[1], length = 0, 
       col = cols[2])
arrows(x0 = high[5:8,1], x1 = midsx[3],
       y0 = high[5:8,2], y1 = midsy[3], length = 0, 
       col = cols[3])
arrows(x0 = high[9:12,1], x1 = midsx[2],
       y0 = high[9:12,2], y1 = midsy[2], length = 0, 
       col = cols[4])

polygon(x = coords[,1], y = coords[,2], 
        col = cols.v.light[1], 
        border = cols.light[1], lwd =2)

plot(seq(10), rep(5,10), xlab = '', ylab='', 
     xaxt = 'n', yaxt = 'n', type = 'l', ylim = c(0,10), xlim = c(0,10))
#points(seq(10), seq(10), type = 'l')
#points(seq(10), rev(seq(10)), type = 'l')
mtext('Species Richness', side = 1, line =1, cex = 0.9)
mtext('DD DC AGB', side = 2, line = 1, cex = 0.9)
par(xpd = TRUE)
arrows(9.5,-1.5,10.5,-1.5, length = 0.1)
arrows(-1, 8.5, -1, 9.5, length = 0.1)
par(xpd = FALSE)
mtext('(f)', side = 3, line = 0, cex = 0.9, adj = 0)

dev.off()

