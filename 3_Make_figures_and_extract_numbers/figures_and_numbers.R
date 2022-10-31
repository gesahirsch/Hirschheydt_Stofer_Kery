load('2_Fit_linear_models_to_SEs/estimates.RData')
setwd('3_Make_figures_and_extract_numbers/')

library(rasterImage)
library(grDevices)
library(scales)

simnames <- c('Case2x150_CovNull','Case2x300','Case4x150','CovOcc','CovDet','CovBoth')


#  Figure 1  ####

# viridis palette
mapPalette <- colorRampPalette(rev(c("#FDE725FF","#DCE319FF","#B8DE29FF","#73D055FF","#55C667FF","#20A387FF",
                                     "#287D8EFF","#39568CFF","#482677FF","#440154FF")))

# Export graphics
# dev.off()
tiff(paste('Fig.1_Heatmaps_slopeSE_Simulation1.tif'), width=2286, height=1440, units='px', compression='lzw')
lay <- layout(matrix(c(1:12), ncol=4, nrow=3, byrow=T), widths=c(1,3,3,3), heights=c(0.3,3,3), respect=T)
# layout.show(lay)
par(mar=c(0,0,0,0))
plot(c(0,1), c(0,1),
     type = 'n', frame=F, axes = F, xlab = '', ylab = '', main = '')
for(s in 1:3){
  plot(c(0,1), c(0,1),
       type = 'n', frame=F, axes = F, xlab = '', ylab = '', main = '')
  text(x=0.5, y=0.4, labels=ifelse(s==1,'Case2x150',simnames[s]), adj=0.5, cex=4, font=2)
}


# SE(psi)
# add legend
par(mar=c(1,2,1,1))
zlim <- range(pretty(range(lmerSE[1:3,'SE psi logit','slope',,])))
plot(c(-0.8,4), c(-0.3,1.3),
     type = 'n', frame=F, axes = F, xlab = '', ylab = '', main = '')
legend_image <- as.raster(matrix(rev(mapPalette(150)), ncol=1))
rasterImage(legend_image, 0, -0.05, 1, 0.95)
rect(0, -0.05, 1, 0.95, col=NULL, border=TRUE)
for(y in seq(-0.05, 0.95, l=length(pretty(zlim)))){
  lines(x=c(1, 1.2), y=c(y,y), lwd=2)
}
text(x=1.5, y=seq(-0.05, 0.95, l=length(pretty(zlim))), labels=pretty(zlim), adj=0, cex=3.5)
text(x=-0.1, y=1.15, labels=expression(paste('Slope ',widehat(gamma[1]),' of')), pos=4, cex=4)
text(x=-0.1, y=1.05, labels=expression(paste(widehat(SE),' of ',widehat(Psi))), pos=4, cex=4)
text(x=-1, y=1.28, labels=expression(paste(bold('(A)'))), pos=4, cex=4.5)

# heatmaps
par(mar=c(14,5,2,3))
for(s in 1:3){
  simtitle <- simnames[s]
  yy <- lmerSE[s,'SE psi logit','slope',,]
  
  image(x=pvals, y=pvals, z=yy,
        asp=1, las=1, col=mapPalette(150), axes=F,
        xlab="", ylab="",
        xlim=c(0.08,0.92), ylim=c(0.08,0.92), zlim=zlim)
  for(h in seq(0.1,0.9,by=0.1)){
    lines(x=c(0.08,0.92), y=c(h,h), col=alpha("black",0.5), lty=3, lwd=3)
    lines(x=c(h,h), y=c(0.08,0.92), col=alpha("black",0.5), lty=3, lwd=3)
  }
  axis(side=1, at=seq(0.1,0.9,by=0.2), tick=F, pos=0.05, cex.axis=3) # tcl=-0.5 (tick-length)
  axis(side=2, at=seq(0.1,0.9,by=0.2), tick=F, pos=0.08, cex.axis=3, las=1)
  mtext(text=expression(paste("Occupancy probability ",Psi)), side=1, line=8.5, cex=2.5)
  mtext(text=expression(paste("Detection probability  ", italic(p))), side=2, line=2, cex=2.5)
}


# SE(p)
# add legend
par(mar=c(1,2,1,1))
zlim <- range(pretty(range(lmerSE[1:3,'SE p logit','slope',,])))
plot(c(-0.8,4), c(-0.3,1.3),
     type = 'n', frame=F, axes = F, xlab = '', ylab = '', main = '')
legend_image <- as.raster(matrix(rev(mapPalette(150)), ncol=1))
rasterImage(legend_image, 0, -0.05, 1, 0.95)
rect(0, -0.05, 1, 0.95, col=NULL, border=TRUE)
for(y in seq(-0.05, 0.95, l=length(pretty(zlim)))){
  lines(x=c(1, 1.2), y=c(y,y), lwd=2)
}
text(x=1.5, y=seq(-0.05, 0.95, l=length(pretty(zlim))), labels=pretty(zlim), adj=0, cex=3.5)
text(x=-0.1, y=1.15, labels=expression(paste('Slope ',widehat(gamma[1]),' of')), pos=4, cex=4)
text(x=-0.1, y=1.05, labels=expression(paste(widehat(SE),' of ',widehat(italic(p)))), pos=4, cex=4)
text(x=-1, y=1.28, labels=expression(paste(bold('(B)'))), pos=4, cex=4.5)

# heatmaps
par(mar=c(14,5,2,3))
for(s in 1:3){
  simtitle <- simnames[s]
  yy <- lmerSE[s,'SE p logit','slope',,]
  
  image(x=pvals, y=pvals, z=yy, # xlim=c(0,1), ylim=c(0,1), zlim=zlim,
        asp=1, las=1, col=mapPalette(150), axes=F,
        xlab="", ylab="",
        xlim=c(0.08,0.92), ylim=c(0.08,0.92), zlim=zlim)
  for(h in seq(0.1,0.9,by=0.1)){
    lines(x=c(0.08,0.92), y=c(h,h), col=alpha("black",0.5), lty=3, lwd=3)
    lines(x=c(h,h), y=c(0.08,0.92), col=alpha("black",0.5), lty=3, lwd=3)
  }
  axis(side=1, at=seq(0.1,0.9,by=0.2), tick=F, pos=0.05, cex.axis=3) # tcl=-0.5 (tick-length)
  axis(side=2, at=seq(0.1,0.9,by=0.2), tick=F, pos=0.08, cex.axis=3, las=1)
  mtext(text=expression(paste("Occupancy probability ",Psi)), side=1, line=8.5, cex=2.5)
  mtext(text=expression(paste("Detection probability  ", italic(p))), side=2, line=2, cex=2.5)
}

dev.off()



#  Figure 2  ####

farben <- c("black","orange1","blue")

p2 <- which(round(pvals,2)==round(0.2,2))
p5 <- which(round(pvals,2)==round(0.5,2))
p8 <- which(round(pvals,2)==round(0.8,2))
combis <- matrix(c(p5,p2,
                   p5,p5,
                   p5,p8), nrow=3, byrow=T)

# dev.off()
tiff(paste('Fig.2_Lineplots_SE_Simulation1.tif'), width=2400, height=2100, units='px', compression='lzw')
lay <- layout(matrix(c(1:12), ncol=4, nrow=3, byrow=T), widths=c(0.5,3,3,3), heights=c(0.5,4,4), respect=T)

par(mar=c(0,0,0,0))
plot(c(0,1), c(0,1),
     type = 'n', frame=F, axes = F, xlab = '', ylab = '', main = '')
for(c in 1:3){
  ipsi <- combis[c,1]
  ip <- combis[c,2]
  
  plot(c(0,1), c(0,1),
       type = 'n', frame=F, axes = F, xlab = '', ylab = '', main = '')
  text(x=0.5, y=0.4, labels=paste('Detection probability =',pvals[ip], sep=' '),
       adj=0.5, cex=4, font=2) # cex=4
}

# add legend
par(mar=c(0,3,0,3))
plot(c(-10,10), c(-10,10),
     type = 'n', frame=F, axes = F, xlab = '', ylab = '', main = '')
text(x=-7, y=10, labels=expression(paste(bold('(A)'))), pos=4, cex=4.5)

# SE psi
par(mar=c(12,6,2,1))

for(c in 1:3){
  ipsi <- combis[c,1]
  ip <- combis[c,2]
  
  for(s in 1:3){
    if(s==1){
      plot(x=1:5, y=estimean[s,'SE psi logit',,ipsi,ip],
           type="o", col=farben[s], pch=c(16,15,17)[s], cex=7, lwd=5,
           xlim=c(0.8,5.5), ylim=c(0, max(estimean[1:3,'SE psi logit',,ipsi,ip])),
           xlab="", ylab="", axes=F, frame=F)
      axis(side=1, at=1:5, labels=c(S[1:5]), pos=0, cex.axis=3.5, padj=1.3)
      axis(side=2, at=seq(0,1.3,by=0.05), labels=seq(0,1.3,by=0.05), hadj=1.1, pos=0.8, las=1, cex.axis=3.5)
      if(c==1){
        mtext(expression(paste('Estimated standard error of  ',widehat(Psi))), side=2, line=8, cex=3)
      }
    }
    if(s>1){
      lines(x=1:5, y=estimean[s,'SE psi logit',,ipsi,ip], type='o',
            col=farben[s], pch=c(16,15,17)[s], lwd=5, cex=7)
    }
  }
}

par(mar=c(0,0,0,3))
plot(c(-10,10), c(-10,10),
     type = 'n', frame=F, axes = F, xlab = '', ylab = '', main = '')
text(x=-1, y=10, labels=expression(paste(bold('(B)'))), pos=4, cex=4.5)


# SE p
par(mar=c(12,6,2,1))

for(c in 1:3){
  ipsi <- combis[c,1]
  ip <- combis[c,2]
  
  for(s in 1:3){
    if(s==1){
      plot(x=1:5, y=estimean[s,'SE p logit',,ipsi,ip],
           type="o", col=farben[s], pch=c(16,15,17)[s], cex=7, lwd=5,
           xlim=c(0.8,5.5), ylim=c(0, max(estimean[1:3,'SE p logit',,ipsi,ip])),
           xlab="", ylab="", axes=F, frame=F)
      axis(side=1, at=1:5, labels=c(S[1:5]), pos=0, cex.axis=3.5, padj=1.3)
      axis(side=2, at=seq(0,1.3,by=0.05), labels=seq(0,1.3,by=0.05), hadj=1.1, pos=0.8, las=1, cex.axis=3.5)
      if(c==1){
        mtext(expression(paste('Estimated standard error of  ',widehat(italic(p)))), side=2, line=8, cex=3)
      }
    }
    if(s>1){
      lines(x=1:5, y=estimean[s,'SE p logit',,ipsi,ip], type='o',
            col=farben[s], pch=c(16,15,17)[s], lwd=5, cex=7)
    }
    if(c==2){
      mtext('Number of single-visit sites', side=1, line=9, cex=3)
    }
    if(c==3){
      legend(x=2.4, y=0.09, col=farben[1:3], pch=c(16,15,17), lwd=3,
             bty='n', legend=c('Case2x150',simnames[2:3]), cex=5, title='Repeated sampling')
      
    }
  }
}

dev.off()



#  Figure 3  ####

# viridis palette
mapPalette <- colorRampPalette(rev(c("#FDE725FF","#DCE319FF","#B8DE29FF","#73D055FF","#55C667FF","#20A387FF",
                                     "#287D8EFF","#39568CFF","#482677FF","#440154FF")))

# Export graphics
# dev.off()
tiff(paste('Fig.3_Heatmaps_slopeSE_Simulation2.tif'), width=2500, height=1200, units='px', compression='lzw')
lay <- layout(matrix(c(1:15), ncol=5, nrow=3, byrow=T), widths=c(1.2,3,3,3,3), heights=c(0.3,3,3), respect=T)

par(mar=c(0,0,0,0))
plot(c(0,1), c(0,1),
     type = 'n', frame=F, axes = F, xlab = '', ylab = '', main = '')
for(s in c(1,4:6)){
  plot(c(0,1), c(0,1),
       type = 'n', frame=F, axes = F, xlab = '', ylab = '', main = '')
  text(x=0.5, y=0.4, labels=ifelse(s==1,'CovNull',simnames[s]), adj=0.5, cex=4, font=2)
}


# SE(psi)
# add legend
par(mar=c(1,2,1,1))
# range(lmerSE[c(1,4:6),'SE psi logit','slope',,])
zlim <- range(pretty(c(-0.13, 0)))
plot(c(-0.8,4), c(-0.3,1.3),
     type = 'n', frame=F, axes = F, xlab = '', ylab = '', main = '')
legend_image <- as.raster(matrix(rev(mapPalette(150)), ncol=1))
rasterImage(legend_image, 0, -0.1, 1, 0.9)
rect(0, -0.1, 1, 0.9, col=NULL, border=TRUE)
for(y in seq(-0.1, 0.9, l=length(pretty(zlim)))){
  lines(x=c(1, 1.2), y=c(y,y), lwd=2)
}
text(x=1.5, y=seq(-0.1, 0.9, l=length(pretty(zlim))), labels=pretty(zlim), adj=0, cex=3.5)
text(x=-0.1, y=1.13, labels=expression(paste('Slope ',widehat(gamma[1]),' of')), pos=4, cex=4)
text(x=-0.1, y=1.01, labels=expression(paste(widehat(SE),' of ',widehat(Psi))), pos=4, cex=4)
text(x=-1, y=1.28, labels=expression(paste(bold('(A)'))), pos=4, cex=4.5)

# heatmaps
par(mar=c(14,5,2,3))
zlim <- range(pretty(range(lmerSE[c(1,4:6),'SE psi logit','slope',,],na.rm=T)))
for(s in c(1,4:6)){
  simtitle <- simnames[s]
  yy <- lmerSE[s,'SE psi logit','slope',,]
  
  image(x=pvals, y=pvals, z=yy, # xlim=c(0,1), ylim=c(0,1), zlim=zlim,
        asp=1, las=1, col=mapPalette(150), axes=F,
        xlab="", ylab="",
        xlim=c(0.08,0.92), ylim=c(0.08,0.92), zlim=zlim)
  for(h in seq(0.1,0.9,by=0.1)){
    lines(x=c(0.08,0.92), y=c(h,h), col=alpha("black",0.5), lty=3, lwd=3)
    lines(x=c(h,h), y=c(0.08,0.92), col=alpha("black",0.5), lty=3, lwd=3)
  }
  axis(side=1, at=seq(0.1,0.9,by=0.2), tick=F, pos=0.05, cex.axis=3) # tcl=-0.5 (tick-length)
  axis(side=2, at=seq(0.1,0.9,by=0.2), tick=F, pos=0.08, cex.axis=3, las=1)
  mtext(text=expression(paste("Occupancy probability ",Psi)), side=1, line=8.5, cex=2.5)
  mtext(text=expression(paste("Detection probability  ", italic(p))), side=2, line=2, cex=2.5)
}


# SE(p)
# add legend
par(mar=c(1,2,1,1))
zlim <- range(pretty(range(lmerSE[c(1,4:6),'SE p logit','slope',,],na.rm=T)))
plot(c(-0.8,4), c(-0.3,1.3),
     type = 'n', frame=F, axes = F, xlab = '', ylab = '', main = '')
legend_image <- as.raster(matrix(rev(mapPalette(150)), ncol=1))
rasterImage(legend_image, 0, -0.1, 1, 0.9)
rect(0, -0.1, 1, 0.9, col=NULL, border=TRUE)
for(y in seq(-0.1, 0.9, l=length(pretty(zlim)))){
  lines(x=c(1, 1.2), y=c(y,y), lwd=2)
}
text(x=1.5, y=seq(-0.1, 0.9, l=length(pretty(zlim))), labels=pretty(zlim), adj=0, cex=3.5)
text(x=-0.1, y=1.13, labels=expression(paste('Slope ',widehat(gamma[1]),' of')), pos=4, cex=4)
text(x=-0.1, y=1.01, labels=expression(paste(widehat(SE),' of ',widehat(italic(p)))), pos=4, cex=4)
text(x=-1, y=1.28, labels=expression(paste(bold('(B)'))), pos=4, cex=4.5)

# heatmaps
par(mar=c(14,5,2,3))
for(s in c(1,4:6)){
  simtitle <- simnames[s]
  yy <- lmerSE[s,'SE p logit','slope',,]
  
  image(x=pvals, y=pvals, z=yy, # xlim=c(0,1), ylim=c(0,1), zlim=zlim,
        asp=1, las=1, col=mapPalette(150), axes=F,
        xlab="", ylab="",
        xlim=c(0.08,0.92), ylim=c(0.08,0.92), zlim=zlim)
  for(h in seq(0.1,0.9,by=0.1)){
    lines(x=c(0.08,0.92), y=c(h,h), col=alpha("black",0.5), lty=3, lwd=3)
    lines(x=c(h,h), y=c(0.08,0.92), col=alpha("black",0.5), lty=3, lwd=3)
  }
  axis(side=1, at=seq(0.1,0.9,by=0.2), tick=F, pos=0.05, cex.axis=3) # tcl=-0.5 (tick-length)
  axis(side=2, at=seq(0.1,0.9,by=0.2), tick=F, pos=0.08, cex.axis=3, las=1)
  mtext(text=expression(paste("Occupancy probability ",Psi)), side=1, line=8.5, cex=2.5)
  mtext(text=expression(paste("Detection probability  ", italic(p))), side=2, line=2, cex=2.5)
}

dev.off()



#  Figure 4  ####

farben <- c("black","#4dac26","goldenrod","#d01c8b")

p2 <- which(round(pvals,2)==round(0.2,2))
p5 <- which(round(pvals,2)==round(0.5,2))
p8 <- which(round(pvals,2)==round(0.8,2))
combis <- matrix(c(p5,p2,
                   p5,p5,
                   p5,p8), nrow=3, byrow=T)

# dev.off()
tiff(paste('Fig.4_Lineplots_SE_Simulation2.tif'), width=2400, height=2100, units='px', compression='lzw')
lay <- layout(matrix(c(1:12), ncol=4, nrow=3, byrow=T), widths=c(0.5,3,3,3), heights=c(0.5,4,4), respect=T)

par(mar=c(0,0,0,0))
plot(c(0,1), c(0,1),
     type = 'n', frame=F, axes = F, xlab = '', ylab = '', main = '')
for(c in 1:3){
  ipsi <- combis[c,1]
  ip <- combis[c,2]
  
  plot(c(0,1), c(0,1),
       type = 'n', frame=F, axes = F, xlab = '', ylab = '', main = '')
  text(x=0.5, y=0.4, labels=paste('Detection probability =',pvals[ip], sep=' '),
       adj=0.5, cex=4, font=2) # cex=4
}

# add legend
par(mar=c(0,3,0,3))
plot(c(-10,10), c(-10,10),
     type = 'n', frame=F, axes = F, xlab = '', ylab = '', main = '')
text(x=-7, y=10, labels=expression(paste(bold('(A)'))), pos=4, cex=4.5)

# SE psi
par(mar=c(12,6,2,1))

for(c in 1:3){
  ipsi <- combis[c,1]
  ip <- combis[c,2]
  
  for(s in c(1,4:6)){
    if(s==1){
      plot(x=1:5, y=estimean[s,'SE psi logit',,ipsi,ip],
           type="o", col=farben[s], pch=c(16,18,15,17)[s], cex=7, lwd=5,
           xlim=c(0.8,5.5), ylim=c(0, max(estimean[c(1,4:6),'SE psi logit',,ipsi,ip])),
           xlab="", ylab="", axes=F, frame=F)
      axis(side=1, at=1:5, labels=c(S[1:5]), pos=0, cex.axis=3.5, padj=1.3)
      axis(side=2, at=seq(0,1.3,by=0.05), labels=seq(0,1.3,by=0.05), hadj=1.1, pos=0.8, las=1, cex.axis=3.5)
      if(c==1){
        mtext(expression(paste('Estimated standard error of  ',widehat(Psi))), side=2, line=8, cex=3)
      }
    }
    if(s>1){
      lines(x=1:5, y=estimean[s,'SE psi logit',,ipsi,ip], type='o',
            col=farben[s-2], pch=c(16,18,15,17)[s-2], lwd=5, cex=7)
    }
  }
}

par(mar=c(0,0,0,3))
plot(c(-10,10), c(-10,10),
     type = 'n', frame=F, axes = F, xlab = '', ylab = '', main = '')
text(x=-1, y=10, labels=expression(paste(bold('(B)'))), pos=4, cex=4.5)


# SE p
par(mar=c(12,6,2,1))

for(c in 1:3){
  ipsi <- combis[c,1]
  ip <- combis[c,2]
  
  for(s in c(1,4:6)){
    if(s==1){
      plot(x=1:5, y=estimean[s,'SE p logit',,ipsi,ip],
           type="o", col=farben[s], pch=c(16,18,15,17)[s], cex=7, lwd=5,
           xlim=c(0.8,5.5), ylim=c(0, max(estimean[c(1,4:6),'SE p logit',,ipsi,ip])),
           xlab="", ylab="", axes=F, frame=F)
      axis(side=1, at=1:5, labels=c(S[1:5]), pos=0, cex.axis=3.5, padj=1.3)
      axis(side=2, at=seq(0,1.3,by=0.05), labels=seq(0,1.3,by=0.05), hadj=1.1, pos=0.8, las=1, cex.axis=3.5)
      if(c==1){
        mtext(expression(paste('Estimated standard error of  ',widehat(italic(p)))), side=2, line=8, cex=3)
      }
    }
    if(s>1){
      lines(x=1:5, y=estimean[s,'SE p logit',,ipsi,ip], type='o',
            col=farben[s-2], pch=c(16,18,15,17)[s-2], lwd=5, cex=7)
    }
    if(c==2){
      mtext('Number of single-visit sites', side=1, line=9, cex=3)
    }
    if(c==3){
      legend(x=2.4, y=0.11, col=farben[1:4], pch=c(16,18,15,17), lwd=3,
             bty='n', legend=c('CovNull',simnames[4:6]), cex=5, title='Covariate settings')
    }
  }
}

dev.off()



#  Get numbers for manuscript  ####

str(lmerSE)
str(estimean)
dimnames(estimean)

# Maximum proportion of numerical failures = non-valid estimates
str(validity)
(1000 - min(validity[,'valid',,,])) / 1000 # 89 %

# Mean and range of proportion of valid estimates for Simulation 1
mean(validity[1:3,'valid',,,])  # 907
range(validity[1:3,'valid',,,]) # 111, 1000

# Mean and range of proportion of valid estimates for Simulation 2
mean(validity[c(1,4:6),'valid',,,])  # 905
range(validity[c(1,4:6),'valid',,,]) # 111, 1000


# SE of psi

# Case2x150, 0 SV sites, psi=0.5, p=0.8
estimean[1,'SE psi logit','SV_0','psi_0.5','p_0.8']   # 0.175

# Case2x150, 500 SV sites, psi=0.5, p=0.8
estimean[1,'SE psi logit','SV_500','psi_0.5','p_0.8'] # 0.121

# % reduction from the former to the latter
(estimean[1,'SE psi logit','SV_0','psi_0.5','p_0.8'] -
    estimean[1,'SE psi logit','SV_500','psi_0.5','p_0.8']) /
  estimean[1,'SE psi logit','SV_0','psi_0.5','p_0.8'] * 100
(0.175-0.121)/0.175 * 100 # rounded


# SE of p

# Case2x150, 0 SV sites, psi=0.5, p=0.8
estimean[1,'SE p logit','SV_0','psi_0.5','p_0.8']    # 0.253

# Case2x150, 500 SV sites, psi=0.5, p=0.8
estimean[1,'SE p logit','SV_500','psi_0.5','p_0.8']  # 0.240

# % reduction from the former to the latter
(estimean[1,'SE p logit','SV_0','psi_0.5','p_0.8'] -
    estimean[1,'SE p logit','SV_500','psi_0.5','p_0.8']) /
  estimean[1,'SE p logit','SV_0','psi_0.5','p_0.8'] * 100
(0.253-0.240)/0.253 * 100 # rounded


# SE of p when SV have the greatest impact

# find the maximum proportional contribution of SV sites on SE of p
max((estimean[1,'SE p logit','SV_0',,] -
       estimean[1,'SE p logit','SV_500',,]) /
      estimean[1,'SE p logit','SV_0',,] * 100)
