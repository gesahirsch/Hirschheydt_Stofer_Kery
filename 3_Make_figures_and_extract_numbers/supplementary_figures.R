load('2_Fit_linear_models_to_SEs/estimates.RData')
setwd('3_Make_figures_and_extract_numbers/')

library(rasterImage)
library(grDevices)
library(scales)

simnames <- c('Case2x150_CovNull','Case2x300','Case4x150','CovOcc','CovDet','CovBoth')


#  Figure S1  ####

mapPalette <- colorRampPalette(c("black","red3","goldenrod2","grey88"))
zlim <- c(0,1)

tiff(paste('Suppl_Fig.1_Validity.tif'), width=3636, height=1528, units='px', compression='lzw')
lay <- layout(matrix(c(1:24), ncol=6, nrow=4, byrow=T), widths=c(1,3,3,3,3,3), heights=c(0.36,3,0.36,3), respect=T)

#  Part A - Case2x150 with all levels of S  ##

par(mar=c(0,0,0,0))
plot(c(0,1), c(0,1),
     type = 'n', frame=F, axes = F, xlab = '', ylab = '', main = '')
for(s in 1:length(S)){
  plot(c(0,1), c(0,1),
       type = 'n', frame=F, axes = F, xlab = '', ylab = '', main = '')
  text(x=0.5, y=0.6, labels=paste(S[s],'SV sites',sep=' '), adj=0.5, cex=4, font=2)
}

# add legend
par(mar=c(1,2,1,1))
plot(c(-0.8,4), c(-0.3,1.3),
     type = 'n', frame=F, axes = F, xlab = '', ylab = '', main = '')
legend_image <- as.raster(matrix(rev(mapPalette(150)), ncol=1))
rasterImage(legend_image, 0, -0.05, 1, 0.95)
rect(0, -0.05, 1, 0.95, col=NULL, border=TRUE)
for(y in seq(-0.05, 0.95, l=length(pretty(zlim)))){
  lines(x=c(1, 1.2), y=c(y,y), lwd=2)
}
text(x=1.5, y=seq(-0.05, 0.95, l=length(pretty(zlim))), labels=pretty(zlim), adj=0, cex=3.5)
text(x=-0.1, y=1.15, labels='Valid', pos=4, cex=4)
text(x=-0.1, y=1.05, labels='proportion', pos=4, cex=4)
text(x=-1, y=1.28, labels=expression(paste(bold('(A)'))), pos=4, cex=4.5)

# heatmaps
par(mar=c(14,5,2,3))
for(s in 1:length(S)){
  yy <- validity[1,'valid',s,,]/simrep
  image(x=pvals, y=pvals, z=yy,
        asp=1, las=1, col=mapPalette(150), axes=F,
        xlab="", ylab="",
        xlim=c(0.08,0.92), ylim=c(0.08,0.92), zlim=zlim)
  mtext(text=paste("min:",round(min(validity[1,'valid',s,,]/simrep),2),
                   " max:",max(validity[1,'valid',s,,]/simrep),
                   " mean:",round(mean(validity[1,'valid',s,,]/simrep),2)),
        side=3, line=1, cex=2.5)
  for(h in seq(0.1,0.9,by=0.1)){
    lines(x=c(0.08,0.92), y=c(h,h), col=alpha("black",0.5), lty=3, lwd=3)
    lines(x=c(h,h), y=c(0.08,0.92), col=alpha("black",0.5), lty=3, lwd=3)
  }
  axis(side=1, at=seq(0.1,0.9,by=0.2), tick=F, pos=0.05, cex.axis=3) # tcl=-0.5 (tick-length)
  axis(side=2, at=seq(0.1,0.9,by=0.2), tick=F, pos=0.08, cex.axis=3, las=1)
  mtext(text=expression(paste("Occupancy probability ",Psi)), side=1, line=8.5, cex=2.5)
  mtext(text=expression(paste("Detection probability  ", italic(p))), side=2, line=2, cex=2.5)
}

#  Part B - all other cases/scenarios  ##

par(mar=c(0,0,0,0))
plot(c(0,1), c(0,1),
     type = 'n', frame=F, axes = F, xlab = '', ylab = '', main = '')
for(c in 2:length(simnames)){
  plot(c(0,1), c(0,1),
       type = 'n', frame=F, axes = F, xlab = '', ylab = '', main = '')
  text(x=0.5, y=0.6, labels=simnames[c], adj=0.5, cex=4, font=2)
}

# add legend
par(mar=c(1,2,1,1))
plot(c(-0.8,4), c(-0.3,1.3),
     type = 'n', frame=F, axes = F, xlab = '', ylab = '', main = '')
legend_image <- as.raster(matrix(rev(mapPalette(150)), ncol=1))
rasterImage(legend_image, 0, -0.05, 1, 0.95)
rect(0, -0.05, 1, 0.95, col=NULL, border=TRUE)
for(y in seq(-0.05, 0.95, l=length(pretty(zlim)))){
  lines(x=c(1, 1.2), y=c(y,y), lwd=2)
}
text(x=1.5, y=seq(-0.05, 0.95, l=length(pretty(zlim))), labels=pretty(zlim), adj=0, cex=3.5)
text(x=-0.1, y=1.15, labels='Valid', pos=4, cex=4)
text(x=-0.1, y=1.05, labels='proportion', pos=4, cex=4)
text(x=-1, y=1.28, labels=expression(paste(bold('(B)'))), pos=4, cex=4.5)

# heatmaps
par(mar=c(14,5,2,3))
for(c in 2:length(simnames)){
  yy <- validity[c,'valid',1,,]/simrep
  image(x=pvals, y=pvals, z=yy,
        asp=1, las=1, col=mapPalette(150), axes=F,
        xlab="", ylab="",
        xlim=c(0.08,0.92), ylim=c(0.08,0.92), zlim=zlim)
  mtext(text=paste("min:",round(min(validity[c,'valid',1,,]/simrep),2),
                   " max:",max(validity[c,'valid',1,,]/simrep),
                   " mean:",round(mean(validity[c,'valid',1,,]/simrep),2)),
        side=3, line=1, cex=2.5)
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



#  Figure S2  ####

# viridis palette
mapPalette <- colorRampPalette(rev(c("#FDE725FF","#DCE319FF","#B8DE29FF","#73D055FF","#55C667FF","#20A387FF",
                                     "#287D8EFF","#39568CFF","#482677FF","#440154FF")))

# SE(psi)
tiff(paste('Suppl_Fig.2_Heatmaps_slopeSEcovariates.tif'), width=2000, height=1821, units='px', compression='lzw')
lay <- layout(matrix(c(1:12), ncol=3, nrow=4, byrow=T), widths=c(1.25,3,3), heights=c(0.3,3,0.3,3), respect=T)
layout.show(lay)
par(mar=c(0,0,0,0))
plot(c(0,1), c(0,1),
     type = 'n', frame=F, axes = F, xlab = '', ylab = '', main = '')
for(s in c(4,6)){
  plot(c(0,1), c(0,1),
       type = 'n', frame=F, axes = F, xlab = '', ylab = '', main = '')
  text(x=0.5, y=0.4, labels=ifelse(s==1,'CovNull',simnames[s]), adj=0.5, cex=4, font=2)
}

# add legend
par(mar=c(1,2,1,1))
# range(lmerSE[c(4,6),'SE slope(OccCov)','slope',,])
zlim <- c(-0.25, 0)
plot(c(-0.8,4), c(-0.3,1.3),
     type = 'n', frame=F, axes = F, xlab = '', ylab = '', main = '')
legend_image <- as.raster(matrix(rev(mapPalette(150)), ncol=1))

rasterImage(legend_image, 0, -0.1, 1, 0.9)
rect(0, -0.1, 1, 0.9, col=NULL, border=TRUE)
for(y in seq(-0.1, 0.9, l=length(pretty(zlim)))){
  lines(x=c(1, 1.2), y=c(y,y), lwd=2)
}
text(x=1.5, y=seq(-0.1, 0.9, l=length(pretty(zlim))), labels=pretty(zlim), adj=0, cex=3.5)
text(x=-0.1, y=1.13, labels=expression(paste('Regr. slope of')), pos=4, cex=4)
text(x=-0.1, y=1.01, labels=expression(paste(widehat(SE),' coef(OccCov)')), pos=4, cex=4)
text(x=-1, y=1.28, labels=expression(paste(bold('(A)'))), pos=4, cex=4.5)

# heatmaps
par(mar=c(14,5,2,3)) # , oma=c(0,0,4,0)
for(s in c(4,6)){
  simtitle <- simnames[s]
  yy <- lmerSE[s,'SE slope(OccCov)','slope',,]
  
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
par(mar=c(0,0,0,0))
plot(c(0,1), c(0,1),
     type = 'n', frame=F, axes = F, xlab = '', ylab = '', main = '')
for(s in 5:6){
  plot(c(0,1), c(0,1),
       type = 'n', frame=F, axes = F, xlab = '', ylab = '', main = '')
  text(x=0.5, y=0.4, labels=ifelse(s==1,'CovNull',simnames[s]), adj=0.5, cex=4, font=2)
}

# add legend
par(mar=c(1,2,1,1))
# range(lmerSE[5:6,'SE slope(DetCov)','slope',,])
zlim <- c(-0.28, 0)
plot(c(-0.8,4), c(-0.3,1.3),
     type = 'n', frame=F, axes = F, xlab = '', ylab = '', main = '')
legend_image <- as.raster(matrix(rev(mapPalette(150)), ncol=1))

rasterImage(legend_image, 0, -0.1, 1, 0.9)
rect(0, -0.1, 1, 0.9, col=NULL, border=TRUE)
for(y in seq(-0.1, 0.9, l=length(pretty(zlim)))){
  lines(x=c(1, 1.2), y=c(y,y), lwd=2)
}
text(x=1.5, y=seq(-0.1, 0.9, l=length(pretty(zlim))), labels=pretty(zlim), adj=0, cex=3.5)
text(x=-0.1, y=1.13, labels=expression(paste('Regr. slope of')), pos=4, cex=4)
text(x=-0.1, y=1.01, labels=expression(paste(widehat(SE),' coef(DetCov)')), pos=4, cex=4)
text(x=-1, y=1.28, labels=expression(paste(bold('(B)'))), pos=4, cex=4.5)

# heatmaps
par(mar=c(14,5,2,3))
for(s in 5:6){
  simtitle <- simnames[s]
  yy <- lmerSE[s,'SE slope(DetCov)','slope',,]
  
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



#  Figure S3  ####

farben <- c("black","#4dac26","goldenrod","#d01c8b")

p2 <- which(round(pvals,2)==round(0.2,2))
p5 <- which(round(pvals,2)==round(0.5,2))
p8 <- which(round(pvals,2)==round(0.8,2))
combis <- matrix(c(p5,p2,
                   p5,p5,
                   p5,p8), nrow=3, byrow=T)

tiff(paste('Suppl_Fig.3_Lineplots_SEcovariates.tif'), width=2400, height=2147, units='px', compression='lzw')
lay <- layout(matrix(c(1:12), ncol=4, nrow=3, byrow=T), widths=c(0.5,3,3,3), heights=c(0.5,4,4), respect=T)
layout.show(lay)

par(mar=c(0,0,0,0))
plot(c(0,1), c(0,1),
     type = 'n', frame=F, axes = F, xlab = '', ylab = '', main = '')
for(c in 1:3){
  ipsi <- combis[c,1]
  ip <- combis[c,2]
  
  plot(c(0,1), c(0,1),
       type = 'n', frame=F, axes = F, xlab = '', ylab = '', main = '')
  text(x=0.5, y=0.4, labels=paste('Detection probability =',pvals[ip], sep=' '),
       adj=0.5, cex=4, font=2)
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
  
  for(s in c(4,6)){
    if(s==4){
      plot(x=1:5, y=estimean[s,'SE slope(OccCov)',,ipsi,ip],
           type="o", col=farben[s-2], pch=c(16,18,15,17)[s-2], cex=7, lwd=5,
           xlim=c(0.8,5.5), ylim=c(0, max(estimean[c(4,6),'SE slope(OccCov)',,ipsi,ip])),
           xlab="", ylab="", axes=F, frame=F)
      axis(side=1, at=1:5, labels=c(S[1:5]), pos=0, cex.axis=3.5, padj=1.3)
      axis(side=2, at=seq(0,1.3,by=0.05), labels=seq(0,1.3,by=0.05), hadj=1.1, pos=0.8, las=1, cex.axis=3.5)
      if(c==1){
        mtext(expression(paste(widehat(SE),' coef(OccCov)')), side=2, line=8, cex=3)
      }
      if(c==1){
        legend(x=1.2, y=0.3, col=farben[c(2,4)], pch=c(18,17), lwd=3,
               bty='n', legend=simnames[c(4,6)], cex=5, title='Covariate settings')
      }
    }
    if(s==6){
      lines(x=1:5, y=estimean[s,'SE slope(OccCov)',,ipsi,ip], type='o',
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
  
  for(s in c(5,6)){
    if(s==5){
      plot(x=1:5, y=estimean[s,'SE slope(DetCov)',,ipsi,ip],
           type="o", col=farben[s-2], pch=c(16,18,15,17)[s-2], cex=7, lwd=5,
           xlim=c(0.8,5.5), ylim=c(0, max(estimean[c(5,6),'SE slope(DetCov)',,ipsi,ip])),
           xlab="", ylab="", axes=F, frame=F)
      axis(side=1, at=1:5, labels=c(S[1:5]), pos=0, cex.axis=3.5, padj=1.3)
      axis(side=2, at=seq(0,1.3,by=0.05), labels=seq(0,1.3,by=0.05), hadj=1.1, pos=0.8, las=1, cex.axis=3.5)
      if(c==1){
        mtext(expression(paste(widehat(SE),' coef(DetCov)')), side=2, line=8, cex=3)
      }
    }
    if(s==6){
      lines(x=1:5, y=estimean[s,'SE slope(DetCov)',,ipsi,ip], type='o',
            col=farben[s-2], pch=c(16,18,15,17)[s-2], lwd=5, cex=7)
    }
    if(c==2){
      mtext('Number of single-visit sites', side=1, line=9, cex=3)
    }
    if(c==1){
      legend(x=1.2, y=0.15, col=farben[3:4], pch=c(15,17), lwd=3,
             bty='n', legend=simnames[5:6], cex=5, title='Covariate settings')
    }
  }
}

dev.off()



#  Figure S4 and S5  ####

mapPalette <- colorRampPalette(rev(c("red","orange","white","skyblue","blue2")))

index <- list(c(1,2), # indicates which of the 4 variables/estimates exist for the given case
              c(1,2),
              c(1,2),
              c(1:3),
              c(1,2,4),
              c(1:4))

for(a in c(1,6)){ # 1:length(simnames) if you want to generate figures for all scenarios
  print(paste('Doing',simnames[a]))
  
  # dev.off()
  if(a==1){
    tiff(paste('Suppl_Fig.4_Bias_',simnames[a],'.tif',sep=''), width=2000, height=2023, units='px', compression='lzw') # 2500 + 2529
    lay <- layout(matrix(c(1,2,2,3:14), ncol=3, nrow=5, byrow=T), widths=c(1,3,3), heights=c(0.36,0.36,3,0.36,3), respect=T)
  }
  if(a==6){
    tiff(paste('Suppl_Fig.5_Bias_',simnames[a],'.tif',sep=''), width=2000, height=3943, units='px', compression='lzw') # 2500 + 2529
    lay <- layout(matrix(c(1,2,2,3:26), ncol=3, nrow=9, byrow=T), widths=c(1,3,3), heights=c(0.36,0.36,3,0.36,3,0.36,3,0.36,3), respect=T)
  }
  # If you want to generate one figure per scenario: comment out the 8 lines above
  #    and uncomment the 10 lines below:
  # if(a <= 3){
  #   tiff(paste('Suppl_Fig.X_bias_',a,'.tif',sep=''), width=2000, height=2023, units='px', compression='lzw') # 2500 + 2529
  #   lay <- layout(matrix(c(1,2,2,3:14), ncol=3, nrow=5, byrow=T), widths=c(1,3,3), heights=c(0.36,0.36,3,0.36,3), respect=T)
  # } else if (a %in% 4:5) {
  #   tiff(paste('Suppl_FigureSX_bias_',a,'.tif',sep=''), width=2000, height=2983, units='px', compression='lzw') # 2500 + 2529
  #   lay <- layout(matrix(c(1,2,2,3:20), ncol=3, nrow=7, byrow=T), widths=c(1,3,3), heights=c(0.36,0.36,3,0.36,3,0.36,3), respect=T)
  # } else {
  #   tiff(paste('Suppl_FigureSX_bias_',a,'.tif',sep=''), width=2000, height=3943, units='px', compression='lzw') # 2500 + 2529
  #   lay <- layout(matrix(c(1,2,2,3:26), ncol=3, nrow=9, byrow=T), widths=c(1,3,3), heights=c(0.36,0.36,3,0.36,3,0.36,3,0.36,3), respect=T)
  # }
  
  par(mar=c(0,0,0,0))
  plot(c(0,1), c(0,1),
       type = 'n', frame=F, axes = F, xlab = '', ylab = '', main = '')
  
  plot(c(0,1), c(0,1),
       type = 'n', frame=F, axes = F, xlab = '', ylab = '', main = '')
  if(a==1){
    text(x=0.5, y=0.6, labels='Case2x150', adj=0.5, cex=4, font=2)
  } else {
    text(x=0.5, y=0.6, labels=simnames[a], adj=0.5, cex=4, font=2)
  }
  
  for(v in dimnames(bias)[[2]][index[[a]]]){
    
    # range(biasmedian[,v,,,])
    if(v %in% c('slope(OccCov)','slope(DetCov)')){
      zlim <- c(-0.9,0.9)
      # zlim <- c(-1.1, 1.)
    }else{
      zlim <- c(-0.4,0.4)
      # zlim <- c(-1.2, 1.2)
    }
    
    par(mar=c(0,0,0,0))
    plot(c(0,1), c(0,1),
         type = 'n', frame=F, axes = F, xlab = '', ylab = '', main = '')
    for(s in c(1,5)){ # only 0 and 5000 SV sites
      plot(c(0,1), c(0,1),
           type = 'n', frame=F, axes = F, xlab = '', ylab = '', main = '')
      text(x=0.5, y=0.6, labels=paste(S[s],'SV sites',sep=' '), adj=0.5, cex=4, font=2)
    }
    
    # add legend
    par(mar=c(1,2,1,1))
    plot(c(-0.8,4), c(-0.3,1.3),
         type = 'n', frame=F, axes = F, xlab = '', ylab = '', main = '')
    legend_image <- as.raster(matrix(rev(mapPalette(150)), ncol=1))
    rasterImage(legend_image, 0, -0.05, 1, 0.95)
    rect(0, -0.05, 1, 0.95, col=NULL, border=TRUE)
    for(y in seq(-0.05, 0.95, l=length(pretty(zlim)))){
      lines(x=c(1, 1.2), y=c(y,y), lwd=2)
    }
    text(x=1.5, y=seq(-0.05, 0.95, l=length(pretty(zlim))), labels=pretty(zlim), adj=0, cex=3.5)
    text(x=-0.1, y=1.15, labels='bias of', pos=4, cex=4)
    if(v=='psi'){
      text(x=-0.1, y=1.05, labels=expression(paste(widehat(Psi))), pos=4, cex=4)
      text(x=-1, y=1.28, labels=expression(paste(bold('(A)'))), pos=4, cex=4.5)
    }
    if(v=='p'){
      text(x=-0.1, y=1.05, labels=expression(paste(widehat(italic(p)))), pos=4, cex=4)
      text(x=-1, y=1.28, labels=expression(paste(bold('(B)'))), pos=4, cex=4.5)
    }
    if(v=='slope(OccCov)'){
      text(x=-0.1, y=1.05, labels='coef(OccCov)', pos=4, cex=4)
      text(x=-1, y=1.28, labels=expression(paste(bold('(C)'))), pos=4, cex=4.5)
    }
    if(v=='slope(DetCov)'){
      text(x=-0.1, y=1.05, labels='coef(DetCov)', pos=4, cex=4)
      text(x=-1, y=1.28, labels=expression(paste(bold('(D)'))), pos=4, cex=4.5)
    }
    
    # heatmaps
    par(mar=c(14,5,2,3))
    for(s in c(1,5)){
      yy <- biasmedian[a,v,s,,]
      image(x=pvals, y=pvals, z=yy,
            asp=1, las=1, col=mapPalette(150), axes=F,
            xlab="", ylab="",
            xlim=c(0.08,0.92), ylim=c(0.08,0.92), zlim=zlim)
      mtext(text=paste("min:",round(min(biasmedian[a,v,s,,]),2),
                       " max:",round(max(biasmedian[a,v,s,,]),2),
                       " mean:",round(mean(biasmedian[a,v,s,,]),2)),
            side=3, line=1, cex=2.5)
      for(h in seq(0.1,0.9,by=0.1)){
        lines(x=c(0.08,0.92), y=c(h,h), col=alpha("black",0.5), lty=3, lwd=3)
        lines(x=c(h,h), y=c(0.08,0.92), col=alpha("black",0.5), lty=3, lwd=3)
      }
      axis(side=1, at=seq(0.1,0.9,by=0.2), tick=F, pos=0.05, cex.axis=3) # tcl=-0.5 (tick-length)
      axis(side=2, at=seq(0.1,0.9,by=0.2), tick=F, pos=0.08, cex.axis=3, las=1)
      mtext(text=expression(paste("Occupancy probability ",Psi)), side=1, line=8.5, cex=2.5)
      mtext(text=expression(paste("Detection probability  ", italic(p))), side=2, line=2, cex=2.5)
    }
  }
  dev.off()
}
