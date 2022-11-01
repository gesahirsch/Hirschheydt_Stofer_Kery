simulated_data <- paste(getwd(), '/1_Simulate_data_and_fit_occupancy_models/', sep='')
setwd('2_Fit_linear_models_to_SEs/')

library(scales) # for alpha()
library(lme4) # for lmer()
library(tidyr) # for gather()

simrep <- 1000  # number of simulation replicates
S <- c(0, 150, 500, 1000, 5000) # levels of single-visit sites S
pvals <- seq(0.1, 0.9, 0.02) # levels for psi/p treatments (psi for occupancy, p for detection)

simnames <- c("Case2x150_CovNull","Case2x300","Case4x150","CovOcc","CovDet","CovBoth")


#  Prepare empty arrays to store estimates  ####

esti <- array(NA, dim=c(length(simnames), 10, length(S), length(pvals), length(pvals), simrep))
dimnames(esti) <- list(simnames,
                       c('psi','SE psi logit','SE psi','p','SE p logit','SE p',
                         'slope(OccCov)','SE slope(OccCov)','slope(DetCov)','SE slope(DetCov)'),
                       paste('SV', S, sep="_"),
                       paste('psi', pvals, sep='_'),
                       paste('p', pvals, sep='_'),
                       NULL)

validity <- array(NA, dim=c(length(simnames), 3, length(S), length(pvals), length(pvals)))
dimnames(validity) <- list(simnames,
                           c('valid','non-valid','NA'),
                           c('SV_0','SV_50','SV_250','SV_1000','SV_10,000'),
                           paste('psi', pvals, sep='_'),
                           paste('p', pvals, sep='_'))

SEs <- c('SE psi logit','SE p logit','SE slope(OccCov)','SE slope(DetCov)')
lmerSE <- array(NA, dim = c(length(simnames), length(SEs), 2, length(pvals), length(pvals)))
dimnames(lmerSE) <- list(simnames,
                         SEs,
                         c('intercept','slope'),
                         dimnames(esti)[[4]],
                         dimnames(esti)[[5]])



#  Start loop: Transfer estimates...  ####

for(a in 1:length(simnames)){
  cat(paste('\n*** Starting', simnames[a], '***'))
  load(paste(simulated_data, simnames[a],".RData", sep=""))

  # Take estimates from individual objects (esti1) into the common object esti
  esti[a,1:6,,,,] <- esti1[1:6,,,,] # parameters common to all scenarios
  if(a %in% c(4,6)){                # parameters only for scenarios with occupancy covariate
    esti[a,c('slope(OccCov)','SE slope(OccCov)'),,,,] <- esti1[c('slope(OccCov)','SE slope(OccCov)'),,,,]
  }
  if(a %in% c(5,6)){                # parameters only for scenarios with detection covariate
    esti[a,c('slope(DetCov)','SE slope(DetCov)'),,,,] <- esti1[c('slope(DetCov)','SE slope(DetCov)'),,,,]
  }

  
  #  ...Calculate validity...  ####
  # count model fitting failures and replace non-valid estimates with NA
  cat(paste('\n*** Counting valid cases in', simnames[a], '***'))
  for(i in 1:length(pvals)){ #
    for(j in 1:length(pvals)){ #
      for(s in 1:length(S)){
        index.na <- which(apply(esti1[c('SE psi logit','SE p logit'),s,i,j,],2,anyNA))
        index.valid <- c()
        for(k in 1:simrep){
          # index simulations where all estimates are valid
          if(sum(is.na(esti1[,s,i,j,k]))==0 &  # is not NA
             esti1['SE psi logit',s,i,j,k]<3 & esti1['SE p logit',s,i,j,k]<3 &         # SEs are not unreasonably high
             esti1['SE psi logit',s,i,j,k]>0.005 & esti1['SE p logit',s,i,j,k]>0.005){ # SEs are not unreasonably low
            index.valid <- c(index.valid, k)
          }
          if(a %in% c(4,6) & k %in% index.valid){  # criterion only for scenarios with occupancy covariate & if otherwise valid
            if(!is.na(esti1['SE slope(OccCov)',s,i,j,k]) & esti1['SE slope(OccCov)',s,i,j,k]>3 |     # SE unreasonably high
               !is.na(esti1['SE slope(OccCov)',s,i,j,k]) & esti1['SE slope(OccCov)',s,i,j,k]<0.005){ # SE unreasonably low
              index.valid <- index.valid[-which(index.valid==k)]
            }
          }
          if(a %in% c(5,6) & k %in% index.valid){  # criterion only for scenarios with detection covariate & if otherwise valid
            if(!is.na(esti1['SE slope(DetCov)',s,i,j,k]) & esti1['SE slope(DetCov)',s,i,j,k]>3 |     # SE unreasonably high
               !is.na(esti1['SE slope(DetCov)',s,i,j,k]) & esti1['SE slope(DetCov)',s,i,j,k]<0.005){ # SE unreasonably low
              index.valid <- index.valid[-which(index.valid==k)]
            }
          }
        }
        esti[a,,s,i,j,-index.valid] <- NA # replace all non-valid cases with NA
        # store counts of NAs, non-valid cases and valid cases
        validity[a,'NA',s,i,j] <- length(index.na)
        validity[a,'non-valid',s,i,j] <- simrep - length(index.na) - length(index.valid)
        validity[a,'valid',s,i,j] <- length(index.valid)
      }
    }
  }

  rm(case,esti1,index.valid,index.na)
  save.image('estimates.RData')

  # # Save plots with proportions of valid estimates
  # mapPalette <- colorRampPalette(c("black","red3","goldenrod2","grey88"))
  # png(paste('validity_',simnames[a],'.png',sep=''), width=800, height=550)
  # par(mfrow=c(2,3), mar=c(3,3,4,2))
  # zlim <- c(0,1)
  # for(s in 1:length(S)){
  #   zvalue <- validity[a,'valid',s,,]/simrep
  #   image(x=pvals, y=pvals, z=zvalue, col=mapPalette(100), axes=T,
  #         xlab="psi", ylab="p", asp=1,
  #         xlim=c(0,1), ylim=c(0,1), zlim=zlim)
  #   title(main=paste("Valid cases with",S[s],"SV \nmin:",round(min(validity[a,'valid',s,,]/simrep),2),
  #                    " max:",max(validity[a,'valid',s,,]/simrep),
  #                    " mean:",round(mean(validity[a,'valid',s,,]/simrep),2)))
  # }
  # plot(x=c(-1,1), y=c(-1,1), axes=F, pch=20, col="white")
  # legend(x=0, y=0, xjust=0.5, yjust=0.5, pch=15, cex=1.6, bty="n",
  #        col=rev(c("black","red3","goldenrod2","forestgreen")), legend=rev(paste(round(seq(zlim[1],zlim[2],length.out=4),2))))
  # dev.off()
  # rm(zvalue,zlim)

  
  # ...Fit linear mixed models to SEs  ####
  cat(paste('\n*** Fitting linear models for', simnames[a], '***'))
  for(i in 1:length(pvals)){
    for(j in 1:length(pvals)){
      # print(paste('Doing psi =',pvals[i],'and p =',pvals[j],'for',simnames[a]))
      if (a==4) {
        sel.SEs <- SEs[1:3]       # incl. SE slope(OccCov)
      } else if (a==5) {
        sel.SEs <- SEs[c(1:2,4)]  # incl. SE slope(DetCov)
      } else if (a==6){
        sel.SEs <- SEs            # all SEs
      } else {
        sel.SEs <- SEs[1:2]       # only SE(psi) and SE(p)
      }
      for(p in sel.SEs){
        # compile data for regression model
        yy <- gather(as.data.frame(esti[a,p,,i,j,]),
                     key="simrepl",     # simulation replicate (1:1000; will be fitted as random effect)
                     value="estimate",  # estimate of the SE
                     1:simrep, factor_key=T)
        yy$InSV <- rep(0:4, times=simrep) # Index value representing the fixed factor level of SV-sites
        # fit regression model
        fm <- try(lmer(estimate ~ InSV + (1|simrepl), data=yy))
        if(isTRUE(class(fm)=="try-error")){
          next # if this model cannot be fitted, the respective lmerSE[] remains NA
        }else{
          lmerSE[a,p,'intercept',i,j] <- fixef(fm)[1]  # intercept of SE
          lmerSE[a,p,'slope',i,j] <- fixef(fm)[2]      # slope of SE
        }
      }
    }

    rm(fm,yy,sel.SEs)
    if(i %in% c(10, 20, 30, 41)){
      cat(paste('\n*** Saving', simnames[a], '***'))
      save.image('estimates.RData')
    }
  }
}



#  Calculate means across simulations  ####

estimean <- array(data=NA, dim=c(6,4,5,41,41))
dimnames(estimean) <- list(simnames,
                           c('SE psi logit','SE p logit','SE slope(OccCov)','SE slope(DetCov)'),
                           dimnames(esti)[[3]],
                           dimnames(esti)[[4]],
                           dimnames(esti)[[5]])
dim(estimean)

for(s in 1:length(simnames)){
  estimean[s,1:2,,,] <- apply(esti[s,c('SE psi logit','SE p logit'),,,,], c(1,2,3,4), mean, na.rm=T)
}
for(s in c(4,6)){ # CovOcc & CovBoth
  estimean[s,'SE slope(OccCov)',,,] <- apply(esti[s,'SE slope(OccCov)',,,,], c(1,2,3), mean, na.rm=T)
}
for(s in c(5,6)){ # CovDet & CovBoth
  estimean[s,'SE slope(DetCov)',,,] <- apply(esti[s,'SE slope(DetCov)',,,,], c(1,2,3), mean, na.rm=T)
}



#  Calculate bias  ####

cat(paste('\n*** Calculating bias ***'))

bias <- array(NA, dim=c(length(simnames), 4, length(S), length(pvals), length(pvals), simrep))
dimnames(bias) <- list(simnames,
                       c('psi','p','slope(OccCov)','slope(DetCov)'),
                       paste('SV', S, sep="_"),
                       paste('psi', pvals, sep='_'),
                       paste('p', pvals, sep='_'),
                       NULL)

biasmedian <- array(data=NA, dim=c(6,4,5,41,41))
dimnames(biasmedian) <- list(simnames,
                             c('psi','p','slope(OccCov)','slope(DetCov)'),
                             dimnames(bias)[[3]],
                             dimnames(bias)[[4]],
                             dimnames(bias)[[5]])
dim(biasmedian)

for(a in 1:length(simnames)){ #
  cat(paste('\n*** Calculating bias for', simnames[a],'***'))
  for(i in 1:length(pvals)){ #
    for(j in 1:length(pvals)){ #
      for(s in 1:length(S)){
        for(k in 1:simrep){
          bias[a,'psi',s,i,j,k] <- esti[a,'psi',s,i,j,k] - pvals[i]
          bias[a,'p',s,i,j,k] <- esti[a,'p',s,i,j,k] - pvals[j]
          if(a %in% c(4,6)){
            bias[a,'slope(OccCov)',s,i,j,k] <- esti[a,'slope(OccCov)',s,i,j,k] - (-1) # -1 is the true (simulated) coefficient
          }
          if(a %in% 5:6){
            bias[a,'slope(DetCov)',s,i,j,k] <- esti[a,'slope(DetCov)',s,i,j,k] - 1 # 1 is the true (simulated) coefficient
          }
        }
      }
    }
  }
}

# Calculate medians across all simulations
# Note: here we use the median (instead of the mean) because the estimates are
#       on the asymmetric probability scale
for(s in 1:length(simnames)){
  biasmedian[s,,,,] <- apply(bias[s,,,,,], c(1,2,3,4), median, na.rm=T)
}



#  Clean and save workspace  ####

cat(paste('\n*** Saving workspace ***'))
rm(a,i,j,k,p,s,J,M,N,mapPalette,SEs,simulated_data,simnames)
save.image('estimates.RData')
