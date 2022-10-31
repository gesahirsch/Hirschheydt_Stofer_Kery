setwd('1_Simulate_data_and_fit_occupancy_models/')

case <- 'CovBoth'

library(AHMbook)
library(unmarked)

#  Set variables  ####
simrep <- 1000  # number of simulation replicates per parameter combination

S <- c(0, 150, 500, 1000, 5000) # levels for the number of single-visit sites S
M <- 150        # number of multi-visit sites
N <- M + S      # resulting total number of sites in each treatment

J <- 2          # number of visits to multi-visit sites

pvals <- seq(0.1, 0.9, 0.02) # levels for psi/p treatments (psi for occupancy, p for detection)
length(pvals)



#  Prepare array to save estimation results  ####
esti1 <- array(NA, dim=c(10, length(S), length(pvals), length(pvals), simrep))
dimnames(esti1) <- list(c('psi','SE psi logit','SE psi','p','SE p logit','SE p',
                          'slope(OccCov)','SE slope(OccCov)','slope(DetCov)','SE slope(DetCov)'),
                        paste('SV', S, sep="_"),
                        paste('psi', pvals, sep='_'), 
                        paste('p', pvals, sep='_'), 
                        NULL)
str(esti1)



#  Run simulation - 9 days 23 hours 52 minutes  ####

seed <- 1 # start a seed for reproducibility

for(k in 1:simrep){
  cat(paste('\n*** Doing simrep', k, '***'))
  for(i in 1:length(pvals)){     # Loop over levels of psi
    for(j in 1:length(pvals)){   # Loop over levels of p
      
      repeat{ # to ensure that data set has at least 1 detection & at least 1 non-detection
        
        # simulate data set
        set.seed(seed)
        dat <- simOcc(M=max(N), J=J, mean.occupancy=pvals[i], 
                      beta1=-1, beta2=0, beta3=0,                               # beta1=-1 coefficient of occupancy covariate OccCov
                      mean.detection=pvals[j], time.effects=c(0, 0), 
                      alpha1=0, alpha2=1, alpha3=0, sd.lp=0, b=0, show.plots=F) # alpha2=1 coefficient of detection covariate DetCov
        
        if(sum(dat$y[1:M,]) > 0 & sum(dat$y[1:M,]) < J*M){
          break
        }
        seed <- seed + 1
      }
      
      # turn additional sites into SV sites
      dat$y[(M+1):nrow(dat$y),2:J] <- NA
      dat$wind[(M+1):nrow(dat$y),2:J] <- NA # DetCov is called 'wind' in simOcc() function
      
      for(s in 1:length(S)){  # loop over levels of S
        
        # grab detection histories
        y <- dat$y[1:N[s],]
        OccCov <- dat$elev[1:N[s]]   # OccCov is called 'elev' in simOcc() function
        DetCov <- dat$wind[1:N[s],]  # DetCov is called 'wind' in simOcc() function
        
        # create unmarked dataframe
        umf <- unmarkedFrameOccu(y=y, siteCovs=data.frame(OccCov=OccCov),
                                 obsCovs=list(DetCov=DetCov))
        
        # provide starting values
        set.seed(seed)
        int <- c(qlogis(runif(1)), rnorm(1), # occupancy parameters (intercept, slope)
                 qlogis(runif(1)), rnorm(1)) # detection parameters (intercept, slope)
        seed <- seed + 1
        
        # fit model
        fm <- occu(~DetCov ~OccCov, data=umf, starts=int)
        
        # save estimates
        esti1[c('psi','SE psi'),s,i,j,k] <- unlist(predict(fm, "state",  newdata=data.frame(OccCov=0)))[1:2] # estimate of psi & SE(psi) on transformed/probability scale
        esti1[c('p','SE p'),s,i,j,k] <- unlist(predict(fm, "det",  newdata=data.frame(DetCov=0)))[1:2]       # estimate of p & SE(p) on transformed/probability scale
        if(sum(is.na(esti1[c('psi','SE psi','p','SE p'),s,i,j,k]))==0){           # only works when no NAs among estimates
          esti1['SE psi logit',s,i,j,k] <- sqrt(vcov(fm)['psi(Int)','psi(Int)'])  # SE(psi) on logit scale
          esti1['SE p logit',s,i,j,k] <- sqrt(vcov(fm)['p(Int)','p(Int)'])        # SE(p) on logit scale
          esti1[c('slope(OccCov)','slope(DetCov)'),s,i,j,k] <- fm@opt$par[c(2,4)] # slopes/coefficients of the two covariates (logit)
          esti1['SE slope(OccCov)',s,i,j,k] <- sqrt(vcov(fm)['psi(OccCov)','psi(OccCov)'])  # SE of slope(OccCov) (logit)
          esti1['SE slope(DetCov)',s,i,j,k] <- sqrt(vcov(fm)['p(DetCov)','p(DetCov)'])      # SE of slope(DetCov) (logit)
        }
      }
    }
  }
  SS <- 20 # define number of simreps before new output is saved
  if(k %% SS == 0){
    save.image(paste(case, '.RData', sep=''))
  }
}



#  Save cleaned RData  ####

rm(i,j,k,s,dat,y,OccCov,DetCov,umf,int,fm,SS,seed,simOcc)
save.image(paste(case, ".RData", sep=""))
