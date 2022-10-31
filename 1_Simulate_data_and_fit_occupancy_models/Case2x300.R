setwd('1_Simulate_data_and_fit_occupancy_models/')

case <- 'Case2x300'

library(AHMbook)
library(unmarked)


#  Set variables  ####
simrep <- 1000  # number of simulation replicates per parameter combination

S <- c(0, 150, 500, 1000, 5000) # levels for the number of single-visit sites S
M <- 300        # number of multi-visit sites
N <- M + S      # resulting total number of sites in each treatment

J <- 2          # number of visits to multi-visit sites

pvals <- seq(0.1, 0.9, 0.02) # levels for psi/p treatments (psi for occupancy, p for detection)
length(pvals)



#  Prepare array to save estimation results  ####
esti1 <- array(NA, dim=c(6, length(S), length(pvals), length(pvals), simrep))
dimnames(esti1) <- list(c('psi','SE psi logit','SE psi','p','SE p logit','SE p'),
                        paste('SV', S, sep="_"),
                        paste('psi', pvals, sep='_'), 
                        paste('p', pvals, sep='_'), 
                        NULL)
str(esti1)



#  Run simulation - 6 days 14 hours 33 minutes  ####

seed <- 1 # start a seed for reproducibility

# for(k in 1:simrep){
#   cat(paste('\n*** Doing simrep', k, '***'))
#   for(i in 1:length(pvals)){     # Loop over levels of psi
#     for(j in 1:length(pvals)){   # Loop over levels of p
#       
#       repeat{ # to ensure that data set has at least 1 detection & at least 1 non-detection
#         
#         # simulate data set
#         set.seed(seed)
#         dat <- simOcc(M=max(N), J=J, mean.occupancy=pvals[i], 
#                       beta1=0, beta2=0, beta3=0, 
#                       mean.detection=pvals[j], time.effects=c(0, 0), 
#                       alpha1=0, alpha2=0, alpha3=0, sd.lp=0, b=0, show.plots=F)
#         
#         if(sum(dat$y[1:M,]) > 0 & sum(dat$y[1:M,]) < J*M){
#           break
#         }
#         seed <- seed + 1
#       }
#       
#       # turn additional sites into SV sites
#       dat$y[(M+1):nrow(dat$y),2:J] <- NA
#       
#       for(s in 1:length(S)){  # loop over levels of S
#         
#         # grab detection histories
#         y <- dat$y[1:N[s],]
#         
#         # create unmarked dataframe
#         umf <- unmarkedFrameOccu(y=y)
#         
#         # provide starting values
#         set.seed(seed)
#         int <- c(qlogis(runif(2))) # two intercepts
#         seed <- seed + 1
#         
#         # fit model
#         fm <- occu(~1 ~1, data=umf, starts=int)
#         
#         # save estimates
#         esti1[c('psi','SE psi'),s,i,j,k] <- unlist(predict(fm, "state",  newdata=data.frame(1)))[1:2] # estimate of psi & SE(psi) on transformed/probability scale
#         esti1[c('p','SE p'),s,i,j,k] <- unlist(predict(fm, "det",  newdata=data.frame(1)))[1:2]       # estimate of p & SE(p) on transformed/probability scale
#         if(sum(is.na(esti1[c('psi','SE psi','p','SE p'),s,i,j,k]))==0){ # only works when no NAs among estimates
#           esti1['SE psi logit',s,i,j,k] <- sqrt(vcov(fm)[1,1]) # SE(psi) on logit scale
#           esti1['SE p logit',s,i,j,k] <- sqrt(vcov(fm)[2,2])   # SE(p) on logit scale
#         }
#       }
#     }
#   }
#   SS <- 20 # define number of simreps before new output is saved
#   if(k %% SS == 0){
#     save.image(paste(case, '.RData', sep=''))
#   }
# }



#  Save cleaned RData  ####

rm(i,j,k,s,dat,y,umf,int,fm,SS,seed,simOcc)
save.image(paste(case, ".RData", sep=""))
