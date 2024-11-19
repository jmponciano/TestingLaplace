rm(list=ls()) # remove everything currently held in the R memory
graphics.off() # close all open graphics window


setwd("/Users/constantingeorgeglen/Documents/Databases/Analyses/Program/TMB/code/LaplaceApproximation/Examples")
source("TestingLaplaceToolkit4.0.R")
source("TestingLaplaceTMB.functions7.0.R") 

setwd("/Users/constantingeorgeglen/Documents/Databases/Analyses/Program/TMB/code/LaplaceApproximation/Examples")
source("TestingLaplaceTMB.models2.0.R") 


### repeat until 500 models have converged
NSIMS      <- 500
counter    <- 0
NN.SIMS    <- list()
tmb.opts   <- list()
# parmsNN    <- list(beta=-1, df1=2, sd1=3, sd2=2) # df=2
parmsNN    <- list(beta=-1, df1=3, sd1=3, sd2=2) # df=2
parms4tmb  <- list(beta=parmsNN$beta, df1=parmsNN$df1, 
                   lnsd1 = log(parmsNN$sd1), lnsd2 = log(parmsNN$sd2))

repeat{
  repeat{
    NN.SIMS.temp <- sim.example(n=50, parm.list=parmsNN, model="Bivar.Norm")
    
    if(var(NN.SIMS.temp[,2])<10000) break
  }
  tempNN <- tryCatch(run.tmb(profile.par.name="beta", 
                             profile.par=seq(-2,2,by=0.1), 
                             known.parm=parms4tmb, single.par.init=0,
                             sim.data=NN.SIMS.temp, model="Bivar.Norm"),
                        error = function(e){ NULL }
                        )
  convergence <- tempNN$MLres$res$pdHess
  # print(tempProbN$MLres$res)
  print( c( counter, is.null(tempNN) )  )
  
  if( !is.null(tempNN) && convergence==T ){
    counter <- counter+1
    NN.SIMS[[counter]]  <- NN.SIMS.temp
    tmb.opts[[counter]] <- tempNN
  }
  
  if(counter==NSIMS){
    break
  }
}
setwd("~/Documents/Databases/Analyses/Program/TMB/code/LaplaceApproximation/Examples/coverageProb")
save(tmb.opts, file="tmboptsNNdf3.RData")
save(parmsNN, NN.SIMS, file="NN500SIMSdf3.RData")
any(unlist( lapply(1:NSIMS, function(x) tmb.opts[[x]]$MLres$res$pdHess ) ) == F)








