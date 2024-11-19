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
PN.SIMS    <- list()
tmb.opts   <- list()
# parmsPN    <- list(beta=-1, df1=2, sd1=3) # df=2
parmsPN    <- list(beta=-1, df1=3, sd1=3) # df=3
parms4tmb  <- list(beta=parmsPN$beta, df1=parmsPN$df1, 
                   lnsd1 = log(parmsPN$sd1))

repeat{
  repeat{
    PN.SIMS.temp <- sim.example(n=50, parm.list=parmsPN, model="Pois.Norm")
    
    if(var(PN.SIMS.temp[,2])<10000) break
  }
  
  tempPN <- tryCatch(run.tmb(profile.par.name="beta", 
                             profile.par=seq(-2,0,by=1), 
                             known.parm=parms4tmb, single.par.init = 0.1,
                             sim.data=PN.SIMS.temp, model="Pois.Norm"),
                        error = function(e){ NULL }
                        )
  convergence <- tempPN$MLres$res$pdHess
  # print(tempProbN$MLres$res)
  print( c( counter, is.null(tempPN) )  )
  
  if( !is.null(tempPN) && convergence==T ){
    counter <- counter+1
    PN.SIMS[[counter]]  <- PN.SIMS.temp
    tmb.opts[[counter]] <- tempPN
  }
  
  if(counter==NSIMS){
    break
  }
}
setwd("~/Documents/Databases/Analyses/Program/TMB/code/LaplaceApproximation/Examples/coverageProb/")
save(tmb.opts, file="tmboptsPNdf3.RData")
save(parmsPN, PN.SIMS, file="PN500SIMSdf3.RData")
any(unlist( lapply(1:NSIMS, function(x) tmb.opts[[x]]$MLres$res$pdHess ) ) == F)








