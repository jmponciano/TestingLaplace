setwd("/Users/constantingeorgeglen/Documents/Databases/Analyses/Program/TMB/code/LaplaceApproximation/Examples")
source("TestingLaplaceToolkit4.0.R")
source("TestingLaplaceTMB.functions7.0.R") 

setwd("/Users/constantingeorgeglen/Documents/Databases/Analyses/Program/TMB/code/LaplaceApproximation/Examples")
source("TestingLaplaceTMB.models2.0.R") 


### repeat until 500 models have converged
NSIMS      <- 500
counter    <- 0
LN.SIMS <- list()
tmb.opts   <- list()
parmsLN <- list(beta=-1, df1=10, sd1=3)
parms4tmb  <- list(beta=parmsLN$beta, df1=parmsLN$df1, 
                   lnsd1 = log(parmsLN$sd1))

repeat{
  repeat{
    # variance doesnt matter since Y is 0 and 1 so every iteration is accepted
    LN.SIMS.temp <- sim.example(n=50, parm.list=parmsLN, model="Probit.Norm")
    # print(var(LN.SIMS.temp[,2]))
    if(var(LN.SIMS.temp[,2])<10000) break
  }
  tempLN <- tryCatch(run.tmb(profile.par.name="beta", profile.par=seq(-2,2,by=0.1), 
                             known.parm=parms4tmb, single.par.init = 0,
                             sim.data=LN.SIMS.temp, model="Probit.Norm"),
                     error = function(e){ NULL }
  )
  convergence <- tempLN$MLres$res$pdHess
  # print(tempLN$MLres$res)
  print( c( counter, is.null(tempLN) )  )
  
  if( !is.null(tempLN) && convergence==T ){
    counter <- counter+1
    LN.SIMS[[counter]] <- LN.SIMS.temp
    tmb.opts[[counter]]   <- tempLN
  }
  
  if(counter==NSIMS){
    break
  }
}
setwd("~/Documents/Databases/Analyses/Program/TMB/code/LaplaceApproximation/Examples/coverageProb")
save(tmb.opts, file="tmboptsLN.RData")
save(parmsLN, LN.SIMS, file="LN500SIMS.RData")
any(unlist( lapply(1:NSIMS, function(x) tmb.opts[[x]]$MLres$res$pdHess ) ) == F)








