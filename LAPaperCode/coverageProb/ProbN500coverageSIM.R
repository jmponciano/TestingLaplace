setwd("/Users/constantingeorgeglen/Documents/Databases/Analyses/Program/TMB/code/LaplaceApproximation/Examples")
source("TestingLaplaceToolkit4.0.R")
source("TestingLaplaceTMB.functions7.0.R") 

setwd("/Users/constantingeorgeglen/Documents/Databases/Analyses/Program/TMB/code/LaplaceApproximation/Examples")
source("TestingLaplaceTMB.models2.0.R") 


### repeat until 500 models have converged
NSIMS      <- 500
counter    <- 0
ProbN.SIMS <- list()
tmb.opts   <- list()
parmsProbN <- list(beta=-1, df1=10, sd1=3)
parms4tmb  <- list(beta=parmsProbN$beta, df1=parmsProbN$df1, 
                   lnsd1 = log(parmsProbN$sd1))

repeat{
  repeat{
    # variance doesnt matter since Y is 0 and 1 so every iteration is accepted
    ProbN.SIMS.temp <- sim.example(n=50, parm.list=parmsProbN, model="Probit.Norm")
    if(var(ProbN.SIMS.temp[,2])<10000) break
  }
  tempProbN <- tryCatch(run.tmb(profile.par.name="beta", profile.par=seq(-2,2,by=0.1), 
                                known.parm=parms4tmb, single.par.init = 0,
                                sim.data=ProbN.SIMS.temp, model="Probit.Norm"),
                        error = function(e){ NULL }
  )
  convergence <- tempProbN$MLres$res$pdHess
  # print(tempProbN$MLres$res)
  print( c( counter, is.null(tempProbN) )  )
  
  if( !is.null(tempProbN) && convergence==T ){
    counter <- counter+1
    ProbN.SIMS[[counter]] <- ProbN.SIMS.temp
    tmb.opts[[counter]]   <- tempProbN
  }
  
  if(counter==NSIMS){
    break
  }
}
setwd("~/Documents/Databases/Analyses/Program/TMB/code/LaplaceApproximation/Examples/coverageProb")
save(tmb.opts, file="tmboptsProbN.RData")
save(parmsProbN, ProbN.SIMS, file="ProbN500SIMS.RData")
any(unlist( lapply(1:NSIMS, function(x) tmb.opts[[x]]$MLres$res$pdHess ) ) == F)








