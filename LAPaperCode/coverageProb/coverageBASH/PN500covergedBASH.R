##### ---------------------------------------------------------------
##### Coverage probability for Poisson-Normal model
##### Author: José Miguel Ponciano and C. George Glen
##### Date: 12/05/2024
##### ---------------------------------------------------------------
rm(list=ls()) # remove everything currently held in the R memory
graphics.off() # close all open graphics window

wd.path.name  <- "/home/george.glen/cmd/TMBPaper_MCcoverage"
# wd.path.name  <- "~/Documents/Databases/Analyses/Program/TMB/code/LaplaceApproximation/Examples/coverageProb/"
setwd(wd.path.name)

load("PN500SIMSdf3.RData")

library(doParallel)
library(parallel)
library(foreach)


### negative log-likelihood function
negll <- function(betaguess,W, Y, sd1=sd1, sd2=sd2, X.B=X.B){
  
  beta <- betaguess
  n <- length(Y)
  log.mean.like.vec <- rep(0,n)
  
  for(i in 1:n){
    
    W.like <- dnorm(rep(W[i],B), mean=X.B,sd=sd1)
    Y.like <- dpois(rep(Y[i],B), lambda=exp(beta*X.B))
    log.mean.like.vec[i] <- log(mean(Y.like*W.like))
  }
  nloglike <- -sum(log.mean.like.vec)
  return(nloglike)
  
}
write.csv(matrix(data=NA),  file=paste0(wd.path.name,"/testPN.csv"))

### Procedure section
cl <- makeCluster(15)
registerDoParallel(cl)
len.PN <- length(PN.SIMS)
res <- foreach(i=1:len.PN, .export = ls(.GlobalEnv)) %dopar%{
  
  temp.data <- PN.SIMS[[i]]
  
  # First we simulate a large number of samples from the t distrib
  parms.PN <- parmsPN
  B       <- 50000
  X.B <- rt(n=B*length(temp.data[,2]), df=parms.PN$df1)
  out <- optim(par=-1, fn=negll, 
               # method="Nelder-Mead", 
               method="BFGS", 
               W=temp.data[,1], Y=temp.data[,2], 
               sd1=parms.PN$sd1, X.B=X.B, hessian=TRUE)

  # save everything into an object
  resObj <- list( obj = out )
  
  return(resObj)
} 
stopCluster(cl = cl)


### SAVE the RDATA
save.image(file = paste0(wd.path.name,"/optimBFGS_PN500covergedBASHdf3.RData"))
# save.image(file = "PN500covergedBASH.RData")

