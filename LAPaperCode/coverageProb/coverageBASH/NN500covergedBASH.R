##### ---------------------------------------------------------------
##### Coverage probability for Normal-Normal model
##### Author: José Miguel Ponciano and C. George Glen
##### Date: 12/05/2024
##### ---------------------------------------------------------------
wd.path.name  <- "/home/george.glen/cmd/TMBPaper_MCcoverage"
setwd(wd.path.name)

load("NN500SIMSdf3.RData")

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
    Y.like <- dnorm(rep(Y[i],B), mean=beta*X.B, sd=sd2)
    log.mean.like.vec[i] <- log(mean(Y.like*W.like))
  }
  nloglike <- -sum(log.mean.like.vec)
  return(nloglike)
  
}
write.csv(matrix(data=NA),  file=paste0(wd.path.name,"/testNN.csv"))

### Procedure section
cl <- makeCluster(11)
registerDoParallel(cl)
res <- foreach(i=1:length(NN.SIMS), .export = ls(.GlobalEnv)) %dopar%{
  
  temp.data <- NN.SIMS[[i]]
  
  # First we simulate a large number of samples from the t distrib
  parms.NN <- parmsNN
  B       <- 50000
  X.B <- rt(n=B*length(temp.data[,2]), df=parms.NN$df1)
  out <- optim(par=-1, fn=negll, 
               # method="Nelder-Mead", 
               method="BFGS", 
               W=temp.data[,1], Y=temp.data[,2], 
               sd1=parms.NN$sd1, sd2=parms.NN$sd2, X.B=X.B, hessian=TRUE)

  # save everything into an object
  resObj <- list( obj = out )
  
  return(resObj)
} 
stopCluster(cl = cl)


### SAVE the RDATA
save.image(file = paste0(wd.path.name,"/optimBFGS_NN500covergedBASHdf3.RData"))
# save.image(file = "NN500covergedBASH.RData")

