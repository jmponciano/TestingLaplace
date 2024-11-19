##### ---------------------------------------------------------------
##### Coverage probability for Probit-Normal model
##### Author: José Miguel Ponciano and C. George Glen
##### Date: 12/05/2024
##### ---------------------------------------------------------------
wd.path.name  <- "/home/george.glen/cmd/TMBPaper_MCcoverage"
setwd(wd.path.name)

load("ProbN500SIMS.RData")

library(doParallel)
library(parallel)
library(foreach)

### negative log-likelihood function
negll <- function(betaguess,W,Y, sd1=sd1,X.B=X.B){
  
  beta <- betaguess
  n <- length(Y)
  log.mean.like.vec <- rep(0,n)
  
  for(i in 1:n){
    
    W.like <- dnorm(rep(W[i],B), mean=X.B,sd=sd1)
    Y.like <- dbinom(rep(Y[i],B),size=1,prob=pnorm(q=beta*X.B))
    log.mean.like.vec[i] <- log(mean(Y.like*W.like))
  }
  nloglike <- -sum(log.mean.like.vec)
  return(nloglike)
  
}
write.csv(matrix(data=NA),  file=paste0(wd.path.name,"/testProbN.csv"))

### Procedure section
cl <- makeCluster(15)
registerDoParallel(cl)
len.ProbN <- length(ProbN.SIMS)
res <- foreach(i=1:len.ProbN, .export = ls(.GlobalEnv)) %dopar%{
  
  temp.data <- ProbN.SIMS[[i]]
  
  # First we simulate a large number of samples from the t distrib
  parms.ProbN <- list(beta=-1, df1=10, sd1=3)
  B           <- 500000
  X.B <- rt(n=B*length(temp.data[,2]), df=parms.ProbN$df1)
  out <- optim(par=-1, fn=negll, 
               method="Nelder-Mead", 
               # method="BFGS", 
               W=temp.data[,1], Y=temp.data[,2], 
               sd1=parms.ProbN$sd1, X.B=X.B, hessian=TRUE)

  # save everything into an object
  resObj <- list( obj = out )
  
  return(resObj)
} 
stopCluster(cl = cl)


### SAVE the RDATA
save.image(file = paste0(wd.path.name,"/ProbN500covergedBASH.RData"))

