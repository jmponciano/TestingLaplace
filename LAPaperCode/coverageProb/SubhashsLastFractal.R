source("TestingLaplaceToolkit4.0.R")
source("TestingLaplaceTMB.functions7.0.R") 
source("TestingLaplaceTMB.models2.0.R") 


# Loading the same simulated data using the same parameter values:

parms.err.varsProbN <- list(beta=-1, df1=2, sd1=3)
out.simProbN <- sim.example(n=50, parm.list=parms.err.varsProbN, model="Probit.Norm")

ProbN.SIMS <- list()
NSIMS <- 500
for(i in 1:NSIMS){
  ProbN.SIMS[[i]] <- sim.example(n=50, parm.list=parms.err.varsProbN, model="Probit.Norm")
}
#save(ProbN.SIMS, file="ProbN500SIMS.RData")

# Getting the estimates from TMB for all data sets
tmb.opts <- list()

parms4tmb <- list(beta=beta, df1=df1,
                  lnsd1 = log(sd1))
for(i in 1:NSIMS){
  testProbN <- run.tmb(profile.par.name="beta", profile.par=seq(-2,2,by=0.1), 
                     known.parm=parms4tmb, sim.data=ProbN.SIMS[[i]], model="Probit.Norm")

  tmb.opts[[i]] <- testProbN$MLres$opt.ML
}
#save(tmb.opts, file="tmboptsProbN.RData")




#load("~/Dropbox/GeorgeGlen/TMB/WhereTMBFails/Examples/PaperDatasets/ProbitNdata.RData")

beta <- parms.err.varsProbN$beta
sd1 <- parms.err.varsProbN$sd1
df1 <- parms.err.varsProbN$df1

W <- out.simProbN[,1]
Y <- out.simProbN[,2]

# First we simulate a large number of samples from the t distrib
B<- 50000
X.B <- rt(n=B*length(Y), df=df1)

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

negll(betaguess=-1, W=W, Y=Y, sd1=sd1, X.B=X.B)

install.packages("tictoc")
library(tictoc)
tic()
out <- optim(par=-1, fn=negll, method="BFGS", W=W, Y=Y, sd1=sd1, X.B=X.B, hessian=TRUE)
toc()

#### Now computing TMB MLE:

parms4tmb <- list(beta=beta, df1=df1,
                  lnsd1 = log(sd1))

testProbN <- run.tmb(profile.par.name="beta", profile.par=seq(-2,2,by=0.1), 
                     known.parm=parms4tmb, sim.data=out.simProbN, model="Probit.Norm")

testProbN$MLres$res
