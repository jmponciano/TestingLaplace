source("JAGS4Profiles.R")
source("simulate.data.R")
library("doParallel")




# Simulating from Model 1: Latent variable is X ~ T_2(0), 
# a reasonably flat distribution 
# W = X + eta where eta ~ N(0,sigma1)
# Y ~ pois( lam= beta*X)
# Response = c(W,Y), a Normal-Poisson given X=x
parms.err.varsProbN <- list(beta=-1, df1=10, sd1=3 )
out.simProbN <- sim.example(n=100, parm.list=parms.err.varsProbN, model="Probit.Norm")
sim.data <- out.simProbN

n.iter <- 10000
n.chains <- 5
n.adapt <- 1000
n.update <- 1000
thin <- 10


n <- nrow(sim.data)
W <- sim.data[,1]
Y <- sim.data[,2]
Y

cl <- makePSOCKcluster(n.chains)
datalist <- list(K=1, n=n, df1=parms.err.varsProbN$df1, 
                 Y=dcdim(data.matrix(Y)), W=dcdim(data.matrix(W)))
dcrunProbN  <- dc.parfit(cl=cl, datalist, c("beta", "lsd1"), ProbitNormDC.MLE, 
                         n.clones=c(1,2,4,8), inits = list("beta"=-0.1, "lsd1"=1),
                         multiply="K", unchanged=c("n","df1"),
                         n.chains = n.chains, n.adapt=n.adapt, n.update=n.update, 
                         n.iter = n.iter, thin=thin)
### had to use initial to specify where to start or estimates returned are wildly off

ModProbN.mles <- as.list(summary(dcrunProbN)[[1]][,1])
stopCluster(cl)

dcdiag(dcrunProbN)

confints.mls <- confint(dcrunProbN)
confints.mls




# MLES to profile over beta
beta.mle <- ModProbN.mles$beta
lsd1.mle <- ModProbN.mles$lsd1



bracket   <- 3
prec.prof <- 0.1
beta.profs <- seq(from=beta.mle-bracket, to=beta.mle+bracket, by=prec.prof)
len.prof <- length(beta.profs)

cl <- makeCluster(parallel::detectCores(),outfile="")
registerDoParallel(cl)
res <- foreach(i=1:len.prof, .combine=rbind)%dopar%{
  
  outlist.lsd1 <- list()
  ith.beta     <- beta.profs[i]
  #cl <- makePSOCKcluster(n.chains)
  ith.datalist <- list(K=1, n=n, beta=ith.beta, df1=parms.err.varsProbN$df1, 
                       Y=dclone::dcdim(data.matrix(Y)), 
                       W=dclone::dcdim(data.matrix(W)))
  clones.seq <- 8
  ith.dcrun <- dclone::dc.fit(ith.datalist, c("lsd1"), 
                              ProbitNorm.betaprof, n.clones=clones.seq,
                              multiply="K", unchanged=c("n","beta", "df1"),
                              n.chains = n.chains, n.adapt=n.adapt, n.update=n.update, 
                              n.iter = n.iter, thin=thin)
  profhatsProbN <- as.list(summary(ith.dcrun)[[1]][,1])
  #stopCluster(cl)
  outlist.lsd1 <- rbind(outlist.lsd1, profhatsProbN[1])
  return(list(lsd1=outlist.lsd1))
}
stopCluster(cl)



# This is the vector of mles of lsd1 and lsd2 for every fixed value of beta
lsdhats4prof <- cbind(lsd1=unlist(res[,1])) 





# Now prepare for the GT algorithm by sampling from g(x|y)
# Sampling from f(X|(W,Y)) 
n.iter <- 10000
n.chains <- 4
n.adapt <- 1000
n.update <- 1000
thin <- 2
datalist <- list(W=W, Y=Y, n=n, 
                 beta=ModProbN.mles$beta, 
                 df1=parms.err.varsProbN$df1, lsd1=ModProbN.mles$lsd1)
out.parms <- c("X")
XpostSample <- dc.fit(datalist, params=out.parms, model=ProbitNorm.hXgY, 
                      multiply=NULL, unchanged=c("n", "df1","beta","lsd1"), 
                      n.chains=n.chains,
                      n.adapt=n.adapt, n.update=n.update, 
                      n.iter = n.iter, thin=thin,n.clones=1) 
XpostBinned <- do.call(rbind, XpostSample)





# Now the GT algorithm itself:

B <- nrow(XpostBinned)
prof.vec <- rep(0,len.prof)

for(i in 1:len.prof){
  
  ith.beta <- beta.profs[i]
  ith.sd1 <- exp(lsdhats4prof[i,1])
  
  GT.logdiff <- rep(0,B)
  for(j in 1:B){
    X.star     <- as.vector(XpostBinned[j,])
    pvec.num   <- pnorm(q=ith.beta*X.star, mean=0, sd=1)
    pvec.denom <- pnorm(q=ModProbN.mles$beta*X.star, mean=0,sd=1)
    pvec.num[pvec.num==1]<- 0.99
    pvec.num[pvec.num==0] <- .Machine$double.xmin
    pvec.denom[pvec.denom==1]<- 0.99
    pvec.denom[pvec.denom==0] <- .Machine$double.xmin
    
    GT.lnumer.persamp <- sum(Y*log(pvec.num) + (1-Y)*log(1-pvec.num)
                             + dnorm(x=W, mean=X.star,sd=ith.sd1, log=TRUE) +
                               dt(x=X.star, df=parms.err.varsProbN$df1, ncp=0,log=TRUE)) 
    GT.ldenom.persamp <- sum(Y*log(pvec.denom) + (1-Y)*log(1-pvec.denom)
                             + dnorm(x=W, mean=X.star,sd=exp(ModProbN.mles$lsd1), log=TRUE) +
                               dt(x=X.star, df=parms.err.varsProbN$df1, ncp=0,log=TRUE))
    
    GT.logdiff[j] <- GT.lnumer.persamp-GT.ldenom.persamp
  }
  
  ratio.prof.vec.star <- exp(GT.logdiff)
  prof.vec[i] <- mean(ratio.prof.vec.star)
  #print(c(  mean(ratio.prof.vec.star) ))
  
}

### Model using TMB
source("TestingLaplaceTMB.functions7.0.R")

ModProbN.TMB <- run.tmb(profile.par.name="beta", search.initial=F, n.search=200,
                        sim.data=sim.data, model="Probit.Norm", method.ci="useconfint",
                        ytol=5, init = list(par=1, hidden=0.1, df1=parms.err.varsProbN$df1))
ModProbN.TMB$model.obj
ModProbN.TMB$model.obj$par


ModProbN.TMB$model.obj$par
check.model <- sdreport(ModProbN.TMB$fun.obj)
check.model$pdHess
ModProbN.TMB$convergence
sqrt(diag(check.model$cov.fixed))
check.model
sqrt(diag(solve(ModProbN.TMB$model.obj$hessian)))


# manual computing sd
import::from("numDeriv", "hessian")
(std <- sqrt(diag(solve(numDeriv::hessian(ModProbN.TMB$fun.obj$fn, ModProbN.TMB$fun.obj$env$last.par.best[1:2])))))
(std <- sqrt(diag(solve(numDeriv::jacobian(ModProbN.TMB$fun.obj$gr,ModProbN.TMB$fun.obj$env$last.par.best[1:2])))))
par.fixed <-ModProbN.TMB$fun.obj$env$last.par.best[1:2]
hessian.fixed <- optimHess(par.fixed, ModProbN.TMB$fun.obj$fn,
                           ModProbN.TMB$fun.obj$gr) ## Marginal precision of theta.
chol(hessian.fixed) # pdHess
chol(ModPN.TMB$model.obj$hessian)




### Model for DD
source("TestingLaplaceDD.functions2.0.R")
source("TestingLaplaceDD.models2.0.R")

MCMCset=list('n.iter.dc'=n.iter,'n.chains.dc'=3,
             'n.iter.dd'=n.iter,'n.chains.dd'=3,
             'n.burn'=n.update,'n.thin'=thin,'n.adapt'=n.adapt)

test.DD.multi.par <- DD_wrapper( profile.par = "beta", model = "Probit.Norm", 
                                 known.params = list(df1=parms.err.varsProbN$df1),
                                 df=sim.data, psi.grid = beta.profs, 
                                 MCMCset=MCMCset, 
                                 CIalpha=0.95, n.clones=c(16))

# Saving objects for joint plot
MP.ProbNtmb      <- data.frame(betavals=ModProbN.TMB$ci.obj$data$beta, tmbrel.like = 
                              ModProbN.TMB$ci.obj$data$rel.like)
MP.ProbNgt       <- data.frame(beta.profs=beta.profs, prof.vec=prof.vec)

MP.ProbNdd <- data.frame(DDxs = test.DD.multi.par$DD.profile.q$Profile.out$psi.grid,
                      DDys = exp(test.DD.multi.par$DD.profile.q$Profile.out$profile.like))

MP.ProbNDCwaldCI <- data.frame(lower=confints.mls[,1], mles=unlist(ModProbN.mles), upper=confints.mls[,2])


### Plot
par(oma=c(1,1,1,1), mar=c(4,5,2,2))
plot(MP.ProbNgt$beta.profs, MP.ProbNgt$prof.vec, type="l", ylim=c(0,1),
     xlab=bquote("Values of" ~ beta), ylab="Profile Likelihood", bty="l",
     cex.lab=1.5, xlim=c(-2, 2), lty=3, lwd=1.5)
#abline(v=-1); abline(h=1)
abline(v=ModProbN.mles$beta, lty=2, lwd=1)
#abline(v=parms.err.varsNN$beta, lty=4)

### plot lines for TMB
lines(MP.ProbNtmb$betavals, MP.ProbNtmb$tmbrel.like, 
      type="l",col="black", lty=2, lwd=1.5);


### plot lines for alg 1 and 2
lines(MP.ProbNdd$DDxs,MP.ProbNdd$DDys,lwd=1.5,lty=1,col="black",type="l")
#lines(test.DD.multi.par$DD.profile.psi$Profile.out$psi.grid, ### ALG 1
#      exp(test.DD.multi.par$DD.profile.psi$Profile.out$profile.like), 
#      lwd=2,lty=1,col="blue",type="l")

### Wald-Type CI and legend
#abline(h=exp(-qchisq(p=0.95,df=1)/2),lty=2, col="gray")

arrows(x0=confints.mls[1,1],x1=confints.mls[1,2],y0=0,y1=0, angle=90, code=3, length=0.06, lwd=2)
#points(x=confints.mls[1,],y=c(0,0), col="black", type="l" ,lty=1,lwd=2)

legend("topleft", 
       legend=c("DD-2","TMB", "GT"), 
       col=c("black", "black","black"),
       lwd=rep(1.5,3), 
       lty=c(1:3), bty="n", cex=0.9)


#save.image("ProbNbetaprof.RData")


