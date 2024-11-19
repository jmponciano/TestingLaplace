rm(list=ls()) # remove everything currently held in the R memory
graphics.off() # close all open graphics windows 


source("JAGS4Profiles.R")
source("simulate.data.R")
library("doParallel")




# Simulating from Model 1: Latent variable is X ~ T_2(0), 
# a reasonably flat distribution 
# W = X + eta where eta ~ N(0,sigma1)
# Y ~ pois( lam= beta*X)
# Response = c(W,Y), a Normal-Poisson given X=x
parms.err.varsPN <- list(beta=-1, df1=2, sd1=3 )
out.simPN <- sim.example(n=100, parm.list=parms.err.varsPN, model="Pois.Norm")
sim.data <- out.simPN

n.iter <- 10000
n.chains <- 5
n.adapt <- 1000
n.update <- 1000
thin <- 2


n <- nrow(sim.data)
W <- sim.data[,1]
Y <- sim.data[,2]
Y

cl <- makePSOCKcluster(n.chains)
datalist <- list(K=1, n=n, df1=parms.err.varsPN$df1, Y=dcdim(data.matrix(Y)), W=dcdim(data.matrix(W)))
dcrunPN  <- dc.parfit(cl=cl, datalist, c("beta", "lsd1"), PoisNormDC.MLE, 
                      n.clones=c(1,2,4,8,16),
                      multiply="K", unchanged=c("n","df1"),
                      n.chains = n.chains, n.adapt=n.adapt, n.update=n.update, 
                      n.iter = n.iter, thin=thin)


ModPN.mles <- as.list(summary(dcrunPN)[[1]][,1])
stopCluster(cl)
#Mod1.mles
#dcdiag(Mod1.dc.all)

confints.mls <- confint(dcrunPN)
confints.mls



# MLES to profile over beta
beta.mle <- ModPN.mles$beta
lsd1.mle <- ModPN.mles$lsd1



bracket <- 3
prec.prof <- 0.1
beta.profs <- seq(from=beta.mle-bracket, to=beta.mle+bracket, by=prec.prof)
len.prof <- length(beta.profs)

cl <- makeCluster(parallel::detectCores(),outfile="")
registerDoParallel(cl)
res <- foreach(i=1:len.prof, .combine=rbind)%dopar%{
  
  outlist.lsd1 <- list()
  outlist.lsd2 <- list()
  ith.beta     <- beta.profs[i]
  #cl <- makePSOCKcluster(n.chains)
  ith.datalist <- list(K=1, n=n, beta=ith.beta, 
                       df1=parms.err.varsPN$df1, 
                       Y=dclone::dcdim(data.matrix(Y)), 
                       W=dclone::dcdim(data.matrix(W)))
  clones.seq <- 16
  ith.dcrun <- dclone::dc.fit(ith.datalist, c("lsd1"), PoisNorm.betaprof, n.clones=clones.seq,
                              multiply="K", unchanged=c("n","beta", "df1"),
                              n.chains = n.chains, n.adapt=n.adapt, n.update=n.update, 
                              n.iter = n.iter, thin=thin)
  profhatsPN <- as.list(summary(ith.dcrun)[[1]][,1])
  #stopCluster(cl)
  outlist.lsd1 <- rbind(outlist.lsd1, profhatsPN[1])
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
                 beta=ModPN.mles$beta, df1=parms.err.varsPN$df1, lsd1=ModPN.mles$lsd1)
out.parms <- c("X")
XpostSample <- dc.fit(datalist, params=out.parms, model=PoisNorm.hXgY, 
                      multiply=NULL, unchanged=c("n", "df1","beta","lsd1"), n.chains=n.chains,
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
    X.star <- XpostBinned[j,]
    GT.lnumer.persamp <- sum(dpois(x=Y, lambda=exp(ith.beta*X.star), log=TRUE) 
                             + dnorm(x=W, mean=X.star,sd=ith.sd1, log=TRUE) +
                               dt(x=X.star, df=parms.err.varsPN$df1, ncp=0, log=TRUE))  #--> this part should cancel! 
    
    GT.ldenom.persamp <- sum(dpois(x=Y, lambda=exp(ModPN.mles$beta*X.star), log=TRUE)
                             + dnorm(x=W, mean=X.star, sd=exp(ModPN.mles$lsd1), log=TRUE) +
                               dt(x=X.star, df=parms.err.varsPN$df1, ncp=0, log=TRUE))  #--> this part should cancel! 
    
    GT.logdiff[j] <- GT.lnumer.persamp-GT.ldenom.persamp
  }
  
  ratio.prof.vec.star <- exp(GT.logdiff)
  prof.vec[i] <- mean(ratio.prof.vec.star)
  # print(c(  mean(ratio.prof.vec.star) ))
  
}

plot(beta.profs, prof.vec, type="l", xlim=c(-2, 2), ylim=c(0,1))
#abline(v=-1); abline(h=1)
abline(v=ModPN.mles$beta, lty=2)




### Model using TMB
source("TestingLaplaceTMB.functions7.0.R")
source("TestingLaplaceTMB.models2.0.R")

ModPN.TMB <- run.tmb(profile.par.name="beta", search.initial=F, n.search=1000,
                     sim.data=sim.data, model="Pois.Norm", method.ci="useconfint",
                     ytol=3, init = list(par=-0.01, hidden=0.1, df1=parms.err.varsPN$df1))

ModPN.TMB$model.obj$par
check.model <- sdreport(ModPN.TMB$fun.obj)
check.model$pdHess
ModPN.TMB$convergence
sqrt(diag(check.model$cov.fixed))
check.model
sqrt(diag(solve(ModPN.TMB$model.obj$hessian)))

plot(ModPN.TMB$ci.obj$data$beta, ModPN.TMB$ci.obj$data$rel.like, type="l", xlim=c(-2,2))
abline(v=parms.err.varsPN$beta)

# manual computing sd
import::from("numDeriv", "hessian")
(std <- sqrt(diag(solve(numDeriv::hessian(ModPN.TMB$fun.obj$fn, ModPN.TMB$fun.obj$env$last.par.best[1:2])))))
(std <- sqrt(diag(solve(numDeriv::jacobian(ModPN.TMB$fun.obj$gr,ModPN.TMB$fun.obj$env$last.par.best[1:2])))))
par.fixed <-ModPN.TMB$fun.obj$env$last.par.best[1:2]
hessian.fixed <- optimHess(par.fixed, ModPN.TMB$fun.obj$fn,
                           ModPN.TMB$fun.obj$gr) ## Marginal precision of theta.
chol(hessian.fixed) # pdHess
ModPN.TMB$model.obj$hessian

# ModPN.TMB$model.obj$par
# lapply( parms.err.varsPN[3], log)
# plot(as.data.frame(ModPN.TMB$ci.obj$data[,c(1,3)]), type="l")
#lines(ModPN.TMB$ci.obj$data$beta, ModPN.TMB$ci.obj$data$rel.like, type="l", lty=2);




### Model for DD
source("TestingLaplaceDD.functions2.0.R")
source("TestingLaplaceDD.models2.0.R")

MCMCset=list('n.iter.dc'=n.iter,'n.chains.dc'=3,
             'n.iter.dd'=n.iter,'n.chains.dd'=3,
             'n.burn'=n.update,'n.thin'=thin,'n.adapt'=n.adapt)

test.DD.multi.par <- DD_wrapper( profile.par = "beta", model = "Pois.Norm", 
                                 known.params = list(df1=parms.err.varsPN$df1),
                                 df=sim.data, psi.grid = seq(-3, 3, 0.01), MCMCset=MCMCset, 
                                 CIalpha=0.95, n.clones=c(16))

### plot lines for alg 1 and 2
# lines(test.DD.multi.par$DD.profile.q$Profile.out$psi.grid, ### ALG 1
#       exp(test.DD.multi.par$DD.profile.q$Profile.out$profile.like), 
#       lty=4, type = "l", lwd=3, col="darkgreen")
lines(test.DD.multi.par$DD.profile.psi$Profile.out$psi.grid, ### ALG 2
      exp(test.DD.multi.par$DD.profile.psi$Profile.out$profile.like), 
      lty=3, type = "l", lwd=3, col="darkred")


points(x=confints.mls[1,],y=c(0,0), col="blue", type="l" ,lty=1,lwd=4)



# Saving objects for joint plot
MP.PNtmb      <- data.frame(betavals=ModPN.TMB$ci.obj$data$beta, tmbrel.like = 
                              ModPN.TMB$ci.obj$data$rel.like)
MP.PNgt       <- data.frame(beta.profs=beta.profs, prof.vec=prof.vec)

MP.PNdd <- data.frame(DDxs = test.DD.multi.par$DD.profile.q$Profile.out$psi.grid,
                      DDys = exp(test.DD.multi.par$DD.profile.q$Profile.out$profile.like))

MP.PNDCwaldCI <- data.frame(lower=confints.mls[,1], mles=unlist(ModPN.mles), upper=confints.mls[,2])



### Plot
par(oma=c(1,1,1,1), mar=c(4,5,2,2))
plot(MP.PNgt$beta.profs, MP.PNgt$prof.vec, type="l", ylim=c(0,1),
     xlab=bquote("Values of" ~ beta), ylab="Profile Likelihood", bty="l",
     cex.lab=1.5, xlim=c(-2, 2), lty=3, lwd=1.5)
abline(v=ModPN.mles$beta, lty=2, lwd=1)
#abline(v=parms.err.varsNN$beta, lty=4)

### plot lines for TMB
lines(MP.PNtmb$betavals, MP.PNtmb$tmbrel.like, 
      type="l",col="black", lty=2, lwd=1.5);


### plot lines for alg 1 and 2
lines(MP.PNdd$DDxs,MP.PNdd$DDys,lwd=1.5,lty=1,col="black",type="l")

### Wald-Type CI and legend
#abline(h=exp(-qchisq(p=0.95,df=1)/2),lty=2, col="gray")

arrows(x0=confints.mls[1,1],x1=confints.mls[1,2],y0=0,y1=0, angle=90, code=3, length=0.06, lwd=2)

legend("topleft", 
       legend=c("DD-2","TMB", "GT"), 
       col=c("black", "black","black"),
       lwd=rep(1.5,3), 
       lty=c(1:3), bty="n", cex=0.9)


#save.image("PNbetaprof.RData")

