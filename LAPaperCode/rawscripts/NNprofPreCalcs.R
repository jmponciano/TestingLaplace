source("JAGS4Profiles.R")
source("simulate.data.R")
library("doParallel")

# Simulating from Model 1: Latent variable is X ~ T_2(0), 
# a reasonably flat distribution 
# W = X + eta where eta ~ N(0,sigma1)
# Y = beta*X + epsilon where epsilon ~ N(0,sigma2)
# Response = c(W,Y), a bivariate Normal given X=x
parms.err.varsNN <- list(beta=-1, df1=2, sd1=3,sd2=2)
out.simNN <- sim.example(n=50, parm.list=parms.err.varsNN, model="Bivar.Norm")
sim.data <- out.simNN

n.iter <- 10000
n.chains <- 5
n.adapt <- 1000
n.update <- 1000
thin <- 2


n <- nrow(sim.data)
W <- sim.data[,1]
Y <- sim.data[,2]


cl <- makePSOCKcluster(n.chains)
datalist <- list(K=1, n=n, df1=parms.err.varsNN$df1, Y=dcdim(data.matrix(Y)), W=dcdim(data.matrix(W)))
dcrunNN <- dc.parfit(cl=cl, datalist, c("beta", "lsd1","lsd2"), BivarNormDC.MLE, n.clones=c(1,2,4,8,16),
                multiply="K", unchanged=c("n","df1"),
                n.chains = n.chains, n.adapt=n.adapt, n.update=n.update, 
                n.iter = n.iter, thin=thin)


ModNN.mles <- as.list(summary(dcrunNN)[[1]][,1])
stopCluster(cl)
#ModNN.mles
#dcdiag(dcrunNN)





confints.mls <- confint(dcrunNN)

# MLES to profile over beta

beta.mle <- ModNN.mles$beta
lsd1.mle <- ModNN.mles$lsd1
lsd2.mle <- ModNN.mles$lsd2


bracket <- 1.5
prec.prof <- 0.01
beta.profs <- seq(from=beta.mle-bracket,to=beta.mle+bracket, by=prec.prof)
len.prof <- length(beta.profs)
len.prof

cl <- makeCluster(8,outfile="")
registerDoParallel(cl)
res <- foreach(i=1:len.prof, .combine=rbind)%dopar%{

  outlist.lsd1 <- list()
  outlist.lsd2 <- list()
  ith.beta <- beta.profs[i]
  #cl <- makePSOCKcluster(n.chains)
  ith.datalist <- list(K=1, n=n, beta=ith.beta, df1=parms.err.varsNN$df1, Y=dclone::dcdim(data.matrix(Y)), 
                       W=dclone::dcdim(data.matrix(W)))
  clones.seq <- c(16)
  
  ith.dcrun <- dclone::dc.fit(ith.datalist, c("lsd1","lsd2"),BivarNorm.betaprof, n.clones=clones.seq,
                  multiply="K", unchanged=c("n","beta", "df1"),
                  n.chains = n.chains, n.adapt=n.adapt, n.update=n.update, 
                  n.iter = n.iter, thin=thin)
  profhatsNN <- as.list(summary(ith.dcrun)[[1]][,1:2])
  #stopCluster(cl)
  outlist.lsd1 <- rbind(outlist.lsd1, profhatsNN[1])
  outlist.lsd2 <- rbind(outlist.lsd2, profhatsNN[2])
  return(list(lsd1=outlist.lsd1,lsd2=outlist.lsd2))
}
stopCluster(cl)

# This is the vector of mles of lsd1 and lsd2 for every fixed value of beta
lsdhats4prof <- cbind(lsd1=unlist(res[,1]),lsd2=unlist(res[,2])) 


# Now prepare for the GT algorithm by sampling from g(x|y)
# Sampling from f(X|(W,Y)) 
n.iter <- 100000
n.chains <- 4
n.adapt <- 1000
n.update <- 1000
thin <- 2
datalist <- list(W=W, Y=Y, beta=ModNN.mles$beta, df1=parms.err.varsNN$df1,lsd1=ModNN.mles$lsd1,lsd2 =ModNN.mles$lsd2, n=n)
out.parms <- c("X")
XpostSample <- dc.fit(datalist, params=out.parms, model=BivarNorm.hXgY,
                 multiply=NULL, unchanged=c("n", "df1", "beta","lsd1","lsd2"), n.chains=n.chains,
                 n.adapt=n.adapt, n.update=n.update, 
                 n.iter = n.iter, thin=thin,n.clones=1) 
XpostBinned <- do.call(rbind, XpostSample)


# Now the GT algorithm itself:

B <- nrow(XpostBinned)
prof.vec <- rep(0,len.prof)


for(i in 1:len.prof){
  
  ith.beta <- beta.profs[i]
  ith.sd1 <- exp(lsdhats4prof[i,1])
  ith.sd2 <- exp(lsdhats4prof[i,2])  

  GT.logdiff <- rep(0,B)
  for(j in 1:B){
    X.star <- XpostBinned[j,]
    GT.lnumer.persamp<- sum(dnorm(x=Y, mean=ith.beta*X.star, sd=ith.sd2, log=TRUE) 
    + dnorm(x=W, mean=X.star,sd=ith.sd1, log=TRUE) +
    dt(x=X.star, df=parms.err.varsNN$df1, ncp=0,log=TRUE))  #--> this part should cancel! 
    
    GT.ldenom.persamp <- sum(dnorm(x=Y, mean=ModNN.mles$beta*X.star, sd=exp(ModNN.mles$lsd2), log=TRUE)
    + dnorm(x=W, mean=X.star,sd=exp(ModNN.mles$lsd1), log=TRUE) +
      dt(x=X.star, df=parms.err.varsNN$df1, ncp=0,log=TRUE))  #--> this part should cancel! 

    GT.logdiff[j] <- GT.lnumer.persamp-GT.ldenom.persamp
  }
  
  ratio.prof.vec.star <- exp(GT.logdiff)
  prof.vec[i] <- mean(ratio.prof.vec.star)
  
}


### Model using TMB
source("TestingLaplaceTMB.functions7.0.R")
source("TestingLaplaceTMB.models2.0.R")

ModNN.TMB <- run.tmb(profile.par.name="beta", search.initial=F, n.search=5000,
                     sim.data=sim.data, model="Bivar.Norm", method.ci="useconfint",
                     ytol=3, init = list(par=-0.1, hidden=0.1, df1=parms.err.varsNN$df1))
ModNN.TMB$model.obj$par
ModNN.TMB$convergence
ModNN.TMB$results

check.model <- sdreport(ModNN.TMB$fun.obj)
check.model$pdHess
ModNN.TMB$convergence
sqrt(diag(check.model$cov.fixed))
check.model
sqrt(diag(solve(ModNN.TMB$model.obj$hessian)))


# manual computing sd
import::from("numDeriv", "hessian")
(std <- sqrt(diag(solve(numDeriv::hessian(ModNN.TMB$fun.obj$fn, ModNN.TMB$fun.obj$env$last.par.best[1:3])))))
(std <- sqrt(diag(solve(numDeriv::jacobian(ModNN.TMB$fun.obj$gr,ModNN.TMB$fun.obj$env$last.par.best[1:3])))))
par.fixed <-ModNN.TMB$fun.obj$env$last.par.best[1:3]
hessian.fixed <- optimHess(par.fixed, ModNN.TMB$fun.obj$fn,
                           ModNN.TMB$fun.obj$gr) ## Marginal precision of theta.
chol(hessian.fixed) # pdHess
# chol(ModNN.TMB$model.obj$hessian)



### NOTE: need to include better diagnostics into TMB model. Maybe use TMB
if(!check.model$pdHess) {
  ## double-check (slower, more accurate hessian)
  env <- environment(ModNN.TMB$fun.obj$fn)
  par <- env$last.par.best
  if (!is.null(rr <- env$random)) {
    par <- par[-rr]
  }
  h <- numDeriv::jacobian(ModNN.TMB$fun.obj$gr, par)
  
  ## fall back to solve(optimHess(par, obj$fn, obj$gr)) ? 
  h <- .5 * (h + t(h))  ## symmetrize
  if (!any(is.na(h))) {
    ev <- try( e_complex_check(eigen(h)$values), silent = TRUE)
    if (!inherits(ev, "try-error") && min(ev)>.Machine$double.eps) {
      ## apparently fit is OK after all ...
      sdr$pdHess <- TRUE
      Vtheta <- try(solve(h), silent=TRUE)
      if (!inherits(Vtheta,"try-error")) {
        check.model$cov.fixed[] <- Vtheta
      } else {
        warning("failed to invert Hessian from numDeriv::jacobian(), falling back to internal vcov estimate")
      }
    } ## eig check OK
  }
} ## !any(is.na(h))
h <- numDeriv::jacobian(ModNN.TMB$fun.obj$gr, ModNN.TMB$fun.obj$env$last.par.best[1:3])
sqrt(diag(solve(h)))

lapply( parms.err.varsNN[3:4], log)
plot(as.data.frame(ModNN.TMB$ci.obj$data[,c(1,3)]), type="l")


### Model for DD
source("TestingLaplaceDD.functions2.0.R")
source("TestingLaplaceDD.models2.0.R")

MCMCset=list('n.iter.dc'=50000,'n.chains.dc'=3,
             'n.iter.dd'=20000,'n.chains.dd'=15,
             'n.burn'=5000,'n.thin'=thin,'n.adapt'=n.adapt)

test.DD.multi.par <- DD_wrapper( profile.par = "beta", model = "Bivar.Norm",
                                 #known.mle = list(MLE=ModNN.mles, FI=solve(vcov(dcrunNN))),
                                 known.params = list(df1=parms.err.varsNN$df1),
                                 df=sim.data, psi.grid = ModNN.TMB$ci.obj$data$beta, 
                                 MCMCset=MCMCset, CIalpha=0.95, n.clones=c(16))



MP.NNtmb      <- data.frame(betavals=ModNN.TMB$ci.obj$data$beta, tmbrel.like = 
                              ModNN.TMB$ci.obj$data$rel.like)
MP.NNgt       <- data.frame(beta.profs=beta.profs, prof.vec=prof.vec)

MP.NNdd <- data.frame(DDxs = test.DD.multi.par$DD.profile.q$Profile.out$psi.grid,
                      DDys = exp(test.DD.multi.par$DD.profile.q$Profile.out$profile.like))

MP.NNDCwaldCI <- data.frame(lower=confints.mls[,1], mles=unlist(ModNN.mles), upper=confints.mls[,2])

### Plot
par(oma=c(1,1,1,1), mar=c(4,5,2,2))
plot(MP.NNgt$beta.profs, MP.NNgt$prof.vec, type="l", ylim=c(0,1),
     xlab=bquote("Values of" ~ beta), ylab="Profile Likelihood", bty="l",
     cex.lab=1.5, xlim=c(-3.25, 3), lty=3, lwd=1.5)
#abline(v=-1); abline(h=1)
abline(v=ModNN.mles$beta, lty=2, lwd=1)
#abline(v=parms.err.varsNN$beta, lty=4)

### plot lines for TMB
lines(MP.NNtmb$betavals, MP.NNtmb$tmbrel.like, 
      type="l",col="black", lty=2, lwd=1.5);


### plot lines for alg 1 and 2
lines(MP.NNdd$DDxs,MP.NNdd$DDys,lwd=1.5,lty=1,col="black",type="l")
#lines(test.DD.multi.par$DD.profile.psi$Profile.out$psi.grid, ### ALG 1
#      exp(test.DD.multi.par$DD.profile.psi$Profile.out$profile.like), 
#      lwd=2,lty=1,col="blue",type="l")
MP.NNdd[which(MP.NNdd$DDys==max(MP.NNdd$DDys)),]
ModNN.mles$beta

# with 5 chains
#         DDxs DDys
# 56 -1.496407    1

## MLE is -1.482772


### Wald-Type CI and legend
#abline(h=exp(-qchisq(p=0.95,df=1)/2),lty=2, col="gray")

arrows(x0=confints.mls[1,1],x1=confints.mls[1,2],y0=0,y1=0, angle=90, code=3, length=0.06, lwd=2)
#points(x=confints.mls[1,],y=c(0,0), col="black", type="l" ,lty=1,lwd=2)

legend("topleft", 
       legend=c("DD-2","TMB", "GT"), 
       col=c("black", "black","black"),
       lwd=rep(1.5,3), 
       lty=c(1:3), bty="n", cex=0.9)


# save.image("NNbetaprof.RData")

