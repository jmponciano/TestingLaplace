rm(list=ls())

##### --------------------------------------------------------------------------
source("TestingLaplaceToolkit4.0.R")
source("TestingLaplaceTMB.functions7.0.R") 
source("TestingLaplaceTMB.models2.0.R") 


##### --------------------------------------------------------------------------
#### Parallel libraries
##### --------------------------------------------------------------------------
library(magrittr)
library(foreach)
library(doParallel)
library(ggplot2); library(patchwork)



##### --------------------------------------------------------------------------
#### kde1d objective function
##### --------------------------------------------------------------------------
library("kde1d")
obj.ftn <- function(xval, Xspost){
  
  Xs.smooth <- kde1d::kde1d(x=Xspost)
  out <- log(kde1d::dkde1d(x=xval,obj=Xs.smooth))
  return(out)
}



##### --------------------------------------------------------------------------
#### MCMC settings
##### --------------------------------------------------------------------------
n.chains <- 5
n.adapt  <- 1000
n.update <- 10000
n.iter   <- 50000
thin     <- 5


##### --------------------------------------------------------------------------
#### Bivar.Norm
##### --------------------------------------------------------------------------


##### --------------------------------------------------------------------------
##### 
# set.seed(156)
theta <- list(beta=-1, sd1=3, sd2=2, df1=2)
out.simNN <- sim.example(n=100, parm.list=theta, model="Bivar.Norm")
orig.data <- as.data.frame( out.simNN )
n <- nrow(orig.data)
W <- orig.data$W
Y <- orig.data$Y

datalist <- list(n=n, df1=theta$df1, Y=Y, W=W,
                 sd1=theta$sd1, sd2=theta$sd2, beta=theta$beta)



##### --------------------------------------------------------------------------
##### 
dcrun <- jags.fit(datalist, c("X[1]"), BivarNorm.Upost, n.chains = n.chains, 
                  n.adapt=n.adapt, n.update=n.update, 
                  n.iter = n.iter, thin=thin)


# do a loop over several theta values
n.thetas  <- 500
theta.mat <- matrix( data=cbind( beta = runif( n.thetas, -3, 0 ),
                                 sd1  = rep(3, n.thetas),
                                 sd2  = rep(2, n.thetas)), 
                     ncol=3, nrow = n.thetas )
# rep(2, n.thetas)

cl <- makeCluster(11)
registerDoParallel(cl)
res <- foreach(j=1:n.thetas, .combine=rbind) %dopar%{
  theta.vec <- numeric(0)
  
  datalist <- list(n=n, df1=theta$df1, Y=Y, W=W,
                   beta=theta.mat[j,1], sd1=theta.mat[j,2], 
                   sd2=theta.mat[j,3])
  
  dcrun <- dclone::jags.fit(datalist, c("X[1]"), 
                            BivarNorm.Upost, n.chains = n.chains, 
                            n.adapt=n.adapt, n.update=n.update, 
                            n.iter = n.iter, thin=thin)
  
  Xs.post <- do.call(rbind,dcrun)
  x4trial <- mean(Xs.post)
  
  mu.opt <- optim(par=x4trial, fn=obj.ftn, method="CG", 
                  hessian=TRUE, Xspost=Xs.post, control=list(fnscale=-1))
  
  theta.vec <- rbind(theta.vec, mu.opt$value - 0.5 * log(-mu.opt$hessian))
  return(theta.vec)
}
stopCluster(cl = cl)

##### --------------------------------------------------------------------------
out.simNN <- orig.data
parms.err.varsNN <- theta


##### --------------------------------------------------------------------------
##### 
MCprof.ModNN <- MCprofile1d(B=50000, parm.list=parms.err.varsNN, sim.data=out.simNN,
                            bracket=6,prec.prof=0.01,model="Bivar.Norm",
                            plot.it=FALSE,overlay=FALSE)
# write.csv(MCprof.ModNN, "MCprof.ModNN.csv")




##### --------------------------------------------------------------------------
# Now computing the TMB profile:
out.simNN <- orig.data
parms.err.varsNN <- theta
parms4tmb <- list(beta=parms.err.varsNN$beta, df1=parms.err.varsNN$df1,
                  lnsd1 = log(parms.err.varsNN$sd1), lnsd2=log(parms.err.varsNN$sd2))
NNtmbMod <- run.tmb(profile.par.name="beta", profile.par=MCprof.ModNN[,1], 
                    known.parm=parms4tmb, sim.data=out.simNN, model="Bivar.Norm")
NNtmbMod$MLres$res
NNtmbMod$MLres$converged


# Plot the MC profile with the TMB profile
theme_set(theme_classic())
theme_update(text = element_text(size=15),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             strip.background = element_blank(),
             axis.ticks.length=unit(0.2, "cm"), 
             legend.key.width = unit(4, "line"),
             legend.position = "top"
)

jth.par <- 1
xrange  <- c(-3,0.05)# range(theta.mat[,1])
NNplot1 <- ggplot(data.frame(MCprof.ModNN), aes(Beta.values, Rel.Prof.Like)) + 
  geom_line(aes( linetype="mc.line"), col="black", lwd=1) + 
  geom_line(data = NNtmbMod$prof, aes( x=param, y=rel.like, linetype="tmb.line"), col="black", lwd=1)+
  labs( x=bquote("Values of" ~ beta), y="Profile Likelihood", title = "A" ) + lims(x=xrange) +
  scale_linetype_manual(name = "", values = c(tmb.line = 4, mc.line = 3),
                        labels = c(tmb.line = "TMB", mc.line = "MC"),
                        limits = c("tmb.line", "mc.line"),
                        guide = guide_legend(override.aes = list(color = c("black", "black"), 
                                                                 linetype = c(4,3), linewidth = c(1, 1) )) ) 
NNplot2 <- ggplot(data.frame(MCprof.ModNN), 
                  aes(Beta.values, NNtmbMod$prof$rel.like/Rel.Prof.Like)) + 
  geom_line(lwd=1, lty=1) + xlim( xrange ) + ylim( c(0,40) ) + 
  labs(  x=bquote("Values of" ~ beta), y="Profile Likelihood Ratio") +
  theme( axis.title.y = element_text(size=12) )
NNplot3 <- ggplot(data.frame(theta=theta.mat[,jth.par], res=exp(res)), 
                  aes(theta, res)) + 
  geom_point(size=2, pch=16, alpha=0.2) + 
  geom_smooth(method = "loess", se=F, lty=1, lwd=1, col="black") +
  xlim( xrange ) +
  labs( x=bquote("Values of" ~ beta), y=bquote("f("~ mu ~ "| y) x det(H)"^-0.5))
NNplot <- NNplot1/(NNplot2|NNplot3) + plot_layout(guides = "collect")
NNplot




### Save files for plots
OP.NNmc  <- data.frame(MCprof.ModNN)
OP.NNtmb <- NNtmbMod
OP.NNc   <- list(theta.mat=theta.mat, res=res, expres=exp(res))


##### --------------------------------------------------------------------------
# ggsave("UhatCalcNNPlot.png", width = 12, height = 10)
# save(parms.err.varsNN, out.simNN, OP.NNmc, OP.NNtmb, OP.NNc, 
#      file="UhatCalcNN-OnePar_plotdata.RData")



##### --------------------------------------------------------------------------
#### Pois.Norm
##### --------------------------------------------------------------------------

##### --------------------------------------------------------------------------
# set.seed(135)
theta <- list(beta=-1, df1=2, sd1=3 )
orig.data <- sim.example(n=100, parm.list=theta, model="Pois.Norm")
orig.data <- as.data.frame( orig.data )

n <- nrow(orig.data)
W <- orig.data[,1]
Y <- orig.data[,2]

datalist <- list(n=n, df1=theta$df1, Y=Y, W=W, sd1=theta$sd1, beta=theta$beta)



##### --------------------------------------------------------------------------
# do a loop over several theta values
n.thetas  <- n.thetas
theta.mat <- matrix( data=cbind( beta = runif( n.thetas, -1.5, -0.2 ),
                                 sd1  = rep(theta$sd1, n.thetas)), 
                     ncol=2, nrow = n.thetas )
# rep(2, n.thetas)

cl <- makeCluster(11)
registerDoParallel(cl)
res <- foreach(j=1:n.thetas, .combine=rbind) %dopar%{
  theta.vec <- numeric(0)
  
  datalist <- list(n=n, df1=theta$df1, Y=Y, W=W,
                   beta=theta.mat[j,1], sd1=theta.mat[j,2])
  
  dcrun <- dclone::jags.fit(datalist, c("X[1]"), 
                            PoisNorm.Upost, n.chains = n.chains, 
                            n.adapt=n.adapt, n.update=n.update, 
                            n.iter = n.iter, thin=thin)
  
  Xs.post <- do.call(rbind,dcrun)
  x4trial <- mean(Xs.post)
  
  mu.opt <- optim(par=x4trial, fn=obj.ftn, method="CG", 
                  hessian=TRUE, Xspost=Xs.post, control=list(fnscale=-1))
  
  theta.vec <- rbind(theta.vec, mu.opt$value - 0.5 * log(-mu.opt$hessian))
  return(theta.vec)
}
stopCluster(cl = cl)

out.simPN <- orig.data
parms.err.varsPN <- theta



##### --------------------------------------------------------------------------
##### 
MCprof.ModPN <- MCprofile1d(B=50000, parm.list=parms.err.varsPN, sim.data=out.simPN,
                            bracket=2,prec.prof=0.01, model="Pois.Norm",
                            plot.it=FALSE,overlay=FALSE)
# write.csv(MCprof.ModPN, "MCprof.ModPN.csv")




##### --------------------------------------------------------------------------
# Now computing the TMB profile:
parms4tmb <- list(beta=parms.err.varsPN$beta, df1=parms.err.varsPN$df1,
                  lnsd1 = log(parms.err.varsPN$sd1))
PNtmbMod <- run.tmb(profile.par.name="beta", profile.par=MCprof.ModPN[,1], 
                    known.parm=parms4tmb, sim.data=out.simPN, model="Pois.Norm")
PNtmbMod$MLres$res
PNtmbMod$MLres$converged




##### --------------------------------------------------------------------------
# Plot the MC profile with the TMB profile
jth.par <- 1
xrange  <- range(theta.mat[,1]) # c(-2,-0.5)
PNplot1 <- ggplot(data.frame(MCprof.ModPN), aes(Beta.values, Rel.Prof.Like)) + 
  geom_line(aes( linetype="mc.line"), col="black", lwd=1) + 
  geom_line(data = PNtmbMod$prof, aes( x=param, y=rel.like, linetype="tmb.line"), col="black", lwd=1)+
  labs( x=bquote("Values of" ~ beta), y="Profile Likelihood", title = "B" ) + lims(x=xrange) +
  scale_linetype_manual(name = "", values = c(tmb.line = 4, mc.line = 3),
                        labels = c(tmb.line = "TMB", mc.line = "MC"),
                        limits = c("tmb.line", "mc.line"),
                        guide = guide_legend(override.aes = list(color = c("black", "black"), 
                                                                 linetype = c(4,3), linewidth = c(1, 1) )) ) 
PNplot2 <- ggplot(data.frame(MCprof.ModPN), 
                  aes(Beta.values, PNtmbMod$prof$rel.like/Rel.Prof.Like)) + 
  geom_line(lwd=1, lty=1) + xlim( c(-1.4, -0.68) ) + 
  labs(  x=bquote("Values of" ~ beta), y="Profile Likelihood Ratio") +
  theme( axis.title.y = element_text(size=12) )
PNplot3 <- ggplot(data.frame(theta=theta.mat[,jth.par], res=exp(res)), 
                  aes(theta, res)) + 
  geom_point(size=2, pch=16, alpha=0.2) + 
  geom_smooth(method = "loess", se=F, lty=1, lwd=1, col="black") +
  xlim( c(-1.4, -0.68) ) +
  labs( x=bquote("Values of" ~ beta), y=bquote("f("~ mu ~ "| y) x det(H)"^-0.5))
PNplot <- PNplot1/(PNplot2|PNplot3) + plot_layout(guides = "collect")
PNplot





##### --------------------------------------------------------------------------
### Save files for plots
OP.PNmc  <- data.frame(MCprof.ModPN)
OP.PNtmb <- PNtmbMod
OP.PNc   <- list(theta.mat=theta.mat, res=res, expres=exp(res))
# save(parms.err.varsPN, out.simPN, OP.PNmc, OP.PNtmb, OP.PNc, 
#      file="UhatCalcPN-OnePar_plotdata.RData")
# ggsave("UhatCalcPNPlot.png", width = 12, height = 10)









##### --------------------------------------------------------------------------
#### Logit.Norm
##### --------------------------------------------------------------------------


##### --------------------------------------------------------------------------
# set.seed(3)
theta <- list(beta=-1, df1=2, sd1=3 )
orig.data <- sim.example(n=50, parm.list=theta, model="Logit.Norm")
n <- nrow(orig.data)
W <- orig.data[,1]
Y <- orig.data[,2]

datalist <- list(n=n, df1=theta$df1, Y=Y, W=W, sd1=theta$sd1, beta=theta$beta)




##### --------------------------------------------------------------------------
# do a loop over several theta values
n.thetas  <- n.thetas
theta.mat <- matrix( data=cbind( beta = runif( n.thetas, -3, 3 ),
                                 sd1  = rep(theta$sd1, n.thetas)), 
                     ncol=2, nrow = n.thetas )
# rep(2, n.thetas)

cl <- makeCluster(11)
registerDoParallel(cl)
res <- foreach(j=1:n.thetas, .combine=rbind) %dopar%{
  theta.vec <- numeric(0)
  
  datalist <- list(n=n, df1=theta$df1, Y=Y, W=W,
                   beta=theta.mat[j,1], sd1=theta.mat[j,2])
  
  dcrun <- dclone::jags.fit(datalist, c("X[1]"), 
                            LogitNorm.Upost, n.chains = n.chains, 
                            n.adapt=n.adapt, n.update=n.update, 
                            n.iter = n.iter, thin=thin)
  
  Xs.post <- do.call(rbind,dcrun)
  x4trial <- mean(Xs.post)
  
  mu.opt <- optim(par=x4trial, fn=obj.ftn, method="CG", 
                  hessian=TRUE, Xspost=Xs.post, control=list(fnscale=-1))
  
  theta.vec <- rbind(theta.vec, mu.opt$value - 0.5 * log(-mu.opt$hessian))
  return(theta.vec)
}
stopCluster(cl = cl)

out.simLN <- orig.data
parms.err.varsLN <- theta



##### --------------------------------------------------------------------------
MCprof.ModLN <- MCprofile1d(B=50000, parm.list=parms.err.varsLN, sim.data=out.simLN,
                            bracket=4,prec.prof=0.01,model="Logit.Norm",
                            plot.it=FALSE,overlay=FALSE)
# write.csv(MCprof.ModPN, "MCprof.ModLN.csv")




##### --------------------------------------------------------------------------
# Now computing the TMB profile:
parms4tmb <- list(beta=parms.err.varsLN$beta, df1=parms.err.varsLN$df1,
                  lnsd1 = log(parms.err.varsLN$sd1))
LNtmbMod <- run.tmb(profile.par.name="beta", profile.par=MCprof.ModLN[,1], 
                    known.parm=parms4tmb, sim.data=out.simLN, model="Logit.Norm")

LNtmbMod$MLres$res
LNtmbMod$MLres$converged



##### --------------------------------------------------------------------------
# Plot the MC profile with the TMB profile
jth.par <- 1
xrange  <- range(theta.mat[,1]) # c(-4,1)
LNplot1 <- ggplot(data.frame(MCprof.ModLN), aes(Beta.values, Rel.Prof.Like)) + 
  geom_line(aes( linetype="mc.line"), col="black", lwd=1) + 
  geom_line(data = LNtmbMod$prof, aes( x=param, y=rel.like, linetype="tmb.line"), col="black", lwd=1)+
  labs( x=bquote("Values of" ~ beta), y="Profile Likelihood", title = "C" ) + lims(x=xrange) +
  scale_linetype_manual(name = "", values = c(tmb.line = 4, mc.line = 3),
                        labels = c(tmb.line = "TMB", mc.line = "MC"),
                        limits = c("tmb.line", "mc.line"),
                        guide = guide_legend(override.aes = list(color = c("black", "black"), 
                                                                 linetype = c(4,3), linewidth = c(1, 1) )) ) 
LNplot2 <- ggplot(data.frame(MCprof.ModLN), 
                  aes(Beta.values, LNtmbMod$prof$rel.like/Rel.Prof.Like)) + 
  geom_line(lwd=1, lty=1) + xlim( xrange ) + ylim( c(0, 3) ) + 
  labs(  x=bquote("Values of" ~ beta), y="Profile Likelihood Ratio") +
  theme( axis.title.y = element_text(size=12) )
LNplot3 <- ggplot(data.frame(theta=theta.mat[,jth.par], res=exp(res)), 
                  aes(theta, res)) + 
  geom_point(size=3, pch=16, alpha=0.2) + 
  geom_smooth(method = "loess", se=F, lty=1, lwd=1, col="black") +
  xlim( xrange ) +
  labs( x=bquote("Values of" ~ beta), y=bquote("f("~ mu ~ "| y) x det(H)"^-0.5))
LNplot <- LNplot1/(LNplot2|LNplot3) + plot_layout(guides = "collect")
LNplot






##### --------------------------------------------------------------------------
### Save files for plots
OP.LNmc  <- data.frame(MCprof.ModLN)
OP.LNtmb <- LNtmbMod
OP.LNc   <- list(theta.mat=theta.mat, res=res, expres=exp(res))
# save(parms.err.varsLN, out.simLN, OP.LNmc, OP.LNtmb, OP.LNc, 
#      file="UhatCalcLN-OnePar_plotdata.RData")
# ggsave("UhatCalcLNPlot.png", width = 12, height = 10)





#### -----------------------------------------------------------------
#### Probit.Norm
#### -----------------------------------------------------------------



##### --------------------------------------------------------------------------
# set.seed(32)
theta <- list(beta=-1, df1=2, sd1=3)
orig.data <- sim.example(n=50, parm.list=theta, model="Probit.Norm")
n <- nrow(orig.data)
W <- orig.data[,1]
Y <- orig.data[,2]

datalist <- list(n=n, df1=theta$df1, Y=Y, W=W, sd1=theta$sd1, beta=theta$beta)





##### --------------------------------------------------------------------------
# do a loop over several theta values
n.thetas  <- n.thetas
theta.mat <- matrix( data=cbind( beta = runif( n.thetas, -3, 1.5 ),
                                 sd1  = rep(theta$sd1, n.thetas)), 
                     ncol=2, nrow = n.thetas )
# rep(2, n.thetas)

cl <- makeCluster(11)
registerDoParallel(cl)
res <- foreach(j=1:n.thetas, .combine=rbind, .export = ls(globalenv())) %dopar%{
  theta.vec <- numeric(0)
  
  datalist <- list(n=n, df1=theta$df1, Y=Y, W=W,
                   beta=theta.mat[j,1], sd1=theta.mat[j,2])
  
  dcrun <- dclone::jags.fit(datalist, c("X[1]"), 
                            ProbitNorm.Upost, n.chains = n.chains, 
                            n.adapt=n.adapt, n.update=n.update, 
                            n.iter = n.iter, thin=thin)
  
  Xs.post <- do.call(rbind,dcrun)
  x4trial <- mean(Xs.post)
  
  mu.opt <- optim(par=x4trial, fn=obj.ftn, method="CG", 
                  hessian=TRUE, Xspost=Xs.post, control=list(fnscale=-1))
  
  theta.vec <- rbind(theta.vec, mu.opt$value - 0.5 * log(-mu.opt$hessian))
  return(theta.vec)
}
stopCluster(cl = cl)


out.simProbitN <- orig.data
parms.err.varsProbitN <- theta



##### --------------------------------------------------------------------------
MCprof.ModProbitN <- MCprofile1d(B=50000, parm.list=parms.err.varsProbitN,
                                 sim.data=out.simProbitN,
                                 bracket=4,prec.prof=0.01,model="Probit.Norm",
                                 plot.it=FALSE,overlay=FALSE)






##### --------------------------------------------------------------------------
# Now computing the TMB profile:
parms4tmb <- list(beta=parms.err.varsProbitN$beta, 
                  df1=parms.err.varsProbitN$df1,
                  lnsd1 = log(parms.err.varsProbitN$sd1))
ProbitNtmbMod <- run.tmb(profile.par.name="beta", profile.par=MCprof.ModProbitN[,1], 
                         known.parm=parms4tmb, sim.data=out.simProbitN, model="Probit.Norm")
ProbitNtmbMod$MLres$res
ProbitNtmbMod$MLres$converged





##### --------------------------------------------------------------------------
# Plot the MC profile with the TMB profile
jth.par <- 1
xrange  <- range(theta.mat[,1]) # c(-4,1)
ProbitNplot1 <- ggplot(data.frame(MCprof.ModProbitN), aes(Beta.values, Rel.Prof.Like)) + 
  geom_line(aes( linetype="mc.line"), col="black", lwd=1) + 
  geom_line(data = ProbitNtmbMod$prof, aes( x=param, y=rel.like, linetype="tmb.line"), col="black", lwd=1)+
  labs( x=bquote("Values of" ~ beta), y="Profile Likelihood", title = "D" ) + lims(x=xrange) +
  scale_linetype_manual(name = "", values = c(tmb.line = 4, mc.line = 3),
                        labels = c(tmb.line = "TMB", mc.line = "MC"),
                        limits = c("tmb.line", "mc.line"),
                        guide = guide_legend(override.aes = list(color = c("black", "black"), 
                                                                 linetype = c(4,3), linewidth = c(1, 1) )) ) 
ProbitNplot2 <- ggplot(data.frame(MCprof.ModProbitN), 
                  aes(Beta.values, ProbitNtmbMod$prof$rel.like/Rel.Prof.Like)) + 
  geom_line(lwd=1, lty=1) + xlim( xrange ) + ylim(c(0,1.3)) +
  labs(  x=bquote("Values of" ~ beta), y="Profile Likelihood Ratio") +
  theme( axis.title.y = element_text(size=12) )
ProbitNplot3 <- ggplot(data.frame(theta=theta.mat[,jth.par], res=exp(res)), 
                  aes(theta, res)) + 
  geom_point(size=3, pch=16, alpha=0.2) + 
  geom_smooth(method = "loess", se=F, lty=1, lwd=1, col="black") +
  xlim( xrange ) +
  labs( x=bquote("Values of" ~ beta), y=bquote("f("~ mu ~ "| y) x det(H)"^-0.5))
ProbitNplot <- ProbitNplot1/(ProbitNplot2|ProbitNplot3) + plot_layout(guides = "collect")
ProbitNplot






##### --------------------------------------------------------------------------
### Save files for plots

OP.ProbitNmc  <- data.frame(MCprof.ModProbitN)
OP.ProbitNtmb <- ProbitNtmbMod
OP.ProbitNc   <- list(theta.mat=theta.mat, res=res, expres=exp(res))
# ggsave("UhatCalcProbitNPlot.png", width = 12, height = 10)
# save(parms.err.varsProbitN, out.simProbitN, OP.ProbitNmc, OP.ProbitNtmb, OP.ProbitNc, 
#      file="UhatCalcProbitN-OnePar_plotdata.RData")

##### --------------------------------------------------------------------------




### final plot
(NNplot|PNplot)/(LNplot|ProbitNplot) #+ plot_layout(widths = c(4, -1.1 ,4.5))
# ggsave("figure4.png", width = 13, height = 11)
 









#### ---------------------------------------------------------------------------
dyn.unload(dynlib(model_name))

if (length(list.files(pattern = "\\.o$", full.names = TRUE))>0) {
  file.remove(list.files(pattern = "\\.o$", full.names = TRUE))
}
if (length(list.files(pattern = "\\.so$", full.names = TRUE))>0) {
  file.remove(list.files(pattern = "\\.so$", full.names = TRUE))
}
if (length(list.files(pattern = "\\.cpp$", full.names = TRUE))>0) {
  file.remove(list.files(pattern = "\\.cpp$", full.names = TRUE))
}
#



