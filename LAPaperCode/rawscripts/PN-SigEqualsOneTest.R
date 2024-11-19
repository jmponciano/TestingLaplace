
source("TestingLaplaceToolkit4.0.R")
source("TestingLaplaceTMB.functions6.0.R") 
source("TestingLaplaceTMB.models2.0.R") 
parms.err.varsPN <- list(beta=-1, df1=2, sd1=1)
out.simPN <- sim.example(n=100, parm.list=parms.err.varsPN, model="Pois.Norm")

ModPN.dc.beta <- onedcfitrun(model="Pois.Norm", true.parms=parms.err.varsPN,
                             sim.data = out.simPN,n.iter=50000,n.adapt=1000,
                             n.update=1000,thin=5, n.chains=3,clones.seq=c(32))

dcdiag(ModPN.dc.beta)

Beta.hat <- summary(ModPN.dc.beta)[[1]][1]
print(Beta.hat)
# Wald-tyep CI's 
Beta.hatPN.CI <- confint(ModPN.dc.beta)


true.parms<-parms.err.varsPN
n <- nrow(out.simPN)
df1 <- true.parms$df1
sd1 <- true.parms$sd1


W <- out.simPN[,1]
Y <- out.simPN[,2]
datalist <- list(K=1, n=n, df1=df1,sd1=sd1, Y=dcdim(data.matrix(Y)), W=dcdim(data.matrix(W)))

#cl<- makePSOCKcluster(5)
dcrun <- dc.fit(datalist, c("beta"), PoisNormDC.beta, n.clones=c(32),
                multiply="K", unchanged=c("n", "df1", "sd1"),
                n.chains = 3, n.adapt=2000, n.update=10000, 
                n.iter = 50000, thin=10, inits=list(beta=-1))
#stopCluster(cl)
dcdiag(dcrun)

Beta.hat <- summary(dcrun)[[1]][1]
print(Beta.hat)
# Wald-tyep CI's 
Beta.hatPN.CI <- confint(dcrun)
Beta.hatPN.CI




# True (assuming all other parms known) profile likelihood for beta using Monte Carlo
MCprof.ModPN <- MCprofile1d(B=300000, parm.list=parms.err.varsPN,sim.data=out.simPN, 
                            bracket=2,prec.prof=0.01, model="Pois.Norm",
                            plot.it=FALSE,overlay=FALSE)

# LOESS smoothing of MC profile
psi.grid <-  seq(Beta.hat-1.9,Beta.hat+1.8,0.01)
MCdataPN <- data.frame(beta_values=MCprof.ModPN[,1],RelPL=MCprof.ModPN[,2])
tmp <- loess(RelPL ~ beta_values,data=MCdataPN, span=0.08)
newdata.Mod <- data.frame(beta_values=psi.grid)
MCprof.smooth <- predict(tmp,newdata=newdata.Mod)
where.nas <- which(is.na(MCprof.smooth), arr.ind=TRUE)
if(length(where.nas)>0){
  MCprof4plot <- MCprof.smooth[-where.nas]/max(MCprof.smooth[-where.nas])
  psi.grid4plot <- psi.grid[-where.nas]
  #MCprof4plot <- MCprof.smooth/max(MCprof.smooth)
}else{
  MCprof4plot <- MCprof.smooth/max(MCprof.smooth)
  psi.grid4plot <- psi.grid
}

# TMB profile likelihood
parms4tmb <- list(beta=parms.err.varsPN$beta, df1=parms.err.varsPN$df1,
                  lnsd1 = log(parms.err.varsPN$sd1))
testPN <- run.tmb(profile.par.name="beta", profile.par=MCprof.ModPN[,1], 
                  known.parm=parms4tmb, sim.data=out.simPN, model="Pois.Norm")


# 	a numeric or complex vector or matrix giving the right-hand side(s) of the linear system. 
# 	If missing, b is taken to be an identity matrix and solve will return the inverse of a
solve(testPN$MLest$opt.ML$hessian)
1/testPN$MLest$opt.ML$hessian


priorvar4dc <- 3*(1/testPN$MLest$opt.ML$hessian)

1/priorvar4dc

# Model 2.a:  estimating only beta using data cloning
PoisNormDC.beta <- function(){
  
  # Prior for beta
  beta~dnorm(0,1/priorvar4dc)  
  for(k in 1:K){
    # Hidden process
    for(i in 1:n){
      X[i,k] ~ dt(0,1,df1)
    }
    # Observation process
    for(i in 1:n){
      
      W[i,k]~dnorm(X[i,k], 1/pow(sd1,2))
      Y[i,k]~dpois(exp(beta*X[i,k]))
    }    
  }
  
}







import::from("numDeriv", "hessian")
library("numDeriv")
tmb.std <- sqrt(diag(solve(hessian(testPN$fun.obj$fn, testPN$model.obj$par))))
std

# Now computing the TMB profile:

# Plot the MC profile with the TMB profile
par(oma=c(1,1,1,1), mar=c(4,5,2,2))
plot(psi.grid4plot,MCprof4plot, type="l", col="black", lwd=2, xlim=c(-1.5,-0.5),
     xlab=bquote("Values of" ~ beta), ylab="Relative Likelihood", bty="l", 
     cex.lab=1.5)
points(testPN$prof$param, testPN$prof$rel.like, type="l",col="red", lty=2, lwd=2);
points(x=Beta.hatPN.CI,y=c(0,0), col="blue", type="l" ,lty=1,lwd=4)
abline(v=Beta.hat, lwd=2,lty=2)
abline(h=exp(-qchisq(p=0.95,df=1)/2),lty=2, col="gray")
legend("topright", legend=c("TMB", "MC", "MLE","DC-CI",bquote("exp(-0.5"~ chi[0.95]^2~")")),
       col=c("red", "black","black", "blue","gray"), lty=c(2,1,2,1,2), lwd=c(2,2,2,4,1), bty="n", cex=1.15)

## saving the ggplot object:

gridmcdf <- data.frame(Beta.values=psi.grid4plot, Rel.Prof.Like=MCprof4plot)





new.par.listPN <- parms.err.varsPN
new.par.listPN$beta <- Beta.hat

# Sampling from f(X|(Y,W))
Xpost.ModPN <- Xpost.samp(model="Pois.Norm", parm.list=new.par.listPN,
                          sim.data=out.simPN, n.iter=51000,n.adapt=1000,
                          n.update=1000,thin=5,n.chains=5)

Xpost.binned <- do.call(rbind,Xpost.ModPN)

# Running the GT algorithm
ModelPN.GTprof <- GT.prof(model="Pois.Norm", parm.list=new.par.listPN,sim.data=out.simPN,
                          XpostSample=Xpost.binned, bracket=2, prec.prof=0.01,
                          plot.it=FALSE, overlay=TRUE,mycol="blue",
                          dotted=1)




Beta.hatPN.CI

### Save files for plots
OP.PNmc       <- gridmcdf
OP.PNtmb      <- testPN$prof
OP.PNgt       <- data.frame(ModelPN.GTprof)
OP.PNDCwaldCI <- c(as.vector(Beta.hatPN.CI)[1],Beta.hat,as.vector(Beta.hatPN.CI)[2])
save(parms.err.varsPN, new.par.listPN, out.simPN,OP.PNmc, OP.PNtmb, OP.PNgt, OP.PNgt, OP.PNDCwaldCI, 
     file="PN-SigOne4plot.RData")
save.image( file="PN-SigOne.RData" )






# Plotting the MC profile with the GT profile
par(oma=c(1,1,1,1), mar=c(4,5,2,2))
plot(psi.grid4plot,MCprof4plot, type="l", col="black", lwd=2, xlim=c(-1.5,-0.5),
     xlab=bquote("Values of" ~ beta), ylab="Relative Likelihood", bty="l", 
     cex.lab=1.5)
points(ModelPN.GTprof[,1],ModelPN.GTprof[,2], type="l", lty=2, lwd=2, col="red")
points(x=Beta.hatPN.CI,y=c(0,0), col="blue", type="l" ,lty=1,lwd=4)
abline(v=Beta.hat, lwd=2,lty=2)
abline(h=exp(-qchisq(p=0.95,df=1)/2),lty=2, col="gray")
legend("topright", legend=c("GT", "MC", "MLE","DC-CI",bquote("exp(-0.5"~ chi[0.95]^2~")")), col=c("red", "black","black", "blue","gray"), lty=c(2,1,2,1,2), lwd=c(2,2,2,4,1), bty="n", cex=1.15)

  




library(ggplot2); library(patchwork)


PN.all <- ggplot(gridmcdf, aes(Beta.values, Rel.Prof.Like) ) + 
  geom_line(linewidth=1, aes(linetype ="black.line"), color="black") + xlim( c(-2,-0.25) ) +
  labs( x=bquote("Values of" ~ beta), y="Relative Likelihood", title = "") +
  geom_line(data = testPN, aes( x=param, y=rel.like, linetype="red.line"), color = "black", linewidth=1) +
  geom_line(data = data.frame(ModelPN.GTprof), aes( x=beta.profs, y=rel.prof, linetype="anothered.line"), color = "black", linewidth=1) +
  
  geom_line(data=data.frame(x=as.vector(Beta.hatPN.CI), y=c(0,0)), aes(x=x,y=y, linetype = "blue.line"), color = "black", linewidth=2)+
  geom_vline(linetype=5,linewidth=0.95, mapping = aes(linetype = "black.line2", xintercept = -1.031482), color = "black")+
  
  geom_hline(linetype=6,linewidth=0.9, mapping = aes(linetype="gray.line", yintercept = exp(-qchisq(p=0.95,df=1)/2)), color = "black")+ 
  scale_linetype_manual(name = "", values = c(red.line = 4, anothered.line=2, black.line = 3, blue.line = 1, black.line2 = 5, gray.line =6), 
                        labels = c(red.line = "TMB", anothered.line="GT", black.line = "MC", blue.line = "DC-CI", black.line2 = "MLE", gray.line = bquote("exp(-0.5"~ chi[0.95]^2~")")), 
                        limits = c("red.line","anothered.line" ,"black.line", "blue.line", "black.line2", "gray.line"), 
                        guide = guide_legend(override.aes = list(color = rep("black",6), linetype = c(4,2,3,1,5,6), linewidth = c(1,1,1, 2,1,1) )) ) +
  theme_classic() + 
  theme(text = element_text(size=15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.key.width = unit(4, "line"),        
        axis.ticks.length=unit(0.2, "cm"), legend.position = "bottom")   


save(PN.all,file="PNall-SigOne.RData")  





