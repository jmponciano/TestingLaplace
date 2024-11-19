rm(list=ls()) # remove everything currently held in the R memory
graphics.off() # close all open graphics window

load("tmboptsPNdf3.RData")
load("PN500SIMSdf3.RData")

library(magrittr)

### Estimate coverage prob
parmsPN    <- list( beta=-1, df1=2, sd1=3 )
NSIMS      <- 500
CI.mat.TMB <- matrix(data=NA, ncol=4, nrow = NSIMS)
colnames(CI.mat.TMB) <- c("Truth","MLE","lCI","uCI")

for(i in 1:NSIMS){
  ci.obj <- cbind(tmb.opts[[i]]$MLres$opt.ML$par - 1.96*sqrt(solve(tmb.opts[[i]]$MLres$opt.ML$hessian)),
                  tmb.opts[[i]]$MLres$opt.ML$par + 1.96*sqrt(solve(tmb.opts[[i]]$MLres$opt.ML$hessian)))
  CI.mat.TMB[i,] <- c( parmsPN$beta, tmb.opts[[i]]$MLres$opt.ML$par, ci.obj )
}

# remove the NAs (if any)
CI.mat.TMB   <- na.omit(CI.mat.TMB)


## coverage estimates for the true perimeters
coverage.TMB <- sum(CI.mat.TMB[,1] >= CI.mat.TMB[,3] & CI.mat.TMB[,1] <= CI.mat.TMB[,4])/dim(CI.mat.TMB)[1]
coverage.TMB
# 0.384: df=2
# 0.408: df=3


## Number of cases where the optimum was the opposite side of the true value
sum(1 >= CI.mat.TMB[,3] & 1 <= CI.mat.TMB[,4])/dim(CI.mat.TMB)[1]
# 0.354: df=2
# 0.432: df=3


### identify which rows contain the true parameter values
contained.TMB <- CI.mat.TMB[,1] >= CI.mat.TMB[,3] & CI.mat.TMB[,1] <= CI.mat.TMB[,4]

# plot(density(CI.mat.TMB[,2]), xlim=c(-10,3))
# abline(v=parmsNN$beta)
# abline(v=-parmsNN$beta)




#### MC results
load(paste0("/Users/constantingeorgeglen/Documents/Databases/Analyses/Program/TMB/code/LaplaceApproximation/Examples/coverageProb/",
            "optimBFGS_PN500covergedBASHdf3.RData"))

### Estimate coverage prob
parmsPN    <- list( beta=-1, df1=2, sd1=3 )
NSIMS      <- 500
CI.mat.MC  <- matrix(data=NA, ncol=4, nrow = NSIMS)
colnames(CI.mat.MC) <- c("Truth","MLE","lCI","uCI")

for(i in 1:NSIMS){
  ci.obj <- cbind(res[[i]]$obj$par - 1.96*sqrt(solve(res[[i]]$obj$hessian)),
                  res[[i]]$obj$par + 1.96*sqrt(solve(res[[i]]$obj$hessian)))
  CI.mat.MC[i,] <- c( parmsPN$beta, res[[i]]$obj$par, ci.obj )
}

# remove the NAs (if any)
CI.mat.MC   <- na.omit(CI.mat.MC)
coverage.MC <- sum(CI.mat.MC[,1] >= CI.mat.MC[,3] & CI.mat.MC[,1] <= CI.mat.MC[,4])/dim(CI.mat.MC)[1]
coverage.MC
# 0.836: df=2
# 0.812: df=3

contained.MC <- CI.mat.MC[,1] >= CI.mat.MC[,3] & CI.mat.MC[,1] <= CI.mat.MC[,4]

# plot(density(CI.mat.TMB[,2]), xlim=c(-10,3))
# abline(v=parmsNN$beta)
# abline(v=-parmsNN$beta)

### par est, hessian, and inverse hessian
tmbparhess <- as.vector(unlist(lapply(1:NSIMS, function(x) tmb.opts[[x]]$MLres$opt.ML$hessian)))
mcparhess  <- unlist(lapply(1:NSIMS, function(x) res[[x]]$obj$hessian))
tmbparinvhess <- as.vector(unlist(lapply(1:NSIMS, function(x) solve(tmb.opts[[x]]$MLres$opt.ML$hessian))))
mcparinvhess  <- unlist(lapply(1:NSIMS, function(x) solve(res[[x]]$obj$hessian)))
tmbparest <- as.vector(unlist(lapply(1:NSIMS, function(x) tmb.opts[[x]]$MLres$opt.ML$par)))
mcparest  <- unlist(lapply(1:NSIMS, function(x) res[[x]]$obj$par))

### make a full data set
full.dataset <- data.frame(tmbparest=tmbparest,
                           tmbparhess=tmbparhess,
                           tmbparinvhess=tmbparinvhess,
                           mcparest=mcparest,
                           mcparhess=mcparhess,
                           mcparinvhess=mcparinvhess)




## MC and TMB hess regress
summary(lm(tmbparhess~mcparhess))


### outlier
# tmb.opts[[379.2075135]]$MLres$opt.ML
# res[[379.2075135]]$obj

png( filename = "PNregress_hess.png", 
     width = 12, height = 8, units = "in", res = 1200, pointsize = 10
)
par(mar=c(5,5,2,3), mfrow=c(1,1))
plot(mcparhess, tmbparhess, pch=16, cex=1.1, 
     xlim=c(0,230), ylim=c(0,230),
     cex.axis=1.5, cex.lab=1.5, cex.main=2,
     xlab=expression(MC ~ H(beta)), ylab=expression(TMB ~ H(beta)), bty="l",
     col="grey30", main = "Poisson-Normal")
abline(reg = c(0,1), lty="dashed", lwd=1.5)
abline(reg = lm(tmbparhess ~ mcparhess, data = full.dataset %>% 
                  dplyr::filter( tmbparhess<=230, mcparhess<=230 )), 
       lty="solid", lwd=1.5)

# points(parmsNN$beta,parmsNN$beta, pch=15, cex=1.2)
dev.off()



## MC and TMB inv hess regress
summary(lm(tmbparinvhess~mcparinvhess))


png( filename = "PNregress_invhess.png", 
     width = 12, height = 8, units = "in", res = 1200, pointsize = 10
)
par(mar=c(5,5,2,3), mfrow=c(1,1))
plot(mcparinvhess, tmbparinvhess, pch=16, 
     xlim=c(0,0.2), ylim=c(0,0.2),
     cex=1.1, cex.axis=1.5, cex.lab=1.5, cex.main=2,
     xlab=expression(MC ~ H(beta)^-1), ylab=expression(TMB ~ H(beta)^-1), 
     bty="l",
     col="grey30", main = "Poisson-Normal")
abline(reg = c(0,1), lty="dashed", lwd=1.5)
abline(reg = lm(tmbparinvhess ~ mcparinvhess, data = full.dataset %>% 
                  dplyr::filter( tmbparinvhess<=0.2, mcparinvhess<=0.2 )), 
       lty="solid", lwd=1.5)
# abline(reg = lm(tmbparinvhess ~ mcparinvhess), lty="solid", lwd=1.5)
# points(parmsNN$beta,parmsNN$beta, pch=15, cex=1.2)
dev.off()






## MC and TMB est regress
summary(lm(tmbparest~mcparest))


png( filename = "PNregress_est.png", 
     width = 12, height = 8, units = "in", res = 1200, pointsize = 10
)
par(mar=c(5,5,2,3), mfrow=c(1,1))
plot(mcparest, tmbparest, pch=16, cex=1.1, 
     xlim=c(-2,1), ylim=c(-2,3),
     cex.axis=1.5, cex.lab=1.5, cex.main=2,
     xlab=expression(MC ~ beta), ylab=expression(TMB ~ beta), bty="l",
     col="grey30", main = "Poisson-Normal")
abline(reg = c(0,1), lty="dashed", lwd=1.5)
abline(reg = lm(tmbparest ~ mcparest, data = full.dataset), 
       lty="solid", lwd=1.5)
# points(parmsNN$beta,parmsNN$beta, pch=15, cex=1.2)
dev.off()




### Coverage plot
png( filename = "PNcoverage.png",
     width = 12, height = 8, units = "in", res = 1200, pointsize = 10
)

par(mar=c(5,5,3,3), mfrow=c(1,1))
cols <- c("grey30", "grey50")
par(mfrow=c(1,2))
plot( c(0,0), type="n", xlim=c(range(CI.mat.TMB[,3:4])), ylim=c(1,NSIMS),
      xlab=expression(beta), ylab="Iteration", bty="l", 
      cex.axis=1.5, cex.lab=1.5, cex.main=1.5,
      main=paste0("Poisson-Normal TMB coverage: ",coverage.TMB) )
for(i in 1:NSIMS){
  
  segments(x0=CI.mat.TMB[i,3],x1=CI.mat.TMB[i,4],y0=i,y1=i, 
           col=ifelse(contained.TMB[i]==T, 
                      scales::alpha(cols[1],0.1), 
                      scales::alpha(cols[2], 0.1)), 
           lwd=2)
  points( x=CI.mat.TMB[i,2], y=i,
          col=ifelse(contained.TMB[i]==T, cols[1], cols[2]), pch=16)
  
}
abline(v=parmsPN$beta, lwd=1, lty="dashed", col="black")
# abline(v=-parmsNN$beta, lwd=2, col="darkred")


plot( c(0,0), type="n", xlim=c(range(CI.mat.MC[,3:4])), ylim=c(1,NSIMS),
      xlab=expression(beta), ylab="", bty="l", 
      cex.axis=1.5, cex.lab=1.5, cex.main=1.5,
      main=paste0("Poisson-Normal MC coverage: ",coverage.MC))
for(i in 1:NSIMS){
  
  segments(x0=CI.mat.MC[i,3],x1=CI.mat.MC[i,4],y0=i,y1=i, 
           col=ifelse(contained.MC[i]==T, 
                      scales::alpha(cols[1],0.1), 
                      scales::alpha(cols[2], 0.1)), 
           lwd=2)
  points( x=CI.mat.MC[i,2], y=i,
          col=ifelse(contained.MC[i]==T, cols[1], cols[2]), pch=16)
  
}
abline(v=parmsPN$beta, lwd=1, lty="dashed", col="black")
# abline(v=-parmsNN$beta, lwd=2, col="darkred")
par(mfrow=c(1,1))
dev.off()




