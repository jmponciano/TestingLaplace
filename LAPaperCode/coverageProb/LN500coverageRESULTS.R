rm(list=ls()) # remove everything currently held in the R memory
graphics.off() # close all open graphics window


load("tmboptsLN.RData")
load("LN500SIMS.RData")


### Estimate coverage prob
parmsLN <- list(beta=-1, df1=10, sd1=3)
NSIMS      <- 500
CI.mat.TMB <- matrix(data=NA, ncol=4, nrow = NSIMS)
colnames(CI.mat.TMB) <- c("Truth","MLE","lCI","uCI")

for(i in 1:NSIMS){
  ci.obj <- cbind(tmb.opts[[i]]$MLres$opt.ML$par - 1.96*sqrt(solve(tmb.opts[[i]]$MLres$opt.ML$hessian)),
                  tmb.opts[[i]]$MLres$opt.ML$par + 1.96*sqrt(solve(tmb.opts[[i]]$MLres$opt.ML$hessian)))
  CI.mat.TMB[i,] <- c( parmsLN$beta, tmb.opts[[i]]$MLres$opt.ML$par, ci.obj )
}

# remove the NAs (if any)
CI.mat.TMB   <- na.omit(CI.mat.TMB)
# CI.mat.TMB <- CI.mat.TMB[CI.mat.TMB[,2]>-2,]
coverage.TMB <- sum(CI.mat.TMB[,1] >= CI.mat.TMB[,3] & CI.mat.TMB[,1] <= CI.mat.TMB[,4])/dim(CI.mat.TMB)[1]
coverage.TMB
# 0.88

contained <- CI.mat.TMB[,1] >= CI.mat.TMB[,3] & CI.mat.TMB[,1] <= CI.mat.TMB[,4]

# plot(density(CI.mat.TMB[,2]), xlim=c(-10,3))
# abline(v=parmsNN$beta)
# abline(v=-parmsNN$beta)




#### MC results
library(Hmisc)
load("/Users/constantingeorgeglen/Documents/Databases/Analyses/Program/TMB/code/LaplaceApproximation/Examples/coverageProb/optimBFGS_LN500covergedBASH.RData")
# load("/Users/constantingeorgeglen/Documents/Databases/Analyses/Program/TMB/code/LaplaceApproximation/Examples/coverageProb/optimNelderM_LN500covergedBASH.RData")


### Estimate coverage prob
### remove models where the hessian is problematic
NSIMS    <- 500
NSIMSvec <- 1:500
obj.rm <- which(sapply(1:NSIMS, function(x) res[[x]]$obj$hessian == 0 ))
if(length(obj.rm)>0){
  res <- lapply(NSIMSvec[NSIMSvec %nin% obj.rm], function(x) res[[x]] )
}else{
  res <- res
}
length(res)
# for( x in 1:NSIMS) solve(res[[x]]$obj$hessian)

### Estimate coverage prob
parmsLN <- list(beta=-1, df1=10, sd1=3)
NSIMS      <- length(res)
CI.mat.MC  <- matrix(data=NA, ncol=4, nrow = NSIMS)
colnames(CI.mat.MC) <- c("Truth","MLE","lCI","uCI")

for(i in 1:NSIMS){
  ci.obj <-tryCatch(cbind(res[[i]]$obj$par - 1.96*sqrt(solve(res[[i]]$obj$hessian)),
                          res[[i]]$obj$par + 1.96*sqrt(solve(res[[i]]$obj$hessian))),
                    error=function(x) NULL)
  CI.mat.MC[i,] <- c( parmsLN$beta, res[[i]]$obj$par, ci.obj )
}


# remove the NAs (if any)
CI.mat.MC   <- na.omit(CI.mat.MC)
NSIMS <- dim(CI.mat.MC)[1]
# CI.mat.MC   <- CI.mat.MC[CI.mat.MC[,4]<20,]

coverage.MC <- sum(CI.mat.MC[,1] >= CI.mat.MC[,3] & CI.mat.MC[,1] <= CI.mat.MC[,4])/dim(CI.mat.MC)[1]
coverage.MC
# 0.9950739

contained <- CI.mat.MC[,1] >= CI.mat.MC[,3] & CI.mat.MC[,1] <= CI.mat.MC[,4]

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


### average bias
median(CI.mat.TMB[,1] - CI.mat.TMB[,2]) # 
median(CI.mat.MC[,1] - CI.mat.MC[,2])   # 


## MC and TMB hess regress
summary(lm(tmbparhess~mcparhess))

png( filename = "LNregress_hess.png", 
     width = 12, height = 8, units = "in", res = 1200, pointsize = 10
)
par(mar=c(5,5,2,3), mfrow=c(1,1))
plot(mcparhess, tmbparhess, pch=16, cex=1.1, 
     cex.axis=1.5, cex.lab=1.5, cex.main=2,
     xlab=expression(MC ~ H(beta)), ylab=expression(TMB ~ H(beta)), bty="l",
     col="grey30", main = "Logit-Normal")
# abline(reg = c(0,1), lty="dashed", lwd=1.5)
# points(parmsNN$beta,parmsNN$beta, pch=15, cex=1.2)
dev.off()



## MC and TMB inv hess regress
summary(lm(tmbparinvhess~mcparinvhess))


png( filename = "LNregress_invhess.png", 
     width = 12, height = 8, units = "in", res = 1200, pointsize = 10
)
par(mar=c(5,5,2,3), mfrow=c(1,1))
plot(mcparinvhess, tmbparinvhess, pch=16, ylim=c(0,1),
     cex=1.1, cex.axis=1.5, cex.lab=1.5, cex.main=2,
     xlab=expression(MC ~ H(beta)^-1), ylab=expression(TMB ~ H(beta)^-1), bty="l",
     col="grey30", main = "Logit-Normal")
# abline(reg = c(0,1), lty="dashed", lwd=1.5)
# points(parmsNN$beta,parmsNN$beta, pch=15, cex=1.2)
dev.off()






## MC and TMB est regress
summary(lm(tmbparest~mcparest))


png( filename = "LNregress_est.png", 
     width = 12, height = 8, units = "in", res = 1200, pointsize = 10
)
par(mar=c(5,5,2,3), mfrow=c(1,1))
plot(mcparest, tmbparest, pch=16, cex=1.1, 
     xlim=c(-3,3), 
     ylim=c(-2,1),
     cex.axis=1.5, cex.lab=1.5, cex.main=2,
     xlab=expression(MC ~ beta), ylab=expression(TMB ~ beta), bty="l",
     col="grey30", main = "Logit-Normal")
# abline(reg = c(0,1), lty="dashed", lwd=1.5)
# abline(reg = lm(tmbparest~mcparest), lty="solid", lwd=1.5)
# points(parmsNN$beta,parmsNN$beta, pch=15, cex=1.2)
dev.off()




### Coverage plot
png( filename = "LNcoverage.png",
     width = 12, height = 8, units = "in", res = 1200, pointsize = 10
)

par(mar=c(5,5,3,3), mfrow=c(1,1))
cols <- c("grey30", "grey50")
par(mfrow=c(1,2))
plot( c(0,0), type="n", xlim=c(range(CI.mat.TMB[,3:4])), ylim=c(1,NSIMS),
      xlab=expression(beta), ylab="Iteration", bty="l", 
      cex.axis=1.5, cex.lab=1.5, cex.main=1.5,
      main=paste0("Logit-Normal TMB coverage: ",coverage.TMB) )
for(i in 1:NSIMS){
  
  segments(x0=CI.mat.TMB[i,3],x1=CI.mat.TMB[i,4],y0=i,y1=i, 
           col=ifelse(contained[i]==T, 
                      scales::alpha(cols[1],0.1), 
                      scales::alpha(cols[2], 0.1)), 
           lwd=2)
  points( x=CI.mat.TMB[i,2], y=i,
          col=ifelse(contained[i]==T, cols[1], cols[2]), pch=16)
  
}
abline(v=parmsLN$beta, lwd=1, lty="dashed", col="black")
# abline(v=-parmsNN$beta, lwd=2, col="darkred")


plot( c(0,0), type="n", 
      # xlim=c(-150,15),
      xlim=c(range(CI.mat.MC[,3:4])), 
      ylim=c(1,NSIMS), 
      xlab=expression(beta), ylab="", bty="l", 
      cex.axis=1.5, cex.lab=1.5, cex.main=1.5,
      main=paste0("Logit-Normal MC coverage: ",coverage.MC))
for(i in 1:NSIMS){
  
  segments(x0=CI.mat.MC[i,3],x1=CI.mat.MC[i,4],y0=i,y1=i, 
           col=ifelse(contained[i]==T, 
                      scales::alpha(cols[1],0.1), 
                      scales::alpha(cols[2], 0.1)), 
           lwd=2)
  points( x=CI.mat.MC[i,2], y=i,
          col=ifelse(contained[i]==T, cols[1], cols[2]), pch=16)
  
}
abline(v=parmsLN$beta, lwd=1, lty="dashed", col="black")
# abline(v=-parmsNN$beta, lwd=2, col="darkred")
par(mfrow=c(1,1))
dev.off()


plot(density(CI.mat.MC[,2]))

