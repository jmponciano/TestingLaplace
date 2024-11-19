# Auxiliary functions for the simulation tests
library(dclone); library(MASS); library(mcmcplots)

# Data simulator function:

# Function to simulate 3 models: Latent variable is ALWAYS called X
#
# Model 1: Latent variable is X ~ T_2(0), a reasonably flat distribution
# W = X + eta where eta ~ N(0,sigma1)
# Y = beta*X + epsilon where epsilon ~ N(0,sigma2)
# Response = c(W,Y), a bivariate Normal given X=x
#
# Model 2: Latent variable  X ~ T_2(0), a reasonably flat distribution
# W = X + eta where eta ~ N(0,sigma1)
# Y ~ Poisson(lambda=exp(beta*X))
# Response = c(W,Y), a bivariate Poisson+Normal given X=x
#
# Model 3: Latent: X ~ Normal(mu, sd), U = X^2 (~Chi-square(1))
# Observation model: Y ~ Poisson(lambda=U)
#
# parm.list = list of parameters needed for simulation, they change according to
# the model
# model options : "Bivar.Norm","Pois.Norm", "Normsq.Pois
sim.example <- function(n=50, parm.list, model="Bivar.Norm"){
  
  switch(model,
         "Bivar.Norm"={
           beta <- parm.list$beta
           sd1  <- parm.list$sd1
           sd2  <- parm.list$sd2
           df1  <- parm.list$df1
           
           X <- rt(n=n, df=df1,ncp=0)
           W <- X + rnorm(n=n, mean=0, sd=sd1)
           Y <- rnorm(n=n,mean=beta*X, sd=sd2)
           out <- cbind(W,Y)
           colnames(out) <- c("W","Y")
         },
         "Pois.Norm"={
           beta <- parm.list$beta
           sd1  <- parm.list$sd1
           df1  <- parm.list$df1
           
           X <- rt(n=n, df=df1,ncp=0)
           W <- X + rnorm(n=n, mean=0, sd=sd1)
           Y <- rpois(n=n,lambda=exp(beta*X))
           out <- cbind(W,Y)
           colnames(out) <- c("W","Y")
           
         },
         "Normsq.Pois"={
           mu <- parm.list$mu
           sd <- parm.list$sd1
           X <- rnorm(n=n, mean=mu, sd=sd)
           Y <- rpois(n=n, lambda=(X^2))
           out <- as.matrix(Y, nrow=n,ncol=1)
           colnames(out) <- "Y"
         },
         "Logit.Norm"={
           require(boot)
           beta <- parm.list$beta
           sd1  <- parm.list$sd1
           df1  <- parm.list$df1
           
           X <- rt(n=n, df=df1, ncp=0)
           W <- X + rnorm(n=n, mean=0, sd=sd1)
           Y <- rbinom(n=n,size=1,prob=inv.logit(beta*X))
           out <- cbind(W,Y)
           colnames(out) <- c("W","Y")
           
         },
         "Probit.Norm"={
           beta <- parm.list$beta
           sd1  <- parm.list$sd1
           df1  <- parm.list$df1
           
           X <- rt(n=n, df=df1,ncp=0)
           W <- X + rnorm(n=n, mean=0, sd=sd1)
           Y <- rbinom(n=n,size=1,prob=pnorm(beta*X))
           out <- cbind(W,Y)
           colnames(out) <- c("W","Y")
           
         },
         "Norm.nonident.var"={
           mu1  <- parm.list$mu1
           sd1  <- parm.list$sd1
           tau1 <- parm.list$tau1
           X1    <- rnorm(n=n, mean=mu1, sd=sd1)
           Y     <- rnorm(n=n, mean=X1, sd=tau1)
           out   <- cbind(X1,Y)
           colnames(out) <- c("X1","Y")
         }
  )
  return(out)
}





###### 3 Model functions for JAGS, called upon using pacakge DCLONE #########

##### Model 1: Bivariate Normal Error in variables model: 

# Model 1.a:  estimating only beta using data cloning
BivarNormDC.beta <- function(){
  
  # Prior for beta
  beta~dnorm(0,1/50)  
  for(k in 1:K){
    # Hidden process
    for(i in 1:n){
      X[i,k] ~ dt(0,1,df1)
    }
    # Observation process
    for(i in 1:n){
      
      W[i,k]~dnorm(X[i,k], 1/pow(sd1,2))
      Y[i,k]~dnorm(beta*X[i,k], 1/pow(sd2,2))
    }    
  }
  
}


# Model 1.b: estimating all parameters from model 1 using data cloning
BivarNormDC.all <- function(){
  
  # likelihood
  for(k in 1:K){
    # Hidden process
    for(i in 1:n){
      X[i,k] ~ dt(0,1,df1)
    }
    # Observation process
    for(i in 1:n){
      
      W[i,k]~dnorm(X[i,k], 1/pow(sd1,2))
      Y[i,k]~dnorm(beta*X[i,k], 1/pow(sd2,2))
    }    
  }
  
  # Prior for beta
  beta  ~ dnorm(0,1/50)
  #sd1   ~ dlnorm(0,1)
  #sd2   ~ dlnorm(0,1)
  lsd1   ~ dnorm(0,1/50)
  lsd2   ~ dnorm(0,1/50)
  sd1 <- exp(lsd1)
  sd2 <- exp(lsd2)
  
  
  #df1  ~ dpois(1) T(1,)
  
}





# Model 1.c:  model to sample from the posterior of X|Y by passing to it all 
# model parameter values

BivarNorm.Upost <- function(){
  
  # Hidden process
  for(i in 1:n){
    X[i] ~ dt(0,1,df1)
  }
  # Observation process
  for(i in 1:n){
    
    W[i]~dnorm(X[i], 1/pow(sd1,2))
    Y[i]~dnorm(beta*X[i], 1/pow(sd2,2))
  }    
}



#############  Poisson-Normal Error in variables model

# Model 2.a:  estimating only beta using data cloning
PoisNormDC.beta <- function(){
  
  # Prior for beta
  beta~dnorm(0,1/25)  
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


# Model 2.b: estimating all parameters from model 1 using data cloning
PoisNormDC.all <- function(){
  
  # likelihood
  for(k in 1:K){
    # Hidden process
    for(i in 1:n){
      X[i,k] ~ dt(0,1,df1)
    }
    # Observation process
    for(i in 1:n){
      
      W[i,k] ~ dnorm(X[i,k], 1/pow(sd1,2))
      Y[i,k] ~ dpois(exp(beta*X[i,k]))
    }    
  }
  
  # Prior for beta
  beta  ~ dnorm(0,1/50)
  #sd1   ~ dlnorm(0,1)
  lsd1   ~ dnorm(0,1/50)
  sd1 <- exp(lsd1)
  #df1  ~ dpois(1) T(1,)
  
  
}


# Model 2.c:  model to sample from the posterior of X|Y by passing to it all 
# model parameter values

PoisNorm.Upost <- function(){
  
  # Hidden process
  for(i in 1:n){
    X[i] ~ dt(0,1,df1)
  }
  # Observation process
  for(i in 1:n){
    
    W[i]~dnorm(X[i], 1/pow(sd1,2))
    Y[i]~dpois(exp(beta*X[i]))
  }    
}


# Model 3.a:  estimating only beta using data cloning
NormsqPoisDC2.mu <- function(){
  
  # Prior for beta 
  mu~dnorm(2,1/50)
  
  for(k in 1:K){
    # Hidden process
    for(i in 1:n){
      X[i,k] ~ dnorm(mu,1/pow(sd1,2))
    }
    # Observation process
    for(i in 1:n){
      
      Usq[i,k] = pow(X[i,k],2)
      Y[i,k]~dpois(Usq[i,k])
    }    
  }
  
}


NormsqPoisDC.mu <- function(){
  
  # Prior for beta 
  mu~dnorm(0,1/50)
  
  for(k in 1:K){
    # Hidden process
    for(i in 1:n){
      X[i,k] ~ dnorm(mu,1/pow(sd1,2))
    }
    # Observation process
    for(i in 1:n){
      
      Usq[i,k] = pow(X[i,k],2)
      Y[i,k]~dpois(Usq[i,k])
    }    
  }
  
}


# Model 3.b: estimating all parameters from model 1 using data cloning
NormsqPoisDC.all <- function(){
  
  # Prior for beta and sd1
  mu~dnorm(0,1/50)
  #sd1~dlnorm(0,1/50)
  lsd1   ~ dnorm(0,1/50)
  sd1 <- exp(lsd1)

  for(k in 1:K){
    # Hidden process
    for(i in 1:n){
      X[i,k] ~ dnorm(mu,1/pow(sd1,2))
    }
    # Observation process
    for(i in 1:n){
      
      Usq[i,k] = pow(X[i,k],2)
      Y[i,k]~dpois(Usq[i,k])
    }    
  }
  
}


# Model 3.c:  model to sample from the posterior of X|Y by passing to it all 
# model parameter values

NormsqPois.Upost <- function(){
  
  # Hidden process
  for(i in 1:n){
    X[i] ~ dnorm(mu,1/pow(sd1,2))
  }
  # Observation process
  for(i in 1:n){
    
    Usq[i] = pow(X[i],2) 
    Y[i]~dpois(Usq[i])
  }    
}



# Model 4.a:  estimating only beta using data cloning
LogitNormDC.beta <- function(){
  
  # Prior for beta
  beta~dnorm(0,1/50)  
  for(k in 1:K){
    # Hidden process
    for(i in 1:n){
      X[i,k] ~ dt(0,1,df1)
    }
    # Observation process
    for(i in 1:n){
      
      W[i,k]~dnorm(X[i,k], 1/pow(sd1,2))
      Y[i,k]~dbin(1/(1+exp(-beta*X[i,k])),1)
    }    
  }
  
}



# Model 4.b: estimating all parameters from model 1 using data cloning
LogitNormDC.all <- function(){
  
  # likelihood
  for(k in 1:K){
    # Hidden process
    for(i in 1:n){
      X[i,k] ~ dt(0,1,df1)
    }
    # Observation process
    for(i in 1:n){
      
      W[i,k] ~ dnorm(X[i,k], 1/pow(sd1,2))
      Y[i,k] ~ dbin(1/(1+exp(-beta*X[i,k])),1)
      
    }    
  }
  # Prior for beta
  beta  ~ dnorm(0,1/50)
  #sd1   ~ dlnorm(0,1)
  lsd1   ~ dnorm(0,1/50)
  sd1 <- exp(lsd1)
  #df1  ~ dpois(1) T(1,)
  
}


# Model 4.c:  model to sample from the posterior of X|Y by passing to it all 
# model parameter values

LogitNorm.Upost <- function(){
  
  # Hidden process
  for(i in 1:n){
    X[i] ~ dt(0,1,df1)
  }
  # Observation process
  for(i in 1:n){
    
    W[i]~dnorm(X[i], 1/pow(sd1,2))
    Y[i]~dbin(1/(1+exp(-beta*X[i])),1)
  }    
}



# Model 5.a:  estimating only beta using data cloning
ProbitNormDC.beta <- function(){
  
  # Prior for beta
  beta~dnorm(-1,2)  
  for(k in 1:K){
    # Hidden process
    for(i in 1:n){
      X[i,k] ~ dt(0,1,df1)
    }
    # Observation process
    for(i in 1:n){
      
      W[i,k]~dnorm(X[i,k], 1/pow(sd1,2))
      Y[i,k]~dbin(pnorm((beta*X[i,k]),0,1) , 1)
    }    
  }
  
}


# Model 5.b: estimating all parameters from model 1 using data cloning
ProbitNormDC.all <- function(){
  
  # likelihood
  for(k in 1:K){
    # Hidden process
    for(i in 1:n){
      X[i,k] ~ dt(0,1,df1)
    }
    # Observation process
    for(i in 1:n){
      
      W[i,k] ~ dnorm(X[i,k], 1/pow(sd1,2))
      Y[i,k] ~ dbin(pnorm((beta*X[i,k]),0,1) , 1)
      
    }    
  }
  # Prior for beta
  beta  ~ dnorm(0,1/50)
  #sd1   ~ dlnorm(0,1)
  lsd1   ~ dnorm(0,1/50)
  sd1 <- exp(lsd1)
  #df1  ~ dpois(1) T(1,)
  
}


# Model 5.c:  model to sample from the posterior of X|Y by passing to it all 
# model parameter values

ProbitNorm.Upost <- function(){
  
  # Hidden process
  for(i in 1:n){
    X[i] ~ dt(0,1,df1)
  }
  # Observation process
  for(i in 1:n){
    
    W[i]~dnorm(X[i], 1/pow(sd1,2))
    Y[i]~dbin(pnorm((beta*X[i]),0,1) , 1)
  }    
}





onedcfitrun <- function(model="Bivar.Norm", true.parms=NULL,sim.data,
                        n.iter=1000,n.adapt=100,n.update=100,thin=10,
                        n.chains=3,clones.seq=c(1,2),known.dfs=2,parallel=TRUE){
  switch(model,
         "Bivar.Norm"={
           
           if(is.null(true.parms)){
             
             n <- nrow(sim.data)
             W <- sim.data[,1]
             Y <- sim.data[,2]
             
             datalist <- list(K=1, n=n, df1=known.dfs, Y=dcdim(data.matrix(Y)), W=dcdim(data.matrix(W)))
             
             dcrun <- dc.fit(datalist, c("beta",  "sd1","sd2"), BivarNormDC.all, n.clones=clones.seq,
                             multiply="K", unchanged=c("n","df1"),
                             n.chains = n.chains, n.adapt=n.adapt, n.update=n.update, 
                             n.iter = n.iter, thin=thin)
             
           }else{
             
             n <- nrow(sim.data)
             df1 <- true.parms$df1
             sd1 <- true.parms$sd1
             sd2 <- true.parms$sd2
             
             W <- sim.data[,1]
             Y <- sim.data[,2]
             datalist <- list(K=1, n=n, df1=df1,sd1=sd1,sd2=sd2, Y=dcdim(data.matrix(Y)), 
                              W=dcdim(data.matrix(W)))
             
             dcrun <- dc.fit(datalist, c("beta"), BivarNormDC.beta, n.clones=clones.seq,
                             multiply="K", unchanged=c("n", "df1", "sd1", "sd2"),
                             n.chains = n.chains, n.adapt=n.adapt, n.update=n.update, 
                             n.iter = n.iter, thin=thin)
             
           }
           
           
         },
         "Pois.Norm"={
           
           if(is.null(true.parms)){
             
             n <- nrow(sim.data)
             W <- sim.data[,1]
             Y <- sim.data[,2]
             datalist <- list(K=1, n=n, df1=known.dfs, Y=dcdim(data.matrix(Y)), W=dcdim(data.matrix(W)))
             
             dcrun <- dc.fit(datalist, c("beta", "sd1"), PoisNormDC.all, n.clones=clones.seq,
                             multiply="K", unchanged=c("n","df1"),
                             n.chains = n.chains, n.adapt=n.adapt, n.update=n.update, 
                             n.iter = n.iter, thin=thin)
             
           }else{
             
             n <- nrow(sim.data)
             df1 <- true.parms$df1
             sd1 <- true.parms$sd1
             
             W <- sim.data[,1]
             Y <- sim.data[,2]
             datalist <- list(K=1, n=n, df1=df1,sd1=sd1, Y=dcdim(data.matrix(Y)), W=dcdim(data.matrix(W)))
             
             dcrun <- dc.fit(datalist, c("beta"), PoisNormDC.beta, n.clones=clones.seq,
                             multiply="K", unchanged=c("n", "df1", "sd1"),
                             n.chains = n.chains, n.adapt=n.adapt, n.update=n.update, 
                             n.iter = n.iter, thin=thin)
             
           }
           
         },
         
         "Normsq.Pois"={
           if(is.null(true.parms)){
             
             n <- nrow(sim.data)
             Y <- sim.data[,1]
             datalist <- list(K=1, n=n, Y=dcdim(data.matrix(Y)))
             
             dcrun <- dc.fit(datalist, c("mu","sd1"), NormsqPoisDC.all, n.clones=clones.seq,
                             multiply="K", unchanged="n", n.chains = n.chains, n.adapt=n.adapt, 
                             n.update=n.update, n.iter = n.iter, thin=thin, inits=list(mu=-1,sd1=1))
             
           }else{
             
             n <- nrow(sim.data)
             sd1 <- true.parms$sd1
             
             Y <- sim.data[,1]
             datalist <- list(K=1, n=n, sd1=sd1, Y=dcdim(data.matrix(Y)))
             
             dcrun <- dc.fit(datalist, c("mu"), NormsqPoisDC2.mu, n.clones=clones.seq,
                             multiply="K", unchanged=c("n", "sd1"), inits=list(mu=-1),
                             n.chains = n.chains, n.adapt=n.adapt, n.update=n.update, 
                             n.iter = n.iter, thin=thin)
             
             
           }
           
         },
         
         "Logit.Norm"={
           if(is.null(true.parms)){
             
             n <- nrow(sim.data)
             W <- sim.data[,1]  
             Y <- sim.data[,2]
             datalist <- list(K=1, n=n,df1=known.dfs, W=dcdim(data.matrix(W)), Y=dcdim(data.matrix(Y)))
             
             dcrun <- dc.fit(datalist, c("beta", "sd1"), LogitNormDC.all, n.clones=clones.seq,
                             multiply="K", unchanged=c("n","df1"),n.chains = n.chains, n.adapt=n.adapt, 
                             n.update=n.update, n.iter = n.iter, thin=thin)
             
           }else{
             
             n <- nrow(sim.data)
             df1 <- true.parms$df1
             sd1 <- true.parms$sd1
             
             W <- sim.data[,1]  
             Y <- sim.data[,2]
             datalist <- list(K=1, n=n, sd1=sd1,df1=df1, W=dcdim(data.matrix(W)), Y=dcdim(data.matrix(Y)))
             
             dcrun <- dc.fit(datalist, c("beta"), LogitNormDC.beta, n.clones=clones.seq,
                             multiply="K", unchanged=c("n", "df1" , "sd1"), inits=list(beta=-1),
                             n.chains = n.chains, n.adapt=n.adapt, n.update=n.update, 
                             n.iter = n.iter, thin=thin)
             
             
           }
           
         },
         
         "Probit.Norm"={
           if(is.null(true.parms)){
             
             n <- nrow(sim.data)
             W <- sim.data[,1]  
             Y <- sim.data[,2]
             datalist <- list(K=1, n=n, df1=known.dfs, W=dcdim(data.matrix(W)), Y=dcdim(data.matrix(Y)))
             
             dcrun <- dc.fit(datalist, c("beta", "sd1"), ProbitNormDC.all, n.clones=clones.seq,
                             multiply="K", unchanged=c("n","df1"),n.chains = n.chains, n.adapt=n.adapt, 
                             n.update=n.update, n.iter = n.iter, thin=thin)
             
           }else{
             
             n <- nrow(sim.data)
             df1 <- true.parms$df1
             sd1 <- true.parms$sd1
             
             W <- sim.data[,1]  
             Y <- sim.data[,2]
             
             if(parallel==TRUE){
             
               cl <- makePSOCKcluster(parallel::detectCores())
               datalist <- list(K=1, n=n, sd1=sd1,df1=df1, W=dcdim(data.matrix(W)), Y=dcdim(data.matrix(Y)))
               
               dcrun <- dc.parfit(cl=cl, datalist, c("beta"), ProbitNormDC.beta, n.clones=clones.seq,
                               multiply="K", unchanged=c("n", "df1" , "sd1"), inits=list(beta=-1),
                               n.chains = n.chains, n.adapt=n.adapt, n.update=n.update, 
                               n.iter = n.iter, thin=thin)
               stopCluster(cl)
                                
             }else{
             
              datalist <- list(K=1, n=n, sd1=sd1,df1=df1, W=dcdim(data.matrix(W)), Y=dcdim(data.matrix(Y)))
             
              dcrun <- dc.fit(datalist, c("beta"), ProbitNormDC.beta, n.clones=clones.seq,
                             multiply="K", unchanged=c("n", "df1" , "sd1"), inits=list(beta=-1),
                             n.chains = n.chains, n.adapt=n.adapt, n.update=n.update, 
                             n.iter = n.iter, thin=thin)
             }
             
           }#end else
           
         }#end probit switch
         
         
  )#end Switch statement
  
  return(dcrun)
  
}  

Xpost.samp <- function(model="Bivar.Norm", parm.list,sim.data,
                       n.iter=1000,n.adapt=100,n.update=100,thin=10,
                       n.chains=3){
  switch(model,
         "Bivar.Norm"={
           
           beta <- parm.list$beta
           sd1  <- parm.list$sd1
           sd2  <- parm.list$sd2
           df1  <- parm.list$df1
           
           n <- nrow(sim.data)
           W <- sim.data[,1]
           Y <- sim.data[,2]
           
           # Sampling from f(X|(W,Y))   
           datalist <- list(W=W, Y=Y, beta=beta, df1=df1,sd1=sd1,sd2 =sd2, n=n)
           out.parms <- c("X")
           outfit <- dc.fit(datalist,params=out.parms, model=BivarNorm.Upost,
                            multiply=NULL, unchanged=c("n","beta", "df1", "sd1", "sd2"), n.chains=n.chains,
                            n.adapt=n.adapt, n.update=n.update, 
                            n.iter = n.iter, thin=thin,n.clones=1) 
           
         },
         
         "Pois.Norm"={
           beta <- parm.list$beta
           sd1  <- parm.list$sd1
           df1  <- parm.list$df1
           
           n <- nrow(sim.data)
           W <- sim.data[,1]
           Y <- sim.data[,2]
           
           
           # Sampling from f(X|(W,Y))   
           datalist <- list(W=W, Y=Y, beta=beta, df1=df1,sd1=sd1, n=n)
           out.parms <- c("X")
           outfit <- dc.fit(datalist,params=out.parms, model=PoisNorm.Upost,
                            multiply=NULL, unchanged=NULL, n.chains=n.chains,
                            n.adapt=n.adapt, n.update=n.update, 
                            n.iter = n.iter, thin=thin,n.clones=1) 
           
         },
         
         "Normsq.Pois"={
           mu <- parm.list$mu
           sd1   <- parm.list$sd1
           
           n <- nrow(sim.data)
           Y <- sim.data[,1]
           
           # Sampling from f(X|(W,Y))   
           datalist <- list(Y=Y, mu=mu, sd1=sd1, n=n)
           out.parms <- c("X")
           outfit <- dc.fit(datalist,params=out.parms, model=NormsqPois.Upost,
                            multiply="n", unchanged=c("mu", "sd1"), 
                            n.chains=n.chains,n.adapt=n.adapt,n.update=n.update, 
                            n.iter = n.iter, thin=thin,n.clones=1) 
           
         },
         "Logit.Norm"={
           require(boot)
           beta <- parm.list$beta
           sd1  <- parm.list$sd1
           df1  <- parm.list$df1
           
           n <- nrow(sim.data)
           W <- sim.data[,1]
           Y <- sim.data[,2]
           
           # Sampling from f(X|(W,Y))   
           datalist <- list(W=W, Y=Y, beta=beta, df1=df1,sd1=sd1, n=n)
           out.parms <- c("X")
           outfit <- dc.fit(datalist,params=out.parms, model=LogitNorm.Upost,
                            multiply=NULL, unchanged=NULL, n.chains=n.chains,
                            n.adapt=n.adapt, n.update=n.update, 
                            n.iter = n.iter, thin=thin,n.clones=1) 
           
         },
         
         "Probit.Norm"={
           beta <- parm.list$beta
           sd1  <- parm.list$sd1
           df1  <- parm.list$df1
           
           n <- nrow(sim.data)
           W <- sim.data[,1]
           Y <- sim.data[,2]
           
           # Sampling from f(X|(W,Y))   
           datalist <- list(W=W, Y=Y, beta=beta, df1=df1,sd1=sd1, n=n)
           out.parms <- c("X")
           outfit <- dc.fit(datalist,params=out.parms, model=ProbitNorm.Upost,
                            multiply=NULL, unchanged=NULL, n.chains=n.chains,
                            n.adapt=n.adapt, n.update=n.update, 
                            n.iter = n.iter, thin=thin,n.clones=1) 
         }         
  )
  return(outfit)
}



###### This function is under construction
GT.prof <- function(model="Bivar.Norm", parm.list,sim.data,
                    XpostSample, bracket=2, prec.prof=0.01,
                    plot.it=TRUE, overlay=TRUE,mycol="blue",
                    dotted=1){
  
  switch(model,
         "Bivar.Norm"={
           
           df1 <- parm.list$df1
           sd1 <- parm.list$sd1
           sd2 <- parm.list$sd2
           beta<- parm.list$beta
           n <- nrow(sim.data)
           W <- sim.data[,1]
           Y <- sim.data[,2]
           
           B <- nrow(XpostSample)
           
           
           # Computing GT's denominator for every sampled X and every value of the
           # parameter 'beta' we are profiling over
           
           beta.profs <- seq(from=beta-bracket,to=beta+bracket, by=prec.prof)
           len.prof <- length(beta.profs)
           prof.vec <- rep(0,len.prof)
           
           for(i in 1:len.prof){
             
             ith.beta <- beta.profs[i]
             #GT.lnumer.vec <- rep(0,B)
             #GT.ldenom.vec <- rep(0,B)
             GT.logdiff <- rep(0,B)
             for(j in 1:B){
               X.star <- XpostSample[j,]
               GT.lnumer.persamp<- sum(dnorm(x=Y, mean=ith.beta*X.star, sd=sd2, log=TRUE)) 
               #+ dnorm(x=W, mean=X.star,sd=sd1, log=TRUE) +
               #dt(x=X.star, df=df1, ncp=0,log=TRUE))  --> this part cancels! 
               
               GT.ldenom.persamp <- sum(dnorm(x=Y, mean=beta*X.star, sd=sd2, log=TRUE))
               
               GT.logdiff[j] <- GT.lnumer.persamp-GT.ldenom.persamp
             }
             
             ratio.prof.vec.star <- exp(GT.logdiff)
             prof.vec[i] <- mean(ratio.prof.vec.star)
             
           }
           
           rel.prof <- prof.vec #/max(prof.vec)  
           if((plot.it==TRUE)&(overlay==FALSE)){plot(beta.profs,rel.prof, type="l", col=mycol, 
                                                     xlab="Values of Beta",
                                                     ylab="Relative Likelihood profile", 
                                                     lty=dotted)
           }else if((plot.it==TRUE)&(overlay==TRUE)){
             points(beta.profs,rel.prof, type="l", col=mycol, 
                    xlab="",
                    ylab="", lty=dotted)
           }
           
           out <- cbind(beta.profs,rel.prof,prof.vec)
           
         },
         
         "Pois.Norm"={
           
           df1 <- parm.list$df1
           sd1 <- parm.list$sd1
           beta<- parm.list$beta
           n <- nrow(sim.data)
           W <- sim.data[,1]
           Y <- sim.data[,2]
           
           B <- nrow(XpostSample)
           
           
           # Computing GT's denominator for every sampled X and every value of the
           # parameter 'beta' we are profiling over
           
           beta.profs <- seq(from=beta-bracket,to=beta+bracket, by=prec.prof)
           len.prof <- length(beta.profs)
           prof.vec <- rep(0,len.prof)
           
           for(i in 1:len.prof){
             
             ith.beta <- beta.profs[i]
             #GT.lnumer.vec <- rep(0,B)
             #GT.ldenom.vec <- rep(0,B)
             GT.logdiff <- rep(0,B)
             for(j in 1:B){
               X.star <- XpostSample[j,]
               GT.lnumer.persamp<- sum(dpois(x=Y, lambda=exp(ith.beta*X.star), log=TRUE)) 
               #+ sum(dnorm(x=W, mean=X.star,sd=sd1, log=TRUE)) +
               #sum(dt(x=X.star, df=df1, ncp=0,log=TRUE))

                 
               GT.ldenom.persamp <-  sum(dpois(x=Y, lambda=exp(beta*X.star), log=TRUE))
               #+ sum(dnorm(x=W, mean=X.star,sd=sd1, log=TRUE)) +
               #  sum(dt(x=X.star, df=df1, ncp=0,log=TRUE))
               
               
               GT.logdiff[j] <- GT.lnumer.persamp-GT.ldenom.persamp
             }
             
             ratio.prof.vec.star <- exp(GT.logdiff) #/max(exp(GT.logdiff))
             
             prof.vec[i] <- mean(ratio.prof.vec.star)
             
           }
           
           rel.prof <- prof.vec
           
           if((plot.it==TRUE)&(overlay==FALSE)){plot(beta.profs,rel.prof, type="l", col=mycol, 
                                                     xlab="Values of Beta",
                                                     ylab="Relative Likelihood profile", lty=dotted)
           }else if((plot.it==TRUE)&(overlay==TRUE)){
             points(beta.profs,rel.prof, type="l",col=mycol,
                    xlab="",
                    ylab="", lty=dotted)
           }
           
           out <- cbind(beta.profs,rel.prof,prof.vec)
           
           
           
         },
         "Normsq.Pois"={
           
           sd1 <- parm.list$sd1
           beta<- parm.list$mu
           n <- nrow(sim.data)
           Y <- sim.data[,1]
           
           B <- nrow(XpostSample)
           
           # Computing GT's denominator for every sampled X and every value of the
           # parameter 'beta' we are profiling over
           
           beta.profs <- seq(from=beta-bracket,to=beta+bracket, by=prec.prof)
           len.prof <- length(beta.profs)
           prof.vec <- rep(0,len.prof)
           
           for(i in 1:len.prof){
             
             ith.beta <- beta.profs[i]
             GT.logdiff <- rep(0,B)
             
             for(j in 1:B){
               X.star <- XpostSample[j,]
               GT.lnumer.persamp<- dnorm(x=X.star, mean=ith.beta, sd=sd1, log=TRUE)+
                 dchisq(x=(X.star^2)/(sd1^2), df=1,ncp=(ith.beta^2), log=TRUE)
               # + dpois(x=Y, lambda=Usq.star, log=TRUE) --> it cancels
               
               GT.ldenom.persamp <-  dnorm(x=X.star, mean=beta, sd=sd1, log=TRUE)+
                 dchisq(x=(X.star^2)/(sd1^2), df=1,ncp=(beta^2), log=TRUE)
               # + dpois(x=Y, lambda=Usq.star, log=TRUE) --> it cancels
               
               GT.logdiff[j] <- sum(GT.lnumer.persamp-GT.ldenom.persamp)
             }
             
             ratio.prof.vec.star <- exp(GT.logdiff)
             prof.vec[i] <- mean(ratio.prof.vec.star)
             
           }
           
           
           rel.prof <- prof.vec/max(prof.vec)  
           if((plot.it==TRUE)&(overlay==FALSE)){plot(beta.profs,rel.prof, type="l", col=mycol, 
                                                     xlab="Values of Beta",
                                                     ylab="Relative Likelihood profile", lty=dotted)
           }else if((plot.it==TRUE)&(overlay==TRUE)){
             points(beta.profs,rel.prof, type="l", col=mycol, 
                    xlab="",
                    ylab="", lty=dotted)
           }
           out <- cbind(beta.profs,rel.prof,prof.vec)
           
         },
         "Logit.Norm"={

           df1 <- parm.list$df1
           sd1 <- parm.list$sd1
           beta<- parm.list$beta
           n <- nrow(sim.data)
           W <- sim.data[,1]
           Y <- sim.data[,2]
           
           B <- nrow(XpostSample)
           # Computing GT's denominator for every sampled X and every value of the
           # parameter 'beta' we are profiling over
           
           beta.profs <- seq(from=beta-bracket,to=beta+bracket, by=prec.prof)
           len.prof <- length(beta.profs)
           prof.vec <- rep(0,len.prof)
           
           for(i in 1:len.prof){
             
             ith.beta <- beta.profs[i]
             #GT.lnumer.vec <- rep(0,B)
             #GT.ldenom.vec <- rep(0,B)
             GT.logdiff <- rep(0,B)
             for(j in 1:B){
               X.star <- XpostSample[j,]
               GT.lnumer.persamp<- dbinom(x=Y, size=1, prob=1/(1+exp(-ith.beta*X.star)), log=TRUE) 
               #+ dnorm(x=W, mean=X.star,sd=sd1, log=TRUE) +
               #dt(x=X.star, df=df1, ncp=0,log=TRUE))  --> this part cancels!
               
               GT.ldenom.persamp <-  dbinom(x=Y, size=1, prob=1/(1+exp(-beta*X.star)), log=TRUE)
               
               GT.logdiff[j] <- sum(GT.lnumer.persamp-GT.ldenom.persamp)
             }
             
             ratio.prof.vec.star <- exp(GT.logdiff)
             prof.vec[i] <- mean(ratio.prof.vec.star)
             
           }
           
           rel.prof <- prof.vec#/max(prof.vec)  
           if((plot.it==TRUE)&(overlay==FALSE)){plot(beta.profs,rel.prof, type="l", col=mycol, 
                                                     xlab="Values of Beta",
                                                     ylab="Relative Likelihood", lty=dotted)
           }else if((plot.it==TRUE)&(overlay==TRUE)){
             points(beta.profs,rel.prof, type="l",col=mycol,
                    xlab="",
                    ylab="", lty=dotted)
           }
           
           out <- cbind(beta.profs,rel.prof,prof.vec)
           
           
         },
         "Probit.Norm"={
           
           #Y[i,k]~dbin(pnorm((beta*X[i,k]),0,1) , 1) 
           df1 <- parm.list$df1
           sd1 <- parm.list$sd1
           beta<- parm.list$beta
           n <- nrow(sim.data)
           W <- sim.data[,1]
           Y <- sim.data[,2]
           
           B <- nrow(XpostSample)
           # Computing GT's denominator for every sampled X and every value of the
           # parameter 'beta' we are profiling over
           
           beta.profs <- seq(from=beta-bracket,to=beta+bracket, by=prec.prof)
           len.prof <- length(beta.profs)
           prof.vec <- rep(0,len.prof)
           
           for(i in 1:len.prof){
             
             ith.beta <- beta.profs[i]
             #GT.lnumer.vec <- rep(0,B)
             #GT.ldenom.vec <- rep(0,B)
             GT.logdiff <- rep(0,B)
             for(j in 1:B){
               X.star <- XpostSample[j,]
               pvec.num   <- pnorm(q=ith.beta*X.star, mean=0, sd=1)
               pvec.denom <- pnorm(q=beta*X.star, mean=0,sd=1)
               pvec.num[pvec.num==1]<- 0.99
               pvec.num[pvec.num==0] <- .Machine$double.xmin
               pvec.denom[pvec.denom==1]<- 0.99
               pvec.denom[pvec.denom==0] <- .Machine$double.xmin

               GT.lnumer.persamp <- sum(Y*log(pvec.num) + (1-Y)*log(1-pvec.num)
                                        + dnorm(x=W, mean=X.star,sd=sd1, log=TRUE) +
                                          dt(x=X.star, df=df1, ncp=0,log=TRUE)) 
               GT.ldenom.persamp <- sum(Y*log(pvec.denom) + (1-Y)*log(1-pvec.denom)
                                        + dnorm(x=W, mean=X.star,sd=sd1, log=TRUE) +
                                          dt(x=X.star, df=df1, ncp=0,log=TRUE))
               
               #+ dnorm(x=W, mean=X.star,sd=sd1, log=TRUE) +
               #dt(x=X.star, df=df1, ncp=0,log=TRUE))  --> this part cancels! 
               
               #GT.lnumer.persamp<-   dbinom(x=Y, size=1, prob=pnorm(q=ith.beta*X.star, mean=0, sd=1), log=TRUE) 
               #+ dnorm(x=W, mean=X.star,sd=sd1, log=TRUE) +
               #dt(x=X.star, df=df1, ncp=0,log=TRUE))  --> this part cancels!
               #GT.ldenom.persamp <-  dbinom(x=Y, size=1, prob=pnorm(q=beta*X.star, mean=0,sd=1), log=TRUE)
               
               GT.logdiff[j] <- GT.lnumer.persamp-GT.ldenom.persamp
             }
             
             ratio.prof.vec.star <- exp(GT.logdiff)
             prof.vec[i] <- mean(ratio.prof.vec.star)
             
           }
           
           # plot(beta.profs, prof.vec, type="l"); abline(v=-1);abline(h=1)
           # 
           # beta.max <- beta.profs[prof.vec==max(prof.vec)]
           # 
           # beta <- beta.max
           # 
           # prof.vec <- rep(0,len.prof)
           # 
           # for(i in 1:len.prof){
           #   
           #   ith.beta <- beta.profs[i]
           #   #GT.lnumer.vec <- rep(0,B)
           #   #GT.ldenom.vec <- rep(0,B)
           #   GT.logdiff <- rep(0,B)
           #   for(j in 1:B){
           #     X.star <- XpostSample[j,]
           #     pvec.num   <- pnorm(q=ith.beta*X.star, mean=0, sd=1)
           #     pvec.denom <- pnorm(q=beta*X.star, mean=0,sd=1)
           #     pvec.num[pvec.num==1]<- 0.99
           #     pvec.num[pvec.num==0] <- .Machine$double.xmin
           #     pvec.denom[pvec.denom==1]<- 0.99
           #     pvec.denom[pvec.denom==0] <- .Machine$double.xmin
           #     
           #     GT.lnumer.persamp <- sum(Y*log(pvec.num) + (1-Y)*log(1-pvec.num)
           #                              + dnorm(x=W, mean=X.star,sd=sd1, log=TRUE) +
           #                                dt(x=X.star, df=df1, ncp=0,log=TRUE)) 
           #     GT.ldenom.persamp <- sum(Y*log(pvec.denom) + (1-Y)*log(1-pvec.denom)
           #                              + dnorm(x=W, mean=X.star,sd=sd1, log=TRUE) +
           #                                dt(x=X.star, df=df1, ncp=0,log=TRUE))
           #     
           #     #+ dnorm(x=W, mean=X.star,sd=sd1, log=TRUE) +
           #     #dt(x=X.star, df=df1, ncp=0,log=TRUE))  --> this part cancels! 
           #     
           #     #GT.lnumer.persamp<-   dbinom(x=Y, size=1, prob=pnorm(q=ith.beta*X.star, mean=0, sd=1), log=TRUE) 
           #     #+ dnorm(x=W, mean=X.star,sd=sd1, log=TRUE) +
           #     #dt(x=X.star, df=df1, ncp=0,log=TRUE))  --> this part cancels!
           #     #GT.ldenom.persamp <-  dbinom(x=Y, size=1, prob=pnorm(q=beta*X.star, mean=0,sd=1), log=TRUE)
           #     
           #     GT.logdiff[j] <- GT.lnumer.persamp-GT.ldenom.persamp
           #   }
           #   
           #   ratio.prof.vec.star <- exp(GT.logdiff)
           #   prof.vec[i] <- mean(ratio.prof.vec.star)
           #   
           # }           
           # 
           # plot(beta.profs, prof.vec, type="l"); abline(v=-1);abline(h=1)
                    
           #rel.prof <- prof.vec/max(prof.vec)  
           if((plot.it==TRUE)&(overlay==FALSE)){plot(beta.profs,rel.prof, type="l", col=mycol, 
                                                     xlab="Values of Beta",
                                                     ylab="Relative Likelihood", lty=dotted)
           }else if((plot.it==TRUE)&(overlay==TRUE)){
             points(beta.profs,rel.prof, type="l",col=mycol,
                    xlab="",
                    ylab="", lty=dotted)
           }
           
           out <- cbind(beta.profs,rel.prof=prof.vec/max(prof.vec),prof.vec)
           
           
         }
  )# End switch
  return(out)
  
}




MCprofile1d <- function(B= 50000, parm.list,sim.data, bracket=2, prec.prof=0.01, 
                        model="Bivar.Norm", plot.it=TRUE,overlay=TRUE){
  
  switch(model,
         "Bivar.Norm"={
           
           beta <- parm.list$beta
           sd1  <- parm.list$sd1
           sd2  <- parm.list$sd2
           df1  <- parm.list$df1
           
           W <- sim.data[,1]
           Y <- sim.data[,2]
           n <- nrow(sim.data)
           X.B <- rt(n=B,df=df1,ncp=0)
           beta.profs <- seq(from=beta-bracket,to=beta+bracket, by=prec.prof)
           log.prof.vec <- rep(0,length(beta.profs))
           for(j in 1:length(beta.profs)){
             
             jth.beta <- beta.profs[j] 
             log.mean.like.vec <- rep(0,n)
             #log.1ton.llike.mat <- matrix(0,nrow=B,ncol=n)
             
             for(i in 1:n){
               #W.llike <- dnorm(rep(W[i],B), mean=X.B,sd=sd1,log=TRUE)
               #Y.llike <- dnorm(rep(Y[i],B),mean=jth.beta*X.B, sd=sd2, log=TRUE)
               #log.1ton.llike.mat[,i] <- W.llike + Y.llike
               
               W.like <- dnorm(rep(W[i],B), mean=X.B,sd=sd1)
               Y.like <- dnorm(rep(Y[i],B),mean=jth.beta*X.B, sd=sd2)
               log.mean.like.vec[i] <- log(mean(Y.like*W.like))
             }
             log.prof.vec[j] <- sum(log.mean.like.vec)
             # sum.across.n <- apply(log.1ton.llike.mat,1,sum)
             # if(sum(is.infinite(sum.across.n))>0){
             #   
             #   which.inf <- which(is.infinite(sum.across.n))
             #   sum.across.n2 <- sum.across.n[-which.inf]
             #   Blikes <- exp(sum.across.n2)
             #   which.inf2 <- which(is.infinite(Blikes))
             #   if(sum(which.inf2)>0){Blikes <- Blikes[-which.inf2]}
             #   Blikes[Blikes==0] <- .Machine$double.xmin
             #   log.prof.vec[j] <- log(mean(Blikes))                            
             #   
             # }else{
             #   Blikes <- exp(sum.across.n)
             #   which.inf2 <- which(is.infinite(Blikes))
             #   if(sum(which.inf2)>0){Blikes <- Blikes[-which.inf2]}
             #   Blikes[Blikes==0] <- .Machine$double.xmin
             #   log.prof.vec[j] <- log(mean(Blikes))             
             # }
             
           }
           rel.prof <- exp(log.prof.vec-max(log.prof.vec))
           out <- cbind(beta.profs,rel.prof)
           colnames(out) <- c("Beta values", "Rel Prof Like")
           if(((plot.it==TRUE)&(overlay==FALSE))){
             plot(beta.profs,rel.prof, type="l", col="red", 
                  xlab= bquote("Values of" ~ beta),lwd=2,
                  ylab="Relative Likelihood profile",bty="l",
                  ylim=c(0,1))
           }else if(((plot.it==TRUE)&(overlay==TRUE))){
             points(beta.profs,rel.prof,lwd=2, type="l", col="blue", lty=2)
           }
         },
         "Pois.Norm"={
           beta <- parm.list$beta
           sd1  <- parm.list$sd1
           df1  <- parm.list$df1
           
           W <- sim.data[,1]
           Y <- sim.data[,2]
           n <- nrow(sim.data)
           X.B <- rt(n=B,df=df1,ncp=0)
           beta.profs <- seq(from=beta-bracket,to=beta+bracket, by=prec.prof)
           log.prof.vec <- rep(0,length(beta.profs))
           for(j in 1:length(beta.profs)){
             
             jth.beta <- beta.profs[j] 
             log.mean.like.vec <- rep(0,n)
             #log.1ton.llike.mat <- matrix(0,nrow=B,ncol=n)
             
             for(i in 1:n){
               #W.llike <- dnorm(rep(W[i],B), mean=X.B,sd=sd1,log=TRUE)
               #Y.llike <- dpois(rep(Y[i],B),lambda=exp(jth.beta*X.B), log=TRUE)
               #log.1ton.llike.mat[,i] <- W.llike+Y.llike
               
               W.like <- dnorm(rep(W[i],B), mean=X.B,sd=sd1)
               Y.like <- dpois(rep(Y[i],B),lambda=exp(jth.beta*X.B))
               log.mean.like.vec[i] <- log(mean(Y.like*W.like))
             }
             log.prof.vec[j] <- sum(log.mean.like.vec)
             # sum.across.n <- apply(log.1ton.llike.mat,1,sum)             
             # if(sum(is.infinite(sum.across.n))>0){
             #   
             #   which.inf <- which(is.infinite(sum.across.n))
             #   sum.across.n2 <- sum.across.n[-which.inf]
             #   Blikes <- exp(sum.across.n2)
             #   which.inf2 <- which(is.infinite(Blikes))
             #   if(sum(which.inf2)>0){Blikes <- Blikes[-which.inf2]} 
             #   Blikes[Blikes==0] <- .Machine$double.xmin
             #   log.prof.vec[j] <- log(mean(Blikes))
             # }else{
             #   Blikes <- exp(sum.across.n)
             #   which.inf2 <- which(is.infinite(Blikes))
             #   if(sum(which.inf2)>0){Blikes <- Blikes[-which.inf2]}
             #   Blikes[Blikes==0] <- .Machine$double.xmin
             #   log.prof.vec[j] <- log(mean(Blikes))             
             # }
             
           }
           
           rel.prof <- exp(log.prof.vec-max(log.prof.vec))
           out <- cbind(beta.profs,rel.prof)
           colnames(out) <- c("Beta values", "Rel Prof Like")
           if(((plot.it==TRUE)&(overlay==FALSE))){
             plot(beta.profs,rel.prof, type="l", col="red", 
                  xlab= bquote("Values of" ~ beta), lwd=2,
                  ylab="Relative Likelihood profile",bty="l",
                  ylim=c(0,1))
           }else if(((plot.it==TRUE)&(overlay==TRUE))){
             points(beta.profs,rel.prof, lwd=2, type="l", col="blue", lty=2)
           }
           
         },
         
         "Normsq.Pois"={
           beta <- parm.list$mu
           sd1 <- parm.list$sd1
           Y <- sim.data[,1]
           n <- nrow(sim.data)
           
           beta.profs <- seq(from=beta-bracket,to=beta+bracket, by=prec.prof)
           log.prof.vec <- rep(0,length(beta.profs))
           
           for(j in 1:length(beta.profs)){
             
             jth.beta <- beta.profs[j] 
             X.B <- rnorm(n=B, mean=jth.beta,sd=sd1)
             log.mean.like.vec <- rep(0,n)
             #log.1ton.llike.mat <- matrix(0,nrow=B,ncol=n)
             for(i in 1:n){
               #Y.llike <- dpois(rep(Y[i],B), lambda=(X.B^2),log=TRUE)
               #log.1ton.llike.mat[,i] <- Y.llike
               Y.like <- dpois(rep(Y[i],B), lambda=(X.B^2))
               log.mean.like.vec[i] <- log(mean(Y.like))
             }
             log.prof.vec[j] <- sum(log.mean.like.vec)
             
             # sum.across.n <- apply(log.1ton.llike.mat,1,sum)
             # if(sum(is.infinite(sum.across.n))>0){
             #   
             #   which.inf <- which(is.infinite(sum.across.n))
             #   sum.across.n2 <- sum.across.n[-which.inf]
             #   Blikes <- exp(sum.across.n2)
             #   which.inf2 <- which(is.infinite(Blikes))
             #   if(sum(which.inf2)>0){Blikes <- Blikes[-which.inf2]} 
             #   Blikes[Blikes==0] <- .Machine$double.xmin
             #   log.prof.vec[j] <- log(mean(Blikes))
             # }else{
             #   Blikes <- exp(sum.across.n)
             #   which.inf2 <- which(is.infinite(Blikes))
             #   if(sum(which.inf2)>0){Blikes <- Blikes[-which.inf2]}
             #   Blikes[Blikes==0] <- .Machine$double.xmin
             #   log.prof.vec[j] <- log(mean(Blikes))             
             # }
             
           }
           rel.prof <- exp(log.prof.vec-max(log.prof.vec))
           out <- cbind(beta.profs,rel.prof)
           colnames(out) <- c("mu values", "Rel Prof Like")
           if(((plot.it==TRUE)&(overlay==FALSE))){
             plot(beta.profs,rel.prof, type="l", col="red", 
                  xlab= bquote("Values of" ~ mu),lwd=2,
                  ylab="Relative Likelihood profile",bty="l",
                  ylim=c(0,1))
           }else if(((plot.it==TRUE)&(overlay==TRUE))){
             points(beta.profs,rel.prof, type="l",lwd=2, col="blue", lty=2)
           }
           
         },
         "Logit.Norm"={
           beta <- parm.list$beta
           sd1  <- parm.list$sd1
           df1  <- parm.list$df1
           
           W <- sim.data[,1]
           Y <- sim.data[,2]
           n <- nrow(sim.data)
           #X.B <- rt(n=B,df=df1,ncp=0)
           beta.profs <- seq(from=beta-bracket,to=beta+bracket, by=prec.prof)
           log.prof.vec <- rep(0,length(beta.profs))
           for(j in 1:length(beta.profs)){
             X.B <- rt(n=B,df=df1,ncp=0)
             jth.beta <- beta.profs[j] 
             log.mean.like.vec <- rep(0,n)
             #log.1ton.llike.mat <- matrix(0,nrow=B,ncol=n)
             
             for(i in 1:n){
               #W.llike <- dnorm(rep(W[i],B), mean=X.B,sd=sd1, log=TRUE)
               #Y.llike <- dbinom(rep(Y[i],B),size=1, prob=inv.logit(jth.beta*X.B), log=TRUE)
               #log.1ton.llike.mat[,i] <- W.llike+Y.llike # This is a vector of length B
               W.like <- dnorm(rep(W[i],B), mean=X.B,sd=sd1)
               Y.like <- dbinom(rep(Y[i],B),size=1, prob=inv.logit(jth.beta*X.B)) #prob=1/(1+exp(-jth.beta*X.B)))
               log.mean.like.vec[i] <- log(mean(Y.like*W.like))
             }
             
             log.prof.vec[j] <- sum(log.mean.like.vec)
             
             #now exponentiate the sum across obs 1 to n
             # sum.across.n <- apply(log.1ton.llike.mat,1,sum)
             # if(sum(is.infinite(sum.across.n))>0){
             # 
             #   which.inf <- which(is.infinite(sum.across.n))
             #   sum.across.n2 <- sum.across.n[-which.inf]
             #   Blikes <- exp(sum.across.n2)
             #   which.inf2 <- which(is.infinite(Blikes))
             #   if(sum(which.inf2)>0){Blikes <- Blikes[-which.inf2]}
             #   Blikes[Blikes==0] <- .Machine$double.xmin
             #   log.prof.vec[j] <- log(mean(Blikes))
             # }else{
             #   Blikes <- exp(sum.across.n)
             #   which.inf2 <- which(is.infinite(Blikes))
             #   if(sum(which.inf2)>0){Blikes <- Blikes[-which.inf2]}
             #   Blikes[Blikes==0] <- .Machine$double.xmin
             #   log.prof.vec[j] <- log(mean(Blikes))
             # }
             
             
           }
           rel.prof <- exp(log.prof.vec-max(log.prof.vec))
           out <- cbind(beta.profs,rel.prof)
           colnames(out) <- c("Beta values", "Rel Prof Like")
           if(((plot.it==TRUE)&(overlay==FALSE))){
             plot(beta.profs,rel.prof, type="l", col="red", 
                  xlab= bquote("Values of" ~ beta),lwd=2,
                  ylab="Relative Likelihood profile",bty="l",
                  ylim=c(0,1))
           }else if(((plot.it==TRUE)&(overlay==TRUE))){
             points(beta.profs,rel.prof, type="l",lwd=2, col="blue", lty=2)
           }
           
           
         },
         "Probit.Norm"={
           beta <- parm.list$beta
           sd1  <- parm.list$sd1
           df1  <- parm.list$df1
           
           W <- sim.data[,1]
           Y <- sim.data[,2]
           n <- nrow(sim.data)
           
           beta.profs <- seq(from=beta-bracket,to=beta+bracket, by=prec.prof)
           log.prof.vec <- rep(0,length(beta.profs))
           for(j in 1:length(beta.profs)){
             
             jth.beta <- beta.profs[j] 
             log.mean.like.vec <- rep(0,n)

             for(i in 1:n){
               
               X.B <- rt(n=B,df=df1)  
               W.like <- dnorm(rep(W[i],B), mean=X.B,sd=sd1)
               Y.like <- dbinom(rep(Y[i],B),size=1,prob=pnorm(q=jth.beta*X.B))
               log.mean.like.vec[i] <- log(mean(Y.like*W.like))
             }
             log.prof.vec[j] <- sum(log.mean.like.vec)

           }
           rel.prof <- exp(log.prof.vec-max(log.prof.vec))
           #rel.prof2 <- exp(log.prof.vec)/max(exp(log.prof.vec))
           beta.profs[rel.prof==max(rel.prof)]
           out <- cbind(beta.profs,rel.prof)
           colnames(out) <- c("Beta values", "Rel Prof Like")
           if(((plot.it==TRUE)&(overlay==FALSE))){
             plot(beta.profs,rel.prof, type="l", col="red", 
                  xlab= bquote("Values of" ~ beta),lwd=2,
                  ylab="Relative Likelihood profile",bty="l",
                  ylim=c(0,1))
           }else if(((plot.it==TRUE)&(overlay==TRUE))){
             points(beta.profs,rel.prof, type="l",lwd=2, col="blue", lty=2)
           }
           
           
         }
  )
  return(out)
}







# 
# GT.denom.vec <- rep(0,B)
# 
# for(j in 1:B){
#   X.star <- Xpost.samp[j,]
#   GT.denom.vec[j] <- sum(dnorm(x=X.star, mean=0, sd=sd1, log=TRUE))
#   + sum(dpois(x=Y, lambda=(beta + X.star)^2, log=TRUE))
# }
# 
# 
# for(i in 1:len.prof){
#   
#   ith.beta <- beta.profs[i]
#   GT.lnumer.vec <- rep(0,B)
#   GT.logratio.vec <- rep(0,B)
#   for(j in 1:B){
#     X.star <- Xpost.samp[j,]
#     GT.lnumer.vec[j] <- sum(dnorm(x=X.star, mean=0, sd=sd1, log=TRUE))
#     + sum(dpois(x=Y, lambda=(ith.beta+X.star)^2, log=TRUE))
#     
#     GT.logratio.vec[j] <- GT.lnumer.vec[j] - GT.denom.vec[j]
#   }
#   
#   ratio.prof.vec.star <- exp(GT.logratio.vec)
#   prof.vec[i] <- mean(ratio.prof.vec.star)
#   
# }



# Function to do a one dimensional Profile Likelihood via simple 
# Monte Carlo averaging for parameter 'beta' in the same 3 models 
# in function "sim.example()" above, and for a given simulated data set

# Latent variable is ALWAYS called X
# Model 1: Latent variable is X ~ T_2(0), a reasonably flat distribution 
# W = X + eta where eta ~ N(0,sigma1)
# Y = beta*X + epsilon where epsilon ~ N(0,sigma2)
# Response = c(W,Y), a bivariate Normal given X=x
#
# Model 2: Latent variable  X ~ T_2(0), a reasonably flat distribution 
# W = X + eta where eta ~ N(0,sigma1)
# Y ~ Poisson(lambda=exp(beta*X))
# Response = c(W,Y), a bivariate Poisson+Normal given X=x
#
# Model 3: Latent: X ~ Normal(mu, sd), U = X^2 (~Chi-square(1))
# Observation model: Y ~ Poisson(lambda=U)
#
# parm.list = list of parameters needed for simulation, they change according to 
# the model

# model options : "Bivar.Norm","Pois.Norm", "Normsq.Pois
# bracket = the bracket around the true value of the model parameter beta
# around which to compute the MC profile likelihood
# prec.prof = the interval in the sequence of the profile likelihood parameter
# values resulting from the bracketing








