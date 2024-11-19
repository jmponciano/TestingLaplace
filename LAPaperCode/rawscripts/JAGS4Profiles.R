# Profiling functions for JAGS, all models

#### ----------------------------------------------------------------
#### Normal-Normal model
#### ----------------------------------------------------------------

BivarNormDC.MLE <- function(){
  
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

# Profiling over beta, estimating the rest:
BivarNorm.betaprof <- function(){
  
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
  #beta  ~ dnorm(0,1/50)
  #sd1   ~ dlnorm(0,1)
  #sd2   ~ dlnorm(0,1)
  lsd1   ~ dnorm(0,1/50)
  lsd2   ~ dnorm(0,1/50)
  sd1 <- exp(lsd1)
  sd2 <- exp(lsd2)
  
  
  #df1  ~ dpois(1) T(1,)
  
}

# Profiling over lsd1, estimating the rest
BivarNorm.lsd1prof <- function(){
  
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
  #lsd1   ~ dnorm(0,1/50)
  lsd2   ~ dnorm(0,1/50)
  sd1 <- exp(lsd1)
  sd2 <- exp(lsd2)
  
  
  #df1  ~ dpois(1) T(1,)
  
}

# Profiling over lsd2, estimating the rest
BivarNorm.lsd2prof <- function(){
  
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
  #lsd2   ~ dnorm(0,1/50)
  sd1 <- exp(lsd1)
  sd2 <- exp(lsd2)
  
  #df1  ~ dpois(1) T(1,)
  
}

# sample the hidden states for X
BivarNorm.hXgY <- function(){
  
  sd1 <- exp(lsd1)
  sd2 <- exp(lsd2)
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





#### ----------------------------------------------------------------
#### Poisson-Normal model
#### ----------------------------------------------------------------
# Model 2.b: estimating all parameters from model 1 using data cloning
PoisNormDC.MLE <- function(){
  
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

# Profiling over beta, estimating the rest:
PoisNorm.betaprof <- function(){
  
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
  #beta  ~ dnorm(0,1/50)
  #sd1   ~ dlnorm(0,1)
  lsd1   ~ dnorm(0,1/50)
  sd1 <- exp(lsd1)
  #df1  ~ dpois(1) T(1,)
  
  
}

# Profiling over lsd1, estimating beta
PoisNorm.lsd1prof <- function(){
  
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
  #lsd1   ~ dnorm(0,1/50)
  sd1 <- exp(lsd1)
  #df1  ~ dpois(1) T(1,)
  
}

# sample the hidden states for X
PoisNorm.hXgY <- function(){
  
  sd1 <- exp(lsd1)
  # Hidden process
  for(i in 1:n){
    X[i] ~ dt(0,1, df1)
  }
  # Observation process
  for(i in 1:n){
    
    W[i] ~ dnorm(X[i], 1/pow(sd1,2))
    Y[i] ~ dpois(exp(beta*X[i]))
    
  }    
}



#### ----------------------------------------------------------------
#### Logit-Normal model
#### ----------------------------------------------------------------
LogitNormDC.MLE <- function(){
  
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


# Profiling over beta, estimating the rest:
LogitNorm.betaprof <- function(){
  
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
  #beta  ~ dnorm(0,1/50)
  #sd1   ~ dlnorm(0,1)
  lsd1  ~ dnorm(0,1/50)
  sd1   <- exp(lsd1)
  #df1   ~ dpois(1) T(1,)
  
  
}

# Profiling over lsd1, estimating beta
LogitNorm.lsd1prof <- function(){
  
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
  #lsd1   ~ dnorm(0,1/50)
  sd1 <- exp(lsd1)
  #df1  ~ dpois(1) T(1,)
  
}

# sample the hidden states for X
LogitNorm.hXgY <- function(){
  
  sd1 <- exp(lsd1)
  # Hidden process
  for(i in 1:n){
    X[i] ~ dt(0,1, df1)
  }
  # Observation process
  for(i in 1:n){
    
    W[i] ~ dnorm(X[i], 1/pow(sd1,2))
    Y[i] ~ dbin(1/(1+exp(-beta*X[i])),1)
    
  }    
}




#### ----------------------------------------------------------------
#### Probit-Normal model
#### ----------------------------------------------------------------
ProbitNormDC.MLE <- function(){
  
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


# Profiling over beta, estimating the rest:
ProbitNorm.betaprof <- function(){
  
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
  #beta  ~ dnorm(0,1/50)
  #sd1   ~ dlnorm(0,1)
  lsd1  ~ dnorm(0,1/50)
  sd1   <- exp(lsd1)
  #df1   ~ dpois(1) T(1,)
  
  
}

# Profiling over lsd1, estimating beta
ProbitNorm.lsd1prof <- function(){
  
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
  #lsd1   ~ dnorm(0,1/50)
  sd1 <- exp(lsd1)
  #df1  ~ dpois(1) T(1,)
  
}

# sample the hidden states for X
ProbitNorm.hXgY <- function(){
  
  sd1 <- exp(lsd1)
  # Hidden process
  for(i in 1:n){
    X[i] ~ dt(0,1, df1)
  }
  # Observation process
  for(i in 1:n){
    
    W[i] ~ dnorm(X[i], 1/pow(sd1,2))
    Y[i] ~ dbin(pnorm((beta*X[i]),0,1) , 1)
    
  }    
}




#### ----------------------------------------------------------------
#### Normal squared Poisson
#### ----------------------------------------------------------------
# Model 3.b: estimating all parameters from model 1 using data cloning
NormsqPoisDC.MLE <- function(){
  
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

# Profiling over mu, estimating lsd1:
NormsqPois.muprof <- function(){
  
  #mu~dnorm(0,1/50)
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

# Profiling oover lsd1, estimating mu:
NormsqPois.lsd1prof <- function(){
  
  # Prior for beta 
  mu~dnorm(0,1/50)
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



NormsqPois.Upost <- function(){
  
  sd1 <- exp(lsd1)
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





#### ----------------------------------------------------------------
#### Non-identifiable variance models
#### ----------------------------------------------------------------

# Model 3.b: estimating all parameters from using data cloning
NonIdentVarDC.MLE <- function(){
  
  # Prior for mu, tau1, and sd1
  mu1    ~ dnorm(0,1/50)
  lsd1  ~ dnorm(0,1/50)
  ltau1 ~ dnorm(0,1/50)
  sd1   <- exp(lsd1)
  tau1  <- exp(ltau1)
  
  for(k in 1:K){
    
    # Hidden process
    for(i in 1:n){
      X[i,k] ~ dnorm(mu1, 1/pow(sd1,2))
    }
    
    # Observation process
    for(i in 1:n){
      Y[i,k] ~ dnorm(X[i,k], 1/pow(tau1,2))
    }    
  }
}


# Profiling over mu, estimating ltau1 and lsd1:
NonIdentVar.muprof <- function(){
  
  # Prior for mu, tau1, and sd1
  #mu    ~ dnorm(0,1/50)
  lsd1  ~ dnorm(0,1/50)
  ltau1 ~ dnorm(0,1/50)
  sd1   <- exp(lsd1)
  tau1  <- exp(ltau1)
  
  for(k in 1:K){
    
    # Hidden process
    for(i in 1:n){
      X[i,k] ~ dnorm(mu1, 1/pow(sd1,2))
    }
    
    # Observation process
    for(i in 1:n){
      Y[i,k] ~ dnorm(X[i,k], 1/pow(tau1,2))
    }    
  }
}



# Profiling oover lsd1, estimating mu and ltau1:
NonIdentVar.lsd1prof <- function(){
  
  # Prior for mu, tau1, and sd1
  mu1    ~ dnorm(0,1/50)
  #lsd1  ~ dnorm(0,1/50)
  ltau1 ~ dnorm(0,1/50)
  #sd1   <- exp(lsd1)
  tau1  <- exp(ltau1)
  
  for(k in 1:K){
    
    # Hidden process
    for(i in 1:n){
      X[i,k] ~ dnorm(mu1, 1/pow(sd1,2))
    }
    
    # Observation process
    for(i in 1:n){
      Y[i,k] ~ dnorm(X[i,k], 1/pow(tau1,2))
    }    
  }
}

# Profiling over ltau1, estimating mu and lsd1:
NonIdentVar.ltau1prof <- function(){
  
  # Prior for mu, tau1, and sd1
  mu1    ~ dnorm(0,1/50)
  lsd1  ~ dnorm(0,1/50)
  #ltau1 ~ dnorm(0,1/50)
  sd1   <- exp(lsd1)
  #tau1  <- exp(ltau1)
  
  for(k in 1:K){
    
    # Hidden process
    for(i in 1:n){
      X[i,k] ~ dnorm(mu1, 1/pow(sd1,2))
    }
    
    # Observation process
    for(i in 1:n){
      Y[i,k] ~ dnorm(X[i,k], 1/pow(tau1,2))
    }    
  }
}



NonIdentVar.hXgY <- function(){
  
  tau1 <- exp(ltau1)
  sd1  <- exp(lsd1)
  
  # Hidden process
  for(i in 1:n){
    X[i] ~ dnorm(mu1, 1/pow(sd1,2))
  }
  # Observation process
  for(i in 1:n){
    
    Y[i] ~ dnorm(X[i], 1/pow(tau1,2))
    
  }    
}















