#### r model file for data doubling
library(dclone); library(MASS); library(mcmcplots)



###### Model functions to do data doubling

#############  Bivariate Normal Error in variables model
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
  beta  ~ dnorm(0,1/100)
  #sd1   ~ dlnorm(0,1)
  #sd2   ~ dlnorm(0,1)
  lsd1   ~ dnorm(0,1/100)
  lsd2   ~ dnorm(0,1/100)
  sd1 <- exp(lsd1)
  sd2 <- exp(lsd2)
  #df1  ~ dpois(1) T(1,)
  
}
BivarNorm.alg2 <- function(){
  
  # likelihood
  for(k in 1:K){
    # Hidden process
    for(i in 1:n){
      X[i,k] ~ dt(0,1,df1)
    }
    # Observation process
    for(i in 1:n){
      W[i,k] ~ dnorm(X[i,k], 1/pow(sd1,2))
      Y[i,k] ~ dnorm(beta*X[i,k], 1/pow(sd2,2))
    }    
  }
  # Prior for beta
  parms  ~ dmnorm(MLE,FI)
  beta  <- parms[1]
  #sd1   <- parms[2]
  #sd2   <- parms[3]
  lsd1   <- parms[2]
  lsd2   <- parms[3]
  sd1 <- exp(lsd1)
  sd2 <- exp(lsd2)
  #df1   <- parms[4]
  
}
BivarNorm.alg1 <- function(){
  # likelihood
  for(k in 1:K){
    # Hidden process
    for(i in 1:n){
      X[i,k] ~ dt(0,1, df1)
    }
    
    # Observation process
    for(i in 1:n){
      
      W[i,k] ~ dnorm(X[i,k], 1/pow(sd1,2))
      Y[i,k] ~ dnorm(beta*X[i,k], 1/pow(sd2,2))
    }    
  }
  
  # Priors on model parameters: 
  beta ~ dnorm(0,1/100)
  #sd1  ~ dlnorm(0,1)
  #sd2  ~ dlnorm(0,1)
  lsd1   ~ dnorm(0,1/100)
  lsd2   ~ dnorm(0,1/100)
  sd1 <- exp(lsd1)
  sd2 <- exp(lsd2)
  #df1  ~ dpois(1) T(1,)
  
}




# #############  Poisson-Normal Error in variables model
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
  beta  ~ dnorm(0,1/100)
  #sd1   ~ dlnorm(0,1)
  lsd1   ~ dnorm(0,1/100)
  sd1 <- exp(lsd1)
  #df1  ~ dpois(1) T(1,)
  
  
}
PoisNorm.alg2 <- function(){
  
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
  parms  ~ dmnorm(MLE,FI)
  beta  <- parms[1]
  #sd1   <- parms[2]
  lsd1   <- parms[2]
  sd1 <- exp(lsd1)
  #df1  ~ dpois(1) T(1,)

}
PoisNorm.alg1 <- function(){
  # likelihood
  for(k in 1:K){
    # Hidden process
    for(i in 1:n){
      X[i,k] ~ dt(0,1, df1)
    }
    
    # Observation process
    for(i in 1:n){
      
      W[i,k] ~ dnorm(X[i,k], 1/pow(sd1,2))
      Y[i,k] ~ dpois(exp(beta*X[i,k]))
    }    
  }
  
  # Priors on model parameters: 
  beta ~ dnorm(0,1/100)
  #sd1  ~ dlnorm(0,1)
  lsd1   ~ dnorm(0,1/100)
  sd1 <- exp(lsd1)
  #df1  ~ dpois(1) T(1,)
  

}





# #############  Normal squared variables model
NormsqPois.MLE <- function(){
  
  # likelihood
  for(k in 1:K){
    # Hidden process
    for(i in 1:n){
      X[i,k] ~ dnorm(mu,1/pow(sd1,2))
    }
    # Observation process
    for(i in 1:n){
      
      Usq[i,k] <- pow(X[i,k], 2)
      Y[i,k]   ~ dpois(Usq[i,k])
      
    }    
  }
  # Prior for beta and sd1
  mu  ~ dnorm(1, 1/100)
  # sd1 ~ dlnorm(1, 1/100)
  lsd1   ~ dnorm(0,1/100)
  sd1 <- exp(lsd1)
  
}
NormsqPois.alg2 <- function(){
  
  # likelihood
  for(k in 1:K){
    # Hidden process
    for(i in 1:n){
      X[i,k] ~ dnorm(mu, 1/pow(sd1,2))
    }
    # Observation process
    for(i in 1:n){
      
      Usq[i,k] <- pow(X[i,k],2)
      Y[i,k]   ~ dpois(Usq[i,k])
    }    
  }
  # Prior for beta
  parms  ~ dmnorm(MLE,FI)
  mu  <- parms[1]
  sd1 <- exp(parms[2])
}
NormsqPois.alg1 <- function(){
  
  # likelihood
  for(k in 1:K){
    # Hidden process
    for(i in 1:n){
      X[i,k] ~ dnorm(mu,1/pow(sd1,2))
    }
    # Observation process
    for(i in 1:n){
      
      Usq[i,k] <- pow(X[i,k],2)
      Y[i,k]   ~ dpois(Usq[i,k])
    }    
  }
  # Prior for beta and sd1
  mu  ~ dnorm(1, 1/100)
  # sd1 ~ dlnorm(1, 1/100)
  lsd1 ~ dnorm(1, 1/100)
  sd1 <- exp(lsd1)
  
}




# #############  LOGIT-Normal Error in variables model
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
  beta  ~ dnorm(0,1/100)
  #sd1   ~ dlnorm(0,1)
  lsd1   ~ dnorm(0,1/100)
  sd1 <- exp(lsd1)
  #df1  ~ dpois(1) T(1,)
  
}
LogitNorm.alg2 <- function(){

  
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
  parms  ~ dmnorm(MLE,FI)
  beta  <- parms[1]
  sd1   <- exp(parms[2])
  #df1  ~ dpois(1) T(1,)
  
}
LogitNorm.alg1 <- function(){
  # likelihood
  for(k in 1:K){
    # Hidden process
    for(i in 1:n){
      X[i,k] ~ dt(0,1, df1)
    }
    
    # Observation process
    for(i in 1:n){
      
      W[i,k] ~ dnorm(X[i,k], 1/pow(sd1,2))
      Y[i,k] ~ dbin(1/(1+exp(-beta*X[i,k])),1)
    }    
  }
  
  # Priors on model parameters: 
  beta ~ dnorm(0,1/100)
  #sd1  ~ dlnorm(0,1)
  lsd1   ~ dnorm(0,1/100)
  sd1 <- exp(lsd1)
  #df1  ~ dpois(1) T(1,)
  
}




# #############  PROBIT-Normal Error in variables model
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
  beta  ~ dnorm(0,1/100)
  #sd1   ~ dlnorm(0,1)
  lsd1   ~ dnorm(0,1/100)
  sd1 <- exp(lsd1)
  #df1  ~ dpois(1) T(1,)
  
}
ProbitNorm.alg2 <- function(){
  
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
  parms  ~ dmnorm(MLE,FI)
  beta  <- parms[1]
  sd1   <- exp(parms[2])
  #df1  ~ dpois(1) T(1,)
}
ProbitNorm.alg1 <- function(){
  # likelihood
  for(k in 1:K){
    # Hidden process
    for(i in 1:n){
      X[i,k] ~ dt(0,1, df1)
    }
    
    # Observation process
    for(i in 1:n){
      
      W[i,k] ~ dnorm(X[i,k], 1/pow(sd1,2))
      Y[i,k] ~ dbin(pnorm((beta*X[i,k]),0,1) , 1)
    }    
  }
  
  # Priors on model parameters: 
  beta ~ dnorm(0,1/100)
  # sd1  ~ dlnorm(0,1)
  lsd1   ~ dnorm(0,1/100)
  sd1 <- exp(lsd1)
  #df1  ~ dpois(1) T(1,)
  
}







