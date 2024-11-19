#### ----------------------------------------------------------
#### TMB models
#### ----------------------------------------------------------


### Automatically set wd to look for .cpp files
#setwd( getSrcDirectory( function(){} )[1] ) ## works in r but not R studio
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
require(TMB)



##  Model 1: bivariate normal error variables model
##   - the observation process Y ~ N( mu= B*X, known_sigma1)
##   -                         W = X + eta, where eta ~ N( 0, known_sigma2 ) 
##   - X ~ T_df(0), is the hidden (latent) process
##   * We only observe (Y, W) and we are trying to estimate the hidden state Z
##   * The parameter of interest is beta. The observed variable is a vector. 

EV_bivnorm <- "
  #include <TMB.hpp>
  #include <stdlib.h>
  
  template<class Type>
  Type objective_function<Type>::operator() () {
    // data:
    DATA_VECTOR(y);
    DATA_VECTOR(w);
    int n = y.size(); // get time series length
    
    // parameters:
    PARAMETER(df1);      // df for t
    PARAMETER(lnsd1);      // standard deviation for w
    PARAMETER(lnsd2);      // standard deviation for y
    PARAMETER(beta); 
    PARAMETER_VECTOR(X); // unobserved state vector
    
    // parameter transformations: 
    Type sd1 = exp(lnsd1);
    Type sd2 = exp(lnsd2);

    // * Likelihood setup and computation * 
    Type nll = 0.0; // initialize negative log likelihood
    
    // process model:
    for(int i = 0; i < n; i++){
      nll -= dt( X(i), df1, true); // ncp=0 for TMB dt
    }
    
    // observation model:
    for (int i = 0; i < n; i++) {
      nll -= dnorm(w(i), X(i), sd1, true) + dnorm(y(i), X(i)*beta, sd2, true);
    }
    
    return nll;
}"
write(EV_bivnorm, file = "EV_bivnorm.cpp")
model_name <- "EV_bivnorm"
if(!file.exists(paste0(model_name, ".o"))) {
  compile("EV_bivnorm.cpp")
}
dyn.load(dynlib("EV_bivnorm"))



##  Model 2: Poisson-Normal error variables model
##   - the observation process Y ~ pois( lambda = exp(B*X))
##   -                         W = X + eta, where eta ~ N( 0, known_sigma ) 
##   - X ~ T_df(0), is the hidden (latent) process
##   * We only observe (Y, W) and we are trying to estimate the hidden state Z
##   * The parameter of interest is beta. The observed variable is a vector. 

EV_normpois <- "
  #include <TMB.hpp>
  
  template<class Type>
  Type objective_function<Type>::operator() () {
    // data:
    DATA_VECTOR(y);
    DATA_VECTOR(w);
    int n = y.size(); // get time series length
    
    // parameters:
    PARAMETER(beta); 
    PARAMETER(df1);      // df for t
    PARAMETER(lnsd1);      // standard deviation for y
    PARAMETER_VECTOR(X); // unobserved state vector
    
    // parameter transformations: 
    Type sd1 = exp(lnsd1);
    
    // * Likelihood setup and computation * 
    Type nll = 0.0; // initialize negative log likelihood
    
    // process model:
    for(int i = 0; i < n; i++){
      nll -= dt(X(i), df1, true); // ncp=0 for TMB dt
    }
    
    // poisson observation model:
    for (int i = 0; i < n; i++) {
      nll -= dnorm(w(i), X(i), sd1, true) + dpois(y(i), exp(beta*X(i)), true);
    }
    
    return nll;

}"

### Compile the model
write(EV_normpois, file = "EV_normpois.cpp")
model_name <- "EV_normpois"
if(!file.exists(paste0(model_name, ".o"))) {
  compile("EV_normpois.cpp")
}
dyn.load(dynlib("EV_normpois"))





##  Model 3: Logit-Normal error variables model
##   - the observation process Y ~ Binomial( size=1, Prob = logit(B*X))
##   -                         W = X + eta, where eta ~ N( 0, known_sigma ) 
##   - X ~ T_df(0), is the hidden (latent) process
##   * We only observe (Y, W) and we are trying to estimate the hidden state Z
##   * The parameter of interest is beta. The observed variable is a vector. 

EV_normlogit <- "
  #include <TMB.hpp>
  
  template<class Type>
  Type objective_function<Type>::operator() () {
    // data:
    DATA_VECTOR(y);
    DATA_VECTOR(w);
    int n = y.size();  // get time series length
    
    // parameters:
    PARAMETER(beta); // beta value to profile over
    PARAMETER(df1);  // df for t
    PARAMETER(lnsd1);  // known standard deviation for w
    PARAMETER_VECTOR(X);       // unobserved state vector
    
        
    // parameter transformations: 
    Type sd1 = exp(lnsd1);
    
    // * Likelihood setup and computation * 
    Type nll = 0.0; // initialize negative log likelihood
    
    // process model:
    for(int i = 0; i < n; i++){
      nll -= dt(X[i], df1, true); // ncp=0 for TMB dt
    }
    
    // poisson observation model:
    for (int i = 0; i < n; i++) {
      nll -= dnorm(w[i], X[i], sd1, true) + dbinom(y[i], Type(1.0), invlogit(beta*X[i]), true);
    }
    
    return nll;

}"

### Compile the model
write(EV_normlogit, file = "EV_normlogit.cpp")
model_name <- "EV_normlogit"
if(!file.exists(paste0(model_name, ".o"))) {
  compile("EV_normlogit.cpp")
}
dyn.load(dynlib("EV_normlogit"))



##  Model 4: Probit-Normal error variables model
##   - the observation process Y ~ Binomial( size=1, Prob = pnorm(B*X))
##   -                         W = X + eta, where eta ~ N( 0, known_sigma ) 
##   - X ~ T_df(0), is the hidden (latent) process
##   * We only observe (Y, W) and we are trying to estimate the hidden state Z
##   * The parameter of interest is beta. The observed variable is a vector. 

EV_normprobit <- "
  #include <TMB.hpp>
  
  template<class Type>
  Type objective_function<Type>::operator() () {
    // data:
    DATA_VECTOR(y);
    DATA_VECTOR(w);
    int n = y.size();  // get time series length
    
    // parameters:
    PARAMETER(beta); // beta value to profile over
    PARAMETER(df1);  // df for t
    PARAMETER(lnsd1);  // known standard deviation for w
    PARAMETER_VECTOR(X);       // unobserved state vector
    
        
    // parameter transformations: 
    Type sd1 = exp(lnsd1);
    
    // * Likelihood setup and computation * 
    Type nll = 0.0; // initialize negative log likelihood
    
    // process model:
    for(int i = 0; i < n; i++){
      nll -= dt(X[i], df1, true); // ncp=0 for TMB dt
    }
    
    // poisson observation model:
    for (int i = 0; i < n; i++) {
      nll -= dnorm(w[i], X[i], sd1, true) + dbinom(y[i], Type(1.0), pnorm(beta*X[i]), true);
    }
    
    return nll;

}"

### Compile the model
write(EV_normprobit, file = "EV_normprobit.cpp")
model_name <- "EV_normprobit"
if(!file.exists(paste0(model_name, ".o"))) {
  compile("EV_normprobit.cpp")
}
dyn.load(dynlib("EV_normprobit"))



##  Model 5: Non-identifiable Poisson model
##   - the observation process Y ~ Poisson( lambda=X^2)
##   - the hidden (latent process) for X ~ N( mu, sigma ) 
##   * We only observe Y and we are trying to estimate the hidden state X
##   * The parameter of interest is beta= B. The observed variable is a vector. 

normsq_pois <- "
  #include <TMB.hpp>
  
  template<class Type>
  Type objective_function<Type>::operator() () {
    // data:
    DATA_VECTOR(y);
    int n = y.size(); // get time series length
    
    // parameters:
    PARAMETER(mu);
    PARAMETER(lnsd1);      // standard deviation for y
    PARAMETER_VECTOR(X);       // unobserved state vector
    
        
    // parameter transformations: 
    Type sd1 = exp(lnsd1);
    
    // * Likelihood setup and computation * 
    Type nll = 0.0; // initialize negative log likelihood
    
    // process model:
    for(int i = 1; i < n; i++){
      nll -= dnorm(X(i), mu, sd1, true);
    }
    
    // poisson observation model:
    for (int i = 0; i < n; i++) {
      nll -= dpois(y(i), pow(X(i), 2), true);
    }
    return nll;
}"

write(normsq_pois, file = "normsq_pois.cpp")
model_name <- "normsq_pois"
if(!file.exists(paste0(model_name, ".o"))) {
  compile("normsq_pois.cpp")
}
dyn.load(dynlib("normsq_pois"))




## Model 6: Errors in variables normal-normal linear regression with non-identifiable variance parameters
## Yi|Xi ∼ N(Xi,σ2)
## Xi    ∼ N(μ,τ2)
nonidentifiable_variance <- "
  #include <TMB.hpp>
  #include <stdlib.h>
  
  template<class Type>
  Type objective_function<Type>::operator() () {
    
    // data:
    DATA_VECTOR(y);
    int n = y.size(); // get time series length
    
    // parameters:
    PARAMETER(mu1);         // mean for x
    PARAMETER(lnsd1);         // standard deviation for x
    PARAMETER(lntau1);        // standard deviation for y
    PARAMETER_VECTOR(X1);   // unobserved state vector for X1

    // parameter transformations: 
    Type tau1    = exp(lntau1);
    Type sd1     = exp(lnsd1);
    Type sdtausq = pow(sd1,2) + pow(tau1,2);

    // * Likelihood setup and computation * 
    Type nll = 0.0; // initialize negative log likelihood
    
    // process model:
    for(int i = 0; i < n; i++){ nll -= dnorm( X1(i), mu1, sd1, true); }
    
    // observation model:
    for (int i = 0; i < n; i++) { nll -= dnorm(y(i), X1(i), tau1, true); }
    
    ADREPORT(sdtausq);
    return nll;
}"

write(nonidentifiable_variance, file = "nonidentifiable_variance.cpp")
model_name <- "nonidentifiable_variance"
if(!file.exists(paste0(model_name, ".o"))) {
  TMB::compile("nonidentifiable_variance.cpp")
}
dyn.load(dynlib("nonidentifiable_variance"))

