# Generating data doubling MCMC samples for profile likelihood
#'
#' @param Pmodel.fn This is the same model function as is used for the 
#'                  computation of the MLE except the prior distribution 
#'                  is the asymptotic distribution of the MLE. 
#' @param dat.p Change the input data list to include MLE and 
#'              Fisher Information matrix as MuP and SigmaP.
#' @param params  Parameter vector to be monitored
#' @param ...     Optional arguments n.adapt, n.iter, n.chains,n.update
#' @param unchanged Nodes in the model that remain unchanged. 
#'                  For example, n the sample size.
#'
#' @return list with elements matrix of MCMC outpute M1, M2 
#'         and diagnostics Rhat.out
#' @export
#'
DDSample_fn = function(Pmodel.fn,dat.p,params,unchanged,...){
  
  model1.fit = dclone::dc.fit(dat.p,params,Pmodel.fn,n.clones=1,multiply="K",unchanged=unchanged,...)
  model2.fit = dclone::dc.fit(dat.p,params,Pmodel.fn,n.clones=2,multiply="K",unchanged=unchanged,...)
  
  # Extract the MCMC output
  M1 = as.matrix(model1.fit)
  M2 = as.matrix(model2.fit)
  
  # MCMC convergence diagnostics
  Rhat.out = rbind(dclone::dcdiag(model1.fit),dclone::dcdiag(model2.fit))
  out = list(model1.fit=model1.fit, model2.fit=model2.fit, 
             M1=M1, M2=M2, Rhat.out=Rhat.out)
  return(out)
}



# Profile likelihood function using the Logistic regression method via DD
#'         
#' @param psi.fn This is a user supplied function to compute the parameter of
#'               of interest from the canonical parameters.              
#' @param M1 This is a matrix of MCMC samples from the posterior with 1 clone
#' @param M2 This is a matrix of MCMC samples from the posterior with 2 clones
#' @param grid.length This is the number of grid points at which PL is computed.
#'                    The range is decided based on the range of psi values. 
#'                    Default value is 100.   
#' @param psi.grid These are the grid points at which PL is computed. Default is 
#'                 NULL. User can supply this in place of the grid length.
#' @param method  Choose between 'q'uartic, 'c'oncave or 's'plines as link function.
#'                Concave fit using 'scam' function in the R package 'mgcv'.
#' @param ...     Optional arguments in the psi.fn if needed.
#'
#' @return        list containing two objects. Profile.like is a matrix with first 
#'                column psi.grid, 2nd column is the relative profile loglikelihood value.
#'                Second object is the profile MLE of the parameter of interest and
#'                associated confidence interval at the prescribed alpha value. This
#'                can be used to make sure the procedure worked properly. Profile MLE
#'                should be the same as the one obtained by transforming the canonical
#'                MLE directly. Default for confidence interval is 90%

profile_fn = function(psi.fn,M1,M2,grid.length=100,
                      psi.grid=NULL,method='c',CIalpha=0.90,...){
  
  # Apply the psi function on the MCMC output to get the psi transformation
  psi1 = psi.fn(M1,...)
  psi2 = psi.fn(M2,...)
  M12 = data.frame(Y.m=c(rep(0,length(psi1)),rep(1,length(psi2))),X = c(psi1,psi2))
  
  if(is.null(psi.grid)){
    psi.range = range(M12$X)
    psi.grid = seq(psi.range[1],psi.range[2],length.out=grid.length)}
  
  if (method == 'q'){
    # We fit a fourth degree polynomial (orthogonal) regression
    glm.fit12.2 = glm(Y.m ~ poly(X,2),family=binomial,data=M12)
    glm.fit12.3 = glm(Y.m ~ poly(X,3),family=binomial,data=M12)
    glm.fit12.4 = glm(Y.m ~ poly(X,4),family=binomial,data=M12)
    tmp = AIC(glm.fit12.2,glm.fit12.3,glm.fit12.4)
    
    # Choose the model fit object that has minimum AIC
    if (AIC(glm.fit12.2) == min(tmp[,2])){glm.fit12=glm.fit12.2}
    if (AIC(glm.fit12.3) == min(tmp[,2])){glm.fit12=glm.fit12.3}
    if (AIC(glm.fit12.4) == min(tmp[,2])){glm.fit12=glm.fit12.4}
  }
  
  if (method == 'c'){ 
    # Fit a concave function
    glm.fit12= scam::scam(Y.m ~ s(X,bs="cv"),family=binomial,data=M12)}
  
  if (method == 's'){ 
    # We will fit using general splines
    glm.fit12= mgcv::gam(Y.m ~ s(X),family=binomial,data=M12)
  }
  
  # Predict on the assigned grid values (psi.grid)
  newdata = data.frame(X = psi.grid)
  profile.like = predict(glm.fit12,newdata) 
  profile.out = data.frame(psi.grid=psi.grid,profile.like=profile.like)
  
  # Compute the PL estimate along with its LPL value
  tmp = dplyr::arrange(profile.out,profile.like)
  Profile.MLE = tmp[nrow(tmp),]
  
  # Standardize the logPL by subtracting the value at the maximum
  profile.out$profile.like = profile.out$profile.like - as.numeric(Profile.MLE[2])
  
  # Compute the confidence interval
  cutoff = -qchisq(CIalpha,1)/2
  index = which(profile.out$profile.like >= cutoff)
  L_index = min(index)
  U_index = max(index)
  CI = c(profile.out[L_index,1],profile.out[U_index,1])
  out = list(Profile.out=profile.out, Profile.MLE=Profile.MLE,Profile.CI = c(CIalpha,CI))
  
  return(out)
}



# Function wrapper for profile_fn and DDSample_fn
#' @data
#' 
#' @params
#' 



DD_profile_fun <- function(data=data, known.mle=NULL, 
                           params=params, psi.fn=psi.fn, 
                           profile.par=profile.par,
                           model.mle=model.mle, unchanged.mle=unchanged.mle,
                           n.clones=n.clones, psi.grid=NULL,
                           model.alg2=model.alg2, unchanged.alg2=unchanged.alg2,
                           model.alg1=model.alg1, unchanged.alg1=unchanged.alg1,
                           CIalpha=0.95, MCMCset=MCMCset, ...){
  
  if( is.null(known.mle) ){
    #### Data cloning to compute MLE
    DCfit = dclone::dc.fit(data=data, params=params, model=model.mle,
                           n.clones=n.clones, multiply="K", unchanged= unchanged.mle,
                           n.iter   = MCMCset$n.iter.dc,
                           n.chains = MCMCset$n.chains.dc,
                           n.adapt  = MCMCset$n.adapt,
                           n.update = MCMCset$n.burn, ...) 
    MLE.dc = coef(DCfit)
    FI.dc  = solve(vcov(DCfit))
    DC.MLE = list(MLE=MLE.dc, FI=FI.dc, fitted=DCfit)
  }else{
    DC.MLE=list(MLE=known.mle$MLE, FI=known.mle$FI)
    MLE.dc = known.mle$MLE
    FI.dc  = known.mle$FI
  }
  
  #### For algorithm 2
  ProfileSample_A2 = DDSample_fn(Pmodel.fn=model.alg2, 
                                 dat.p=append(data, 
                                              list(MLE=MLE.dc,FI=FI.dc)),
                                 params=profile.par, 
                                 unchanged=unchanged.alg2,
                                 n.iter   = MCMCset$n.iter.dd,
                                 n.chains = MCMCset$n.chains.dd,
                                 n.adapt  = MCMCset$n.adapt,
                                 n.update = MCMCset$n.burn)
  
  DD.profile.q = profile_fn(psi.fn=psi.fn, 
                            M1=ProfileSample_A2$M1, 
                            M2=ProfileSample_A2$M2, 
                            method="c",CIalpha=CIalpha)
  if(is.null(psi.grid)){psi.grid = DD.profile.q$Profile.out$psi.grid}
  
  
  #### For algorithm 1
  ProfileSample_A1 = DDSample_fn(Pmodel.fn=model.alg1,
                                 dat.p=data, params=profile.par, 
                                 unchanged = unchanged.alg1,
                                 n.iter   = MCMCset$n.iter.dd,
                                 n.chains = MCMCset$n.chains.dd,
                                 n.adapt  = MCMCset$n.adapt,
                                 n.update = MCMCset$n.burn)
  
  DD.profile.psi   = profile_fn(psi.fn=psi.fn,
                                M1=ProfileSample_A1$M1,
                                M2=ProfileSample_A1$M2, 
                                method="c",psi.grid=psi.grid,CIalpha=CIalpha)
  
  #### Return an object
  return( list(DC.MLE=DC.MLE, 
               alg1=DD.profile.psi,
               alg2=DD.profile.q, 
               ProfileSample_A1=ProfileSample_A1, 
               ProfileSample_A2=ProfileSample_A2,
               DD.profile.psi=DD.profile.psi, 
               DD.profile.q=DD.profile.q) )
  
}





# Function wrapper by model for DD_profile_fun
#' @

DD_wrapper <- function(profile.par=profile.par, known.params=known.params, known.mle=NULL,
                       model="Bivar.Norm", sim.data=sim.data, psi.grid=NULL, MCMCset=MCMCset, 
                       CIalpha=0.95, n.clones=n.clones, ...){
  # function that takes in the parameter of interest and returns it
  psi.fn <- function(M){ M }
  
  switch(model,
         "Bivar.Norm"={
           # unchanged
           unchanged.mle  <- c("n","df1")
           unchanged.alg1 <- c("n","df1")
           unchanged.alg2 <- c("n","df1","MLE","FI")
           # model parameters
           # params <- c("beta","sd1","sd2")
           params <- c("beta","lsd1","lsd2")
           # prepare data
           df.mod   <- list(K=1, n=dim(sim.data)[1], df1=known.params[['df1']],
                            W=dclone::dcdim(data.matrix(sim.data[,1])),
                            Y=dclone::dcdim(data.matrix(sim.data[,2])))
           out <- DD_profile_fun(data=df.mod, params=params, known.mle=known.mle, 
                                 psi.fn=psi.fn, profile.par=profile.par,
                                 model.mle=BivarNormDC.MLE, 
                                 model.alg2=BivarNorm.alg2, 
                                 model.alg1=BivarNorm.alg1, 
                                 n.clones=n.clones, psi.grid=psi.grid, MCMCset=MCMCset,
                                 unchanged.mle=unchanged.mle, 
                                 unchanged.alg1=unchanged.alg1,
                                 unchanged.alg2=unchanged.alg2, 
                                 CIalpha=CIalpha)
         },
         "Pois.Norm"={
           # unchanged
           unchanged.mle  <- c("n","df1")
           unchanged.alg1 <- c("n","df1")
           unchanged.alg2 <- c("n","df1","MLE","FI")
           # model parameters
           #params=c("beta","sd1")
           params=c("beta","lsd1")
           # prepare data
           df.mod   <- list(K=1, n=dim(sim.data)[1], df1=known.params[['df1']],
                            W=dclone::dcdim(data.matrix(sim.data[,1])),
                            Y=dclone::dcdim(data.matrix(sim.data[,2])))
           out <- DD_profile_fun(data=df.mod, params=params, psi.fn=psi.fn, 
                                 profile.par=profile.par, known.mle=known.mle, 
                                 model.mle=PoisNormDC.MLE, 
                                 model.alg2=PoisNorm.alg2,
                                 model.alg1=PoisNorm.alg1, 
                                 n.clones=n.clones, psi.grid=psi.grid, MCMCset=MCMCset,
                                 unchanged.mle=unchanged.mle, 
                                 unchanged.alg1=unchanged.alg1,
                                 unchanged.alg2=unchanged.alg2, 
                                 CIalpha=CIalpha)
         },
         "Normsq.Pois"={
           # unchanged
           unchanged.mle  <- c("n")
           unchanged.alg1 <- c("n")
           unchanged.alg2 <- c("n","MLE","FI")
           # model parameters
           params=c("mu","lsd1")
           #params=c("mu","sd1")
           # prepare data
           df.mod   <- list(K=1, n=dim(sim.data)[1],
                            Y=dclone::dcdim(data.matrix(sim.data)))
           out <- DD_profile_fun(data=df.mod, params=params, psi.fn=psi.fn, 
                                 profile.par=profile.par, known.mle=known.mle, 
                                 model.mle=NormsqPois.MLE, 
                                 model.alg2=NormsqPois.alg2, 
                                 model.alg1=NormsqPois.alg1, 
                                 n.clones=n.clones, psi.grid=psi.grid, MCMCset=MCMCset,
                                 unchanged.mle=unchanged.mle, 
                                 unchanged.alg1=unchanged.alg1,
                                 unchanged.alg2=unchanged.alg2, 
                                 CIalpha=CIalpha)
           
         },
         "Logit.Norm"={
           # unchanged
           unchanged.mle  <- c("n","df1")
           unchanged.alg1 <- c("n","df1")
           unchanged.alg2 <- c("n","df1","MLE","FI")
           # model parameters
           #params=c("beta","sd1")
           params=c("beta","lsd1")
           # prepare data
           df.mod   <- list(K=1, n=dim(sim.data)[1], df1=known.params[['df1']],
                            W=dclone::dcdim(data.matrix(sim.data[,1])),
                            Y=dclone::dcdim(data.matrix(sim.data[,2])))
           out <- DD_profile_fun(data=df.mod, params=params, psi.fn=psi.fn, 
                                 profile.par=profile.par, known.mle=known.mle, 
                                 model.mle=LogitNormDC.MLE, 
                                 model.alg2=LogitNorm.alg2, 
                                 model.alg1=LogitNorm.alg1, 
                                 n.clones=n.clones, psi.grid=psi.grid, MCMCset=MCMCset,
                                 unchanged.mle=unchanged.mle, 
                                 unchanged.alg1=unchanged.alg1,
                                 unchanged.alg2=unchanged.alg2, 
                                 CIalpha=CIalpha)
           
         },
         
         "Probit.Norm"={
           # unchanged
           unchanged.mle  <- c("n","df1")
           unchanged.alg1 <- c("n","df1")
           unchanged.alg2 <- c("n","df1","MLE","FI")
           # model parameters
           #params=c("beta","sd1")
           params=c("beta","lsd1")
           # prepare data
           df.mod   <- list(K=1, n=dim(sim.data)[1], df1=known.params[['df1']],
                            W=dclone::dcdim(data.matrix(sim.data[,1])),
                            Y=dclone::dcdim(data.matrix(sim.data[,2])))
           out <- DD_profile_fun(data=df.mod, params=params, psi.fn=psi.fn, 
                                 profile.par=profile.par, known.mle=known.mle, 
                                 model.mle=ProbitNormDC.MLE, 
                                 model.alg2=ProbitNorm.alg2, 
                                 model.alg1=ProbitNorm.alg1, 
                                 n.clones=n.clones, psi.grid=psi.grid, MCMCset=MCMCset,
                                 unchanged.mle=unchanged.mle, 
                                 unchanged.alg1=unchanged.alg1,
                                 unchanged.alg2=unchanged.alg2, 
                                 CIalpha=CIalpha)
         },         
  )
  return(out)
  
}


