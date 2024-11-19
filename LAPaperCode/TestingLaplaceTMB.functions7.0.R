#### ----------------------------------------------------------
#### Single- and multi-parameter TMB models
#### ----------------------------------------------------------

### Automatically set wd to look for .cpp files
#setwd( getSrcDirectory( function(){} )[1] ) ## works in r but not R studio
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
require(TMB)



# Profile likelihood calculator using TMB: three modalities:
#.   "usemanual", "useconfint", "root finding".
#.    usemanual = slowest but manual Chi-square confints for profile likes, stable
#.    useconfint = fast but fails sometimes and who knows why
#.    useroot = fastest but will not compute relative profile likelihoods
tmbprofileCI <- function(obj, par="beta", prec=0.95, ystep=0.01, ytol=3.5, use="useconfint"){
  
  if(use %in% "usemanual"){
    ### Compute the profile CI
    profilecurveTMB <- tmbprofile(obj=obj, name=par, trace=F, ystep=ystep)
    ### Check for any duplicated values
    which.dup       <- which( duplicated(profilecurveTMB) )
    profilecurveTMB <- na.omit( profilecurveTMB[-which.dup,] )
    ### Define the CI to compute
    crit.LR <- exp(-qchisq(prec,1)/2); 
    ### Compute the relative profile CI
    profilecurveTMB$rel.like <- exp(-profilecurveTMB$value-max(-profilecurveTMB$value))
    ### Compute the desires profile CI 
    prof <- profilecurveTMB[profilecurveTMB$rel.like > crit.LR,]
    return(list(data=profilecurveTMB, prof=prof, ci=range(prof[[par]])))
  }
  
  if(use %in% "useconfint"){
    
    ### Compute the profile CI
    profilecurveTMB <- tmbprofile(obj=obj, name=par, trace=F, ytol=ytol)
    #### The ytol value controls the range of the likelihood values 
    #### so increase this by adding 0.0,...,9 to the prec value
    #### We need ytol to be greater than qnorm(0.99)
    
    ### Check for any duplicated values
    which.dup       <- which( duplicated(profilecurveTMB) )
    profilecurveTMB <- na.omit( profilecurveTMB[-which.dup,] )
    
    ### Compute the desires profile CI 
    prof <- confint(profilecurveTMB, level = prec) 
    ### Compute the relative profile CI
    profilecurveTMB$rel.like <- exp(-profilecurveTMB$value-max(-profilecurveTMB$value))
    
    return(list(data=profilecurveTMB, prof=prof, ci=prof))
  }
  
  if(use %in% "useroot"){
    ### Compute the profile CI using root-finding 
    prof <- tmbroot(obj=obj, name=par, trace=F, target = (0.5 * qchisq(prec, sim.data = 1)))
    return(list(ci=prof))
  }
  
}





### function to run the TMB model with single parameters (nothing unknown)
run.tmb <- function( profile.par=NULL, profile.par.name=NULL, 
                     known.parm=NULL, sim.data, model="Bivar.Norm", init=NULL,
                     method.ci="usemanual", ytol=3, use.optim=TRUE, single.par.init=0.01,
                     search.initial=F, n.search=100){
  
  # Ensure to be a data.frame to call objects from using $
  if(!is.data.frame(sim.data)) sim.data <- as.data.frame(sim.data)
  
  # include packages if not loaded
  require(Hmisc); require(TMB);#require(TMBhelper);
  
  
  if( is.null(known.parm) ){
    
    # run cases where there are no knowns for the multi-parameter cases
    switch(model,
           "Bivar.Norm"={
             # function to search parameter space for the best initial values
             if( search.initial == T){
               
               model.param  <- c("lnsd1","lnsd2","beta")
               unknown.pars <- list()
               for( var_name in model.param ){ unknown.pars[[var_name]] <- init[["par"]] }
               
               ### Set up the data and par list to read into TMB
               data <- append(list(y=sim.data$Y, w=sim.data$W),known.parm)
               par  <- append(list(X=rep(init[["hidden"]], dim(sim.data)[1]), 
                                   df1=init[["df1"]]), unknown.pars)
               
               
               ### Search for the "best" initial value of a grid of values
               obj <- MakeADFun(data, par, random = "X", DLL = "EV_bivnorm", silent = T, 
                                map = list(df1=factor(NA)))
               start.beta  <- matrix( data=cbind(lnsd1=runif(n.search, -10, 10 ),
                                                 lnsd2=runif(n.search, -10, 10 ),
                                                 beta=runif(n.search, -10, 10 )),
                                      ncol = length(unknown.pars))
               colnames(start.beta) <- model.param
               
               # if none of the initial starting values work revert back to the initial values in par
               LLvalues    <- na.omit( apply(start.beta, 1, function(i){ obj$fn(x=i)}) )
               if( length(LLvalues)>0 ){
                 best.start  <- as.list(start.beta[which(na.omit(LLvalues) == min(na.omit(LLvalues))),])
               }else{
                 best.start  <- unknown.pars
               }
               
               ### Using min LLvalues fit the TMB model
               par <- append(list(X=rep(init[['hidden']], dim(sim.data)[1]), df1=init[["df1"]]), best.start)
               
               obj <- MakeADFun(data, par, random = "X", DLL = "EV_bivnorm", silent = T, map = list(df1=factor(NA)))
               
             }else{
               
               model.param  <- c("lnsd1","lnsd2","beta")
               unknown.pars <- list()
               for( var_name in model.param ){ unknown.pars[[var_name]] <- init[["par"]] }
               
               ### Set up the data and par list to read into TMB
               data <- append(list(y=sim.data$Y, w=sim.data$W),known.parm)
               par  <- append(list(X=rep(init[["hidden"]], dim(sim.data)[1]), df1=init[["df1"]]),unknown.pars)
               
               ### Fit the TMB model
               obj  <- MakeADFun(data, par, random = "X", DLL = "EV_bivnorm", silent = T, map = list(df1=factor(NA)))
               
             }
           },
           
           "Pois.Norm"={
             # function to search parameter space for the best initial values
             if( search.initial == T){
               
               model.param  <- c("lnsd1","beta")
               unknown.pars <- list()
               for( var_name in model.param ){ unknown.pars[[var_name]] <- init[["par"]] }
               
               ### Set up the data and par list to read into TMB
               data <- append(list(y=sim.data$Y, w=sim.data$W),known.parm)
               par  <- append(list(X=rep(init[["hidden"]], dim(sim.data)[1]), df1=init[["df1"]]),unknown.pars)
               
               
               ### Search for the "best" initial value of a grid of values
               obj <- MakeADFun(data, par, random = "X", DLL = "EV_normpois", silent = T, map = list(df1=factor(NA)))
               start.beta  <- matrix( data=cbind(lnsd1=runif(n.search, -10, 10 ),
                                                 beta=runif(n.search, -10, 10 )),
                                      ncol = length(unknown.pars))
               colnames(start.beta) <- model.param
               
               # if none of the initial starting values work revert back to the initial values in par
               LLvalues    <- na.omit( apply(start.beta, 1, function(i){ obj$fn(x=i)}) )
               if( length(LLvalues)>0 ){
                 best.start  <- as.list(start.beta[which(na.omit(LLvalues) == min(na.omit(LLvalues))),])
               }else{
                 best.start  <- unknown.pars
               }
               
               ### Using min LLvalues fit the TMB model
               par <- append(list(X=rep(init[['hidden']], dim(sim.data)[1]), df1=init[["df1"]]), best.start)
               
               obj <- MakeADFun(data, par, random = "X", DLL = "EV_normpois", silent = T, map = list(df1=factor(NA)))
               
               
             }else{
               model.param  <- c("lnsd1","beta")
               unknown.pars <- list()
               for( var_name in model.param ){ unknown.pars[[var_name]] <- init[["par"]] }
               
               ### Set up the data and par list to read into TMB
               data <- append(list(y=sim.data$Y, w=sim.data$W),known.parm)
               par  <- append(list(X=rep(init[["hidden"]], dim(sim.data)[1]), df1=init[["df1"]]),unknown.pars)
               
               ### Fit the TMB model
               obj  <- MakeADFun(data, par, random = "X", DLL = "EV_normpois", silent = T, map = list(df1=factor(NA)))
             }
           },
           
           "Logit.Norm"={
             if( search.initial == T){
               
               model.param  <- c("lnsd1","beta")
               unknown.pars <- list()
               for( var_name in model.param ){ unknown.pars[[var_name]] <- init[["par"]] }
               
               ### Set up the data and par list to read into TMB
               data <- append(list(y=sim.data$Y, w=sim.data$W), known.parm)
               par  <- append(list(X=rep(init[["hidden"]], dim(sim.data)[1]), df1=init[["df1"]]),unknown.pars)
               
               
               ### Search for the "best" initial value of a grid of values
               obj <- MakeADFun(data, par, random = "X", DLL = "EV_normlogit", silent = T, map = list(df1=factor(NA)))
               start.beta  <- matrix( data=cbind(lnsd1=runif(n.search, -10, 10 ),
                                                 beta=runif(n.search, -10, 10 )),
                                      ncol = length(unknown.pars))
               colnames(start.beta) <- model.param
               
               # if none of the initial starting values work revert back to the initial values in par
               LLvalues    <- na.omit( apply(start.beta, 1, function(i){ obj$fn(x=i)}) )
               if( length(LLvalues)>0 ){
                 best.start  <- as.list(start.beta[which(na.omit(LLvalues) == min(na.omit(LLvalues))),])
               }else{
                 best.start  <- unknown.pars
               }
               
               ### Using min LLvalues fit the TMB model
               par  <- append(list(X=rep(init[['hidden']], dim(sim.data)[1]), df1=init[["df1"]]), best.start)
               
               obj  <- MakeADFun(data, par, random = "X", DLL = "EV_normlogit", silent = T, map = list(df1=factor(NA)))
               
               
             }else{
               model.param  <- c("lnsd1","beta")
               unknown.pars <- list()
               for( var_name in model.param ){ unknown.pars[[var_name]] <- init[["par"]] }
               
               ### Set up the data and par list to read into TMB
               data <- append(list(y=sim.data$Y, w=sim.data$W),known.parm)
               par  <- append(list(X=rep(init[["hidden"]], dim(sim.data)[1]), df1=init[["df1"]]),unknown.pars)
               
               ### Fit the TMB model
               obj  <- MakeADFun(data, par, random = "X", DLL = "EV_normlogit", silent = T, map = list(df1=factor(NA)))
               
               
             }
           },
           
           "Probit.Norm"={
             if( search.initial == T){
               
               model.param  <- c("lnsd1","beta")
               unknown.pars <- list()
               for( var_name in model.param ){ unknown.pars[[var_name]] <- init[["par"]] }
               
               ### Set up the data and par list to read into TMB
               data <- append(list(y=sim.data$Y, w=sim.data$W), known.parm)
               par  <- append(list(X=rep(init[["hidden"]], dim(sim.data)[1]), df1=init[["df1"]]),unknown.pars)
               
               
               ### Search for the "best" initial value of a grid of values
               obj <- MakeADFun(data, par, random = "X", DLL = "EV_normprobit", silent = T, map = list(df1=factor(NA)))
               start.beta  <- matrix( data=cbind(lnsd1=runif(n.search, -10, 10 ),
                                                 beta=runif(n.search, -10, 10 )),
                                      ncol = length(unknown.pars))
               colnames(start.beta) <- model.param
               
               # if none of the initial starting values work revert back to the initial values in par
               LLvalues    <- na.omit( apply(start.beta, 1, function(i){ obj$fn(x=i)}) )
               if( length(LLvalues)>0 ){
                 best.start  <- as.list(start.beta[which(na.omit(LLvalues) == min(na.omit(LLvalues))),])
               }else{
                 best.start  <- unknown.pars
               }
               
               ### Using min LLvalues fit the TMB model
               par  <- append(list(X=rep(init[['hidden']], dim(sim.data)[1]), df1=init[["df1"]]), best.start)
               
               obj  <- MakeADFun(data, par, random = "X", DLL = "EV_normprobit", silent = T, map = list(df1=factor(NA)))
               
             }else{
               model.param  <- c("lnsd1","beta")
               unknown.pars <- list()
               for( var_name in model.param ){ unknown.pars[[var_name]] <- init[["par"]] }
               
               ### Set up the data and par list to read into TMB
               data <- append(list(y=sim.data$Y, w=sim.data$W),known.parm)
               par  <- append(list(X=rep(init[["hidden"]], dim(sim.data)[1]), df1=init[["df1"]]),unknown.pars)
               
               ### Fit the TMB model
               obj  <- MakeADFun(data, par, random = "X", DLL = "EV_normprobit", silent = T, map = list(df1=factor(NA)))
               
             }
           },
           
           "Normsq.Pois"={
             if( search.initial == T){
               
               model.param  <- c("lnsd1","mu")
               unknown.pars <- list()
               for( var_name in model.param ){ unknown.pars[[var_name]] <- init[["par"]] }
               
               ### Set up the data and par list to read into TMB
               data <- append(list(y=sim.data$Y), known.parm)
               par  <- append(list(X=rep(init[["hidden"]], dim(sim.data)[1])),unknown.pars)
               
               
               ### Search for the "best" initial value of a grid of values
               obj <- MakeADFun(data, par, random = "X", DLL = "normsq_pois", silent = T)
               start.beta  <- matrix( data=cbind(lnsd1=runif(n.search, -10, 10 ),
                                                 mu=runif(n.search, 0.001, 10 )),
                                      ncol = length(unknown.pars))
               colnames(start.beta) <- model.param
               
               # if none of the initial starting values work revert back to the initial values in par
               LLvalues    <- na.omit( apply(start.beta, 1, function(i){ obj$fn(x=i)}) )
               if( length(LLvalues)>0 ){
                 best.start  <- as.list(start.beta[which(na.omit(LLvalues) == min(na.omit(LLvalues))),])
               }else{
                 best.start  <- unknown.pars
               }
               
               ### Using min LLvalues fit the TMB model
               par  <- append(list(X=rep(init[['hidden']], dim(sim.data)[1])), best.start)
               
               obj  <- MakeADFun(data, par, random = "X", DLL = "normsq_pois", silent = T)
               
               
             }else{
               model.param  <- c("lnsd1","mu")
               unknown.pars <- list()
               for( var_name in model.param ){ unknown.pars[[var_name]] <- init[["par"]] }
               
               ### Set up the data and par list to read into TMB
               data <- list(y=sim.data$Y)
               par  <- append(list(X=rep(init[["hidden"]], dim(sim.data)[1])),unknown.pars)
               
               ### Fit the TMB model
               obj  <- MakeADFun(data, par, random = "X", 
                                 DLL = "normsq_pois", silent = T)
               
             }
           },
           
           "Norm.nonident.var"={
             if( search.initial == T){
               
               model.param  <- c("mu1","lnsd1","lntau1")
               unknown.pars <- list()
               for( var_name in model.param ){ unknown.pars[[var_name]] <- init[["par"]] }
               
               ### Set up the data and par list to read into TMB
               data <- append(list(y=sim.data$Y), known.parm)
               par  <- append(list(X1=rep(init[["hidden"]], dim(sim.data)[1])),
                              unknown.pars)
               
               
               ### Search for the "best" initial value of a grid of values
               obj <- MakeADFun(data, par, random = c("X1"), DLL = "nonidentifiable_variance", silent = T)
               start.beta  <- matrix( data=cbind(mu1=runif(n.search, 0.1, 10 ),
                                                 lnsd1=runif(n.search, -10, 10 ),
                                                 lntau1=runif(n.search, -10, 10 )),
                                      ncol = length(unknown.pars))
               colnames(start.beta) <- model.param
               
               # if none of the initial starting values work revert back to the initial values in par
               LLvalues    <- na.omit( apply(start.beta, 1, function(i){ obj$fn(x=i)}) )
               if( length(LLvalues)>0 ){
                 best.start  <- as.list(start.beta[which(na.omit(LLvalues) == min(na.omit(LLvalues))),])
               }else{
                 best.start  <- unknown.pars
               }
               
               ### Using min LLvalues fit the TMB model
               par  <- append(list(X1=rep(init[["hidden"]], dim(sim.data)[1])), best.start)
               
               obj  <- MakeADFun(data, par, random = c("X1"), DLL = "nonidentifiable_variance", silent = T)
               
               
             }else{
               model.param  <- c("mu1","lnsd1","lntau1")
               unknown.pars <- list()
               for( var_name in model.param ){ unknown.pars[[var_name]] <- init[["par"]] }
               
               ### Set up the data and par list to read into TMB
               data <- append(list(y=sim.data$Y), known.parm)
               par  <- append(list(X1=rep(init[["hidden"]], dim(sim.data)[1])),
                              unknown.pars)
               
               ### Fit the TMB model
               obj  <- MakeADFun(data, par, random = c("X1"), DLL = "nonidentifiable_variance", silent = T)
               
             }
           }
    )
    
    # use optim or nlminb
    if(use.optim==TRUE){
      opt  <- suppressWarnings( optim(par=obj$par, fn=obj$fn, gr=obj$gr, hessian = TRUE, method = "BFGS") )
    }else{
      opt  <- suppressWarnings( nlminb(start=obj$par, objective=obj$fn, gradient=obj$gr) )
    }
    
    ### Compute profile CI and return model objects
    res.sd <- sdreport(obj, getJointPrecision = T, bias.correct = T)
    internal.checks( obj )
    out.res <- list(fun.obj=obj, model.obj =opt, results = res.sd,
                    convergence = list( optim.check=ifelse(opt$convergence >0, 
                                                           "not converged", 
                                                           "converged"),
                                        pdHess.check = ifelse(!res.sd$pdHess, 
                                                              "not converged: Hessian of fixed effects was not positive definite", 
                                                              "converged")),
                    ci.obj = suppressWarnings( tmbprofileCI(obj, par=profile.par.name, 
                                                            use = method.ci, ystep=0.01, 
                                                            ytol=ytol) ) )
  }else{
    
    # run cases where the only unknown is beta/mu for the single parameter cases
    switch(model,
           "Bivar.Norm"={
             
             ### assign names to the parameters in the model
             model.param  <- c("lnsd1","lnsd2","df1","beta")
             if(is.null(profile.par.name)){
               profile.par.name <- model.param[which(model.param %nin% names(known.parm))]
             }
             true.par.vals <- known.parm
             
             ### Fit the TMB model
             LLvec     <- matrix(data = NA, ncol=3, nrow = length(profile.par))
             LLvec[,1] <- profile.par
             for(i in 1:length(profile.par)){
               
               ## assign the profile parameter the ith value
               known.parm[[profile.par.name]] <- profile.par[i]
               
               ## create the objects for TMB
               data <- list(y=sim.data$Y, w=sim.data$W)
               par  <- append( list(X=rep(1, dim(sim.data)[1])), known.parm)
               obj  <- MakeADFun(data, par, random = "X", 
                                 DLL = "EV_bivnorm", silent = T,
                                 # Parameter entries with NAs in the factor are fixed)
                                 map = list(beta=factor(NA),lnsd1=factor(NA), 
                                            lnsd2=factor(NA),df1=factor(NA))) 
               opt <- tryCatch(
                 optim(obj$par, obj$fn, obj$gr, "BFGS"),
                 error=function(e) e
               )
               if(!inherits(opt, "error")){
                 LLvec[i,2] <- -opt$value
               }else{break}
               
             }
             LLvec <- LLvec[which(!is.na(LLvec[,2])),]
             LLvec[,3] <- exp(LLvec[,2]-max(LLvec[,2]))
             out <- LLvec
             out <- as.data.frame(LLvec)
             colnames(out) <- c("param","loglik","rel.like")
             
             
             
             ### compute the ML estimate
             data.ML <- list(y=sim.data$Y, w=sim.data$W)
             known.parm.ML <- true.par.vals[names(true.par.vals)!=profile.par.name]
             
             fixed.par <- known.parm.ML
             for(param in names(fixed.par)){ fixed.par[[param]] <- factor(NA) }
             
             par.ML  <- append( list(X=rep(1, dim(sim.data)[1])), true.par.vals)
             par.ML[[profile.par.name]] <- single.par.init
             
             obj.ML  <- MakeADFun(data.ML, par.ML, random = "X", 
                                  DLL = "EV_bivnorm", silent = T,
                                  # Parameter entries with NAs in the factor are fixed)
                                  map = fixed.par)
           },
           
           "Pois.Norm"={
             
             ### assign names to the parameters in the model
             model.param  <- c("lnsd1","df1","beta")
             if(is.null(profile.par.name)){
               profile.par.name <- model.param[which(model.param %nin% names(known.parm))]
             }
             true.par.vals <- known.parm
             
             ### Fit the TMB model
             LLvec     <- matrix(data = NA, ncol=3, nrow = length(profile.par))
             LLvec[,1] <- profile.par
             for(i in 1:length(profile.par)){
               
               ## assign the profile parameter the ith value
               known.parm[[profile.par.name]] <- profile.par[i]
               
               ## create the objects for TMB
               data <- list(y=sim.data$Y, w=sim.data$W)
               par  <- append( list(X=rep(1, dim(sim.data)[1])), known.parm)
               
               obj       <- MakeADFun(data, par, random = "X", 
                                      DLL = "EV_normpois", silent = T,
                                      #Parameter entries with NAs in the factor are fixed)
                                      map = list(beta=factor(NA),
                                                 lnsd1=factor(NA),
                                                 df1=factor(NA))) 
               opt <- tryCatch(
                 optim(obj$par, obj$fn, obj$gr, "BFGS"),
                 error=function(e) e
               )
               if(!inherits(opt, "error")){
                 LLvec[i,2] <- -opt$value
               }else{break}
               
             }
             LLvec <- LLvec[which(!is.na(LLvec[,2])),]
             LLvec[,3] <- exp(LLvec[,2]-max(LLvec[,2]))
             out <- LLvec
             out <- as.data.frame(LLvec)
             colnames(out) <- c("param","loglik","rel.like")
             
             ### now compute the ML estimate
             data.ML <- list(y=sim.data$Y, w=sim.data$W)
             known.parm.ML <- true.par.vals[names(true.par.vals)!=profile.par.name]
             
             fixed.par <- known.parm.ML
             for(param in names(fixed.par)){ fixed.par[[param]] <- factor(NA) }
             
             par.ML  <- append( list(X=rep(1, dim(sim.data)[1])), true.par.vals)
             par.ML[[profile.par.name]] <- single.par.init
             
             obj.ML  <- MakeADFun(data.ML, par.ML, random = "X", 
                                  DLL = "EV_normpois", silent = T,
                                  # Parameter entries with NAs in the factor are fixed)
                                  map = fixed.par)
           },
           
           "Logit.Norm"={
             
             ### assign names to the parameters in the model
             model.param  <- c("lnsd1","df1","beta")
             if(is.null(profile.par.name)){
               profile.par.name <- model.param[which(model.param %nin% names(known.parm))]
             }
             true.par.vals <- known.parm
             
             ### Fit the TMB model
             LLvec     <- matrix(data = NA, ncol=3, nrow = length(profile.par))
             LLvec[,1] <- profile.par
             for(i in 1:length(profile.par)){
               
               ## assign the profile parameter the ith value
               known.parm[[profile.par.name]] <- profile.par[i]
               
               ## create the objects for TMB
               data <- list(y=sim.data$Y, w=sim.data$W)
               par  <- append( list(X=rep(1, dim(sim.data)[1])), known.parm)
               
               obj       <- MakeADFun(data, par, random = "X", 
                                      DLL = "EV_normlogit", silent = T,
                                      #Parameter entries with NAs in the factor are fixed)
                                      map = list(beta=factor(NA),
                                                 lnsd1=factor(NA),
                                                 df1=factor(NA))) 
               
               opt <- tryCatch(
                 optim(obj$par, obj$fn, obj$gr, "BFGS"),
                 error=function(e) e
               )
               if(!inherits(opt, "error")){
                 LLvec[i,2] <- -opt$value
               }else{break}
               
             }
             LLvec <- LLvec[which(!is.na(LLvec[,2])),]
             LLvec[,3] <- exp(LLvec[,2]-max(LLvec[,2]))
             out <- LLvec
             out <- as.data.frame(LLvec)
             colnames(out) <- c("param","loglik","rel.like")
             
             ### now compute the ML estimate
             data.ML <- list(y=sim.data$Y, w=sim.data$W)
             known.parm.ML <- true.par.vals[names(true.par.vals)!=profile.par.name]
             
             fixed.par <- known.parm.ML
             for(param in names(fixed.par)){ fixed.par[[param]] <- factor(NA) }
             
             par.ML  <- append( list(X=rep(1, dim(sim.data)[1])), true.par.vals)
             par.ML[[profile.par.name]] <- single.par.init
             
             obj.ML  <- MakeADFun(data.ML, par.ML, random = "X", 
                                  DLL = "EV_normlogit", silent = T,
                                  # Parameter entries with NAs in the factor are fixed)
                                  map = fixed.par)
           },
           
           "Probit.Norm"={
             
             ### assign names to the parameters in the model
             model.param  <- c("lnsd1","df1","beta")
             if(is.null(profile.par.name)){
               profile.par.name <- model.param[which(model.param %nin% names(known.parm))]
             }
             true.par.vals <- known.parm
             
             ### Fit the TMB model
             LLvec     <- matrix(data = NA, ncol=3, nrow = length(profile.par))
             LLvec[,1] <- profile.par
             for(i in 1:length(profile.par)){
               
               ## assign the profile parameter the ith value
               known.parm[[profile.par.name]] <- profile.par[i]
               
               ## create the objects for TMB
               data <- list(y=sim.data$Y, w=sim.data$W)
               par  <- append( list(X=rep(1, dim(sim.data)[1])), known.parm)
               
               obj       <- MakeADFun(data, par, random = "X", 
                                      DLL = "EV_normprobit", silent = T,
                                      #Parameter entries with NAs in the factor are fixed)
                                      map = list(beta=factor(NA),
                                                 lnsd1=factor(NA),
                                                 df1=factor(NA))) 
               opt <- tryCatch(
                 optim(obj$par, obj$fn, obj$gr, "BFGS"),
                 error=function(e) e
               )
               if(!inherits(opt, "error")){
                 LLvec[i,2] <- -opt$value
               }else{break}
               
             }
             LLvec <- LLvec[which(!is.na(LLvec[,2])),]
             
             LLvec[,3] <- exp(LLvec[,2]-max(LLvec[,2]))
             out <- LLvec
             out <- as.data.frame(LLvec)
             colnames(out) <- c("param","loglik","rel.like")
             
             ### now compute the ML estimate
             data.ML <- list(y=sim.data$Y, w=sim.data$W)
             known.parm.ML <- true.par.vals[names(true.par.vals)!=profile.par.name]
             
             fixed.par <- known.parm.ML
             for(param in names(fixed.par)){ fixed.par[[param]] <- factor(NA) }
             
             par.ML  <- append( list(X=rep(1, dim(sim.data)[1])), true.par.vals)
             par.ML[[profile.par.name]] <- single.par.init
             
             obj.ML  <- MakeADFun(data.ML, par.ML, random = "X", 
                                  DLL = "EV_normprobit", silent = T,
                                  # Parameter entries with NAs in the factor are fixed)
                                  map = fixed.par)
           },
           
           "Normsq.Pois"={
             
             
             ### assign names to the parameters in the model
             model.param  <- c("lnsd1","mu")
             if(is.null(profile.par.name)){
               profile.par.name <- model.param[which(model.param %nin% names(known.parm))]
             }
             true.par.vals <- known.parm
             
             ### Fit the TMB model
             LLvec     <- matrix(data = NA, ncol=3, nrow = length(profile.par))
             LLvec[,1] <- profile.par
             for(i in 1:length(profile.par)){
               
               ## assign the profile parameter the ith value
               known.parm[[profile.par.name]] <- profile.par[i]
               
               ## create the objects for TMB
               data <- list(y=sim.data$Y)
               par  <- append( list(X=rep(1, dim(sim.data)[1])), known.parm)
               
               obj       <- MakeADFun(data, par, random = "X", 
                                      DLL = "normsq_pois", silent = T,
                                      #Parameter entries with NAs in the factor are fixed)
                                      map = list(mu=factor(NA),
                                                 lnsd1=factor(NA))) 
               opt <- tryCatch(
                 optim(obj$par, obj$fn, obj$gr, "BFGS"),
                 error=function(e) e
               )
               if(!inherits(opt, "error")){
                 LLvec[i,2] <- -opt$value
               }else{break}
               
             }
             LLvec <- LLvec[which(!is.na(LLvec[,2])),]
             
             LLvec[,3] <- exp(LLvec[,2]-max(LLvec[,2]))
             out <- LLvec
             out <- as.data.frame(LLvec)
             colnames(out) <- c("param","loglik","rel.like")
             
             ### now compute the ML estimate
             data.ML <- list(y=sim.data$Y)
             known.parm.ML <- true.par.vals[names(true.par.vals)!=profile.par.name]
             
             fixed.par <- known.parm.ML
             for(param in names(fixed.par)){ fixed.par[[param]] <- factor(NA) }
             
             par.ML  <- append( list(X=rep(1, dim(sim.data)[1])), true.par.vals)
             par.ML[[profile.par.name]] <- single.par.init
             
             obj.ML  <- MakeADFun(data.ML, par.ML, random = "X", 
                                  DLL = "normsq_pois", silent = T,
                                  # Parameter entries with NAs in the factor are fixed)
                                  map = fixed.par)
           },
    )
    # use optim or nlminb
    if(use.optim==TRUE){
      opt.ML  <- suppressWarnings( optim(par=obj.ML$par, fn=obj.ML$fn, gr=obj.ML$gr, hessian = TRUE, method = "BFGS") )
    }else{
      opt.ML  <- suppressWarnings( nlminb(start=obj.ML$par, objective=obj.ML$fn, gradient=obj.ML$gr) )
    }
    res.sd <- sdreport(obj.ML, getJointPrecision = T, bias.correct = T)
    internal.checks( obj.ML )
    MLres  <- list( res= res.sd, obj.ML=obj.ML, opt.ML=opt.ML,
                    converged = ifelse(!res.sd$pdHess, 
                                       "not converged: Hessian of fixed effects was not positive definite", 
                                       "converged"))
    out.res <- list( prof=out, MLres=MLres )
    
  }
  ### return the ci
  return(out.res)
}

# from https://github.com/glmmTMB/glmmTMB/blob/7a20635602a4aeeea840b8fffdffe2413254a471/glmmTMB/R/glmmTMB.R
internal.checks <- function(obj, h = NULL, data.tmb.old = NULL){
  
  sdr <- sdreport(obj)
  
  ## generic complex eigenvalue checker
  e_complex_check <- function(ev, tol=sqrt(.Machine$double.eps)) {
    if (is.complex(ev)) {
      if ((maxim <- max(abs(Im(ev)))) > tol) {
        ev <- Re(ev)
      } else {
        stop(sprintf("detected complex eigenvalues of covariance matrix (max(abs(Im))=%g: try se=FALSE?",
                     maxim))
      }
    }
    return(ev)
  }
  if(!is.null(sdr$pdHess)) {
    if(!sdr$pdHess) {
      ## double-check (slower, more accurate hessian)
      env <- environment(obj$fn)
      par <- env$last.par.best
      if (!is.null(rr <- env$random)) {
        par <- par[-rr]
      }
      h <- numDeriv::jacobian(obj$gr, par)
      ## fall back to solve(optimHess(par, obj$fn, obj$gr)) ? 
      h <- .5 * (h + t(h))  ## symmetrize
      if (!any(is.na(h))) {
        ev <- try(
          e_complex_check(eigen(h)$values),
          silent = TRUE)
        if (!inherits(ev, "try-error") && min(ev)>.Machine$double.eps) {
          ## apparently fit is OK after all ...
          sdr$pdHess <- TRUE
          Vtheta <- try(solve(h), silent=TRUE)
          if (!inherits(Vtheta,"try-error")) {
            sdr$cov.fixed[] <- Vtheta
          } else {
            warning("failed to invert Hessian from numDeriv::jacobian(), falling back to internal vcov estimate")
          }
        } ## eig check OK
      } ## !any(is.na(h))
    } ## !sdr$pdHess          
    if (!sdr$pdHess) { ## still bad after trying numDeriv ...
      warning(paste0("Model convergence problem; ",
                     "non-positive-definite Hessian matrix. "))
    }
  } else if(length(sdr$cov.fixed)>0) {
    eigval <- try(1/eigen(sdr$cov.fixed)$values, silent=TRUE)
    if( is(eigval, "try-error") || ( min(e_complex_check(eigval)) < .Machine$double.eps*10 ) ) {
      warning(paste0("Model convergence problem; ",
                     "extreme or very small eigenvalues detected. "))
    } ## extreme eigval
  }  ## do eigval check
  
}



# modified from https://github.com/glmmTMB/glmmTMB/blob/master/glmmTMB/R/diagnose.R
diagnose <- function(fit,
                     eval_eps=1e-5, evec_eps=1e-2,
                     check_hessian=TRUE, explain = TRUE) {
  prt_explain <- function(...) {
    if (explain) {
      s <- do.call(paste, c(list(...), list(collapse = "")))
      cat(strwrap(s), "\n", sep = "\n")
    }
    return(invisible(NULL))
  }
  model_OK <- TRUE
  ## pull out the TMB object from the fit
  obj <- fit$fun.obj
  ee  <- obj$env
  
  ## extract parameters
  pp <- ee$last.par.best
  if (!is.null(r <- ee$random)) { pp <- pp[-r] }
  ss <- suppressWarnings(summary( sdreport(obj) ))
  ss <- ss[grepl("^(beta|ln|mu)", rownames(ss)), ]
  
  fit$results$par.fixed
  
  ## easiest way to get names corresponding to all of the parameters
  vv <- solve(fit$model.obj$hessian)
  nn <- make.unique(names(pp))
  rownames(ss) <- names(pp) <- nn
  
  
  if (check_hessian) {
    if (!"results" %in% names(fit)) {
      warning("sdreport was not computed, skipping Hessian check")
    }
    if (!fit$results$pdHess) {
      model_OK <- FALSE
      cat("Non-positive definite (NPD) Hessian\n\n")
      prt_explain("The Hessian matrix represents the curvature of the",
                  "log-likelihood surface at the maximum likelihood estimate (MLE) of the parameters",
                  "(its inverse is the estimate of the parameter covariance matrix). ",
                  "A non-positive-definite Hessian means that the likelihood surface is approximately flat ",
                  "(or upward-curving) at the MLE, which means the model is overfitted or poorly posed ",
                  "in some way.",
                  "NPD Hessians are often associated with extreme parameter estimates."
                  
      )
      nonfinite_sd <- !is.finite(suppressWarnings(sqrt(diag(fit$results$cov.fixed))))
      if (any(nonfinite_sd)) {
        cat("parameters with non-finite standard deviations:\n")
        if (all(nonfinite_sd)) {
          cat("(all of them!)\n\n")
        } else {
          cat(strwrap(paste(nn[nonfinite_sd],
                            collapse=", ")),"\n\n",sep="\n")
        }
      }
      
      ## fit hessian with Richardson extrapolation (more accurate/slower than built-in optimHess)
      cat("recomputing Hessian via Richardson extrapolation.",
          "If this is too slow, consider setting check_hessian = FALSE",
          "\n\n")
      h <- numDeriv::jacobian(fit$fun.obj$gr, pp)
      if (any(is.na(h))) {
        bad <- NA
        eigs <- eigen(fit$results$cov.fixed)
        complex_eigs <- is.complex(eigs)
      } else {
        eigs <- eigen(h)
        ## non-positive definite means some of the eigenvectors are <= 0
        complex_eigs <- is.complex(eigs$values)
        if (!complex_eigs) {
          bad <- which(eigs$values/max(eigs$values) <= eval_eps)
        }
      }
      if (length(bad) == 0 && !complex_eigs) {
        cat("Hessian seems OK\n")
        prt_explain("Internal calculations suggested that the Hessian was bad/non-positive definite;",
                    "however, a slower and more precise calculation suggests that it's actually OK. Your model",
                    "may be somewhat numerically unstable.")
        return(invisible(h)) ## bail out here
      }
      if (complex_eigs) {
        cat("Hessian has complex eigenvalues\n\n")
        prt_explain("We would have used the smallest eigenvalues of the ",
                    "Hessian to determine which components were bad",
                    "but instead we got complex eigenvalues.",
                    "(Not really sure what to do with this ...)")
      } else {
        prt_explain("The next set of diagnostics attempts to determine which elements of the Hessian",
                    "are causing the non-positive-definiteness. ",
                    "Components with very small eigenvalues represent 'flat' directions, ",
                    "i.e., combinations of parameters for which the data may contain very little information. ",
                    sprintf("So-called 'bad elements' represent the dominant components (absolute values >%1.3g)", evec_eps),
                    "of the eigenvectors corresponding to the 'flat' directions")
        cat(sprintf("maximum Hessian eigenvalue = %1.3g",eigs$values[1]),"\n")
        ## there could be more than one 'bad' direction/eigenvector ..
        for (b in bad) {
          cat(sprintf("Hessian eigenvalue %d = %1.3g (relative val = %1.3g)",
                      b, eigs$values[b],eigs$values[b]/eigs$values[1]),"\n")
          bad_vec <- eigs$vectors[,b]
          bad_elements <- which(abs(bad_vec)>evec_eps)
          cat("   bad elements:", nn[bad_elements], "\n")
        }
      } ## hessian has real eigs
    } ## bad hessian
  } ## check hessian
  if (model_OK) cat("model looks OK!\n")
  return(invisible(model_OK))
}
