.inspect.engine <- function(tpts, log_shift, concentrations, rates
                            , nInit=10, nIter=300, na.rm=TRUE, BPPARAM=bpparam() #nCores=2L
                            , verbose=TRUE, estimateRatesWith=c('der', 'int'), nAttempts=1
                            , sigmoidDegradation=FALSE, sigmoidSynthesis=FALSE, sigmoidTotal=FALSE
                            , sigmoidProcessing=FALSE, sigmoidPre=FALSE
                            , testOnSmooth=TRUE, seed=NULL)
{
  
  chisq.test.inspect <- function(D, df)
    pchisq(D, df, lower.tail=TRUE)
  
  fisher.test.inspect <- function(p1, p2)
    pchisq(-2*(log(p1) + log(p2)), 4, lower.tail=FALSE)
  
  optimParamsSimple <- function(interpRates, tpts_exp, alpha_exp, alpha_var
                                , total_exp, total_var, maxit=500
                                , log_shift, .time_transf, .rxnrateSimple, ode, .makeSimpleModel, .logLikelihood
                                , .emptyGene)
  {
    
    if( is.null(interpRates) ) return(.emptyGene())
    
    simpleModelChisq <- function(par, tpts, fun, df, alpha_exp, alpha_var
                                 , total_exp, total_var)
    {
      splitpar <- split(par 
                        , c(rep('alpha',df[1]), rep('beta',df[2]))
      )
      #
      params <- list()
      params$alpha <- function(x) 
        fun$alpha$value(.time_transf(x, log_shift), splitpar$alpha)
      params$beta <- function(x)
        fun$beta$value(.time_transf(x, log_shift), splitpar$beta)
      #
      cinit <- c(params$alpha(tpts[1]) / params$beta(tpts[1]))
      names(cinit) <- 't'
      model <- ode(y=cinit, times=tpts, func=.rxnrateSimple, parms=params)
      #
      alpha_model <- params$alpha(tpts)
      total_model <- model[,'t']
      #
      D <- chisq(alpha_exp, alpha_model, alpha_var) +
        chisq(total_exp, total_model, total_var)
      # df <- length(alpha_exp) + length(total_exp) - sum(df[1:2])
      # return(log(chisq.test.inspect(D, df)))
      return(D)
    }
    
    optOut <- tryCatch(
      optim(
        par=c(interpRates$alpha$par, interpRates$beta$par)
        , fn=simpleModelChisq, tpts=tpts_exp 
        , fun=list(alpha=interpRates$alpha$fun
                   , beta=interpRates$beta$fun) 
        , df=c(interpRates$alpha$df, interpRates$beta$df)
        , alpha_exp=alpha_exp, total_exp=total_exp
        , alpha_var=alpha_var, total_var=total_var
        , control=list(maxit=maxit) #, trace=1)
        , method='Nelder-Mead'
      )
      , error=function(e)
        optim(
          par=c(interpRates$alpha$par, interpRates$beta$par)
          , fn=simpleModelChisq, tpts=tpts_exp 
          , fun=list(alpha=interpRates$alpha$fun
                     , beta=interpRates$beta$fun) 
          , df=c(interpRates$alpha$df, interpRates$beta$df)
          , alpha_exp=alpha_exp, total_exp=total_exp
          , alpha_var=alpha_var, total_var=total_var
          , control=list(maxit=maxit) #, trace=1)
          , method='BFGS'
        )
    )
    #
    splitpar <- split(optOut$par, c(rep('alpha',interpRates$alpha$df)
                                    , rep('beta',interpRates$beta$df)))
    interpRates$alpha$params <- splitpar$alpha
    interpRates$beta$params  <- splitpar$beta
    #
    model <- .makeSimpleModel(tpts=tpts_exp
                              , hyp=list(alpha=interpRates$alpha, beta=interpRates$beta)
                              , log_shift, .time_transf, ode, .rxnrateSimple)
    logLik <- .logLikelihood(alpha_exp, model$alpha, alpha_var) + 
      .logLikelihood(total_exp, model$total, total_var)
    k <- interpRates$alpha$df + interpRates$beta$df
    n <- length(alpha_exp) + length(total_exp)
    chisqTest <- log(pchisq(optOut$value, n-k, lower.tail=TRUE))
    AIC <- 2*k - 2*logLik
    AICc <- AIC + 2*k*(k+1)/(n-k-1)
    return(list(
      alpha=interpRates$alpha
      , beta=interpRates$beta
      , gamma=.emptyGene()$gamma
      , test=chisqTest
      , logLik=logLik
      , AIC=AIC
      , AICc=AICc
      , counts=optOut$counts
      , convergence=optOut$convergence
      , message=optOut$message
    ))
    
  }
  
  optimParams <- function(interpRates, tpts_exp, alpha_exp, alpha_var, total_exp
                          , total_var, preMRNA_exp, preMRNA_var
                          # , test=c('merged', 'combined'), pval=c('lin', 'log')
                          , maxit=500, log_shift, .time_transf, .rxnrate, ode, .makeModel, .logLikelihood)
  {
    modelChisq <- function(par, tpts, fun, df, alpha_exp, alpha_var #, pval
                           , total_exp, total_var, preMRNA_exp, preMRNA_var)
    {
      splitpar <- split(par 
                        , c(rep('alpha',df[1]), rep('beta',df[2]), rep('gamma',df[3])) 
      )
      #
      params <- list()
      params$alpha <- function(x) 
        fun$alpha$value(.time_transf(x, log_shift), splitpar$alpha)
      params$beta  <- function(x)
        fun$beta$value(.time_transf(x, log_shift), splitpar$beta)
      params$gamma <- function(x)
        fun$gamma$value(.time_transf(x, log_shift), splitpar$gamma)
      #
      cinit <- c(params$alpha(tpts[1]) / params$gamma(tpts[1])
                 , params$alpha(tpts[1]) / params$beta(tpts[1]) + 
                   params$alpha(tpts[1]) / params$gamma(tpts[1]))
      names(cinit) <- c('p', 't')
      model <- ode(y=cinit, times=tpts, func=.rxnrate, parms=params)
      #
      alpha_model <- params$alpha(tpts)
      total_model <- model[,'t']
      preMRNA_model <- model[,'p']
      #
      D <- chisq(alpha_exp, alpha_model, alpha_var) +
        chisq(total_exp, total_model, total_var) +
        chisq(preMRNA_exp, preMRNA_model, preMRNA_var)
      df <- length(alpha_exp) + length(total_exp) + length(preMRNA_exp) - sum(df)
      testValue <- chisq.test.inspect(D, df)
      return(testValue)
    }
    optOut <- tryCatch(
      #
      # optimize the model to the experiment using 
      # Nelder-Mead method
      #
      optim(
        par=unlist(sapply(interpRates, '[[', 'params'))
        , fn=modelChisq , tpts=tpts_exp #, pval=pval 
        , fun=sapply(interpRates, '[[', 'fun')
        , df=unlist(sapply(interpRates, '[[', 'df'))
        , alpha_exp=alpha_exp, total_exp=total_exp
        , preMRNA_exp=preMRNA_exp
        , alpha_var=alpha_var, total_var=total_var
        , preMRNA_var=preMRNA_var
        , control = list(maxit = maxit)
        , method = 'Nelder-Mead'
      )
      #
      # in case Nelder-Mead method fails, try with BFGS
      #
      , error=function(e) {
        optim(
          par=unlist(sapply(interpRates, '[[', 'params'))
          , fn=modelChisq, tpts=tpts_exp #, pval=pval
          , fun=sapply(interpRates, '[[', 'fun')
          , df=unlist(sapply(interpRates, '[[', 'df'))
          , alpha_exp=alpha_exp, total_exp=total_exp
          , preMRNA_exp=preMRNA_exp
          , alpha_var=alpha_var, total_var=total_var
          , preMRNA_var=preMRNA_var
          , control = list(maxit = maxit)
          , method = 'BFGS'
        )
      })
    if( !is.null(optOut$error) ) return(optOut$error)
    #
    # assign the parameters belonging to the optimized model
    # to the parameter functions
    #
    splitpar <- split(optOut$par
                      , c(rep('alpha',interpRates$alpha$df)
                          , rep('beta',interpRates$beta$df)
                          , rep('gamma',interpRates$gamma$df))
    )
    interpRates$alpha$params <- splitpar$alpha
    interpRates$beta$params  <- splitpar$beta
    interpRates$gamma$params <- splitpar$gamma
    # even if the minimization used the linear pvalue
    # give back the log one
    # if( pval=='lin' ) 
    optOut$value <- log(optOut$value)
    #
    # return parameter functions and other output
    # from the optimization procedure
    #
    model <- .makeModel(tpts=tpts_exp
                        , hyp=list(alpha=interpRates$alpha, beta=interpRates$beta, 
                                   gamma=interpRates$gamma)
                        , log_shift, .time_transf, ode, .rxnrate)
    logLik <- .logLikelihood(alpha_exp, model$alpha, alpha_var) + 
      .logLikelihood(total_exp, model$total, total_var) +
      .logLikelihood(preMRNA_exp, model$preMRNA, preMRNA_var)
    
    k <- interpRates$alpha$df + interpRates$beta$df + interpRates$gamma$df
    n <- length(alpha_exp) + length(total_exp) + length(preMRNA_exp)
    chisqTest <- optOut$value
    AIC <- 2*k - 2*logLik
    AICc <- AIC + 2*k*(k+1)/(n-k-1)
    return(list(
      alpha=interpRates$alpha
      , beta=interpRates$beta
      , gamma=interpRates$gamma
      , test=chisqTest
      , logLik=logLik
      , AIC=AIC
      , AICc=AICc
      , counts=optOut$counts
      , convergence=optOut$convergence
      , message=optOut$message
    ))
  }
  
  modelOneGene <- function(i, seed=NULL,
                           .chooseModel, .time_transf, .DimpulseModel, .DsigmoidModel, .constantModelP,
                           .emptyGene, .sigmoidModel, .impulseModel, .sigmoidModelP, .impulseModelP,
                           .polynomialModelP, .makeModel, .makeSimpleModel, .logLikelihood, .rxnrate,
                           .rxnrateSimple, optimParams, optimParamsSimple, verbose, nAttempts,
                           concentrations, rates, tpts, log_shift, na.rm, sigmoidTotal,
                           sigmoidSynthesis, nInit, nIter, testOnSmooth, estimateRatesWith) 
  {
    

    
    

    ## set the mode of the gene, "only exons gene" or 
    ## "introns exons gene"
    if( !is.null(seed) ) set.seed(seed)
    if( all(is.na(concentrations$preMRNA[i,])) |  
        all(is.na(rates$gamma[i,])) ) 
    {
      intExMode <- FALSE
    } else {
      intExMode <- TRUE
    }
    
    
    ## start the analysis
    paramAttempts <- sapply(1:nAttempts, function(k)
      # tryCatch(
    {
      modelTotalRNAfun <- tryCatch(.chooseModel(tpts=tpts
                                                , log_shift=log_shift
                                                , experiment=concentrations$total[i,]
                                                , variance=concentrations$total_var[i]
                                                , na.rm=na.rm, sigmoid=sigmoidTotal
                                                , impulse=TRUE, polynomial=FALSE
                                                , nInit=nInit, nIter=nIter
                                                , .time_transf=.time_transf
                                                , .sigmoidModel=.sigmoidModel
                                                , .impulseModel=.impulseModel
                                                , .sigmoidModelP=.sigmoidModelP
                                                , .impulseModelP=.impulseModelP
                                                , .polynomialModelP=.polynomialModelP
      ), error=function(e) return(.emptyGene(e)))
     
      
       modelSynthesisRatefun <- tryCatch(.chooseModel(tpts=tpts
                                                     , log_shift=log_shift
                                                     , experiment=rates$alpha[i,] 
                                                     , variance=rates$alpha_var[i]
                                                     , na.rm=na.rm, sigmoid=sigmoidSynthesis
                                                     , impulse=TRUE, polynomial=FALSE
                                                     , nInit=nInit, nIter=nIter
                                                     , .time_transf=.time_transf
                                                     , .sigmoidModel=.sigmoidModel
                                                     , .impulseModel=.impulseModel
                                                     , .sigmoidModelP=.sigmoidModelP
                                                     , .impulseModelP=.impulseModelP
                                                     , .polynomialModelP=.polynomialModelP
      ), error=function(e) return(.emptyGene(e)))
      
       
       modelTotalRNA <- 
        if( testOnSmooth ) {
          modelTotalRNAfun$fun$value(
            .time_transf(tpts, log_shift)
            , modelTotalRNAfun$params)
        } else { concentrations$total[i,] }
      modelSynthesisRate <- 
        if( testOnSmooth ) {
          modelSynthesisRatefun$fun$value(
            .time_transf(tpts, log_shift)
            , modelSynthesisRatefun$params)	
        } else { rates$alpha[i,] }
      if( intExMode ) {
        modelPreMRNAfun <- tryCatch(.chooseModel(tpts=tpts
                                                 , log_shift=log_shift
                                                 , experiment=concentrations$preMRNA[i,]
                                                 , variance=concentrations$preMRNA_var[i]
                                                 , na.rm=na.rm, sigmoid=sigmoidPre
                                                 , impulse=TRUE, polynomial=FALSE
                                                 , nInit=nInit, nIter=nIter
                                                 , .time_transf=.time_transf
                                                 , .sigmoidModel=.sigmoidModel
                                                 , .impulseModel=.impulseModel
                                                 , .sigmoidModelP=.sigmoidModelP
                                                 , .impulseModelP=.impulseModelP
                                                 , .polynomialModelP=.polynomialModelP
        ), error=function(e) return(.emptyGene(e)))
        modelPreMRNA <- 
          if( testOnSmooth ) 
            modelPreMRNAfun$fun$value(
              .time_transf(tpts, log_shift)
              , modelPreMRNAfun$params) 
        else concentrations$preMRNA[i,]
      }
      #
      if( estimateRatesWith == 'der' ) {
        if( testOnSmooth ) {
          # total RNA derivative
          if( modelTotalRNAfun$type == 'impulse' )
            modelTotalRNAderivative <- .DimpulseModel(
              .time_transf(tpts, log_shift)
              , modelTotalRNAfun$params)
          if( modelTotalRNAfun$type == 'sigmoid' )
            modelTotalRNAderivative <- .DsigmoidModel(
              .time_transf(tpts, log_shift)
              , modelTotalRNAfun$params)
          if( modelTotalRNAfun$type == 'constant' )
            modelTotalRNAderivative <- rep(0, length(tpts))
          if( intExMode ) {
            # pre mRNA derivative
            if( modelPreMRNAfun$type == 'impulse' )
              modelPreMRNAderivative <- .DimpulseModel(
                .time_transf(tpts, log_shift)
                , modelPreMRNAfun$params)
            if( modelPreMRNAfun$type == 'sigmoid' )
              modelPreMRNAderivative <- .DsigmoidModel(
                .time_transf(tpts, log_shift)
                , modelPreMRNAfun$params)
            if( modelPreMRNAfun$type == 'constant' )
              modelPreMRNAderivative <- rep(0, length(tpts))
          }
        } else {
          # total RNA derivative
          spfun <- splinefun(tpts, concentrations$total[i,] )
          modelTotalRNAderivative <- spfun(tpts, deriv=1)
          if( intExMode ) {
            # pre mRNA derivative
            spfun <- splinefun(tpts, concentrations$preMRNA[i,] )
            modelPreMRNAderivative <- spfun(tpts, deriv=1)
          }
        }
        # degradation rate
        modelTotalRNAderivative[1] <- 0
        firstGuessDegrRate <- (modelSynthesisRate - 
                                 modelTotalRNAderivative ) / (modelTotalRNA - modelPreMRNA )
        if( intExMode ) {
          # processing rate
          modelPreMRNAderivative[1] <- 0
          firstGuessProcessingRate <- (modelSynthesisRate - 
                                         modelPreMRNAderivative ) / modelPreMRNA					
        }
      } else {
        # degradation rate
        firstGuessDegrRate <- rates$beta[i,]
        if( intExMode ) {
          # processing rate
          firstGuessProcessingRate <- rates$gamma[i,]
        }
      }
      idx <- firstGuessDegrRate<0 | 
        !is.finite(firstGuessDegrRate)
      firstGuessDegrRate[idx] <- NA
      if( intExMode ) {
        idx <- firstGuessProcessingRate<0 | 
          !is.finite(firstGuessProcessingRate)
        firstGuessProcessingRate[idx] <- NA
      }
      #
      constantRates <- list(
        alpha=list(fun=.constantModelP
                   , type='constant', df=1
                   , params=mean(modelSynthesisRate, na.rm=TRUE))
        , beta=list(fun=.constantModelP
                    , type='constant', df=1
                    , params=mean(firstGuessDegrRate, na.rm=TRUE))
        , gamma=if( intExMode ) {
          list(fun=.constantModelP
               , type='constant', df=1
               , params=mean(firstGuessProcessingRate, na.rm=TRUE))
        } else { NULL }
      )
      varyingRates <- list(
        alpha=modelSynthesisRatefun
        , beta=tryCatch(
          .chooseModel(tpts=tpts
                       , log_shift=log_shift
                       , experiment=firstGuessDegrRate
                       , variance=1
                       , na.rm=na.rm, sigmoid=sigmoidDegradation
                       , impulse=TRUE, polynomial=FALSE
                       , nInit=nInit, nIter=nIter
                       , .time_transf=.time_transf
                       , .sigmoidModel=.sigmoidModel
                       , .impulseModel=.impulseModel
                       , .sigmoidModelP=.sigmoidModelP
                       , .impulseModelP=.impulseModelP
                       , .polynomialModelP=.polynomialModelP
          )
          , error=function(e) return(.emptyGene(e)$beta))
        , gamma=if( intExMode ) { tryCatch(
          .chooseModel(tpts=tpts
                       , log_shift=log_shift
                       , experiment=firstGuessProcessingRate
                       , variance=1
                       , na.rm=na.rm, sigmoid=sigmoidProcessing
                       , impulse=TRUE, polynomial=FALSE
                       , nInit=nInit, nIter=nIter
                       , .time_transf=.time_transf
                       , .sigmoidModel=.sigmoidModel
                       , .impulseModel=.impulseModel
                       , .sigmoidModelP=.sigmoidModelP
                       , .impulseModelP=.impulseModelP
                       , .polynomialModelP=.polynomialModelP
          ), error=function(e) return(.emptyGene(e)$gamma))
        } else { NULL }
      )
      ratesToTest <- list('0'=constantRates
                          , a=list(alpha=varyingRates$alpha, beta=constantRates$beta
                                   , gamma=constantRates$gamma)
                          , b=list(alpha=constantRates$alpha, beta=varyingRates$beta
                                   , gamma=constantRates$gamma)
                          , c=if( intExMode ) 
                            list(alpha=constantRates$alpha, beta=constantRates$beta
                                 , gamma=varyingRates$gamma) else NULL
                          , ab=list(alpha=varyingRates$alpha, beta=varyingRates$beta
                                    , gamma=constantRates$gamma)
                          , bc=if( intExMode ) 
                            list(alpha=constantRates$alpha, beta=varyingRates$beta
                                 , gamma=varyingRates$gamma) else NULL
                          , ac=if( intExMode ) 
                            list(alpha=varyingRates$alpha, beta=constantRates$beta
                                 , gamma=varyingRates$gamma) else NULL
                          , abc=if( intExMode ) varyingRates else NULL
      )
      if( intExMode ) {
        results <- lapply(ratesToTest, function(interpRates) {
          tryCatch(
            optimParams(interpRates
                        , tpts_exp=tpts
                        , alpha_exp=modelSynthesisRate
                        , alpha_var=rates$alpha_var[i]
                        , total_exp=modelTotalRNA
                        , total_var=concentrations$total_var[i]
                        , preMRNA_exp=modelPreMRNA
                        , preMRNA_var=concentrations$preMRNA_var[i]
                        , maxit=nIter
                        , log_shift=log_shift
                        , ode=deSolve::ode
                        , .time_transf=.time_transf
                        , .rxnrate=.rxnrate
                        , .makeModel=.makeModel
                        , .logLikelihood=.logLikelihood
            )
            , error=function(e) .emptyGene(e)
          )
        })
      } else {
        results <- lapply(ratesToTest, function(interpRates) {
          tryCatch(
            optimParamsSimple(interpRates
                              , tpts_exp=tpts
                              , alpha_exp=modelSynthesisRate
                              , alpha_var=rates$alpha_var[i]
                              , total_exp=modelTotalRNA
                              , total_var=concentrations$total_var[i]
                              , maxit=nIter
                              , log_shift=log_shift
                              , ode=deSolve::ode
                              , .time_transf=.time_transf
                              , .rxnrateSimple=.rxnrateSimple
                              , .makeSimpleModel=.makeSimpleModel
                              , .logLikelihood=.logLikelihood
                              , .emptyGene=.emptyGene
            )
            , error=function(e) .emptyGene(e)
          )
        })
      }
      ## sometimes tryCatch is not able to return
      ## the empty gene when an error occurr
      check <- sapply(results, length) == 10
      if( !all(check) ) {
        for( i in which(!check) ) {
          results[[i]] <- .emptyGene('Unknown error')
        }
      }
      return(results)
    } ) 
    
    
    
    
    
    
    
    
    
    
    
    
    
    if( verbose ) {
      if( is.null(rownames(concentrations$total)[i]) )
        message('Gene "no_name" completed.' )
      else message(paste('Gene "', 
                         rownames(concentrations$total)[i],'" completed.', sep=''))
    }
    ## choose the best model for each test, out of the many attempts
    chisqPvals <- sapply(1:nAttempts, function(i) 
      sapply(paramAttempts[,i], '[[', 'test'))
    
    ix <- apply(chisqPvals, 1, function(x) {
      x <- na.omit(x)
      if(length(x)>0) return(which.min(x)) else return(1)
    })
    ## in case the model is simple there is no minumum in all tests
    ## involving 'c', therefore force to assign to them the first attempt
    if( !intExMode ) ix[grep('c', names(ix))] <- 1
    ## extract the selected model, re-assign names and return
    selectedParams <- paramAttempts[cbind(1:length(ix),unlist(ix))]
    names(selectedParams) <- rownames(paramAttempts)
    # paramSpecs<-selectedParams
     return(selectedParams)
  }
  
  ######################
  ## MAIN FUNCTION ###
  ##################
  
  nGenes <- nrow(rates$alpha)
  tpts <- tpts
 
  paramSpecs <- lapply(1, modelOneGene, seed=seed, 
                         .chooseModel=.chooseModel,
                         .time_transf=.time_transf,
                         .DimpulseModel=.DimpulseModel,
                         .DsigmoidModel=.DsigmoidModel,
                         .constantModelP=.constantModelP,
                         .emptyGene=.emptyGene,
                         optimParams=optimParams,
                         optimParamsSimple=optimParamsSimple,
                         verbose=verbose,
                         nAttempts=nAttempts,
                         concentrations=concentrations,
                         rates=rates,
                         tpts=tpts,
                         log_shift=log_shift,
                         na.rm=na.rm,
                         sigmoidTotal=sigmoidTotal,
                         nInit=nInit,
                         nIter=nIter,
                         sigmoidSynthesis=sigmoidSynthesis,
                         testOnSmooth=testOnSmooth,
                         estimateRatesWith=estimateRatesWith,
                         .sigmoidModel=.sigmoidModel,
                         .impulseModel=.impulseModel,
                         .sigmoidModelP=.sigmoidModelP,
                         .impulseModelP=.impulseModelP,
                         .polynomialModelP=.polynomialModelP,
                         .rxnrate=.rxnrate,
                         .rxnrateSimple=.rxnrateSimple,
                         .makeModel=.makeModel,
                         .makeSimpleModel=.makeSimpleModel,
                         .logLikelihood=.logLikelihood
  )
  ratesSpecs<-paramSpecs
  return(paramSpecs)
  
  # i=1:nGenes
  # seed=seed
  # .chooseModel=.chooseModel
  # .time_transf=.time_transf
  # .DimpulseModel=.DimpulseModel
  # .DsigmoidModel=.DsigmoidModel
  # .constantModelP=.constantModelP
  # .emptyGene=.emptyGene
  # optimParams=optimParams
  # optimParamsSimple=optimParamsSimple
  # verbose=verbose
  # nAttempts=nAttempts
  # concentrations=concentrations
  # rates=rates
  # tpts=tpts
  # log_shift=log_shift
  # na.rm=na.rm
  # sigmoidTotal=sigmoidTotal
  # nInit=nInit
  # nIter=nIter
  # sigmoidSynthesis=sigmoidSynthesis
  # testOnSmooth=testOnSmooth
  # estimateRatesWith=estimateRatesWith
  # .sigmoidModel=.sigmoidModel
  # .impulseModel=.impulseModel
  # .sigmoidModelP=.sigmoidModelP
  # .impulseModelP=.impulseModelP
  # .polynomialModelP=.polynomialModelP
  # .rxnrate=.rxnrate
  # .rxnrateSimple=.rxnrateSimple
  # .makeModel=.makeModel
  # .makeSimpleModel=.makeSimpleModel
  # .logLikelihood=.logLikelihood
  
}
