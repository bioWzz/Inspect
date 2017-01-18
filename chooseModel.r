.chooseModel <- function(tpts, log_shift, experiment, variance=NULL, na.rm=TRUE
                         , sigmoid=TRUE, impulse=TRUE, polynomial=TRUE, nInit=10, nIter=500
                         , .time_transf, .impulseModel, .sigmoidModel, .sigmoidModelP, .impulseModelP
                         , .polynomialModelP)
  #### choose a functional form between impulse and sigmoid according 
  #### to the one that has the gratest pvalue in the chi squared test
{
  
  
  
  chisq.test.default <- function(experiment, model, variance=NULL, df)
  {
    if( is.null(variance) ) variance <- stats::var(experiment )
    D = chisq(experiment, model, variance)
    modelDF <- max(0, length(experiment)-df)
    pchisq(D,  modelDF, lower.tail=TRUE)
  }
  
  optimFailOut <- function(e)
    list(par=NA, value=NA, counts=NA, convergence=1, message=e)
  #
  # impulse model functions
  #
  im.parguess <- function(tpts , values ) {
    # values = expressions.avgd(eD)
    # tp = tpts(eD)
    ntp   <- length(tpts)
    peaks <- which(diff(sign(diff(values)))!=0)+1
    if( length(peaks) == 1 ) peak <- peaks
    if( length(peaks)  > 1 ) peak <- sample(peaks, 1)
    if( length(peaks) == 0 ) peak <- round(length(tpts)/2)
    #
    initial_values <- runif( 1, min=min(values[1:3])
                             , max=max(values[1:3]))
    intermediate_values <- values[peak]
    if( intermediate_values==0 ) intermediate_values <- mean(values[seq(peak-1,peak+1)])
    end_values <- runif( 1, min=min(values[(ntp-2):ntp])
                         , max=max(values[(ntp-2):ntp]))
    time_of_first_response  <- tpts[peak-1]
    time_of_second_response <- tpts[peak+1]
    slope_of_response <- diff(range(tpts)) / 
      (time_of_second_response-time_of_first_response)
    #
    return(c(h0=initial_values, h1=intermediate_values
             , h2=end_values, t1=time_of_first_response
             , t2=time_of_second_response, b=slope_of_response))
  }
  #
  im.chisq <- function(par, tpts, experiment, variance=NULL, .impulseModel) 
  {
    model <- .impulseModel(tpts, par)
    chisq(experiment, model, variance)
  }
  #
  im.optim.chisq <- function(tpts, experiment, variance=NULL, ninit=10
                             , maxit=500) 
    sapply(1:ninit, function(x) 
    
        
        optim(
        par=im.parguess(tpts,experiment)
        , fn=im.chisq, tpts=tpts
        , experiment=experiment
        , variance=variance
        , .impulseModel=.impulseModel
        , control=list(maxit=maxit)
      ))
      
      
  #
  # sigmoid model functions
  #
  sm.parguess <- function(tpts , values ) {
    # values = expressions.avgd(eD)
    # tp = tpts(eD)
    
    time_span <- diff(range(tpts))
    # sample the time uniformely
    time_of_response <- runif( 1, min=min(tpts), max=max(tpts))
    # slope of response must be high if the time of response is close to one
    # of the two boundaries
    distance_from_boundary <- min(time_of_response - min(tpts)
                                  , max(tpts) - time_of_response)
    slope_of_response <- time_span / distance_from_boundary
    ntp   <- length(tpts)
    initial_values <- runif( 1, min=min(values[1:3])
                             , max=max(values[1:3]))
    end_values <- runif( 1, min=min(values[(ntp-2):ntp])
                         , max=max(values[(ntp-2):ntp]))
    #
    return(c(h0=initial_values, h1=end_values, t1=time_of_response
             , b=slope_of_response))
  }
  #
  sm.chisq <- function(par, tpts, experiment, variance=NULL, .sigmoidModel) 
  {
    model <- .sigmoidModel(tpts, par)
    chisq(experiment, model, variance)
  }
  #
  sm.optim.chisq <- function(tpts, experiment, variance=NULL, ninit=10
                             , maxit=500) 
    sapply(1:ninit, function(x) 
      tryCatch(optim(
        par=sm.parguess(tpts, experiment)
        , fn=sm.chisq, tpts=tpts
        , experiment=experiment
        , variance=variance
        , .sigmoidModel=.sigmoidModel
        , control=list(maxit=maxit)
      ), error=function(e) optimFailOut(e))) 
  
  pn.optim.aic <- function(tpts , experiment, variance=NULL) {
    if( length(experiment) < 3 ) return(NA)
    polyOrderChisq <- function(i) {
      model <- lm(experiment~poly(tpts, i, raw=TRUE ))
      return(list(par=model$coeff, value=AIC(model)))}
    return(sapply(1:min(7,length(tpts)-2), polyOrderChisq))
  }
  
  # remove missing values
  if( na.rm) {
    idx <- is.finite(experiment)
    tpts <- tpts[idx]
    experiment <- experiment[idx]
  }
  
  ## 
  if( length(experiment)==0 ) {
    stop('.chooseModel: no time points have a finite value.
         Impossible to evaluate any kind of model.')
    return(list(type='constant', fun=.constantModelP
                , params=mean(experiment, na.rm=TRUE), pval=NA, df=1))
  }
  ## 
  if( length(experiment)<=2 ) {
    warning('.chooseModel: less than three time points have a finite value. 
            Impossible evaluate a variable model.
            Returning a constant model.')
    return(list(type='constant', fun=.constantModelP
                , params=mean(experiment, na.rm=TRUE), pval=NA, df=1))
  }
  
  ## re-evaluate flags of function to evaluate according to the lenght 
  ## of the experiment
  sigmoid <- sigmoid
  impulse <- impulse & length(experiment)>2
  polynomial <- polynomial & length(experiment)>2
  
  tptslog <- .time_transf(tpts, log_shift)
  
  # sigmoid
  if( sigmoid ) {
    outSM  <- sm.optim.chisq(tpts=tptslog, experiment=experiment
                             , variance=variance, ninit=nInit, maxit=nIter)
    bestSM <- which.min(unlist(outSM[2,]))
    pvalSM <- chisq.test.default(experiment=experiment
                                 , model=.sigmoidModel(tptslog, outSM[,bestSM]$par)
                                 , variance=variance, df=length(outSM[,bestSM]$par))
    dfSM <- length(outSM[,bestSM]$par)
  } else dfSM <- NA
  # impulse
  if( impulse ) {
    outIM  <- im.optim.chisq(tpts=tptslog, experiment=experiment, 
                             variance=variance, ninit=nInit, maxit=nIter)
    bestIM <- which.min(unlist(outIM[2,]))
    pvalIM <- chisq.test.default(experiment=experiment
                                 , model=.impulseModel(tptslog, outIM[,bestIM]$par) 
                                 , variance=variance, df=length(outIM[,bestIM]$par))
    dfIM <- length(outIM[,bestIM]$par)
  } else dfIM <- NA
  
  # polynomial
  if( polynomial ) {
    outPN  <- pn.optim.aic(tptslog, experiment, variance )
    bestPN <- which.min(unlist(outPN[2,]))
    pvalPN <- chisq.test.default(experiment=experiment
                                 , model=.polynomialModel(tptslog, outPN[,bestPN]$par) 
                                 , variance=variance, df=length(outPN[,bestPN]$par))
    dfPN <- length(outPN[,bestPN]$par)
  } else dfPN <- NA
  
  pvals <- c(
    sigmoid=if( sigmoid ) pvalSM else NA
    , impulse=if( impulse ) pvalIM else NA
    , polynomial=if( polynomial ) pvalPN else NA
  )
  funcs <- c(.sigmoidModelP, .impulseModelP, .polynomialModelP)
  dfs <- c(dfSM, dfIM, dfPN)
  type <- names(pvals)[which.max(pvals)]
  df   <- dfs[which.max(pvals)]
  
  if( type=='sigmoid'    ) params <- outSM[,bestSM]$par
  if( type=='impulse'    ) params <- outIM[,bestIM]$par
  if( type=='polynomial' ) params <- outPN[,bestPN]$par
  
  pval <- pvals[which.max(pvals)]
  fun  <- funcs[[which.max(pvals)]]
  
  return(list(type=type, fun=fun , params=params, pval=pval, df=df))
  
}