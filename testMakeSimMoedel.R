setMethod('makeSimModel', 'INSPEcT', 
          function(object, nGenes, newTpts=NULL, 
                   probs=c(constant=.5,sigmoid=.3,impulse=.2), na.rm=TRUE, seed=NULL) {
            
            tpts <- object@tpts
            concentrations <- list(
              total=ratesFirstGuess(object, 'total')
              , total_var=ratesFirstGuessVar(object, 'total')
              , preMRNA=ratesFirstGuess(object, 'preMRNA')
              , preMRNA_var=ratesFirstGuessVar(object, 'preMRNA')
            )
            rates <- list(
              alpha=ratesFirstGuess(object, 'synthesis')
              , alpha_var=ratesFirstGuessVar(object, 'synthesis')
              , beta=ratesFirstGuess(object, 'degradation')
              , gamma=ratesFirstGuess(object, 'processing')
            )
            #
            out <- .makeSimData(nGenes, tpts, concentrations, rates
                                , newTpts=newTpts, probs=probs, na.rm=na.rm, seed=seed)
            # arrange simdataSpecs form .makeSimData
            simdataSpecs <- out$simdataSpecs
            simdataSpecs <- lapply(simdataSpecs, function(x) list(x))
            #
            newObject <- new('INSPEcT_model')
            newObject@ratesSpecs <- simdataSpecs
            newObject@params$sim$flag <- TRUE
            newObject@params$sim$foldchange <- out$simulatedFC
            newObject@params$sim$noiseVar <- out$noiseVar
            
            return(newObject)
            
          })
