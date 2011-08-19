#####################################################################################
## Author: Daniel Sabanes Bove [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: Bayesian FPs for GLMs
## 
## Time-stamp: <[GlmBayesMfpSamples-methods.R] by DSB Die 03/08/2010 13:49 (CEST)>
##
## Description:
## Additional convenience methods for GlmBayesMfpSamples class objects.
##
## History:
## 03/08/2010   file creation with a subset method
#####################################################################################

##' Subset method for GlmBayesMfpSamples objects
##'
##' Index the samples to select a subset of samples.
##'
##' @name GlmBayesMfpSamples-subsetting
##' @aliases [,GlmBayesMfpSamples,ANY,missing,missing-method
##'
##' @usage \S4method{[}{GlmBayesMfpSamples,ANY,missing,missing}(x, i)
##' @param x the \code{\linkS4class{GlmBayesMfpSamples}} object
##' @param i the vector defining the subset of samples
##' @return The subset of the same class.
##' @note The function call will fail if any of the saved bfpCurves or ucCoefs
##' does not have enough samples to be subset by \code{i} !
##'
##' @seealso \code{\linkS4class{GlmBayesMfpSamples}} 
##' @keywords methods
##' @exportMethod "["
setMethod("[",
          signature=
          signature(x="GlmBayesMfpSamples", i="ANY", j="missing", drop="missing"),
          def=
          function(x, i){
              x@fitted <- x@fitted[, i, drop=FALSE]

              if(length(x@predictions))
              {
                  x@predictions <- x@predictions[, i, drop=FALSE]
              }

              x@fixed <- x@fixed[i]
              x@z <- x@z[i]

              bfpCurves <- x@bfpCurves
              for(p in names(bfpCurves))
              {
                  ## save attributes
                  sg <- attr(bfpCurves[[p]], "scaledGrid")
                  wov <- attr(bfpCurves[[p]], "whereObsVals")
                  
                  bfpCurves[[p]] <- bfpCurves[[p]][, i, drop=FALSE]

                  ## get back the attributes
                  attr(bfpCurves[[p]], "scaledGrid") <- sg
                  attr(bfpCurves[[p]], "whereObsVals") <- wov
              }
              x@bfpCurves <- bfpCurves

              ucCoefs <- x@ucCoefs
              for(p in names(ucCoefs))
              {
                  ucCoefs[[p]] <- ucCoefs[[p]][, i, drop=FALSE]
              }
              x@ucCoefs <- ucCoefs

              x@nSamples <- length(x@fixed)              

              return(x)
          })



