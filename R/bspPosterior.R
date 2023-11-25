#' Update a BSP prior with data and return a BSP posterior
#'
#' @param bspPriorObject The BSP object that defines the BSP prior
#' @param data An nx2 matrix of data the first column are the actual data
#'             the second column contains the censoring variable
#'             0 - if right censored 1 - if fully observed
#' @param calculateMoments A flag to indicate whether second moments should be
#' calculated. If the posterior is a component that will be merged with other components
#' this should be TRUE; otherwise, setting to FALSE can save time.
#' @param genSamples Should the resulting posterior object generate monte carlo samples at each support point for the credible intervals?
#' @param suppress.messages Should messages that warn of infinite precision be suppressed? Default FALSE
#'
#' @return An object representing a Beta-Stacy Process.
#' @export
#'
#' @examples
#' prior<-bsp(c(1,5),c(.1, .9),1)
#' data<-matrix(c(5,1,31,1, 2,0), nrow=3, byrow=TRUE)
#' bspPosterior(prior, data)
#'
bspPosterior <- function(bspPriorObject, data, calculateMoments=TRUE, genSamples=FALSE, suppress.messages=FALSE) {
  if(class(bspPriorObject)!="betaStacyProcess")stop("Prior must be class betaStacyProcess")
  ###CHECKS
  if (!is.matrix(data))stop("data must be an nx2 matrix")
  if (nrow(data) == 0) return(bspPriorObject)

  ####ARRANGING DATA
  data<-data[order(data[,1]),, drop=F]
  censor<-data[,2]
  data<-data[,1]
  if (data[1]<=0)stop("All recorded failure times must be greater than 0")
  prior<-bspPriorObject

  #####COMPUTATION
  prior_ts<-prior$support
  prior_gs<-prior$centeringMeasure

  m<-function(data, t) sapply(t, FUN=function(x) sum(data>=x))
  j<-function(data, censor, t) sapply(t, FUN=function(x) sum(censor*(data==x)))
  prior_g<-function(t)sapply(t, FUN=function(t)c(0,prior_gs)[sum(prior_ts<=t)+1])
  allJumps<-sort(unique(c(prior_ts, data)))
  precAtJumps<-evaluate_precision(prior, allJumps)
  GAtJumps<-prior_g(allJumps)
  GBeforeJumps<-c(0, GAtJumps[-length(GAtJumps)])
  allJs<-j(data, censor, allJumps)
  allMs<-m(data, allJumps)

  cumulativeProduct = cumprod(1-
                                ((precAtJumps)*(GAtJumps-GBeforeJumps)+allJs)/
                                (precAtJumps*(1-GBeforeJumps)+allMs)
  )

  centeringMeasure <- 1-cumulativeProduct
  alpha <- (precAtJumps*(1-GAtJumps)+allMs-allJs)/(1-centeringMeasure)
  isna <- any(is.na(alpha))
  # In case where posterior precision is 1, get alpha from just before last jump 
  if(isna){ 
    first.na <-which(is.na(alpha))[1]
    newPoint<-mean(allJumps[(first.na-1):first.na]) 

    new.alpha<-(evaluate_precision(prior, newPoint)*(1-prior_g(newPoint))+
                  m(data, newPoint)-
                  j(data, censor, newPoint))/
      (1-centeringMeasure[first.na-1]) 
  }
  lastTime <- allJumps[length(allJumps)]+1
  lastAlpha<-(evaluate_precision(prior, lastTime)*(1-prior_gs[length(prior_gs)]))/#+m(data, lastTime)-j(data, censor, lastTime))/
    (1-centeringMeasure[length(centeringMeasure)])

  alphaAfter<-c(alpha[-1], lastAlpha)
  names(alphaAfter) <- allJumps
  if(isna) alphaAfter[first.na-1]<-new.alpha
  inf.precs <- is.nan(alphaAfter)
  if(any(inf.precs)){
    alphaAfter[which(inf.precs)] <- Inf
    if(!suppress.messages) message('At least one precision value is Infinity. Likely due to a prior precision of 0.')
  }
  return(makeBSP(allJumps, centeringMeasure, alphaAfter, calculateMoments, genSamples))

}
