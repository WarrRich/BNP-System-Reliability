#' Define a Beta-Stacy Process
#'
#' @param support A numeric vector indicating the points at which the support for the BSP is defined
#' @param centeringMeasure The mean or centering measure of the BSP at each point on the support
#' or a function that can be evaluated at each point
#' @param precision A constant, vector, or function indicating the precision of the BSP. Precisions cannot be defined at t=0. IMPORTANT NOTE: If heterogeneous precisions are provided, this is defined
#' just after time t because the BSP precision is left continuous. For example, with support 1,2,3, four precision values are needed. 
#' If four precisions are supplied, the first precision value will be the precision just after t=0 and the fourth precision value will correspond to the precision just after t=3. 
#' If only 3 precision values are provided, this function will mark the precision just after t=3 the same as the precision at t=3, unless t=0 is explicitly supplied in the support. 
#' In that case, the precisions will just be defined after each support point. 
#' @param calculateMoments This is used to calculate the second moment while
#' creating the bsp. The moment is used when merging components with other
#' components. If this function is only being used to create a prior, it can
#' normally be set to false.
#' @param genSamples Logical indicating wither or not to generate samples at each given support value. 
#' @param reps Only used if genSamples is true. Number of Monte Carlo samples to generate at each support value. Default is 10000
#'
#' @return An object representing a Beta-Stacy Process.
#' @export
#'
#' @details
#' The most common use case for this function will be defining BSPs to use as priors.
#' The precision reflects the uncertainty in the prior. In general, a one unit increase
#' in the precision is equivalent to one additional observation.
#'
#' @examples
#' b <- bsp(support=c(1,3,5), centeringMeasure=c(0.2,0.4, .9), precision=2)
#'
#'
#'
bsp <- function(support, centeringMeasure, precision, calculateMoments=FALSE, genSamples=FALSE, reps=10000) {

  # Add checks to ensure the input is valid
  # support is a numeric vector (sorted)
  # centeringMeasure is a vector of the same length as support all numbers must be in [0,1]
  #     and nondecreasing
  #     An alternate input is a CDF function that will be evaluated at the
  #     support points
  # precision is a nonnegative vector of numeric values and should be the same length as support
  #     if only one number is supplied it is replicated n times
  #    should these be nonincreasing???????
  if (is.function(centeringMeasure)){
    centeringMeasure<-centeringMeasure(support)
  }

  if (is.unsorted(support)) stop("support should be an increasing series of time points")

  if (length(centeringMeasure)!=length(support))stop("centeringMeasure and support length differ")
  if(any(centeringMeasure>1 |centeringMeasure<0))stop("All centeringMeasure points must be between 0 and 1")

  if (min(support)!=0){
    support<-c(0, support)
    centeringMeasure<-c(0, centeringMeasure)
  }else if (centeringMeasure[1]!=0)stop("Centering measure cannot be greater than 0 when t=0")
  #if (!is.numeric(precision)| length(precision)!=1) 
    #stop("Precision must be a single number")
  #precision<-rep(precision, length(support))

  if (is.function(precision)){
    precision<-precision(support)
    if (length(precision) != length(support))stop("Precision function must evaluate to same length as support")
  }
  if (length(precision)==1){
    precision <- rep(precision, length(support))
  }
  if (length(precision) == (length(support)-1)){
    # This is to make the precision function input friendly to the user. The precision will correspond exactly to the support.
    precision <- c(precision, precision[length(precision)]) # Make the precision after the last point the same as the last point. 
  }
  if (length(precision) != length(support))stop("precision and support length differ. Precision is measured at all points just after t.")

  makeBSP(support, centeringMeasure, precision, calculateMoments, genSamples, reps)

}

#' @keywords internal
makeBSP<-function(support, centeringMeasure, precision, calculateMoments, genSamples=FALSE, reps=10000){
  bsp=structure(list(support=support, centeringMeasure=centeringMeasure,
                 precision=precision), class="betaStacyProcess")
  if (calculateMoments)bsp$E2<-E2(bsp)
  if (genSamples) bsp$Samples<-bspSampling(bsp, reps)
  return(bsp)

}

#' Evaluate centering measure at specific times
#'
#' @param bsp A BSP object
#' @param times An optional vector of times that the centering measure should be determined for. If none provided, the default is over the entire support of the bsp object.
#' @return A named vector of same length as times with the centering measure for each time. The names correspond to the times.
#' @param round Number of decimal places to round the respective times to. A value of -1 means no rounding. Default is 4
#' @param linear.interpolation Should linear interpolation be used to evalue the centering measure between support points? Note: this will not work past the edge of the bsp support. 
#' @export
#'
#'@note For times in between jumps on the support, the centering measure is
#'considered equal to its value after the last jump
#' @examples
#' evaluate_centering_measure(bsp(c(1,2), c(.2,.6), 1), c(.5,1.5,2.5))

evaluate_centering_measure<-function(bsp, times=NULL, round=4, linear.interpolation=FALSE){
  if(inherits(bsp, 'bspPosteriorList')){
    l <- lapply(bsp, function(x) evaluate_centering_measure(x, times=times))
    names(l) <- names(bsp)
    return(l)
  }
  if (is.null(times)){
    res <- bsp$centeringMeasure
    names(res) <- if(round > -1) round(bsp$support,round) else bsp$support
    return(res)
  }
  if (any(times < 0)) stop('Specified times cannot be negative')
  support<-bsp$support
  centeringMeasure<-bsp$centeringMeasure
  indices<-sapply(times, FUN=function(x)sum(support<=x))
  
  if(linear.interpolation){
    delta.x <- diff(support)[indices]
    delta.y <- diff(centeringMeasure)[indices]
    added.vals <- (times - support[indices])*delta.y/delta.x
    interp.vals <- ifelse(is.na(added.vals), 0, added.vals) + centeringMeasure[indices]
    names(interp.vals) <- times
    return(interp.vals)
  }
  
  
  cMeas<-centeringMeasure[indices]
  indices[indices!=0]<-cMeas
  names(indices) <- times
  indices
}

#' Evaluate precision at specific times
#'
#' @param bsp A BSP object
#' @param times An optional vector of times that the precesion should be determined for. If none provided, the default is over the entire support of the bsp object.
#' Default is over the entire support of the bsp.
#' @param round Number of decimal places to round the respective times to. A value of -1 means no rounding. Default is 4
#' @param alphaAfter If TRUE, the returned value will reflect the precision just after the given support times. Default is
#' FALSE. Recall that the precision function is left continous.
#' @return A named vector of same length as times with the precision for each time. The names correspond to the times.
#' @export
#'
#' @note IMPORTANT: The precision function is left continous, unlike the centering measure. For computational efficiency, the precision stored
#' with the bsp class (bsp$precision) gives the precision value epsilon after time t. The evaluate_precision() function gives the precision at exactly time t, unless alphaAfter=TRUE. 
#' @examples
#' evaluate_precision(bsp(c(1,2), c(.2,.6), 1), c(.5,1.5,2.5))

evaluate_precision<-function(bsp, times=NULL, round=4, alphaAfter=FALSE){
  if(inherits(bsp, 'bspPosteriorList')){
    l <- lapply(bsp, function(x) evaluate_precision(x, times=times,round=round, alphaAfter=alphaAfter))
    names(l) <- names(bsp)
    return(l)
  }
  support<-bsp$support
  eps <- 0
  
  if (is.null(times)) times <- support
  
  if(alphaAfter){
    eps <- min(diff(support))/2
    times <- times + eps
  }
  
  support[support==0]=-.1
  precision<-bsp$precision
  s <- sapply(times, FUN=function(t)precision[sum(t>support)][1])
  times <- times - eps
  names(s) <- if(round > -1) round(times,round) else times
  s
}

#' Evaluate second moment at specific times
#'
#' @param bsp A BSP object
#' @param times An optional vector of times that the second moment should be determined for. If none provided, the default is over the entire support of the bsp object.
#' @param round Number of decimal places to round the respective times to. A value of -1 means no rounding. Default is 4
#'
#' @return A named vector of same length as times with the second moment for each time. The names correspond to the times.
#' @export
#'
#' @examples
#' evaluate_precision(bsp(c(1,2), c(.2,.6), 1), c(.5,1.5,2.5))

evaluate_second_moment<-function(bsp, times=NULL, round=4){
  if(inherits(bsp, 'bspPosteriorList')){
    l <- lapply(bsp, function(x) evaluate_second_moment(x, times=times))
    names(l) <- names(bsp)
    return(l)
  }
  if(is.null(bsp$E2)){
    bsp$E2<-E2(bsp)
    warning("Moments for this bsp were not pre-calculated.")
  }
  if (is.null(times)){
    res <- bsp$E2
    names(res) <- if(round > -1) round(bsp$support,round) else bsp$support
    return(res)
  }
  support<-bsp$support
  E2<-bsp$E2
  indices<-sapply(times, FUN=function(x)sum(support<=x))
  sec_mom<-E2[indices]
  indices[indices!=0]<-sec_mom
  indices
}

