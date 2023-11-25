
#' Calculates the second moment of a bsp object at each spot on the support. Internal function
#'
#' @param bsp The bsp object
#'
#' @return Vector of second moments that aligns with bsp$support
#' @keywords internal
#' @examples
#' bsp=bsp(c(1:3), centeringMeasure = c(.1,.9, .98), precision = 2)
#' bsp$E2=BnpSysRel:::E2(bsp)
#' bsp$E2
#'
E2 <- function(bsp) {
  base<-bsp$centeringMeasure
  prec<-evaluate_precision(bsp, bsp$support)
  n<-length(base)
  ## Calculate the 2nd Moment
  temp <- ((1-base[-1])*(prec[-1]*(1-base[-1])+1))/
    ((1-base[-n])*(prec[-1]*(1-base[-n])+1))
  E2 <- cumprod(temp)-1+2*base[-1]
  #if (any(is.nan(E2)))print(paste("E2 had a NAN", E2))
  #E2[is.nan(E2)] <- 0

  return(c(0, E2))

}

#' Creates a bsp with given moments at each spot on the support
#'
#' @param mList A list with three elements, all equal length, containing the
#' support, E1, and E2 for the new bsp
#'
#' @return A bsp object with those moments
#' @export
#'
#' @examples
#' bsp=bspFromMoments(list(E1=c(.2, .4,.8), E2=c(.5, .5, .5), support=1:3))
#' bsp
#'
#'
bspFromMoments <- function(mList) {
  E2<-mList$E2
  E1<-mList$E1
  support<-mList$support

  pE1 <- c(0,E1[-length(support)])
  pE2 <- c(0,E2[-length(support)])

  ## Precision
  prec <- ((pE2-2*pE1+1)*(1-E1)-(E2-2*E1+1)*(1-pE1))/((E2-2*E1+1)*(1-pE1)^2-(pE2-2*pE1+1)*(1-E1)^2)
  ## Return list base, prec, and times
  prec<-c(prec[-1], prec[length(prec)]) # Makes it return the precision just after t
  hasNAs<-is.nan(prec)
  if (any(hasNAs)){
    warning("NAN in bspFromMoments()")
    prec[hasNAs]<-0
  }
  #for (i in 2:length(support)) {
   # if(is.nan(prec[i])) {
    #  prec[i]=0
     # warning("NAN in bspFromMoments()")
      #}
  #}
  prec[prec>999999999]<-999999999
  negative.precs <- prec<0
  if(any(prec<0)) warning('Method of moments produced negative precision(s). Setting precision(s) to 0.')
  prec[negative.precs]<-0
  
  return(makeBSP(support, E1, prec, calculateMoments = TRUE))
}


#b=E2(bsp)
#varx = b$E2 - b$E1^2


#' Calculates the first and second moments of the subsystem if bsp1 and bsp2 are in series
#'
#' @param bsp1 The first bsp object
#' @param bsp2 The second bsp object
#'
#' @return A list with three elements, all equal length, containing the
#' support, E1, and E2 for the new bsp
#' @export
#'

E1E2_series <- function(bsp1, bsp2) {

  times = unique(sort(c(bsp1$support,bsp2$support)))
  n <- length(times)

  new_C1_E1 <- evaluate_centering_measure(bsp1, times)
  new_C1_E2 <- evaluate_second_moment(bsp1, times)

  new_C2_E1 <- evaluate_centering_measure(bsp2, times)
  new_C2_E2 <- evaluate_second_moment(bsp2, times)

  E1 <- 1-(1-new_C1_E1)*(1-new_C2_E1)

  #Why are we multiplying by 1-new_c1_E1, page nine seems to say times by the base
  E2 <- 1-2*(1-new_C1_E1)*(1-new_C2_E1)+(1-2*new_C1_E1+new_C1_E2)*(1-2*new_C2_E1+new_C2_E2)

  dups<-which(duplicated(E1) & duplicated(E2))
  if(length(dups) > 0){
    E1<-E1[-dups]
    E2<-E2[-dups]
    times<-times[-dups]
  }
  m <- which.max(round(E1, 8))
  #print(m)
  ## Return list E1, E2, and support
  return( list("E1" = E1[1:m], "E2" = E2[1:m], "support"= times[1:m]) )
}

#' Calculates the first and second moments of the subsystem if bsp1 and bsp2 are in parallel
#'
#' @param bsp1 The first bsp object
#' @param bsp2 The second bsp object
#'
#' @return A list with three elements, all equal length, containing the
#' support, E1, and E2 for the new bsp
#' @export
#'
E1E2_parallel <- function(bsp1, bsp2) {

  times = unique(sort(c(bsp1$support,bsp2$support)))
  n <- length(times)

  new_C1_E1 <- evaluate_centering_measure(bsp1, times)
  new_C1_E2 <- evaluate_second_moment(bsp1, times)

  new_C2_E1 <- evaluate_centering_measure(bsp2, times)
  new_C2_E2 <- evaluate_second_moment(bsp2, times)

  E1 <- new_C1_E1*new_C2_E1

  E2 <- new_C1_E2*new_C2_E2

  dups<-which(duplicated(E1) & duplicated(E2))
  if(length(dups) > 0){
    E1<-E1[-dups]
    E2<-E2[-dups]
    times<-times[-dups]
  }

  ExtraZeros=which(round(E1, 8)==0)[-1]
  if(length(ExtraZeros)>0){
    E1<-E1[-ExtraZeros]
    E2<-E2[-ExtraZeros]
    times<-times[-ExtraZeros]
  }
  #m <- n-which.min(rev(E1))+1

  #if (E1[m]==0 & E2[m]==0)m=m+1

  ## Return list E1, E2, and times
  return( list("E1" = E1, "E2" = E2, "support"= times) )

}
