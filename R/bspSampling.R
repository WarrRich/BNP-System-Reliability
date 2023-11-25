##BSPSampling.R
#'Creates a matrix where each row represents a sampled CDF from the provided BSP
#'@importFrom stats rbeta
#'@param bsp Any bsp object
#'@param reps Number of samples to draw, default 10000
#'@param seed Seet to set before generating the samples. If NULL, don't set the seed. Deafult 0. 
#'
#'@return A matrix. Each row can be plotted against
#'the BSP's support to view a sampled CDF from the BSP.
#'Column names are the support point for that column. Use samplesAt to
#'get extract sampled distribution at single point
#'@export
#'
#' @examples
#' bsp<-bsp(c(1:3), centeringMeasure = c(.1,.5, .98), precision = 4)
#' samples<-bspSampling(bsp, reps=200)
#' #To visualize
#' matplot(matrix(bsp$support[-1], nrow=3), t(samples), type='l', xlab='support', ylab='Centering Measure')
bspSampling<-function(bsp, reps=10000, seed=0){
  if(reps < 1) stop('Reps must be a positive integer')
  support <- bsp$support[-1]
  draws<-matrix(NA,nrow = reps, ncol=length(support))
  colnames(draws)<-support
  center <- evaluate_centering_measure(bsp, support)
  precs <- evaluate_precision(bsp, support)
  pmf <- diff(bsp$centeringMeasure)
  if(!is.null(seed)) set.seed(seed)
  
  draws[,1]<- rbeta(reps, precs[1]*pmf[1], precs[1]*(1-center[1]))
  for (i in 2:ncol(draws)){
    alphak <-precs[i]*(pmf[i])
    betak <- precs[i]*(1-center[i])
    epsk <- 1-draws[,i-1]
    draws[,i] <- draws[,i-1] + rbeta(reps, alphak, betak)*epsk
  }
  draws
}

#'Finds approximate (simulated) credible intervals for the percent failed at a given time
#'
#'@param bsp The bsp object
#'@param times Place where you would like credible intervals. If none provided, the default is the entire support of the bsp object. 
#'@param cred.level Credible level for intervals (default is 95\%)
#'@param reps Number of samples at each time point used in calculations for the credible interval. Default is 1000. If the bsp object already has samples stored in bsp$Samples, reps is ignored.
#'@param round Number of decimal places to round the respective times to. A value of -1 means no rounding. Default is -1
#'@param cols (Depricated) Output credible intervals on the columns? (result will be a length(times)x3 data frame with times also as a column). 
#'Default is TRUE. To get credible interval results similar to package version 0.1.0, use FALSE. 
#'@note If you intend to call this function multiple times on the same object,
#'changing the genSamples arguement to TRUE in the bsp() object or calling
#'bsp$Samples<-bspSampling(bsp) will improve performance. 
#'
#'@return A length(times)x3 data frame with the support, lower, and upper bounds of the credible interval on the columns.
#'@export
#'
#' @examples
#' bsp=bsp(c(1:3), centeringMeasure = c(.1,.9, .98), precision = 2)
#' bspConfint(bsp, 1:3)
#'
bspConfint<-function(bsp, times=NULL, cred.level=.95, reps=1000, round=-1, cols=TRUE){
  if (cred.level<0 |cred.level>1)stop("cred.level must be between 0 and 1")
  #if (alpha.level>.6)warning(paste0("Alpha level set to ", alpha.level, ". Are you sure you didn't mean ",
                                   #1-alpha.level,"?"))
  if(cred.level < .3) warning(paste0('Credible level set to ', cred.level,', did you mean ',1-cred.level, '?'))
  if(inherits(bsp, 'bspPosteriorList')){
    l <- lapply(bsp, function(x) bspConfint(x, times=times, cred.level=cred.level, reps=reps, cols=cols, round=round))
    names(l) <- names(bsp)
    return(l)
  }
  alpha.level <- 1-cred.level
  if(is.null(bsp$Samples)) bsp$Samples <- bspSampling(bsp, reps=reps)
  times <- if(is.null(times)) bsp$support else times
  columns <- sapply(times, FUN=function(time) sum(time>=bsp$support)-1)
  intervals<- sapply(columns, FUN= function(col) quantile(bsp$Samples[,col ], c(alpha.level/2, 1-alpha.level/2)))
  intervals[is.na(intervals)]<-0
  support <- if(round > -1) round(times,round) else times
  if (cols){
    intervals<- as.data.frame(t(intervals))
    return(cbind.data.frame(support,intervals))
  }
  colnames(intervals) <- support
  intervals
}

#'Extract distribution at select point
#'
#'@param samples Either a BSP object or the samples produced by bspSampling(). If the bsp object doesn't have samples
#'stored, it will generate them.
#'@param time The spot on the support at which you would like to examine the distribution
#'
#'@return A vector of samples drawn from that point on the support
#'
#'@export
#'
#'@examples
#'bsp <- bsp(1:3, c(.2,.5,.8), 3)
#'samples <- bspSampling(bsp)
#'marginal <- samplesAt(samples, 2.5)
samplesAt<-function(samples, time){
  if(inherits(samples, 'betaStacyProcess')){
    samples <- if(is.null(samples$Samples)) bspSampling(samples) else samples$Samples
  }
  times<-as.numeric(colnames(samples))
  index<-sum(times<=time)
  if(index==0)stop("Time less than all of support")
  return(samples[,index])
}


# #'Finds approximate confidence interval by fitting a beta distribution to the BSP moemnts
# #'
# #'
# bspConfint2<-function(bsp,times, alpha.level=.05){
#   E1<-evaluate_centering_measure(bsp, times)
#   E2<-evaluate_second_moment(bsp, times)
#   v=E2-E1^2
#   alpha<-E1*(E1*(1-E1)/v-1)
#   beta<-alpha*(1/E1-1)
#   uppers<-qbeta(1-alpha.level/2, alpha, beta)
#   lowers<-qbeta(alpha.level/2, alpha, beta)
#  return(rbind(lowers,uppers))
#
# }
