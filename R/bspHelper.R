# Helper functions

#' Plots a list of bsp objects in a grid-like pattern
#' @importFrom utils tail
#' @importFrom graphics plot
#' @param l The list of bsp objects to be plotted. These are plotted row-wise.
#' @param nrow (Optional) Number of rows in the plotted grid
#' @param ncol (Optional) Number of columns in the plotted grid
#' @param titles A vector of titles for the plots. Default paste('Comp', 1:length(l))
#' @param xlims list of x-axis limits to use for each respective plot. If an element is NULL, make no changes to the xlimits. 
#' The default is a list of NULLs.
#' @param priors Optional list of priors that are plotted on top of the bsp objects. If provided, it will plot a blue line indicating the prior centering measure. Note:
#' This only works as long as each prior is explicitly defined (including any subsystem components). 
#' @param data Optional list of data matrices who's kaplan-meier curves are plotted on top of the bsp objects. If provided, it will plot dashed black lines for each object.
#' @param prior.col Optional color string for the prior line. Only used if priors is provided.
#' @param ... further arguments to be passed to the plot.betaStacyProcess() function. Allows the user to suppress the legend or add confidence bands, for example.
#' @return A gtable grid of bsp ggplots. 
#' @note This function requires the gridExtra package. If the data is plotted, it requires the Survival package.
#' @export
#'
#' @examples
#' set.seed(0)
#' C.A <- cbind(rweibull(15, 3, 10), rbinom(15, 1, .6))
#' C.B <- cbind(rweibull(10, 4, 10), rbinom(10, 1, .6))
#' C.C <- cbind(rweibull(20, 5, 10), rbinom(20, 1, .6))
#' C <- list(C.A, C.B, C.C)
#' 
#' ps <- lapply(1:length(C), function(x) bsp(c(8,9,10,11), c(.1, .3, .6, .9), 5))
#' sys.names <- paste0(paste0('c',1:3))
#' names(ps) <- names(C) <- sys.names
#' txt <- paste0('S','(',paste0('c',1:3, collapse = ',') ,'):System')
#' post <- estimateSystemReliability(textConnection(txt), priorList = ps, dataList = C)
#' 
#' plotBspGrid(post[1:3], ncol=2, withConfInt=TRUE, withLegend=FALSE, priors=ps, data=C, ylab='G(t)')
#'
plotBspGrid <- function(l, nrow=NULL, ncol=NULL, titles=paste('Comp', 1:length(l)), xlims=rep(NULL, length(l)), priors=NULL, data=NULL, prior.col='blue', ...){
  if(!requireNamespace('gridExtra', quietly=FALSE)) {    
    stop("Package gridExtra needed for this function to work. Please install it.", call. = FALSE)}
  if(!requireNamespace('ggplot2', quietly=FALSE)) {    
    stop("Package ggplot2 needed for this function to work. Please install it.", call. = FALSE)}
  if(length(titles)==1) titles<-rep(titles, length(l))
  if(length(xlims)==1) titles<-rep(xlims, length(l))
  prior.lines <- if(!is.null(priors)){
    if(length(priors) != length(l)) stop('List lengths between the priors and posteriors do not match')
    prior.lines <- lapply(seq_along(l), function(p) 
      ggplot2::geom_step(ggplot2::aes(x=c(min(l[[p]]$support[-1]), priors[[p]]$support[-1], max(l[[p]]$support)*1.04), 
                                      y=c(0, priors[[p]]$centeringMeasure[-1],tail(priors[[p]]$centeringMeasure,1))), colour=prior.col, lwd=.45))
  }
  dat.lines <- if(!is.null(data)){
    if(length(data) != length(l)) stop('List lengths between the data and posteriors do not match')
    if(!requireNamespace('survival', quietly=FALSE)) {    
      stop("Package survival needed for this function to work. Please install it.", call. = FALSE)}
    km.est <- function(d){
      surv <- survival::Surv(d[,1], d[,2], type='right')
      km <- survival::survfit(surv ~ 1, se.fit=FALSE) 
      cbind(km$time,1-km$surv)
    }
    dats <- lapply(data, km.est)
    lapply(seq_along(dats), function(d) 
      ggplot2::geom_step(ggplot2::aes(x=dats[[d]][,1], 
                             y=dats[[d]][,2]), lwd=.6, colour='black', linetype='dashed'))
  }
  
  plts <- lapply(seq_along(l), function(x) plot(l[[x]], title=titles[x], xlim=xlims[[x]], ...) + prior.lines[[x]] + dat.lines[[x]])
  gridExtra::grid.arrange(grobs=plts,nrow=nrow, ncol=ncol)
}


#' Creates a BSP prior from some data. User is able to shift/scale the prior and set its precision.
#' @param data The data which to create the prior from. Needs to be a matrix with 2 columns for the failure times and censored indicator
#' @param EF A modified Effectiveness Factor (more info: http://reliawiki.org/index.php/Crow_Extended). A number between -1 and 1.
#' A value of 0 indicates no shifting in the support. A non-negative value (EF) indicates a system that's 1/(1-EF) as good. A negative value
#' indicates a system that is 1/(1-|EF|) times as worse. Mathematically, this scales the support by 1/(1-EF) if EF is positive
#' or 1-|EF| if EF is negative. Default is 0 (or no scaling)
#' @param centering.max Rescales all centering measure values such that the centering measure at the maximum support is centering.max.
#' Default is NULL, which means no rescaling.
#' @param prior.prec A single number or vector indicating the precision to set for the bsp() object. Recall that this is the precision
#' just after each support point. Default is 1. If NULL, use the precisions generated from a zero precision prior with the data provided.
#' @param suppress.messages If prior.prec is NULL, Sometimes messages are printed due to a zero precision prior. The default is FALSE, 
#' but if this is TRUE, it suppress them. 
#' @return A scaled or shifted bsp object created from provided data. Can be used as a prior. 
#' @export
#'
#' @examples
#' d<-matrix(c(5,1,15, 1, 31,1, 2,0, 20, 0), ncol=2, byrow=TRUE)
#' bspPrior(d, EF=2/3, prior.prec=0.5) # Three times as good
#' 
bspPrior <- function(data, EF=0, centering.max=NULL, prior.prec=1, suppress.messages=FALSE){
  l <- length(prior.prec)
  if(!(l %in% c(0,1,nrow(data)+1))) stop('Prior precision must be either a constant or of length nrow(data)+1')
  if(abs(EF) >= 1) stop("Effectiveness factor should be between -1 and 1")
  if(!is.null(centering.max)) if( centering.max > 1 | centering.max < 0) stop('Maximum centering measure should be between 0 and 1')
  
  if(is.null(prior.prec)){
    p <- bspPosterior(createUninformativePrior(data), data, suppress.messages=suppress.messages)
    prior.prec <- p$precision
  }else{
    p <- suppressMessages(bspPosterior(createUninformativePrior(data), data))
    if(length(prior.prec) == 1){
      prior.prec <- rep(prior.prec, length(p$support))
    }
  }
  cent <- p$centeringMeasure
  if(!is.null(centering.max)){
    cent <- cent/max(cent)*centering.max
  }
  bsp(p$support*(1-abs(EF))^ifelse(EF >= 0, -1, 1), cent, prior.prec)
}