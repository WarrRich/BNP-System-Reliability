#' Converts the bsp output (support, centeringMeasure, precision) from a bsp object to a data frame
#'
#'@param x The bsp object to get the output from
#'@param row.names Optional row names for the new data frame.
#'@param optional Part of as.data.frame() signature. Ignored.
#'@param alphaAfter Should the precision be evaluated just after time t? Default TRUE
#'@param ... further arguments to be passed to or from methods
#' @return A data frame of parameters for the given bsp object
#' @export
#'
#'@details This function Converts the bsp output (support, centeringMeasure, precision) from a bsp object to a data frame.
#'@note The precision returned by this object is the precision JUST AFTER time t. 
as.data.frame.betaStacyProcess <- function(x, row.names=NULL, optional=FALSE, alphaAfter=TRUE, ...){
  precision <- if(alphaAfter) x$precision else evaluate_precision(x)
  m <- cbind.data.frame(x$support, x$centeringMeasure, precision)
  colnames(m) <- c('support', 'centeringMeasure', 'precision')
  row.names(m) <- NULL
  m
}

#' Prints a Beta-Stacy process object as a data frame
#'
#'@param x The bsp object to be printed
#'@param alphaAfter Should the precision be evaluated just after time t? Default TRUE
#'@param ... further arguments to be passed to or from methods
#' @return none
#' @export
#'
#'@note The precision printed by this object is the precision JUST AFTER time t, unless alphaAfter is FALSE. 
#' @examples
#' print(bsp(c(1,5),c(0.2,0.4),1))
#'
#'
print.betaStacyProcess<-function(x, alphaAfter=TRUE, ...){
  print(as.data.frame(x, alphaAfter=alphaAfter))
}

#' Plots a Beta-Stacy process over its support
#' @param x The bsp object to be plotted
#' @param withConfInt If true, will calculate and plot approximate 95\% confidence intervals, default FALSE
#' @param ncdfs Plots ncdfs CDFs in gray. If FALSE or 0, no cdfs will be plotted. Default is FALSE
#' @param cred.level Level for credible interval, if none provided, default is 0.95
#' @param withLegend Only when withConfInt=TRUE. Should the legend be included on the plot?
#' @param reps Only valid when withConfint is TRUE. Generates max(reps, ncdfs) samples. Default is 1000
#' @param title Title for the plot. Default is an empty string
#' @param xlim X-limits of the plot
#' @param xlab label for the x axis. Default is "Time (t)"
#' @param ylab label for the y axis. Defualt is "Centering Measure G(t)"
#' @param ... further arguments to be passed to or from methods
#' @return A ggplot object
#' @export
#'
#' @examples
#' plot(bsp(c(1,5,10),c(0.2,0.4, .8),1))
#'
#'
plot.betaStacyProcess<- function(x, withConfInt=FALSE,
                                 ncdfs=FALSE, cred.level=NULL, withLegend=T, reps=1000, title='', xlim=NULL, xlab="Time (t)", ylab="Centering Measure G(t)", ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package ggplot2 needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  ribbons <- requireNamespace("pammtools", quietly = TRUE)
  if(!ribbons & withConfInt){
    warning('Package pammtools not installed. Will use lines instead of ribbons for credible intervals')
  }
  if(ncdfs < 0) stop("The number of CDFs must be a non-negative integer!")
  if(!is.null(cred.level) & !withConfInt) warning('Credible level not used')
  cred.level <- if(is.null(cred.level)) .95 else cred.level
  if (cred.level<=0 | cred.level>=1)stop("Confidence must be in (0, 1)")
  conf.text <- paste0(round(cred.level*100,1),'% Credible Interval')
  times <- centeringMeasure <- support <- lower <- upper <- NULL # For package checking purposes
  t.after <- diff(range(x$support))*.04 + max(x$support)
  line.support <- c(x$support[-1], t.after)
  
  data <- data.frame(times=line.support, centeringMeasure=evaluate_centering_measure(x, line.support))
  p<-ggplot2::ggplot() +
    ggplot2::geom_step(data=data, ggplot2::aes(x=times, y=centeringMeasure),size=1)+
    ggplot2::geom_point(data=data[-nrow(data),], ggplot2::aes(x=times, y=centeringMeasure))+
    ggplot2::xlab(xlab)+ggplot2::ylab(ylab)
  #add some precision stuff
  if(withConfInt | ncdfs){
    # Use existing samples if already stored
    samples <- if(is.null(x$Samples)){
      bspSampling(x, reps=ifelse(ncdfs & !withConfInt, ncdfs, max(reps, round(ncdfs))))
      }else{
        nsamps <- nrow(x$Samples)
        if(nsamps < ncdfs){
          rbind(x$Samples, bspSampling(x, reps=ncdfs-nsamps))
        }else{x$Samples}
      }
    if(withConfInt){
    x$Samples <- samples # Stores samples inside the object
    ls <- length(line.support)
    cred_int <- data.frame(support=line.support, lower=NA, upper=NA)
    cred_int[1:(ls-1),2:3] <- bspConfint(x, cred.level=cred.level, cols=TRUE)[-1,-1]
    cred_int[ls, 2:3] <- cred_int[ls-1, 2:3]
    
      if(ribbons){
        p <- p + pammtools::geom_stepribbon(data=cred_int,ggplot2::aes(x=support, ymin=lower, ymax=upper, color=conf.text),fill='red',alpha=.2)
      }else{
        p<-p+ggplot2::geom_step(data=cred_int,
                              ggplot2::aes(x=support,y=lower, color=conf.text))+
          ggplot2::geom_step(data=cred_int,
                           ggplot2::aes(x=support,y=upper, color=conf.text))
      }
    }
    if (ncdfs){
      paths=cbind.data.frame(samples, samples[, ncol(samples)])
      
      for (i in 1:round(ncdfs)){
        p <- p + ggplot2::geom_line(data=NULL,
                               ggplot2::aes_string(y=as.numeric(paths[i,]), x=line.support),
                              inherit.aes = F,
                              alpha=.2, show.legend = F)
        # cat('Line ', i, '\n')
      }
    }
  }
  scale <- if(is.null(xlim)){
      ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(.05, 0)))
  }else{
    ggplot2::coord_cartesian(xlim = xlim)
    }
  p <- p + 
    ggplot2::theme_bw()  +
    ggplot2::theme(text = ggplot2::element_text(size=14), legend.text = ggplot2::element_text(size = 11)) +
    ggplot2::ggtitle(title) + scale + 
    if (withLegend){
              ggplot2::theme(legend.title = ggplot2::element_blank(), legend.position='bottom')
            }else ggplot2::theme(legend.position = "none")
  return(p)
}

#' Calculates the expected time to failure for various quantiles
#'
#'@param x The bsp object to calculate from
#'@param probs numeric vector of probabilities
#'@param ... further arguments to be passed to or from methods
#' @return a named numeric vector of times associated with each quantile
#' @export
#'
#'@details This function calculates the expected failure time for a given quantile from
#'the centering measure of the BSP.
#' @examples
#' bsp <- bsp(1:3, c(.25,.5,.75), 1)
#' quantile(bsp, c(.333, .666))
#'
quantile.betaStacyProcess<-function(x, probs, ...){
  qs<-sapply(probs, function(p)x$support[sum(x$centeringMeasure<=p)])
  names(qs)<-paste0(probs*100,"%")
  qs
}

#' Gives the median time to failure of a BSP object
#' @importFrom stats median
#' @param x The bsp object to calculate from
#' @param na.rm Part of median() signature. Ignored.
#' @param ... further arguments to be passed to or from methods
#' @return Median time to Failure
#' @export
#' @details This function calculates the median failure time for a given BSP object
median.betaStacyProcess <- function(x, na.rm=TRUE, ...){
  x$support[which.max(which(x$centeringMeasure<=.5))]
}

#' Gives the mean time to failure of a BSP object
#' @param x The bsp object to calculate from
#' @param ... further arguments to be passed to or from methods
#' @return Mean time to Failure
#' @export
#' @details This function calculates the expected failure time for a given BSP object
mean.betaStacyProcess <- function(x, ...){
  sum(diff(x$centeringMeasure)*x$support[-1])
}

#' Gives the standard deviation of time to failure of a BSP object
#' @param bsp The bsp object to calculate from
#' @return The standard deviation of time to failure
#' @export
#' @details This function calculates the standard deviation of the failure time for a given BSP object
bspSdFailureTime <- function(bsp){
  sqrt(sum(diff(bsp$centeringMeasure)*bsp$support[-1]^2)-mean(bsp)^2)
}

#' Gives summary statistics of time to failure of a BSP object
#' @importFrom stats quantile
#' @param object The bsp object to calculate from
#' @param quants Optional failure time quantiles to include in the summary. 
#' Default is c(.05, .1, .25, .5, .75, .9, .95)
#' @param ... further arguments to be passed to or from methods
#' @return Summary statistics for the time to failure
#' @export
#' @details This function calculates summary statistics of the failure time for a given BSP object
summary.betaStacyProcess<-function(object, quants=c(.05, .1, .25, .5, .75, .9, .95), ...){
  c(quantile.betaStacyProcess(object, quants), 'mean'=mean(object), 'sd'=bspSdFailureTime(object))
  #print(paste("Median time to failure:", round(x$support[which.max(which(x$centeringMeasure<=.5))], 3)))
}

