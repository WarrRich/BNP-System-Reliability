% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bsp.R
\name{bsp}
\alias{bsp}
\title{Define a Beta-Stacy Process}
\usage{
bsp(
  support,
  centeringMeasure,
  precision,
  calculateMoments = FALSE,
  genSamples = FALSE,
  reps = 10000
)
}
\arguments{
\item{support}{A numeric vector indicating the points at which the support for the BSP is defined}

\item{centeringMeasure}{The mean or centering measure of the BSP at each point on the support
or a function that can be evaluated at each point}

\item{precision}{A constant, vector, or function indicating the precision of the BSP. Precisions cannot be defined at t=0. IMPORTANT NOTE: If heterogeneous precisions are provided, this is defined
just after time t because the BSP precision is left continuous. For example, with support 1,2,3, four precision values are needed. 
If four precisions are supplied, the first precision value will be the precision just after t=0 and the fourth precision value will correspond to the precision just after t=3. 
If only 3 precision values are provided, this function will mark the precision just after t=3 the same as the precision at t=3, unless t=0 is explicitly supplied in the support. 
In that case, the precisions will just be defined after each support point.}

\item{calculateMoments}{This is used to calculate the second moment while
creating the bsp. The moment is used when merging components with other
components. If this function is only being used to create a prior, it can
normally be set to false.}

\item{genSamples}{Logical indicating wither or not to generate samples at each given support value.}

\item{reps}{Only used if genSamples is true. Number of Monte Carlo samples to generate at each support value. Default is 10000}
}
\value{
An object representing a Beta-Stacy Process.
}
\description{
Define a Beta-Stacy Process
}
\details{
The most common use case for this function will be defining BSPs to use as priors.
The precision reflects the uncertainty in the prior. In general, a one unit increase
in the precision is equivalent to one additional observation.
}
\examples{
b <- bsp(support=c(1,3,5), centeringMeasure=c(0.2,0.4, .9), precision=2)



}
