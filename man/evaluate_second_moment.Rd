% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bsp.R
\name{evaluate_second_moment}
\alias{evaluate_second_moment}
\title{Evaluate second moment at specific times}
\usage{
evaluate_second_moment(bsp, times = NULL, round = 4)
}
\arguments{
\item{bsp}{A BSP object}

\item{times}{An optional vector of times that the second moment should be determined for. If none provided, the default is over the entire support of the bsp object.}

\item{round}{Number of decimal places to round the respective times to. A value of -1 means no rounding. Default is 4}
}
\value{
A named vector of same length as times with the second moment for each time. The names correspond to the times.
}
\description{
Evaluate second moment at specific times
}
\examples{
evaluate_precision(bsp(c(1,2), c(.2,.6), 1), c(.5,1.5,2.5))
}
