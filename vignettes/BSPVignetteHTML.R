## ---- results=FALSE-----------------------------------------------------------
library(BnpSysRel)
text <- "S(BackTire, FrontTire):Tires
P(BackBrake, FrontBrake):Brakes
S(Tires, Brakes):Bike"

file <- textConnection(text)

## -----------------------------------------------------------------------------
times <- seq(10, 100, length.out = 10)
tirePrior <- bsp(support = times,
               centeringMeasure = pnorm(times, 50, 15),
               precision = 3)

## ----fig1, fig.height = 4.5, fig.width = 6, fig.align = "center"--------------
plot(tirePrior, withConfInt = T)

## -----------------------------------------------------------------------------
brakePrior <- bsp(support=c(22,44,66), centeringMeasure=c(.25,.5,.75), precision = 1)
priorList <- list(FrontTire = tirePrior, BackTire=tirePrior,
               FrontBrake = brakePrior, BackBrake = brakePrior)

## -----------------------------------------------------------------------------
backTireData <- matrix(c(14, 1,
                       29, 1,
                       67, 1,
                       75, 0,
                       75, 0), byrow = T, nrow=5)

## -----------------------------------------------------------------------------
set.seed(20)
frontTireData <- cbind(rnorm(50, 50, 10), 1)
backBrake <- cbind(rnorm(30, 20, 5), 1)
frontBrake <- cbind(rchisq(40, 20), 1)
bike <- cbind(c(40,47,42,46), 1)
dataList <- list(BackTire=backTireData,
              FrontTire=frontTireData,
               BackBrake=backBrake,
               FrontBrake=frontBrake,
               Bike=bike)

## -----------------------------------------------------------------------------
posteriors <- estimateSystemReliability(file=file, priorList=priorList, dataList=dataList)

## ----fig3, fig.height = 5, fig.width = 7, fig.align = "center", message=FALSE----
plotBspGrid(posteriors[1:6], nrow=2, titles=c(names(posteriors)[1:6]), withConfInt=T, withLegend=F)

## ----fig2, fig.height = 4.5, fig.width = 6, fig.align = "center"--------------
plot(posteriors$Bike, withConfInt=TRUE, title='Bike')

## -----------------------------------------------------------------------------
evaluate_centering_measure(posteriors$Bike, times=c(20,50))
bspConfint(posteriors$Bike, times=c(20,50), cred.level = .95)

## -----------------------------------------------------------------------------
quantile(posteriors$Bike, probs=c(.1,.9))

## -----------------------------------------------------------------------------
BackTirePosterior <- bspPosterior(priorList$BackTire, dataList$BackTire)

## ----results='hide'-----------------------------------------------------------
t1 <- print(posteriors$BackTire)
t2 <- print(BackTirePosterior)

## -----------------------------------------------------------------------------
knitr::kable(t1)
knitr::kable(t2)

