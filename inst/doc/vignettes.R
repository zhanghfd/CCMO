## ----eval=FALSE---------------------------------------------------------------
#  library(CCMO)
#  data("SampleData",package="CCMO")
#  dim(SampleData)
#  head(SampleData)

## ----eval=FALSE---------------------------------------------------------------
#  # library(CCMO)
#  Y = SampleData[,1]
#  Gc = SampleData[,2]
#  Gm = SampleData[,12]
#  X = SampleData[,-(1:21)]
#  fit = singleSNP(Y,Gc,Gm,Xo=X,Xc=X,Xm=X,X.Gm=X)

## ----eval=FALSE---------------------------------------------------------------
#  names(fit)

## ----eval=FALSE---------------------------------------------------------------
#  fit$new

## ----eval=FALSE---------------------------------------------------------------
#  fit$cov.new

## ----eval=FALSE---------------------------------------------------------------
#  fit$log
#  fit$cov.log

## ----eval=FALSE---------------------------------------------------------------
#  fit = OmnibusTest(fit,test=7:10)
#  fit$Omnibus

## ----eval=FALSE---------------------------------------------------------------
#  Gc = SampleData[,2:11]
#  Gm = SampleData[,12:21]
#  system.time(fit1 <- multipleSNP(Y,Gc,Gm,Xo=X,Xc=X,Xm=X,X.Gm=X,cores=1))
#  system.time(fit2 <- multipleSNP(Y,Gc,Gm,Xo=X,Xc=X,Xm=X,X.Gm=X,cores=2))

## ----eval=FALSE---------------------------------------------------------------
#  fit2[[2]]

## ----eval=FALSE---------------------------------------------------------------
#  fit <- multipleSNP(Y,Gc,Gm,Xo=X,Xc=X,Xm=X,X.Gm=X,test=7:10)
#  fit[[2]]$Omnibus

## ----eval=FALSE---------------------------------------------------------------
#  # library(CCMO)
#  Y = POESampleData[,1]
#  gmm = POESampleData[,2:6]
#  gcc = POESampleData[,7:11]
#  X = POESampleData[,12]
#  data = MultiLociPOE.input(gmm,gcc,0)
#  gmm = data$gmm
#  gcc = data$gcc
#  hap = data$hap
#  ppi = data$ppi
#  loci = 1
#  f = 0.01
#  fit = MultiLociPOE(Y,gmm,gcc,X,loci,hap,f,ppi)

## ----eval=FALSE---------------------------------------------------------------
#  names(fit)

## ----eval=FALSE---------------------------------------------------------------
#  fit$new

## ----eval=FALSE---------------------------------------------------------------
#  fit$cov.new

## ----eval=FALSE---------------------------------------------------------------
#  fit$log
#  fit$cov.log

