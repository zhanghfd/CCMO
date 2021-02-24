#' @title Association analysis with multiple SNPs
#' @description Multiple SNPs can be analyzed utilizing multiple CPU cores (requires the R package \link{doMC}), which calls \link{singleSNP}. Estimation results include the parameter estimates and their estimated standard errors, the p-values for significance tests (need to call \link{OmnibusTest}), some model information, and omnibus test if required (the argument \code{test} should be specified).
#' @param Y a \code{n}-vector of disease statuses for n offspring (1 for case and 0 for control).
#' @param Gc a \code{n} x \code{m} matrix of genotypes for offspring (n: number of offspring; m: number of SNPs). The possible values should be 0, 1, 2, or 1/2 (see \code{normalized.genotype} for details).
#' @param Gm a \code{n} x \code{m} matrix of genotypes for mothers (n: number of mothers; m: number of SNPs). The possible values should be 0, 1, 2, or 1/2 (see \code{normalized.genotype} for details).
#' @param Xo a \code{n} x \code{n_o} matrix of maternal covariates for main effects (n: number of mothers; \code{n_o}: the numbers of genotypes; default value: NULL)
#' @param Xc a \code{n} x \code{n_c} matrix of maternal covariates for intercation effects with offspring genotypes (\code{n}: number of mothers; \code{n_c}: the numbers of genotypes; default value: NULL)
#' @param Xm a \code{n} x \code{n_m} matrix of maternal covariates for intercation effects with maternal genotypes (\code{n}: number of mothers; \code{n_m}: the numbers of genotypes; default value: NULL)
#' @param X.Gm a \code{n} x \code{n_g} matrix of maternal covariates potentially associated with maternal genotypes (\code{n}: number of mothers; \code{n_g}: the numbers of genotypes; default value: NULL)
#' @param G.main a vector containing the main SNP effects taking values 'Gc' and/or 'Gm' (default value: c('Gc','Gm'))
#' @param G.int an indicator for the presence of Gc x Gm interaction (TRUE for yes and FALSE for no; default value: FALSE)
#' @param mode mode of inheritance ('rec' for recessive, 'add' for additive, 'dom' for dominant; default value: 'add')
#' @param prev specified disease prevalence (default value: 0.01)
#' @param ind a logical variable indicating whether Gm and X.Gm are independent (TRUE for independence and FALSE for dependence, default value: FALSE)
#' @param HWE a logical variable indicating whether the HWE assumption is incorporated (TRUE for incorprated and FALSE otherwise; default value: TRUE)
#' @param normalized.genotype a logical variable indicating whether the genotypes are normalized (TRUE for normalized and FALSE otherwise, default value: FALSE). If FALSE, Gc and Gm should take values 0 (genotype \code{AA}), 1 (genotype \code{AB}), 2 (genotype \code{BB}); otherwise, their values depend on the mode of inheritance ('add': 0 for genotype \code{AA}, 0.5 for genotype \code{AB}, 1 for genotype \code{BB}; 'rec': 0 for genotype \code{AA} or AB, 1 for genotype BB; 'dom': 0 for genotype \code{AA}, 1 for genotype AB or BB).
#' @param test a vector of predictor indices in Omnibus test (default value: NULL; if not NULL, call \link{OmnibusTest})
#' @param cores the number of CPU cores used for parallele execution. If not specified, the system will determine its value.
#' @return a list of length \code{m}, each element is a list with the following elements
#' @return \item{\code{new}}{estimation and significance test results for the new method}
#' @return \item{\code{log}}{estimation and significance test results for the standard logistic regression method}
#' @return \item{\code{cov.new}}{covariance matrix of the estimated parameters by the new method}
#' @return \item{\code{cov.log}}{covariance matrix of the estimated parameters by the standard logistic regression method)}
#' @return \item{\code{penetrance}}{logistic regression model for the penetrance function}
#' @return \item{\code{daLOG}}{daLOG model relating maternal genotype and maternal covariates}
#' @return \item{\code{Omnibus}}{Omnibus test results (test statistic, degrees of freedom, p-value) if \code{test} is not NULL}
#' @seealso \code{\link{singleSNP}} \code{\link{OmnibusTest}}
#' @examples
#' \dontrun{
#' library(foreach)
#' data(SampleData)
#' Y = SampleData[,1]
#' M = 10 # the number of SNPs
#' Gc = SampleData[,1+1:M]
#' Gm = SampleData[,1+M+(1:M)]
#' X = SampleData[,-(1:(1+2*M))]
#' res = multipleSNP(Y,Gc,Gm,Xo=X,Xc=X,Xm=X,X.Gm=X)
#' }
#' @export
multipleSNP <- function(Y,Gc,Gm,Xo=NULL,Xc=NULL,Xm=NULL,X.Gm=NULL,G.main=c('Gc','Gm'),G.int=FALSE,mode='add',prev=0.01,ind=FALSE,HWE=TRUE,normalized.genotype=FALSE,test=NULL,cores=NULL){

  M = ncol(Gc)
  registerDoMC(cores=cores)
  i = 0
  foreach::foreach(i = 1:M) %dopar% {
    d = singleSNP(Y,Gc[,i],Gm[,i],Xo,Xc,Xm,X.Gm,G.main,G.int,mode,prev,ind,HWE,normalized.genotype)
    if(!is.null(test)) d = OmnibusTest(d,test)
    d
  }
}

#' @title A sample dataset
#' @name SampleData
#' @description This sample dataset include the genotypes for 2000 mother-offspring pairs (1000 with case offspring and 1000 with control offspring), together with two covariate values for each of the 2000 mothers.
#' @import aod
#' @import foreach
#' @import doMC
#' @import nloptr
#' @import numDeriv
#' @importFrom genetics genotype allele.names
#' @importFrom haplo.stats haplo.em geno1to2
#' @importFrom stats binomial glm pchisq complete.cases
#' @importFrom MASS ginv
#' @examples
#' \dontrun{
#' data(SampleData)
#' }
NULL

#' @title A sample dataset for POE
#' @name POESampleData
#' @description This sample dataset for POE based five loci.
#' @examples
#' \dontrun{
#' data(POESampleData)
#' }
NULL



#' @title Association analysis with a single SNP
#' @description Association analysis is conducted between a single SNP and the disease of interest. Estimation results include the parameter estimates and their estimated standard errors, the p-values for significance tests, covariance matrices, and some model information.
#' @param Y a \code{n}-vector of disease statuses for n offspring (1 for case and 0 for control)
#' @param Gc a \code{n}-vector of genotypes for offspring, taking values 0, 1, and 2 (n: number of offspring).  The possible values should be 0, 1, 2, or 1/2 (see \code{normalized.genotype} for details).
#' @param Gm a \code{n}-vector of genotypes for mothers, taking values 0, 1, and 2 (n: number of mothers).  The possible values should be 0, 1, 2, or 1/2 (see \code{normalized.genotype} for details).
#' @param Xo a \code{n} x \code{n_o} matrix of maternal covariates for main effects (n: number of mothers; \code{n_o}: the numbers of genotypes; default value: NULL)
#' @param Xc a \code{n} x \code{n_c} matrix of maternal covariates for intercation effects with offspring genotypes (\code{n}: number of mothers; \code{n_c}: the numbers of genotypes; default value: NULL)
#' @param Xm a \code{n} x \code{n_m} matrix of maternal covariates for intercation effects with maternal genotypes (\code{n}: number of mothers; \code{n_m}: the numbers of genotypes; default value: NULL)
#' @param X.Gm a \code{n} x \code{n_g} matrix of maternal covariates potentially associated with maternal genotypes (\code{n}: number of mothers; \code{n_g}: the numbers of genotypes; default value: NULL)
#' @param G.main a vector containing the main SNP effects taking values 'Gc' and/or 'Gm' (default value: c('Gc','Gm'))
#' @param G.int an indicator for the presence of Gc x Gm interaction (TRUE for yes and FALSE for no; default value: FALSE)
#' @param mode mode of inheritance ('rec' for recessive, 'add' for additive, 'dom' for dominant; default value: 'add')
#' @param prev specified disease prevalence (default value: 0.01)
#' @param ind a logical variable indicating whether Gm and X.Gm are independent (TRUE for independence and FALSE for dependence, default value: FALSE)
#' @param HWE a logical variable indicating whether the HWE assumption is incorporated (TRUE for incorprated and FALSE otherwise; default value: TRUE)
#' @param normalized.genotype a logical variable indicating whether the genotypes are normalized (TRUE for normalized and FALSE otherwise, default value: FALSE). If FALSE, Gc and Gm should take values 0 (genotype \code{AA}), 1 (genotype \code{AB}), 2 (genotype \code{BB}); otherwise, their values depend on the mode of inheritance ('add': 0 for genotype \code{AA}, 0.5 for genotype \code{AB}, 1 for genotype \code{BB}; 'rec': 0 for genotype \code{AA} or AB, 1 for genotype BB; 'dom': 0 for genotype \code{AA}, 1 for genotype AB or BB).
#' @return a list with the following elements
#' @return \item{\code{new}}{estimation and significance test results for the new method}
#' @return \item{\code{log}}{estimation and significance test results for the standard logistic regression method}
#' @return \item{\code{cov.new}}{covariance matrix of the estimated parameters by the new method}
#' @return \item{\code{cov.log}}{covariance matrix of the estimated parameters by the standard logistic regression method)}
#' @return \item{\code{penetrance}}{logistic regression model for the penetrance function}
#' @return \item{\code{daLOG}}{daLOG model relating maternal genotype and maternal covariates}
#' @seealso \code{\link{multipleSNP}} \code{\link{OmnibusTest}}
#' @examples
#' \dontrun{
#' data(SampleData)
#' Y = SampleData[,1]
#' M = 10 # the number of SNPs
#' Gc = SampleData[,2]
#' Gm = SampleData[,12]
#' X = SampleData[,-(1:(1+2*M))]
#' res = singleSNP(Y,Gc,Gm,Xo=X,Xc=X,Xm=X,X.Gm=X)
#' }
#' @export
singleSNP <- function(Y,Gc,Gm,Xo=NULL,Xc=NULL,Xm=NULL,X.Gm=NULL,G.main=c('Gc','Gm'),G.int=FALSE,mode='add',prev=0.01,ind=FALSE,HWE=TRUE,normalized.genotype=FALSE){

  n.o = ifelse(is.null(Xo),0,max(1,ncol(Xo)))
  n.c = ifelse(is.null(Xc),0,max(1,ncol(Xc)))
  n.m = ifelse(is.null(Xm),0,max(1,ncol(Xm)))
  n.Gm = ifelse(is.null(X.Gm),0,max(1,ncol(X.Gm)))

  if(!normalized.genotype) Gc[(Gc==0&Gm==2)|(Gc==2&Gm==0)] = NA
  if(normalized.genotype & mode=='add') Gc[(Gc==0&Gm==1)|(Gc==1&Gm==0)] = NA
  ct = cumsum(c(3,n.o,n.c,n.m,n.Gm))
  dat = cbind(Y,Gc,Gm,Xo,Xc,Xm,X.Gm)
  dat = dat[complete.cases(dat),]
  Y = dat[,1]
  Gc = dat[,2]
  Gm = dat[,3]
  if(!is.null(Xo))  Xo = matrix(dat[,4:ct[2]],ncol=n.o)
  if(!is.null(Xc))  Xc = matrix(dat[,(ct[2]+1):ct[3]],ncol=n.c)
  if(!is.null(Xm))  Xm = matrix(dat[,(ct[3]+1):ct[4]],ncol=n.m)
  if(!is.null(X.Gm))  X.Gm = matrix(dat[,(ct[4]+1):ct[5]],ncol=n.Gm)


  I <- numeric(2)
  I[1] <- ifelse('Gc' %in% G.main,1,0)
  I[2] <- ifelse('Gm' %in% G.main,1,0)

  if(!normalized.genotype){
    if(mode=='add'){Gc <- Gc/2; Gm <- Gm/2}
    if(mode=='dom'){Gc <- ifelse(Gc==0,0,1); Gm <- ifelse(Gm==0,0,1)}
    if(mode=='rec'){Gc <- ifelse(Gc==2,1,0); Gm <- ifelse(Gm==2,1,0)}
  }

  if(is.null(X.Gm)){name.X.Gm <- NULL}else{
    # daLOG model
    X.Gm.name = colnames(X.Gm)
    if(!is.null(colnames(X.Gm))){
      name.X.Gm <- paste0('Gm~',X.Gm.name)
      daLOG <- paste0('Gm ~ ',X.Gm.name[1])
      if(length(X.Gm.name)>1) for(k in 2:length(X.Gm.name)) daLOG <- paste(daLOG,X.Gm.name[k],sep=' + ')
    }else{
      n.X.Gm = ncol(X.Gm)
      name.X.Gm <- paste0('Gm~Z',seq(n.X.Gm))
      daLOG <- paste0('Gm ~ ','Z1')
      if(n.X.Gm > 1) for(k in 2:n.X.Gm) daLOG <- paste(daLOG,paste0('Z',k),sep=' + ')
    }
  }
  if(ind) daLOG = NULL else daLOG = daLOG

  if(is.null(Xo)){name.Xo <- NULL}else{
    if(!is.null(colnames(Xo))){
      name.Xo <- colnames(Xo)
    }else{
      name.Xo <- paste0('Xo',seq_len(ncol(Xo)))
    }
  }
  if(is.null(Xc)){name.Xc <- NULL}else{
    if(!is.null(colnames(Xc))){
      name.Xc <- paste0('Gc:',colnames(Xc))
    }else{
      name.Xc <- paste0('Gc:Xc',seq_len(ncol(Xc)))
    }
  }
  if(is.null(Xm)){name.Xm <- NULL}else{
    if(!is.null(colnames(Xm))){
      name.Xm <- paste0('Gm:',colnames(Xm))
    }else{
      name.Xm <- paste0('Gm:Xm',seq_len(ncol(Xm)))
    }
  }
  name.X <- c(name.Xo,name.Xc,name.Xm)

  # penetrance model
  penetrance <- 'Y ~ '
  tmp <- c('Gc','Gm')[which(I==1)]
  if(length(tmp)>=1) tt = tmp[1]
  if(length(tmp)==2) tt = paste(tt,tmp[2],sep=' + ')
  penetrance = paste0(penetrance,tt)
  if(length(name.X)>=1) tt = name.X[1]
  if(length(name.X)>1) for(k in 2:length(name.X))  tt = paste(tt,name.X[k],sep=' + ')
  penetrance = paste(penetrance,tt,sep=' + ')
  if(G.int) penetrance = paste0(penetrance,'Gc:Gm',sep=' + ')

  if(!ind){
    # DEP
    res=solution.dep(Y,Gc,Gm,Xo,Xc,Xm,X.Gm,G.main,G.int,mode,prev,HWE)
    est.new = res$est #estimated values of parameters
    cov.new = res$cov #covariance of parameters
    est.log <- res$est.log
    cov.log <- res$cov.log
    names(est.log) <- colnames(cov.log) <- rownames(cov.log) <- c('Intercept',c('Gc','Gm')[which(I==1)],name.X,'Gc:Gm'[G.int])
    if(HWE){
      names(est.new) <- colnames(cov.new) <- rownames(cov.new) <- c('theta','Intercept',c('Gc','Gm')[which(I==1)],name.X,'Gc:Gm'[G.int],name.X.Gm)
    }else{
      names(est.new) <- colnames(cov.new) <- rownames(cov.new) <- c('theta','F','Intercept',c('Gc','Gm')[which(I==1)],name.X,'Gc:Gm'[G.int],name.X.Gm)
    }
  }else{
    # IND
    res=solution.ind(Y,Gc,Gm,Xo,Xc,Xm,G.main,G.int,mode,prev,HWE)
    est.new = res$est #estimated values of parameters
    cov.new = res$cov #covariance of parameters
    est.log <- res$est.log
    cov.log <- res$cov.log
    names(est.log) <- colnames(cov.log) <- rownames(cov.log) <- c('Intercept',c('Gc','Gm')[which(I==1)],name.X,'Gc:Gm'[G.int])
    if(HWE){
      names(est.new) <- colnames(cov.new) <- rownames(cov.new) <- c('theta','Intercept',c('Gc','Gm')[which(I==1)],name.X,'Gc:Gm'[G.int])
    }else{
      names(est.new) <- colnames(cov.new) <- rownames(cov.new) <- c('theta','F','Intercept',c('Gc','Gm')[which(I==1)],name.X,'Gc:Gm'[G.int])
    }
  }

  se.new <- sqrt(diag(cov.new))
  z.new <- est.new/se.new
  pval.new <- 1-pchisq(z.new^2,1)
  z.new[1] <- NA
  pval.new[1] <- NA

  se.log <- sqrt(diag(cov.log))
  z.log <- est.log/se.log
  z.log <- est.log/sqrt(diag(cov.log))
  pval.log <- 1-pchisq(z.log^2,1)

  list(new=data.frame(est=est.new,se=se.new,`z score`=z.new,`p value`=pval.new),
       log=data.frame(est=est.log,se=se.log,`z score`=z.log,`p.value`=pval.log),
       cov.new=cov.new,cov.log=cov.log,HWE=HWE,IND=ind,penetrance=penetrance,daLOG=daLOG)

}

#' @title Omnibus test for multiple effects
#' @description The interested effects include the genetic effects and/or environment effects and/or gene-environment interaction effects and/or gene-gene interaction effects. Association analysis is conducted between a single SNP and the disease of interest.
#' @param d an object returned by the R function \link{singleSNP}
#' @param test a vector of predictor indices in Omnibus test (default value: NULL)
#' @return a list with an additional element for Omnibus test p-values
#' @return \item{\code{Omnibus}}{Omnibus test results (test statistic, degrees of freedom, p-value) if \code{test} is not NULL}
#' @seealso \code{\link{singleSNP}} \code{\link{multipleSNP}}
#' @examples
#' \dontrun{
#' data(SampleData)
#' Y = SampleData[,1]
#' M = 10 # the number of SNPs
#' Gc = SampleData[,2]
#' Gm = SampleData[,12]
#' X = SampleData[,-(1:(1+2*M))]
#' d = singleSNP(Y,Gc,Gm,Xo=X,Xc=X,Xm=X,X.Gm=X)
#' test = c(3,4,6,7) # main SNP effects and SNP-covariate interaction effects
#' d = OmnibusTest(d,test)
#' }
#' @export
OmnibusTest <- function(d,test){

  p.value <- numeric(2)
  est = d$new$est[test];
  cov = MASS::ginv(d$cov.new[test,test]);

  df = length(test)

  stat.new = as.numeric((est%*%cov)%*%est)
  p.new = 1-pchisq(stat.new,df)
  res.new = c(stat.new,df,p.new)

  if(d$HWE){
    rr = wald.test(b = d$log$est, Sigma = d$cov.log, Terms = test-1)
  }else{
    rr = wald.test(b = d$log$est, Sigma = d$cov.log, Terms = test-2)
  }

  res.log = rr$result$chi2
  Omnibus <- rbind(res.new,res.log)
  rownames(Omnibus) = c('new','log')
  colnames(Omnibus) = c('stat','df','p.value')

  d$Omnibus = Omnibus
  d
}


#' @title POEs analysis using multi-locus genotype data from case-control mother-offspring pairs
#' @description This function implements the statistical method M-HAP for detecting POEs using genotype data from case-control mother-offspring pairs by adjusting for covariates. M-HAP utilizes available information such as Medelian inheritance law, Hardy-Weinberg equilibrium, and conditional independence between offspring genotypes and maternal covariates given maternal genotypes.
#' @param Y a \code{n}-vector of disease statuses for \code{n} offspring (1 for case and 0 for control).
#' @param gcc a \code{n} x \code{m} matrix of genotypes for offspring (\code{n}: number of offspring; \code{m}: number of SNPs). The possible values should be 0, 1, 2.
#' @param gmm a \code{n} x \code{m} matrix of genotypes for mothers (\code{n}: number of mothers; \code{m}: number of SNPs). The possible values should be 0, 1, 2.
#' @param X a \code{n}-vector of maternal covariates for main effects (\code{n}: number of mothers).
#' @param loci an indicator for the test locus. The possible values could be 1,...,\code{m} (\code{m}: number of SNPs).
#' @param hap a \code{l} x \code{m} matrix of possible haplotypes in the population of interest (\code{l}: number of possible haplotypes; \code{m}: number of SNPs).
#' @param f specified disease prevalence.
#' @param ppi a \code{l}-vector of the corresponding haplotype frequencies (\code{m}: number of possible haplotypes).
#' @return a list with the following elements
#' @return \item{\code{new}}{estimation and significance test results for the new method M-HAP}
#' @return \item{\code{log}}{estimation and significance test results for the standard logistic regression method}
#' @return \item{\code{cov.new}}{covariance matrix of the estimated parameters by the new method M-HAP}
#' @return \item{\code{cov.log}}{covariance matrix of the estimated parameters by the standard logistic regression method)}
#' @examples
#' \dontrun{
#' data(POESampleData)
#' Y = POESampleData[,1]
#' gmm = POESampleData[,2:6]
#' gcc = POESampleData[,7:11]
#' X = POESampleData[,12]
#' loci = 1
#' f = 0.01
#' data = MultiLociPOE.input(gmm,gcc,0)
#' gmm = data$gmm
#' gcc = data$gcc
#' hap = data$hap
#' ppi = data$ppi
#' fit = MultiLociPOE(Y,gmm,gcc,X,loci,hap,f,ppi)
#' }
#' @export
MultiLociPOE = function(Y,gmm,gcc,X,loci,hap,f,ppi){

  res = imprinting(Y,gmm,gcc,X,loci,hap,f,ppi)

  est.new = res$est #estimated values of parameters
  cov.new = res$cov #covariance matrix of estimated parameters
  names(est.new) = colnames(cov.new) = rownames(cov.new) = c('theta','Intercept','Gm','Gc','Gmc-Gpc','X','eta')
  se.new = sqrt(diag(cov.new))
  z.new = est.new/se.new
  pval.new = 1-pchisq(z.new^2,1)
  z.new[1] = NA
  pval.new[1] = NA
  new = data.frame(est=est.new,se=se.new,`z score`=z.new,`p value`=pval.new)

  est.log = res$est.log
  cov.log = res$cov.log
  names(est.log) = colnames(cov.log) = rownames(cov.log) = c('Intercept','Gm','Gc','Gmc-Gpc','X')
  se.log = sqrt(diag(cov.log))
  z.log = est.log/se.log
  z.log = est.log/sqrt(diag(cov.log))
  pval.log = 1-pchisq(z.log^2,1)
  log=data.frame(est=est.log,se=se.log,`z score`=z.log,`p.value`=pval.log)

  return(list(new=new,log=log,cov.new=cov.new,cov.log=cov.log))

}

#' @title Inputs for the function MultiLociPOE
#' @description This function can be used to obtain required inputs for the function MultiLociPOE: all possible haplotypes in the population of interest and the corresponding haplotype frequencies.
#' @param gcc a \code{n} x \code{m} matrix of genotypes for offspring (n: number of offspring; m: number of SNPs). The possible values should be 0, 1, 2.
#' @param gmm a \code{n} x \code{m} matrix of genotypes for mothers (n: number of mothers; m: number of SNPs). The possible values should be 0, 1, 2.
#' @param type an indicator for the genotype data format. (0 for genotype coded as 0, 1, 2; 1 for normalized genotype coded as 0, 1/2, 1; 2 for genotype coded like "A/A", "A/B", "B/B")
#' @return a list with the following elements
#' @return \item{\code{gmm}}{matrix of genotypes for mothers (n: number of mothers; m: number of SNPs)}
#' @return \item{\code{gcc}}{matrix of genotypes for offspring (n: number of offspring; m: number of SNPs)}
#' @return \item{\code{hap}}{a \code{l} x \code{m} matrix of possible haplotypes in the cohort (l: number of possible haplotypes; m: number of SNPs)}
#' @return \item{\code{ppi}}{a \code{l}-vector of the corresponding haplotype frequencies (l: number of possible haplotypes)}
#' @export

MultiLociPOE.input = function(gmm,gcc,type){

  Gcc = Gmm = matrix(NA,nrow(gmm),ncol(gmm))

  if(type==0){

    Gmm = gmm
    Gcc = gcc
  }
  else if(type==1){

    Gmm = gmm*2
    Gcc = gcc*2
  }
  else {

    for(i in 1:ncol(gmm)){
      a = allele.names(genotype(gmm[,i]))
      a0 = paste0(a[2], "/", a[2])
      a1 = paste0(a[1], "/", a[2])
      a2 = paste0(a[1], "/", a[1])
      for (j in 1:nrow(gmm) ) {
        if(genotype(gmm[,i])[j]==a0) Gmm[j,i] = 0
        else if(genotype(gmm[,i])[j]==a1) Gmm[j,i] = 1
        else Gmm[j,i] = 2
      }
    }
  }

  K =ncol(Gmm)
  gmm1 = geno1to2(Gmm)
  gcc1 = geno1to2(Gcc)
  haplo_raw = haplo.em(rbind(gmm1, gcc1))
  ppi_raw = haplo_raw$hap.prob
  n_hap_raw = length(ppi_raw)
  hap_raw = matrix(as.numeric(as.matrix(haplo_raw$haplotype)), n_hap_raw, K, byrow = FALSE) - 1

  return(list(gmm=Gmm,gcc=Gcc,hap=hap_raw,ppi=ppi_raw))

}



# Gm and X are dependent
solution.dep=function(Y,Gc,Gm,Xo,Xc,Xm,X.Gm,G.main,G.int,mode,prev,HWE){

  n <- length(Y)
  Z = design.matrix(n,Gc,Gm,Xc,Xm,Xo,G.main,G.int)
  f = prev

  fit <- glm(Y ~ 0 + Z,family = binomial)
  res <- summary(fit)$coef
  est.log <- as.vector(res[,1])
  cov.log = vcov(fit)

  theta = 1-sqrt(mean(Gm==0))
  n.eta = max(1,ncol(X.Gm)) # number of covariates possibly associated with Gm
  n.beta = ncol(Z) # number of all regression parameters in the penetrance model
  theta0 = log(f/(1-f))
  beta = rep(theta0,n.beta)
  beta[-1] = est.log[-1]
  Theta = c(theta,0,beta,rep(0,n.eta))

  P=list()
  if(mode=='add'){
    for(j in 1:3){
      P[[j]]=list()
      for(k in 1:3){
        P[[j]][[k]]=design.matrix(n,(j-1)/2,(k-1)/2,Xc,Xm,Xo,G.main,G.int)
      }
    }
  } else {
    for(j in 1:2){
      P[[j]]=list()
      for(k in 1:2){
        P[[j]][[k]]=design.matrix(n,j-1,k-1,Xc,Xm,Xo,G.main,G.int)
      }
    }
  }

  if(mode=='add'){
    m11 = sum(Gc==0&Gm==0)
    m12 = sum(Gc==0&Gm==0.5)
    m21 = sum(Gc==0.5&Gm==0)
    m22 = sum(Gc==0.5&Gm==0.5)
    m23 = sum(Gc==0.5&Gm==1)
    m32 = sum(Gc==1&Gm==0.5)
    m33 = sum(Gc==1&Gm==1)
  } else{
    m00 = sum(Gc==0&Gm==0)
    m01 = sum(Gc==0&Gm==1)
    m10 = sum(Gc==1&Gm==0)
    m11 = sum(Gc==1&Gm==1)
  }

  n.beta = ncol(Z)

  n = length(Y)
  ns = c(sum(Y==1),sum(Y==0))
  II = list()
  II[[1]] = which(Y==0)
  II[[2]] = which(Y==1)

  lambda = ns[1]/n/f - ns[2]/n/(1-f)
  p = length(Theta) # number of all parameters

  upper = lower = Theta
  upper[1] = 0.99
  lower[1] = 0.01
  if(HWE){lower[2] = upper[2] = 0 }else{
    lower[2] = -0.1
    upper[2] = 0.5
  }
  upper[-(1:2)] = Theta[-(1:2)] + 2
  lower[-(1:2)] = Theta[-(1:2)] - 2


  fn.rec = function(Theta){
    theta=Theta[1];
    F = ifelse(HWE,0,Theta[2])
    beta = Theta[3:(2+n.beta)];
    eta = c(log(xi.rec(theta,F)),Theta[-(1:(2+n.beta))]);

    tr = tr.rec(theta,F);
    h.fn = H.rec(n.beta,Theta,X.Gm,P);

    res = m00*log(tr[1,1])+m10*log(tr[2,1])+m01*log(tr[1,2])+m11*log(tr[2,2]);
    res = res + sum(log(dalog(Gm,X.Gm,eta)));
    res = res + sum(log(penetrance1(Z,Y,beta)));
    res = res - sum(log(1+lambda*(h.fn-f)));
    return(-res);
  }
  fn.dom = function(Theta){
    theta=Theta[1];
    F = ifelse(HWE,0,Theta[2])
    beta = Theta[3:(2+n.beta)];
    eta = c(log(xi.dom(theta,F)),Theta[-(1:(2+n.beta))]);

    tr = tr.dom(theta,F);
    h.fn = H.dom(n.beta,Theta,X.Gm,P);

    res = m00*log(tr[1,1])+m10*log(tr[2,1])+m01*log(tr[1,2])+m11*log(tr[2,2]);
    res = res + sum(log(dalog(Gm,X.Gm,eta)));
    res = res + sum(log(penetrance1(Z,Y,beta)));
    res = res - sum(log(1+lambda*(h.fn-f)));
    return(-res);
  }
  fn.add = function(Theta){
    theta=Theta[1];
    F = ifelse(HWE,0,Theta[2])
    beta = Theta[3:(2+n.beta)];
    eta = Theta[-(1:(2+n.beta))];

    tr = tr.add(theta);
    h.fn = H.add(n.beta,Theta,X.Gm,P);

    res = m11*log(tr[1,1])+m12*log(tr[1,2])+m21*log(tr[2,1])+m22*log(tr[2,2])+m23*log(tr[2,3])+m32*log(tr[3,2])+m33*log(tr[3,3]);
    res = res + sum(log(dalog.add(Gm,X.Gm,theta,F,eta)));
    res = res + sum(log(penetrance1(Z,Y,beta)));
    res = res - sum(log(1+lambda*(h.fn-f)));
    return(-res);
  }

  gr.rec = function(Theta){
    theta=Theta[1];
    F = ifelse(HWE,0,Theta[2])
    beta = Theta[3:(2+n.beta)];

    eta = c(log(xi.rec(theta,F)),Theta[-(1:(2+n.beta))]);

    tr.fn = tr.rec(theta,F)
    tr.gr = tr.gr.rec(theta,F)

    logq.theta = logq.F = rep(0,n);

    logq.gr = tr.gr[[1]]/tr.fn

    # derivative of log(q) w.r.t. theta
    logq.theta[Gc==0&Gm==0] = logq.gr[1,1];
    logq.theta[Gc==1&Gm==0] = logq.gr[2,1];
    logq.theta[Gc==0&Gm==1] = logq.gr[1,2];
    logq.theta[Gc==1&Gm==1] = logq.gr[2,2];

    logq.gr = tr.gr[[2]]/tr.fn

    # derivative of log(q) w.r.t. F
    logq.F[Gc==0&Gm==0] = logq.gr[1,1];
    logq.F[Gc==1&Gm==0] = logq.gr[2,1];
    logq.F[Gc==0&Gm==1] = logq.gr[1,2];
    logq.F[Gc==1&Gm==1] = logq.gr[2,2];

    eta0.gr = logxi.gr.rec(theta,F);
    logr.gr = dalog.gr(Gm,X.Gm,eta)/dalog(Gm,X.Gm,eta);

    # derivative of log(r) w.r.t. theta
    logr.theta = logr.gr*eta0.gr[1];
    # derivative of log(r) w.r.t. F
    logr.F     = logr.gr*eta0.gr[2]
    # derivative of log(r) w.r.t. eta
    logr.eta   = sweep(X.Gm,1,logr.gr,'*');
    logp.gr = penetrance1.gr(Z,Y,beta)/penetrance1(Z,Y,beta);
    # derivative of log(p) w.r.t. beta
    logp.beta = sweep(Z,1,logp.gr,'*');

    h.fn = H.rec(n.beta,Theta,X.Gm,P);
    h.gr = H.gr.rec(n.beta,Theta,X.Gm,P);

    # derivative of H w.r.t. all parameters
    others = sweep(lambda*h.gr,1,1+lambda*(h.fn-f),'/');
    all = cbind(logq.theta+logr.theta,logq.F+logr.F,logp.beta,logr.eta) - others;
    -colSums(all)
  }

  gr.dom = function(Theta){
    theta=Theta[1];
    F = ifelse(HWE,0,Theta[2])
    beta = Theta[3:(2+n.beta)];

    eta = c(log(xi.dom(theta,F)),Theta[-(1:(2+n.beta))]);

    tr.fn = tr.dom(theta,F)
    tr.gr = tr.gr.dom(theta,F)

    logq.theta = logq.F = rep(0,n);

    logq.gr = tr.gr[[1]]/tr.fn

    # derivative of log(q) w.r.t. theta
    logq.theta[Gc==0&Gm==0] = logq.gr[1,1];
    logq.theta[Gc==1&Gm==0] = logq.gr[2,1];
    logq.theta[Gc==0&Gm==1] = logq.gr[1,2];
    logq.theta[Gc==1&Gm==1] = logq.gr[2,2];

    logq.gr = tr.gr[[2]]/tr.fn

    # derivative of log(q) w.r.t. F
    logq.F[Gc==0&Gm==0] = logq.gr[1,1];
    logq.F[Gc==1&Gm==0] = logq.gr[2,1];
    logq.F[Gc==0&Gm==1] = logq.gr[1,2];
    logq.F[Gc==1&Gm==1] = logq.gr[2,2];

    eta0.gr = logxi.gr.dom(theta,F);
    logr.gr = dalog.gr(Gm,X.Gm,eta)/dalog(Gm,X.Gm,eta);

    # derivative of log(r) w.r.t. theta
    logr.theta = logr.gr*eta0.gr[1];
    # derivative of log(r) w.r.t. F
    logr.F     = logr.gr*eta0.gr[2]
    # derivative of log(r) w.r.t. eta
    logr.eta   = sweep(X.Gm,1,logr.gr,'*');
    logp.gr = penetrance1.gr(Z,Y,beta)/penetrance1(Z,Y,beta);
    # derivative of log(p) w.r.t. beta
    logp.beta = sweep(Z,1,logp.gr,'*');

    h.fn = H.dom(n.beta,Theta,X.Gm,P);
    h.gr = H.gr.dom(n.beta,Theta,X.Gm,P);

    # derivative of H w.r.t. all parameters
    others = sweep(lambda*h.gr,1,1+lambda*(h.fn-f),'/');
    all = cbind(logq.theta+logr.theta,logq.F+logr.F,logp.beta,logr.eta) - others;
    -colSums(all)
  }

  gr.add = function(Theta){
    theta=Theta[1];
    F = ifelse(HWE,0,Theta[2])
    beta = Theta[3:(2+n.beta)];
    eta = Theta[-(1:(2+n.beta))];

    logq.gr = tr.gr.add(theta)[[1]]/tr.add(theta)

    logq.theta = rep(NA,n);
    logq.theta[Gc==0&Gm==0] = logq.gr[1,1];
    logq.theta[Gc==0.5&Gm==0] = logq.gr[2,1];
    logq.theta[Gc==0&Gm==0.5] = logq.gr[1,2];
    logq.theta[Gc==0.5&Gm==0.5] = logq.gr[2,2];
    logq.theta[Gc==1&Gm==0.5] = logq.gr[3,2];
    logq.theta[Gc==0.5&Gm==1] = logq.gr[2,3];
    logq.theta[Gc==1&Gm==1] = logq.gr[3,3];


    logr.gr = sweep(dalog.gr.add(Gm,X.Gm,theta,F,eta),1,dalog.add(Gm,X.Gm,theta,F,eta),'/');

    logr.theta = logr.gr[,1];
    logr.F = logr.gr[,2]
    logr.eta = logr.gr[,-(1:2)];

    logp.gr = penetrance1.gr(Z,Y,beta)/penetrance1(Z,Y,beta);
    logp.beta = sweep(Z,1,logp.gr,'*');

    h.fn = H.add(n.beta,Theta,X.Gm,P);
    h.gr = H.gr.add(n.beta,Theta,X.Gm,P);
    others = sweep(lambda*h.gr,1,1+lambda*(h.fn-f),'/');

    all = cbind(logq.theta+logr.theta,logr.F,logp.beta,logr.eta) - others;
    -colSums(all)
  }

  if(mode=='rec'){
    fit = nloptr(x0=Theta,eval_f=fn.rec,eval_grad_f=gr.rec,lb=lower,ub=upper,opts=list("algorithm"="NLOPT_LD_LBFGS","xtol_rel"=1.0e-4) );
    hess = hessian(func=fn.rec,x=fit$solution);
  }
  if(mode=='dom'){
    fit = nloptr(x0=Theta,eval_f=fn.dom,eval_grad_f=gr.dom,lb=lower,ub=upper,opts=list("algorithm"="NLOPT_LD_LBFGS","xtol_rel"=1.0e-4) );
    hess = hessian(func=fn.dom,x=fit$solution);
  }
  if(mode=='add'){
    fit = nloptr(x0=Theta,eval_f=fn.add,eval_grad_f=gr.add,lb=lower,ub=upper,opts=list("algorithm"="NLOPT_LD_LBFGS","xtol_rel"=1.0e-4) );
    hess = hessian(func=fn.add,x=fit$solution);
  }
  est = fit$solution

  score.cov.rec = function(Theta){
    theta=Theta[1];
    F = Theta[2]
    beta = Theta[3:(2+n.beta)];

    eta = c(log(xi.rec(theta,F)),Theta[-(1:(2+n.beta))]);

    tr.fn = tr.rec(theta,F)
    tr.gr = tr.gr.rec(theta,F)

    logq.theta = logq.F = rep(0,n);

    logq.gr = tr.gr[[1]]/tr.fn

    # derivative of log(q) w.r.t. theta
    logq.theta[Gc==0&Gm==0] = logq.gr[1,1];
    logq.theta[Gc==1&Gm==0] = logq.gr[2,1];
    logq.theta[Gc==0&Gm==1] = logq.gr[1,2];
    logq.theta[Gc==1&Gm==1] = logq.gr[2,2];

    logq.gr = tr.gr[[2]]/tr.fn

    # derivative of log(q) w.r.t. F
    logq.F[Gc==0&Gm==0] = logq.gr[1,1];
    logq.F[Gc==1&Gm==0] = logq.gr[2,1];
    logq.F[Gc==0&Gm==1] = logq.gr[1,2];
    logq.F[Gc==1&Gm==1] = logq.gr[2,2];

    eta0.gr = logxi.gr.rec(theta,F);
    logr.gr = dalog.gr(Gm,X.Gm,eta)/dalog(Gm,X.Gm,eta);

    # derivative of log(r) w.r.t. theta
    logr.theta = logr.gr*eta0.gr[1];
    # derivative of log(r) w.r.t. F
    logr.F     = logr.gr*eta0.gr[2]
    # derivative of log(r) w.r.t. eta
    logr.eta   = sweep(X.Gm,1,logr.gr,'*');
    logp.gr = penetrance1.gr(Z,Y,beta)/penetrance1(Z,Y,beta);
    # derivative of log(p) w.r.t. beta
    logp.beta = sweep(Z,1,logp.gr,'*');

    h.fn = H.rec(n.beta,Theta,X.Gm,P);
    h.gr = H.gr.rec(n.beta,Theta,X.Gm,P);

    # derivative of H w.r.t. all parameters
    others = sweep(lambda*h.gr,1,1+lambda*(h.fn-f),'/');
    all = cbind(logq.theta+logr.theta,logq.F+logr.F,logp.beta,logr.eta) - others;

    res = matrix(0,p,p);

    for(k in 1:2){
      tmp = all[II[[k]],];
      tt = as.matrix(scale(tmp,center=TRUE,scale=FALSE))*ns[k]/(ns[k]-1);
      res = res + t(tt)%*%tt;
    }

    return(res);
  }
  score.cov.dom = function(Theta){
    theta=Theta[1];
    F = Theta[2]
    beta = Theta[3:(2+n.beta)];

    eta = c(log(xi.dom(theta,F)),Theta[-(1:(2+n.beta))]);

    tr.fn = tr.dom(theta,F)
    tr.gr = tr.gr.dom(theta,F)

    logq.theta = logq.F = rep(0,n);

    logq.gr = tr.gr[[1]]/tr.fn

    # derivative of log(q) w.r.t. theta
    logq.theta[Gc==0&Gm==0] = logq.gr[1,1];
    logq.theta[Gc==1&Gm==0] = logq.gr[2,1];
    logq.theta[Gc==0&Gm==1] = logq.gr[1,2];
    logq.theta[Gc==1&Gm==1] = logq.gr[2,2];

    logq.gr = tr.gr[[2]]/tr.fn

    # derivative of log(q) w.r.t. F
    logq.F[Gc==0&Gm==0] = logq.gr[1,1];
    logq.F[Gc==1&Gm==0] = logq.gr[2,1];
    logq.F[Gc==0&Gm==1] = logq.gr[1,2];
    logq.F[Gc==1&Gm==1] = logq.gr[2,2];

    eta0.gr = logxi.gr.dom(theta,F);
    logr.gr = dalog.gr(Gm,X.Gm,eta)/dalog(Gm,X.Gm,eta);

    # derivative of log(r) w.r.t. theta
    logr.theta = logr.gr*eta0.gr[1];
    # derivative of log(r) w.r.t. F
    logr.F     = logr.gr*eta0.gr[2]
    # derivative of log(r) w.r.t. eta
    logr.eta   = sweep(X.Gm,1,logr.gr,'*');
    logp.gr = penetrance1.gr(Z,Y,beta)/penetrance1(Z,Y,beta);
    # derivative of log(p) w.r.t. beta
    logp.beta = sweep(Z,1,logp.gr,'*');

    h.fn = H.dom(n.beta,Theta,X.Gm,P);
    h.gr = H.gr.dom(n.beta,Theta,X.Gm,P);

    # derivative of H w.r.t. all parameters
    others = sweep(lambda*h.gr,1,1+lambda*(h.fn-f),'/');
    all = cbind(logq.theta+logr.theta,logq.F+logr.F,logp.beta,logr.eta) - others;

    res = matrix(0,p,p);

    for(k in 1:2){
      tmp = all[II[[k]],];
      tt = as.matrix(scale(tmp,center=TRUE,scale=FALSE))*ns[k]/(ns[k]-1);
      res = res + t(tt)%*%tt;
    }

    return(res);
  }
  score.cov.add = function(Theta){
    beta = Theta[3:(2+n.beta)];

    theta=Theta[1];
    F = Theta[2]
    eta = Theta[-(1:(2+n.beta))];

    logq.gr = tr.gr.add(theta)[[1]]/tr.add(theta)

    logq.theta = rep(NA,n);
    logq.theta[Gc==0&Gm==0] = logq.gr[1,1];
    logq.theta[Gc==0.5&Gm==0] = logq.gr[2,1];
    logq.theta[Gc==0&Gm==0.5] = logq.gr[1,2];
    logq.theta[Gc==0.5&Gm==0.5] = logq.gr[2,2];
    logq.theta[Gc==1&Gm==0.5] = logq.gr[3,2];
    logq.theta[Gc==0.5&Gm==1] = logq.gr[2,3];
    logq.theta[Gc==1&Gm==1] = logq.gr[3,3];


    logr.gr = sweep(dalog.gr.add(Gm,X.Gm,theta,F,eta),1,dalog.add(Gm,X.Gm,theta,F,eta),'/');

    logr.theta = logr.gr[,1];
    logr.F = logr.gr[,2]
    logr.eta = logr.gr[,-(1:2)];

    logp.gr = penetrance1.gr(Z,Y,beta)/penetrance1(Z,Y,beta);
    logp.beta = sweep(Z,1,logp.gr,'*');

    h.fn = H.add(n.beta,Theta,X.Gm,P);
    h.gr = H.gr.add(n.beta,Theta,X.Gm,P);
    others = sweep(lambda*h.gr,1,1+lambda*(h.fn-f),'/');

    all = cbind(logq.theta+logr.theta,logr.F,logp.beta,logr.eta) - others;

    res = matrix(0,p,p);

    for(k in 1:2){
      tmp = all[II[[k]],];
      tt = as.matrix(scale(tmp,center=TRUE,scale=FALSE))*ns[k]/(ns[k]-1);
      res = res + t(tt)%*%tt;
    }

    return(res);
  }

  V = hess;
  if(mode=='rec') U = score.cov.rec(est);
  if(mode=='dom') U = score.cov.dom(est);
  if(mode=='add') U = score.cov.add(est);

  if(HWE){
    U = U[-2,-2]
    V = V[-2,-2]
    est = est[-2]
  }
  V.inv = solve(V);
  cov = (V.inv%*%U)%*%V.inv;

  return(list(est =est,cov =cov,est.log=est.log,cov.log=cov.log ));

}

# Gm and X are independent
solution.ind = function(Y,Gc,Gm,Xo,Xc,Xm,G.main,G.int,mode,prev,HWE){
  n <- length(Y)
  Z = design.matrix(n,Gc,Gm,Xc,Xm,Xo,G.main,G.int)
  f = prev

  fit <- glm(Y ~ 0 + Z,family = binomial)
  res <- summary(fit)$coef
  est.log <- as.vector(res[,1])
  cov.log = vcov(fit)
  theta = 1-sqrt(mean(Gm==0))
  n.beta = ncol(Z) # number of all regression parameters in the penetrance model
  theta0 = log(f/(1-f))
  beta = rep(theta0,n.beta)
  beta[-1] = est.log[-1]
  Theta = c(theta,0,beta)

  P=list()
  if(mode=='add'){
    for(j in 1:3){
      P[[j]]=list()
      for(k in 1:3){
        P[[j]][[k]]=design.matrix(n,(j-1)/2,(k-1)/2,Xc,Xm,Xo,G.main,G.int)
      }
    }
  } else {
    for(j in 1:2){
      P[[j]]=list()
      for(k in 1:2){
        P[[j]][[k]]=design.matrix(n,j-1,k-1,Xc,Xm,Xo,G.main,G.int)
      }
    }
  }

  if(mode=='rec'){
    m00 = sum(Gc==0&Gm==0);
    m01 = sum(Gc==0&Gm==1);
    m10 = sum(Gc==1&Gm==0);
    m11 = sum(Gc==1&Gm==1);
  }
  if(mode=='dom'){
    m00 = sum(Gc==0&Gm==0);
    m01 = sum(Gc==0&Gm==1);
    m10 = sum(Gc==1&Gm==0);
    m11 = sum(Gc==1&Gm==1);
  }
  if(mode=='add'){
    m11 = sum(Gc==0&Gm==0);
    m12 = sum(Gc==0&Gm==0.5);
    m21 = sum(Gc==0.5&Gm==0);
    m22 = sum(Gc==0.5&Gm==0.5);
    m23 = sum(Gc==0.5&Gm==1);
    m32 = sum(Gc==1&Gm==0.5);
    m33 = sum(Gc==1&Gm==1);
  }

  n.beta = ncol(Z)

  n = length(Y);
  ns = c(sum(Y==1),sum(Y==0));



  II = list();
  II[[1]] = which(Y==0);
  II[[2]] = which(Y==1);

  lambda = ns[1]/n/f - ns[2]/n/(1-f);

  p = 2+n.beta; # number of all parameters

  upper = lower = Theta;
  upper[1] = 0.99;
  lower[1] = 0.01;
  if(HWE){lower[2]=upper[2]=0}else{
    lower[2] = -0.1
    upper[2] = 0.5
  }
  upper[-(1:2)] = Theta[-(1:2)] + 2;
  lower[-(1:2)] = Theta[-(1:2)] - 2;

  fn.rec = function(Theta){

    theta=Theta[1];
    F = ifelse(HWE,0,Theta[2])
    beta = Theta[-(1:2)];

    mc = mc.fn.rec(theta,F);
    h.fn = H0.rec(n,Theta,P);

    res = m00*log(mc[1,1])+m10*log(mc[2,1])+m01*log(mc[1,2])+m11*log(mc[2,2]);
    res = res + sum(log(penetrance1(Z,Y,beta)));
    res = res - sum(log(1+lambda*(h.fn-f)));
    -res
  }
  fn.dom = function(Theta){
    theta=Theta[1];
    F = ifelse(HWE,0,Theta[2])
    beta = Theta[-(1:2)];

    mc = mc.fn.dom(theta,F);
    h.fn = H0.dom(n,Theta,P);

    res = m00*log(mc[1,1])+m10*log(mc[2,1])+m01*log(mc[1,2])+m11*log(mc[2,2]);
    res = res + sum(log(penetrance1(Z,Y,beta)));
    res = res - sum(log(1+lambda*(h.fn-f)));
    -res
  }
  fn.add = function(Theta){
    theta=Theta[1];
    F = ifelse(HWE,0,Theta[2])
    beta = Theta[-(1:2)];

    mc = mc.fn.add(theta,F);
    h.fn = H0.add(n,Theta,P);

    res = m11*log(mc[1,1])+m12*log(mc[1,2])+m21*log(mc[2,1])+m22*log(mc[2,2])+m23*log(mc[2,3])+m32*log(mc[3,2])+m33*log(mc[3,3]);
    res = res + sum(log(penetrance1(Z,Y,beta)));
    res = res - sum(log(1+lambda*(h.fn-f)));
    return(-res);
  }
  gr.rec = function(Theta){

    theta=Theta[1];
    F = ifelse(HWE,0,Theta[2])
    beta = Theta[-(1:2)];

    mc.fn = mc.fn.rec(theta,F)
    mc.gr = mc.gr.rec(theta,F)

    logq.theta.F = list(NA,NA)
    logq.theta.F[[1]] = logq.theta.F[[2]] = rep(NA,n);

    for(k in 1:2){
      logq.gr = mc.gr[[k]]/mc.fn
      logq.theta.F[[k]][Gc==0&Gm==0] = logq.gr[1,1];
      logq.theta.F[[k]][Gc==1&Gm==0] = logq.gr[2,1];
      logq.theta.F[[k]][Gc==0&Gm==1] = logq.gr[1,2];
      logq.theta.F[[k]][Gc==1&Gm==1] = logq.gr[2,2];
    }

    logp.gr = penetrance1.gr(Z,Y,beta)/penetrance1(Z,Y,beta);
    logp.beta = sweep(Z,1,logp.gr,'*');

    h.fn = H0.rec(n,Theta,P);
    h.gr = H0.gr.rec(n,Theta,P);
    others = sweep(lambda*h.gr,1,1+lambda*(h.fn-f),'/');

    all = cbind(logq.theta.F[[1]],logq.theta.F[[2]],logp.beta) - others;

    colSums(-all)
  }
  gr.dom = function(Theta){

    theta=Theta[1];
    F = ifelse(HWE,0,Theta[2])
    beta = Theta[-(1:2)];

    mc.fn = mc.fn.dom(theta,F)
    mc.gr = mc.gr.dom(theta,F)

    logq.theta.F = list(NA,NA)
    logq.theta.F[[1]] = logq.theta.F[[2]] = rep(NA,n);

    for(k in 1:2){
      logq.gr = mc.gr[[k]]/mc.fn
      logq.theta.F[[k]][Gc==0&Gm==0] = logq.gr[1,1];
      logq.theta.F[[k]][Gc==1&Gm==0] = logq.gr[2,1];
      logq.theta.F[[k]][Gc==0&Gm==1] = logq.gr[1,2];
      logq.theta.F[[k]][Gc==1&Gm==1] = logq.gr[2,2];
    }

    logp.gr = penetrance1.gr(Z,Y,beta)/penetrance1(Z,Y,beta);
    logp.beta = sweep(Z,1,logp.gr,'*');

    h.fn = H0.dom(n,Theta,P);
    h.gr = H0.gr.dom(n,Theta,P);
    others = sweep(lambda*h.gr,1,1+lambda*(h.fn-f),'/');

    all = cbind(logq.theta.F[[1]],logq.theta.F[[2]],logp.beta) - others;

    colSums(-all)
  }
  gr.add = function(Theta){
    theta=Theta[1];
    F = ifelse(HWE,0,Theta[2])
    beta = Theta[-(1:2)];

    mc.gr = mc.gr.add(theta,F)
    mc.fn = mc.fn.add(theta,F)

    logq.theta.F = list(NA,NA)
    logq.theta.F[[1]] = logq.theta.F[[2]] = rep(NA,n);
    for(k in 1:2){
      logq.gr = mc.gr[[k]]/mc.fn
      logq.theta.F[[k]][Gc==0&Gm==0] = logq.gr[1,1];
      logq.theta.F[[k]][Gc==0.5&Gm==0] = logq.gr[2,1];
      logq.theta.F[[k]][Gc==0&Gm==0.5] = logq.gr[1,2];
      logq.theta.F[[k]][Gc==0.5&Gm==0.5] = logq.gr[2,2];
      logq.theta.F[[k]][Gc==1&Gm==0.5] = logq.gr[3,2];
      logq.theta.F[[k]][Gc==0.5&Gm==1] = logq.gr[2,3];
      logq.theta.F[[k]][Gc==1&Gm==1] = logq.gr[3,3];
    }

    logp.gr = penetrance1.gr(Z,Y,beta)/penetrance1(Z,Y,beta);
    logp.beta = sweep(Z,1,logp.gr,'*');

    h.fn = H0.add(n,Theta,P);
    h.gr = H0.gr.add(n,Theta,P);
    others = sweep(lambda*h.gr,1,1+lambda*(h.fn-f),'/');
    all = cbind(logq.theta.F[[1]],logq.theta.F[[2]],logp.beta) - others;
    colSums(-all)
  }

  if(mode=='rec'){
    fit = nloptr(x0=Theta,eval_f=fn.rec,eval_grad_f=gr.rec,lb=lower,ub=upper,opts=list("algorithm"="NLOPT_LD_LBFGS","xtol_rel"=1.0e-4) );
    hess = hessian(func=fn.rec,x=fit$solution);
  }
  if(mode=='dom'){
    fit = nloptr(x0=Theta,eval_f=fn.dom,eval_grad_f=gr.dom,lb=lower,ub=upper,opts=list("algorithm"="NLOPT_LD_LBFGS","xtol_rel"=1.0e-4) );
    hess = hessian(func=fn.dom,x=fit$solution);
  }
  if(mode=='add'){
    fit = nloptr(x0=Theta,eval_f=fn.add,eval_grad_f=gr.add,lb=lower,ub=upper,opts=list("algorithm"="NLOPT_LD_LBFGS","xtol_rel"=1.0e-4) );
    hess = hessian(func=fn.add,x=fit$solution);
  }
  est = fit$solution;

  score.cov.rec = function(Theta){


    theta=Theta[1];
    F = Theta[2]
    beta = Theta[-(1:2)];

    mc.fn = mc.fn.rec(theta,F)
    mc.gr = mc.gr.rec(theta,F)

    logq.theta.F = list(NA,NA)
    logq.theta.F[[1]] = logq.theta.F[[2]] = rep(NA,n);

    for(k in 1:2){
      logq.gr = mc.gr[[k]]/mc.fn
      logq.theta.F[[k]][Gc==0&Gm==0] = logq.gr[1,1];
      logq.theta.F[[k]][Gc==1&Gm==0] = logq.gr[2,1];
      logq.theta.F[[k]][Gc==0&Gm==1] = logq.gr[1,2];
      logq.theta.F[[k]][Gc==1&Gm==1] = logq.gr[2,2];
    }

    logp.gr = penetrance1.gr(Z,Y,beta)/penetrance1(Z,Y,beta);
    logp.beta = sweep(Z,1,logp.gr,'*');

    h.fn = H0.rec(n,Theta,P);
    h.gr = H0.gr.rec(n,Theta,P);
    others = sweep(lambda*h.gr,1,1+lambda*(h.fn-f),'/');

    all = cbind(logq.theta.F[[1]],logq.theta.F[[2]],logp.beta) - others;

    res = matrix(0,p,p);

    for(k in 1:2){
      tmp = all[II[[k]],];
      tt = as.matrix(scale(tmp,center=TRUE,scale=FALSE))*ns[k]/(ns[k]-1);
      res = res + t(tt)%*%tt;
    }

    return(res);
  }
  score.cov.dom = function(Theta){

    theta=Theta[1];
    F = Theta[2]
    beta = Theta[-(1:2)];

    mc.fn = mc.fn.dom(theta,F)
    mc.gr = mc.gr.dom(theta,F)

    logq.theta.F = list(NA,NA)
    logq.theta.F[[1]] = logq.theta.F[[2]] = rep(NA,n);

    for(k in 1:2){
      logq.gr = mc.gr[[k]]/mc.fn
      logq.theta.F[[k]][Gc==0&Gm==0] = logq.gr[1,1];
      logq.theta.F[[k]][Gc==1&Gm==0] = logq.gr[2,1];
      logq.theta.F[[k]][Gc==0&Gm==1] = logq.gr[1,2];
      logq.theta.F[[k]][Gc==1&Gm==1] = logq.gr[2,2];
    }

    logp.gr = penetrance1.gr(Z,Y,beta)/penetrance1(Z,Y,beta);
    logp.beta = sweep(Z,1,logp.gr,'*');

    h.fn = H0.dom(n,Theta,P);
    h.gr = H0.gr.dom(n,Theta,P);
    others = sweep(lambda*h.gr,1,1+lambda*(h.fn-f),'/');

    all = cbind(logq.theta.F[[1]],logq.theta.F[[2]],logp.beta) - others;

    res = matrix(0,p,p);

    for(k in 1:2){
      tmp = all[II[[k]],];
      tt = as.matrix(scale(tmp,center=TRUE,scale=FALSE))*ns[k]/(ns[k]-1);
      res = res + t(tt)%*%tt;
    }

    return(res);
  }
  score.cov.add = function(Theta){
    theta=Theta[1];
    F = Theta[2]
    beta = Theta[-(1:2)];

    mc.gr = mc.gr.add(theta,F)
    mc.fn = mc.fn.add(theta,F)

    logq.theta.F = list(NA,NA)
    logq.theta.F[[1]] = logq.theta.F[[2]] = rep(NA,n);
    for(k in 1:2){
      logq.gr = mc.gr[[k]]/mc.fn
      logq.theta.F[[k]][Gc==0&Gm==0] = logq.gr[1,1];
      logq.theta.F[[k]][Gc==0.5&Gm==0] = logq.gr[2,1];
      logq.theta.F[[k]][Gc==0&Gm==0.5] = logq.gr[1,2];
      logq.theta.F[[k]][Gc==0.5&Gm==0.5] = logq.gr[2,2];
      logq.theta.F[[k]][Gc==1&Gm==0.5] = logq.gr[3,2];
      logq.theta.F[[k]][Gc==0.5&Gm==1] = logq.gr[2,3];
      logq.theta.F[[k]][Gc==1&Gm==1] = logq.gr[3,3];
    }

    logp.gr = penetrance1.gr(Z,Y,beta)/penetrance1(Z,Y,beta);
    logp.beta = sweep(Z,1,logp.gr,'*');

    h.fn = H0.add(n,Theta,P);
    h.gr = H0.gr.add(n,Theta,P);
    others = sweep(lambda*h.gr,1,1+lambda*(h.fn-f),'/');
    all = cbind(logq.theta.F[[1]],logq.theta.F[[2]],logp.beta) - others;

    res = matrix(0,p,p);

    for(k in 1:2){
      tmp = all[II[[k]],];
      tt = as.matrix(scale(tmp,center=TRUE,scale=FALSE))*ns[k]/(ns[k]-1);
      res = res + t(tt)%*%tt;
    }
    return(res);
  }

  V = hess;
  if(mode=='rec') U = score.cov.rec(est);
  if(mode=='dom') U = score.cov.dom(est);
  if(mode=='add') U = score.cov.add(est);

  if(HWE){
    U = U[-2,-2]
    V = V[-2,-2]
    est = est[-2]
  }

  V.inv = solve(V);
  cov = (V.inv%*%U)%*%V.inv;

  return(list(est =est,cov=cov,est.log=est.log,cov.log=cov.log));

}

# design matrix
design.matrix = function(n,Gc,Gm,Xc,Xm,Xo,G.main,G.int){
  I <- numeric(2)
  I[1] <- ifelse('Gc' %in% G.main,1,0)
  I[2] <- ifelse('Gm' %in% G.main,1,0)
  if(length(Gc)==1){
    l=rep(1,n)
    Gc=Gc*l
    Gm=Gm*l
  }
  cbind(1,cbind(Gc,Gm)[,which(I==1)],Xo,Xc*Gc,Xm*Gm,Gc*Gm[G.int])
}

#logistic  sigmoid pr(Y|Gc,Gm,X;beta)
penetrance1=function(E,Y,beta){
  tmp = exp(as.vector(E%*%beta))
  return(1-Y+(2*Y-1)*tmp/(1+tmp))
}

#logistic  sigmoid pr(Y=1|Gc,Gm,X;beta)
penetrance0=function(E,beta){
  tmp = exp(as.vector(E%*%beta))
  return(tmp/(1+tmp))
}

# d(pr(Y|Gc,Gm,X;beta))/d(beta)
penetrance1.gr = function(E,Y,beta){
  tmp = exp(as.vector(E%*%beta))
  return((2*Y-1)*tmp/(1+tmp)^2)
}

# d(pr(Y=1|Gc,Gm,X;beta))/d(beta)
penetrance0.gr = function(E,beta){
  tmp = exp(as.vector(E%*%beta))
  return(tmp/(1+tmp)^2)
}

#Pr(Gc|Gm;theta)
tr.rec = function(theta,F){
  D = F/(1-F)
  res = matrix(0,2,2);
  res[2,1] = theta^2/(1+theta+D);
  res[2,2] = theta;
  res[1,] = 1-res[2,];
  return(res);
}
tr.dom = function(theta,F){
  D = F/(1-F)
  res = matrix(0,2,2);
  res[1,1] = 1-theta;
  res[1,2] = (1-theta)^2/(2-theta+D);
  res[2,] = 1-res[1,];
  return(res);
}
tr.add = function(theta){
  res = matrix(0,3,3);
  res[1,1] = 1-theta;
  res[2,1] = theta;
  res[1,2] = (1-theta)/2;
  res[2,2] = 1/2;
  res[3,2] = theta/2;
  res[2,3] = 1-theta;
  res[3,3] = theta;
  return(res);
}

# d(Pr(Gc|Gm;theta))/d(theta)
tr.gr.rec = function(theta,F){
  res = list(NULL,NULL)
  res[[1]] = res[[2]] = matrix(0,2,2);
  D = F/(1-F)

  res[[1]][2,1] = 2*theta/(1+theta+D)-theta^2/(1+theta+D)^2;
  res[[1]][2,2] = 1;
  res[[1]][1,] = -res[[1]][2,];

  res[[2]][2,1] = -theta^2/(1+theta+D)^2/(1-F)^2
  res[[2]][1,1] = -res[[2]][2,1]

  return(res);
}
tr.gr.dom = function(theta,F){
  res = list(NULL,NULL)
  res[[1]] = res[[2]] = matrix(0,2,2);
  D = F/(1-F)

  res[[1]][1,1] = -1;
  res[[1]][1,2] = -2*(1-theta)/(2-theta+D)+(1-theta)^2/(2-theta+D)^2;
  res[[1]][2,] = -res[[1]][1,];

  res[[2]][1,2] = -(1-theta)^2/(2-theta+D)^2;
  res[[2]][2,2] = -res[[2]][1,2]

  return(res);
}
tr.gr.add = function(theta){
  res = list(NULL,NULL)
  res[[1]] =res[[2]] = matrix(0,3,3);
  res[[1]][1,1] = -1;
  res[[1]][2,1] = 1;
  res[[1]][1,2] = -1/2;
  res[[1]][3,2] = 1/2;
  res[[1]][2,3] = -1;
  res[[1]][3,3] = 1;
  return(res);
}

# xi(theta,F)
xi.rec = function(theta,F){
  tmp = (1-F)*theta^2+F*theta;
  tmp/(1-tmp)
}
xi.dom = function(theta,F){
  tmp = (1-F)*(1-theta)^2+F*(1-theta);
  (1-tmp)/tmp
}
xi.add= function(theta,F){
  c((1-F)*(1-theta)^2+F*(1-theta),2*(1-F)*theta*(1-theta),(1-F)*theta^2+F*theta)
}

#d(xi(theta))/d(theta)
xi.gr.rec = function(theta,F){
  c(2*(1-F)*theta+F,theta-theta^2)/(1-(1-F)*theta^2-F*theta)^2
}
xi.gr.dom = function(theta,F){
  tmp = 1-theta
  c(2*(1-F)*tmp+F,tmp^2-tmp)/((1-F)*tmp^2+F*tmp)^2
}
xi.gr.add= function(theta,F){
  list(c(-2*(1-F)*(1-theta)-F,2*(1-F)*(1-2*theta),2*(1-F)*theta+F),c(1-theta-(1-theta)^2,-2*theta*(1-theta),theta-theta^2))
}

#d(eta(theta))/d(theta)
logxi.gr.rec = function(theta,F){
  return(xi.gr.rec(theta,F)/xi.rec(theta,F));
}
logxi.gr.dom = function(theta,F){
  return(xi.gr.dom(theta,F)/xi.dom(theta,F));
}

# Pr(Gm|X;eta)
dalog = function(Gm,X,eta){
  # Gm = 0 or 1
  tmp = exp(eta[1]+as.numeric(X%*%eta[-1]));
  return(1-Gm+(2*Gm-1)*tmp/(1+tmp));
}
dalog.add = function(Gm,X,theta,F,eta){
  # additive mode: Gm = 0 or 1/2 or 1
  Gm1 = Gm*2+1;
  tmp = as.numeric(X%*%eta);
  etmp.2 = exp(tmp/2);
  xi = xi.add(theta,F);
  tmp0 = xi[Gm1]*exp(Gm*tmp);
  tmp1 = xi[1]+xi[2]*etmp.2+xi[3]*etmp.2^2;
  return(tmp0/tmp1);
}

# d(Pr(Gm|X;eta)/d(eta)
dalog.gr = function(Gm,X,eta){
  # Gm = 0 or 1
  tmp = exp(eta[1]+as.numeric(X%*%eta[-1]));
  return((2*Gm-1)*tmp/(1+tmp)^2);
}
dalog.gr.add = function(Gm,X,theta,F,eta){
  Gm1 = Gm*2+1;
  tmp = as.numeric(X%*%eta);
  xi = xi.add(theta,F);
  xi.d = xi.gr.add(theta,F);

  etmp.2 = exp(tmp/2);
  etmp.1 = etmp.2^2;
  eGtmp = exp(Gm*tmp);

  tmp0 = xi[1] + xi[2]*etmp.2 + xi[3]*etmp.1;
  tmp11 = xi.d[[1]][Gm1]*eGtmp;
  tmp12 = xi.d[[2]][Gm1]*eGtmp;
  tmp2 = xi[Gm1]*eGtmp;

  tmp31 = tmp11/tmp0 - (xi.d[[1]][1]+xi.d[[1]][2]*etmp.2 + xi.d[[1]][3]*etmp.1) * tmp2/tmp0^2;
  tmp32 = tmp12/tmp0 - (xi.d[[2]][1]+xi.d[[2]][2]*etmp.2 + xi.d[[2]][3]*etmp.1) * tmp2/tmp0^2;
  tmp4 = sweep(X,1,tmp2*Gm/tmp0 - tmp2/tmp0^2*(xi[2]*etmp.2/2 + xi[3]*etmp.1),'*');

  return(cbind(tmp31,tmp32,tmp4));
}

# H_i = \sum_{Gm,Gc} Pr(Y=1|Gm,Gc,X_i)Pr(Gc|Gm)Pr(Gm|X_i)
H.rec = function(n.beta,Theta,X,P){
  theta = Theta[1];
  F = Theta[2]
  beta = Theta[3:(2+n.beta)];

  eta0 = log(xi.rec(theta,F));
  eta = c(eta0,Theta[-(1:(2+n.beta))]);

  n = nrow(X);
  h = rep(0,n);
  gt = tr.rec(theta,F);

  for(k in 1:2){
    for(j in 1:2){
      h = h + penetrance0(P[[j]][[k]],beta)*dalog(k-1,X,eta)*gt[j,k];
    }
  }
  return(h);
}
H.dom = function(n.beta,Theta,X,P){
  theta = Theta[1];
  F = Theta[2]
  beta = Theta[3:(2+n.beta)];

  eta0 = log(xi.dom(theta,F));
  eta = c(eta0,Theta[-(1:(2+n.beta))]);

  n = nrow(X);
  h = rep(0,n);
  gt = tr.dom(theta,F);

  for(k in 1:2){
    for(j in 1:2){
      h = h + penetrance0(P[[j]][[k]],beta)*dalog(k-1,X,eta)*gt[j,k];
    }
  }
  return(h);
}
H.add = function(n.beta,Theta,X,P){

  theta = Theta[1];
  F = Theta[2]
  beta = Theta[3:(2+n.beta)];
  eta = Theta[-(1:(2+n.beta))];

  n = nrow(X);
  h = rep(0,n);
  gt = tr.add(theta);

  for(k in 1:3){
    for(j in 1:3){
      h = h + penetrance0(P[[j]][[k]],beta)*dalog.add((k-1)/2,X,theta,F,eta)*gt[j,k];
    }
  }
  return(h);
}

# derivative of H_i
H.gr.rec = function(n.beta,Theta,X,P){

  theta = Theta[1];
  F = Theta[2]
  beta = Theta[3:(2+n.beta)];
  eta0 = log(xi.rec(theta,F));
  eta = c(eta0,Theta[-(1:(2+n.beta))]);

  n = nrow(X);
  n.X = ncol(X);
  p = length(Theta);

  h.gr = matrix(0,n,p);

  eta0.gr = logxi.gr.rec(theta,F);
  q.fn = tr.rec(theta,F);
  q.gr = tr.gr.rec(theta,F);

  r.fn = cbind(dalog(0,X,eta),dalog(1,X,eta));
  r.gr = cbind(dalog.gr(0,X,eta),dalog.gr(1,X,eta));

  for(k in 1:2){
    for(j in 1:2){
      p.fn = penetrance0(P[[j]][[k]],beta)
      p.gr = penetrance0.gr(P[[j]][[k]],beta)

      h.gr.theta = p.fn*(q.gr[[1]][j,k]*r.fn[,k]+q.fn[j,k]*r.gr[,k]*eta0.gr[1])
      h.gr.F     = p.fn*(q.gr[[2]][j,k]*r.fn[,k]+q.fn[j,k]*r.gr[,k]*eta0.gr[2])
      h.gr.eta = sweep(X,1,p.fn*q.fn[j,k]*r.gr[,k],'*')
      tmp = p.gr*q.fn[j,k]*r.fn[,k]
      h.gr[,1] = h.gr[,1] + h.gr.theta
      h.gr[,2] = h.gr[,2] + h.gr.F
      h.gr[,3:(2+n.beta)] = h.gr[,3:(2+n.beta)] + sweep(P[[j]][[k]],1,tmp,'*')
      h.gr[,-(1:(2+n.beta))] = h.gr[,-(1:(2+n.beta))] + h.gr.eta
    }
  }

  return(h.gr)

}
H.gr.dom = function(n.beta,Theta,X,P){

  theta = Theta[1];
  F = Theta[2]
  beta = Theta[3:(2+n.beta)];
  eta0 = log(xi.dom(theta,F));
  eta = c(eta0,Theta[-(1:(2+n.beta))]);

  n = nrow(X);
  n.X = ncol(X);
  p = length(Theta);

  h.gr = matrix(0,n,p);

  eta0.gr = logxi.gr.dom(theta,F);
  q.fn = tr.dom(theta,F);
  q.gr = tr.gr.dom(theta,F);

  r.fn = cbind(dalog(0,X,eta),dalog(1,X,eta));
  r.gr = cbind(dalog.gr(0,X,eta),dalog.gr(1,X,eta));

  for(k in 1:2){
    for(j in 1:2){
      p.fn = penetrance0(P[[j]][[k]],beta)
      p.gr = penetrance0.gr(P[[j]][[k]],beta)

      h.gr.theta = p.fn*(q.gr[[1]][j,k]*r.fn[,k]+q.fn[j,k]*r.gr[,k]*eta0.gr[1])
      h.gr.F     = p.fn*(q.gr[[2]][j,k]*r.fn[,k]+q.fn[j,k]*r.gr[,k]*eta0.gr[2])
      h.gr.eta = sweep(X,1,p.fn*q.fn[j,k]*r.gr[,k],'*')
      tmp = p.gr*q.fn[j,k]*r.fn[,k]
      h.gr[,1] = h.gr[,1] + h.gr.theta
      h.gr[,2] = h.gr[,2] + h.gr.F
      h.gr[,3:(2+n.beta)] = h.gr[,3:(2+n.beta)] + sweep(P[[j]][[k]],1,tmp,'*')
      h.gr[,-(1:(2+n.beta))] = h.gr[,-(1:(2+n.beta))] + h.gr.eta
    }
  }


  return(h.gr)

}
H.gr.add = function(n.beta,Theta,X,P){

  theta = Theta[1];
  F = Theta[2]
  beta = Theta[3:(2+n.beta)];
  eta = Theta[-(1:(2+n.beta))];

  n = nrow(X);
  p = length(Theta);
  n.X = ncol(X);

  h.gr = matrix(0,n,p);

  q.fn = tr.add(theta);
  q.gr = tr.gr.add(theta)[[1]];

  r.fn = cbind(dalog.add(0,X,theta,F,eta),dalog.add(0.5,X,theta,F,eta),dalog.add(1,X,theta,F,eta));
  r.gr = list(NA,NA,NA);
  r.gr[[1]] = dalog.gr.add(0,X,theta,F,eta); # derivative of pr(Gm|X) w.r.t. (theta,F,eta)
  r.gr[[2]] = dalog.gr.add(0.5,X,theta,F,eta);
  r.gr[[3]] = dalog.gr.add(1,X,theta,F,eta);

  for(k in 1:3){
    for(j in 1:3){
      p.fn = penetrance0(P[[j]][[k]],beta)
      p.gr = penetrance0.gr(P[[j]][[k]],beta)

      h.gr.theta = p.fn*(q.gr[j,k]*r.fn[,k]+q.fn[j,k]*r.gr[[k]][,1]);
      h.gr.F = p.fn*q.fn[j,k]*r.gr[[k]][,2];
      if(n.X==1){
        h.gr.eta = r.gr[[k]][,-(1:2)]*p.fn*q.fn[j,k];
      }else{
        h.gr.eta = sweep(r.gr[[k]][,-(1:2)],1,p.fn*q.fn[j,k],'*');
      }
      tmp = p.gr*q.fn[j,k]*r.fn[,k]

      h.gr[,1] = h.gr[,1] + h.gr.theta
      h.gr[,2] = h.gr[,2] + h.gr.F
      h.gr[,3:(2+n.beta)] = h.gr[,3:(2+n.beta)] + sweep(P[[j]][[k]],1,tmp,'*')
      h.gr[,-(1:(2+n.beta))] = h.gr[,-(1:(2+n.beta))] + h.gr.eta
    }
  }

  return(h.gr)

}

######### IND #########
# Gm and X are independent
mc.fn.add = function(theta,F){
  ps = xi.add(theta,F)
  res = matrix(0,3,3);
  res[1,1] = ps[1]^2+ps[1]*ps[2]/2;
  res[1,2] = ps[1]*ps[2]/2+ps[2]^2/4
  res[2,1] = ps[1]*ps[2]/2+ps[1]*ps[3];
  res[2,2] = (ps[1]*ps[2]+ps[2]^2+ps[2]*ps[3])/2;
  res[2,3] = ps[1]*ps[3]+ps[2]*ps[3]/2
  res[3,2] = ps[3]*ps[2]/2+ps[2]^2/4;
  res[3,3] = ps[3]^2+ps[2]*ps[3]/2
  res
}
mc.gr.add = function(theta,F){
  res = list(NA,NA)
  res[[1]] = res[[2]] = matrix(0,3,3);
  ps = xi.add(theta,F)
  ps.gr = xi.gr.add(theta,F)
  for(i in 1:2){
    res[[i]][1,1] = 2*ps[1]*ps.gr[[i]][1]+(ps.gr[[i]][1]*ps[2]+ps.gr[[i]][2]*ps[1])/2
    res[[i]][1,2] = (ps.gr[[i]][1]*ps[2]+ps.gr[[i]][2]*ps[1]+ps[2]*ps.gr[[i]][2])/2
    res[[i]][2,1] = (ps.gr[[i]][1]*ps[2]+ps.gr[[i]][2]*ps[1])/2+ps.gr[[i]][1]*ps[3]+ps[1]*ps.gr[[i]][3]
    res[[i]][2,2] = (ps.gr[[i]][1]*ps[2]+ps.gr[[i]][2]*ps[1]+ps.gr[[i]][3]*ps[2]+ps.gr[[i]][2]*ps[3])/2+ps[2]*ps.gr[[i]][2]
    res[[i]][2,3] = ps.gr[[i]][1]*ps[3]+ps[1]*ps.gr[[i]][3]+(ps.gr[[i]][3]*ps[2]+ps.gr[[i]][2]*ps[3])/2
    res[[i]][3,2] = ps.gr[[i]][2]*ps[2]/2+(ps.gr[[i]][3]*ps[2]+ps.gr[[i]][2]*ps[3])/2
    res[[i]][3,3] = 2*ps[3]*ps.gr[[i]][3]+(ps.gr[[i]][3]*ps[2]+ps.gr[[i]][2]*ps[3])/2
  }
  res
}
mc.fn.dom = function(theta,F){
  pp = mc.fn.add(theta,F)
  res = matrix(0,2,2);
  res[1,1] = pp[1,1];
  res[1,2] = pp[1,2]+pp[1,3]
  res[2,1] = pp[2,1]+pp[3,1]
  res[2,2] = 1-sum(res)
  res
}
mc.gr.dom = function(theta,F){
  pp = mc.gr.add(theta,F)
  res = list(NA,NA)
  res[[1]]=res[[2]]=matrix(0,2,2);
  for(i in 1:2){
    res[[i]][1,1] = pp[[i]][1,1];
    res[[i]][1,2] = pp[[i]][1,2]+pp[[i]][1,3]
    res[[i]][2,1] = pp[[i]][2,1]+pp[[i]][3,1]
    res[[i]][2,2] = -sum(res[[i]])
  }
  res
}
mc.fn.rec = function(theta,F){
  pp = mc.fn.add(theta,F)
  res = matrix(0,2,2);
  res[2,2] = pp[3,3];
  res[1,2] = pp[1,3]+pp[2,3]
  res[2,1] = pp[3,1]+pp[3,2]
  res[1,1] = 1-sum(res)
  res
}
mc.gr.rec = function(theta,F){
  pp = mc.gr.add(theta,F)
  res = list(NA,NA)
  res[[1]]=res[[2]]=matrix(0,2,2);
  for(i in 1:2){
    res[[i]][2,2] = pp[[i]][3,3];
    res[[i]][1,2] = pp[[i]][1,3]+pp[[i]][2,3]
    res[[i]][2,1] = pp[[i]][3,1]+pp[[i]][3,2]
    res[[i]][1,1] = -sum(res[[i]])
  }
  res
}

# H_i = \sum_{Gm,Gc} Pr(Y=1|Gm,Gc,X_i)Pr(Gc|Gm)Pr(Gm|X_i)
H0.rec = function(n,Theta,P){

  theta = Theta[1];
  F = Theta[2]
  beta = Theta[-(1:2)];

  h = rep(0,n);
  mc = mc.fn.rec(theta,F); # joint genotype distribution

  for(k in 1:2){
    for(j in 1:2){
      h = h + penetrance0(P[[j]][[k]],beta)*mc[j,k];
    }
  }

  return(h);
}
H0.dom = function(n,Theta,P){

  theta = Theta[1];
  F = Theta[2]
  beta = Theta[-(1:2)];

  h = rep(0,n);
  mc = mc.fn.dom(theta,F); # joint genotype distribution

  for(k in 1:2){
    for(j in 1:2){
      h = h + penetrance0(P[[j]][[k]],beta)*mc[j,k];
    }
  }

  return(h);
}
H0.add = function(n,Theta,P){

  theta = Theta[1];
  F = Theta[2]
  beta = Theta[-(1:2)];

  h = rep(0,n);
  mc = mc.fn.add(theta,F); # joint genotype distribution

  for(k in 1:3){
    for(j in 1:3){
      h = h + penetrance0(P[[j]][[k]],beta)*mc[j,k];
    }
  }

  return(h);
}

H0.gr.rec = function(n,Theta,P){

  theta = Theta[1];
  F = Theta[2]
  beta = Theta[-(1:2)];

  p = length(Theta);

  h.gr = matrix(0,n,p);

  q.fn = mc.fn.rec(theta,F);
  q.gr = mc.gr.rec(theta,F);

  for(k in 1:2){
    for(j in 1:2){
      p.fn = penetrance0(P[[j]][[k]],beta)
      p.gr = penetrance0.gr(P[[j]][[k]],beta)

      h.gr[,1] = h.gr[,1] + q.gr[[1]][j,k]*p.fn;
      h.gr[,2] = h.gr[,2] + q.gr[[2]][j,k]*p.fn;
      h.gr[,3:p] = h.gr[,3:p] + sweep(P[[j]][[k]],1,q.fn[j,k]*p.gr,'*')
    }
  }

  return(h.gr);

}
H0.gr.dom = function(n,Theta,P){

  theta = Theta[1];
  F = Theta[2]
  beta = Theta[-(1:2)];

  p = length(Theta);

  h.gr = matrix(0,n,p);

  q.fn = mc.fn.dom(theta,F);
  q.gr = mc.gr.dom(theta,F);

  for(k in 1:2){
    for(j in 1:2){
      p.fn = penetrance0(P[[j]][[k]],beta)
      p.gr = penetrance0.gr(P[[j]][[k]],beta)

      h.gr[,1] = h.gr[,1] + q.gr[[1]][j,k]*p.fn;
      h.gr[,2] = h.gr[,2] + q.gr[[2]][j,k]*p.fn;
      h.gr[,3:p] = h.gr[,3:p] + sweep(P[[j]][[k]],1,q.fn[j,k]*p.gr,'*')
    }
  }

  return(h.gr);

}
H0.gr.add = function(n,Theta,P){

  theta = Theta[1];
  F = Theta[2]
  beta = Theta[-(1:2)];

  p = length(Theta);

  h.gr = matrix(0,n,p);

  q.fn = mc.fn.add(theta,F);
  q.gr = mc.gr.add(theta,F);

  for(k in 1:3){
    for(j in 1:3){
      p.fn = penetrance0(P[[j]][[k]],beta)
      p.gr = penetrance0.gr(P[[j]][[k]],beta)

      h.gr[,1] = h.gr[,1] + q.gr[[1]][j,k]*p.fn;
      h.gr[,2] = h.gr[,2] + q.gr[[2]][j,k]*p.fn
      h.gr[,3:p] = h.gr[,3:p] + sweep(P[[j]][[k]],1,q.fn[j,k]*p.gr,'*')
    }
  }

  return(h.gr);

}

######### Imprinting #########


imprinting = function(Y,gmm,gcc,X,loci,hap,f,ppi){

  F = distinguish(gmm, gcc, loci, hap)
  gm = F$gm
  gc = F$gc
  gcm = F$gcm
  gcp = F$gcp
  phi = F$phi

  gm.log = gm[phi==1]
  gc.log = gc[phi==1]
  gcm.log = gcm[phi==1]
  gcp.log = gcp[phi==1]
  Y.log = Y[phi==1]
  X.log = X[phi==1]

  Z = Im.design.matrix1(gm,gc,gcm,gcp,X)
  if(mean(phi==1)==1){
    Z1 = Z
    Z2 = Z
    Pm = rep(1,length(Y))
    Pp = rep(0,length(Y))
  } else{
    Z1 = Im.design.matrix2_1(gm,gc,gcm,gcp,X,phi)
    Z2 = Im.design.matrix2_2(gm,gc,gcm,gcp,X,phi)
    W = distinguish0(gmm[-which(phi==1),],gcc[-which(phi==1),],loci,hap,ppi)
    Pm = rep(1,length(Y))
    Pp = rep(0,length(Y))
    Pm[-which(phi==1)] = W$Pm
    Pp[-which(phi==1)] = W$Pp
  }
  P1 = Pm/(Pm+Pp)
  P2 = Pp/(Pm+Pp)

  Z.log = Im.design.matrix1(gm.log,gc.log,gcm.log,gcp.log,X.log)
  fit = glm(Y.log ~ 0 + Z.log,family = binomial)
  res = summary(fit)$coef
  est.log = as.vector(res[,1])
  cov.log = vcov(fit)
  pval.log = as.vector(res[,4])

  theta = 1-sqrt(mean(gm==0))
  n.beta = ncol(Z) # number of all regression parameters in the penetrance model
  n.eta = max(1,ncol(X)) # number of covariates possibly associated with Gm
  theta0 = log(f/(1-f))
  beta = rep(theta0,n.beta)
  beta[-1] = est.log[-1]
  Theta = c(theta,beta,rep(0,n.eta))
  n = length(Y)
  ns = c(sum(Y==1),sum(Y==0))
  lambda = ns[1]/n/f - ns[2]/n/(1-f)
  p = length(Theta)

  II = list()
  II[[1]] = which(Y==0)
  II[[2]] = which(Y==1)

  upper = lower = Theta
  upper[1] = 0.99
  lower[1] = 0.01
  upper[-1] = Theta[-1] + 2
  lower[-1] = Theta[-1] - 2

  #number of (gc,gm)
  m11 = sum(gc==0&gm==0)
  m12 = sum(gc==0&gm==1)
  m21 = sum(gc==1&gm==0)
  m22 = sum(gc==1&gm==1)
  m23 = sum(gc==1&gm==2)
  m32 = sum(gc==2&gm==1)
  m33 = sum(gc==2&gm==2)

  #likelihood
  fn.add = function(Theta){
    theta=Theta[1]
    beta = Theta[2:(1+n.beta)]
    eta = Theta[-(1:(1+n.beta))]

    tr = Im.tr.add(theta)
    h.fn = Im.H.add(n.beta,Theta,X)

    res = m11*log(tr[1,1])+m12*log(tr[1,2])+m21*log(tr[2,1])+m22*log(tr[2,2])+m23*log(tr[2,3])+m32*log(tr[3,2])+m33*log(tr[3,3])
    res = res + sum(log(Im.dalog.add(gm,X,theta,eta)))
    res = res + sum(log(P1*Im.penetrance1(Z1,Y,beta)+P2*Im.penetrance1(Z2,Y,beta)))
    res = res - sum(log(1+lambda*(h.fn-f)))
    return(-res);
  }

  #gradient
  gr.add = function(Theta){
    theta=Theta[1]
    beta = Theta[2:(1+n.beta)]
    eta = Theta[-(1:(1+n.beta))]

    logq.gr = Im.tr.gr.add(theta)/Im.tr.add(theta)

    logq.theta = rep(NA,n)
    logq.theta[gc==0&gm==0] = logq.gr[1,1]
    logq.theta[gc==1&gm==0] = logq.gr[2,1]
    logq.theta[gc==0&gm==1] = logq.gr[1,2]
    logq.theta[gc==1&gm==1] = logq.gr[2,2]
    logq.theta[gc==2&gm==1] = logq.gr[3,2]
    logq.theta[gc==1&gm==2] = logq.gr[2,3]
    logq.theta[gc==2&gm==2] = logq.gr[3,3]


    logr.gr = sweep(Im.dalog.gr.add(gm,X,theta,eta),1,Im.dalog.add(gm,X,theta,eta),'/')

    logr.theta = logr.gr[,1]
    logr.eta = logr.gr[,-1]

    logp.gr1 = P1*Im.penetrance1.gr(Z1,Y,beta)/(P1*Im.penetrance1(Z1,Y,beta)+P2*Im.penetrance1(Z2,Y,beta))
    logp.gr2 = P2*Im.penetrance1.gr(Z2,Y,beta)/(P1*Im.penetrance1(Z1,Y,beta)+P2*Im.penetrance1(Z2,Y,beta))
    logp.beta = sweep(Z1,1,logp.gr1,'*') + sweep(Z2,1,logp.gr2,'*')

    h.fn = Im.H.add(n.beta,Theta,X)
    h.gr = Im.H.gr.add(n.beta,Theta,X)
    others = sweep(lambda*h.gr,1,1+lambda*(h.fn-f),'/')

    all = cbind(logq.theta+logr.theta,logp.beta,logr.eta) - others
    -colSums(all)
  }


  fit = nloptr(x0=Theta,eval_f=fn.add,eval_grad_f=gr.add,lb=lower,ub=upper,opts=list("algorithm"="NLOPT_LD_LBFGS","xtol_rel"=1.0e-4) )
  hess = hessian(func=fn.add,x=fit$solution)

  est = fit$solution

  #covariance
  score.cov.add = function(Theta){

    theta=Theta[1]
    beta = Theta[2:(1+n.beta)]
    eta = Theta[-(1:(1+n.beta))]

    logq.gr = Im.tr.gr.add(theta)/Im.tr.add(theta)

    logq.theta = rep(NA,n)
    logq.theta[gc==0&gm==0] = logq.gr[1,1]
    logq.theta[gc==1&gm==0] = logq.gr[2,1]
    logq.theta[gc==0&gm==1] = logq.gr[1,2]
    logq.theta[gc==1&gm==1] = logq.gr[2,2]
    logq.theta[gc==2&gm==1] = logq.gr[3,2]
    logq.theta[gc==1&gm==2] = logq.gr[2,3]
    logq.theta[gc==2&gm==2] = logq.gr[3,3]


    logr.gr = sweep(Im.dalog.gr.add(gm,X,theta,eta),1,Im.dalog.add(gm,X,theta,eta),'/')

    logr.theta = logr.gr[,1]
    logr.eta = logr.gr[,-1]

    logp.gr1 = P1*Im.penetrance1.gr(Z1,Y,beta)/(P1*Im.penetrance1(Z1,Y,beta)+P2*Im.penetrance1(Z2,Y,beta))
    logp.gr2 = P2*Im.penetrance1.gr(Z2,Y,beta)/(P1*Im.penetrance1(Z1,Y,beta)+P2*Im.penetrance1(Z2,Y,beta))
    logp.beta = sweep(Z1,1,logp.gr1,'*') + sweep(Z2,1,logp.gr2,'*')

    h.fn = Im.H.add(n.beta,Theta,X)
    h.gr = Im.H.gr.add(n.beta,Theta,X)
    others = sweep(lambda*h.gr,1,1+lambda*(h.fn-f),'/')

    all = cbind(logq.theta+logr.theta,logp.beta,logr.eta) - others

    res = matrix(0,p,p)

    for(k in 1:2){
      tmp = all[II[[k]],]
      tt = as.matrix(scale(tmp,center=TRUE,scale=FALSE))*ns[k]/(ns[k]-1)
      res = res + t(tt)%*%tt;
    }

    return(res);
  }

  V = hess
  U = score.cov.add(est);

  V.inv = solve(V);
  cov = (V.inv%*%U)%*%V.inv;

  return(list(est =est,cov =cov,est.log=est.log,cov.log=cov.log))

}


#identify the parental origin of off-spring alleles

distinguish = function(gmm,gcc,loci,hap){

  l = nrow(gmm)
  K = nrow(hap)
  gm = gmm[,loci]
  gc = gcc[,loci]
  gcm = rep(3,l)
  gcp = rep(0,l)
  phi = rep(0,l)

  for(u in 1:l){

    if(gc[u]==0){
      gcm[u] = 0
      gcp[u] = 0
      phi[u] = 1
    }
    else if(gm[u]==0&gc[u]==1){
      gcm[u] = 0
      gcp[u] = 1
      phi[u] = 1
    }
    else if(gm[u]==1&gc[u]==1){
      for(k in 1:K){
        #mother
        hm1 = hap[k,]
        hm2 = gmm[u,] - hm1
        r = 0
        for (kk in 1:K){
          r = r + all(hm2==hap[kk,])
        }
        if(r==1){
          #pass hm1 to child
          hcm = hm1
          hcp = gcc[u,] - hcm
          t = 0
          for (kk in 1:K){
            t = t + all(hcp==hap[kk,])
          }
          if(t==1&gcm[u]!=hcm[loci]){
            gcm[u] = hcm[loci]
            gcp[u] = hcp[loci]
            if(phi[u]==0) phi[u] = 1
            else phi[u] = phi[u] + 1
          }
        }
      }
    }
    else if(gm[u]==2&gc[u]==1){
      gcm[u] = 1
      gcp[u] = 0
      phi[u] = 1
    }
    else{
      gcm[u] = 1
      gcp[u] = 1
      phi[u] = 1
    }
  }

  return(list(gm = gm, gc = gc, gcm = gcm, gcp = unlist(gcp), phi = phi))

}


distinguish0 = function(gmm,gcc,loci,hap,ppi){
  l = max(nrow(gmm),1)
  K = nrow(hap)
  gcm = rep(3,l)
  Pm = numeric(l)
  Pp = numeric(l)
  if(l>1){
    for(u in 1:l){
      for(k in 1:K){
        #mother
        hm1 = hap[k,]
        hm2 = gmm[u,] - hm1
        r = 0
        for (kk in 1:K){
          r = r + all(hm2==hap[kk,])
        }
        if(r==1){
          #pass hm1 to child
          hcm = hm1
          hcp = gcc[u,] - hcm
          t = 0
          for (kk in 1:K){
            t = t + all(hcp==hap[kk,])
          }
          if(t==1&gcm[u]!=hcm[loci]){
            gcm[u] = hcm[loci]
            if(gcm[u]==1) Pm[u] = Pm[u] + ppi[whichrow_new(hm1,hap)]*ppi[whichrow_new(hm2,hap)]*ppi[whichrow_new(hcp,hap)]
            else Pp[u] = Pp[u] + ppi[whichrow_new(hm1,hap)]*ppi[whichrow_new(hm2,hap)]*ppi[whichrow_new(hcp,hap)]

          }
        }
      }
    }
  } else{
    for(k in 1:K){
      #mother
      hm1 = hap[k,]
      hm2 = gmm - hm1
      r = 0
      for (kk in 1:K){
        r = r + all(hm2==hap[kk,])
      }
      if(r==1){
        #pass hm1 to child
        hcm = hm1
        hcp = gcc - hcm
        t = 0
        for (kk in 1:K){
          t = t + all(hcp==hap[kk,])
        }
        if(t==1&gcm!=hcm[loci]){
          gcm = hcm[loci]
          if(gcm==1) Pm = Pm + ppi[whichrow_new(hm1,hap)]*ppi[whichrow_new(hm2,hap)]*ppi[whichrow_new(hcp,hap)]
          else Pp = Pp + ppi[whichrow_new(hm1,hap)]*ppi[whichrow_new(hm2,hap)]*ppi[whichrow_new(hcp,hap)]

        }
      }
    }
  }
  return(list(Pm = Pm, Pp = Pp))
}

whichrow_new = function(hcm,hap){
  m = 0
  K = nrow(hap)
  for (kk in 1:K) m = c(m,all(hcm==hap[kk,]))
  return(which(m==1)-1)
}

#Pr(Gc|Gm;theta)

Im.tr.add = function(theta){
  res = matrix(0,3,3);
  res[1,1] = 1-theta;
  res[2,1] = theta;
  res[1,2] = (1-theta)/2;
  res[2,2] = 1/2;
  res[3,2] = theta/2;
  res[2,3] = 1-theta;
  res[3,3] = theta;
  return(res);
}

#d(Pr(Gc|Gm;theta))/d(theta)

Im.tr.gr.add = function(theta){
  res = matrix(0,3,3);
  res[1,1] = -1;
  res[2,1] = 1;
  res[1,2] = -1/2;
  res[3,2] = 1/2;
  res[2,3] = -1;
  res[3,3] = 1;
  return(res);
}


#double additive logistic regression model

#xi(theta)

Im.xi.add = function(theta){
  c((1-theta)^2, 2*theta*(1-theta), theta^2)
}

#d(xi(theta))/d(theta)

Im.xi.gr.add= function(theta){
  c(-2*(1-theta),2*(1-2*theta),2*theta)
}

#Pr(Gm|X;eta)

Im.dalog.add = function(Gm,X,theta,eta){
  # additive mode: Gm = 0 or 1 or 2
  Gm1 = Gm + 1
  Gm = Gm/2
  tmp = as.numeric(X*eta)
  etmp.2 = exp(tmp/2)
  xi = Im.xi.add(theta)
  tmp0 = xi[Gm1]*exp(Gm*tmp)
  tmp1 = xi[1]+xi[2]*etmp.2+xi[3]*etmp.2^2
  return(tmp0/tmp1)
}

#d(Pr(Gm|X;eta)/d(eta)

Im.dalog.gr.add = function(Gm,X,theta,eta){
  Gm1 = Gm + 1;
  Gm = Gm/2
  tmp = as.numeric(X*eta);
  xi = Im.xi.add(theta);
  xi.d = Im.xi.gr.add(theta);

  etmp.2 = exp(tmp/2);
  etmp.1 = etmp.2^2;
  eGtmp = exp(Gm*tmp);

  tmp0 = xi[1] + xi[2]*etmp.2 + xi[3]*etmp.1;
  tmp1 = xi.d[Gm1]*eGtmp;
  tmp2 = xi[Gm1]*eGtmp;

  tmp3 = tmp1/tmp0 - (xi.d[1]+xi.d[2]*etmp.2 + xi.d[3]*etmp.1) * tmp2/tmp0^2;
  tmp4 = X*(tmp2*Gm/tmp0 - tmp2/tmp0^2*(xi[2]*etmp.2/2 + xi[3]*etmp.1))

  return(cbind(tmp3,tmp4));
}


#Pr(Gpc)

gpc.add= function(theta){
  c(1-theta,theta)
}

#d(Pr(Gpc;theta))/d(theta)

gpc.gr.add= function(theta){
  c(-1,1)
}


#H_u = \sum_{Gm,Gmc,Gpc} Pr(Y=1|Gm,Gmc,Gpc,X_i)Pr(Gmc|Gm)Pr(Gpc)Pr(Gm|X_u)

Im.H.add = function(n.beta,Theta,X){

  theta = Theta[1]
  beta = Theta[2:(1+n.beta)]
  eta = Theta[-(1:(1+n.beta))]

  n = length(X)
  h = rep(0,n)

  #Pr(Gmc|Gm)
  gt = matrix(0,2,3)
  gt[1,1] = 1
  gt[1,2] = 1/2
  gt[2,2] = 1/2
  gt[2,3] = 1

  #Pr(Gpc)
  gp = gpc.add(theta)

  for(j in 1:3){
    for(k in 1:2){
      for (l in 1:2) {
        h = h + Im.penetrance0(j-1,k-1,l-1,X,beta)*Im.dalog.add(j-1,X,theta,eta)*gt[k,j]*gp[l]
      }
    }
  }
  return(h)
}

#derivative of H_u

Im.H.gr.add = function(n.beta,Theta,X){

  theta = Theta[1]
  beta = Theta[2:(1+n.beta)]
  eta = Theta[-(1:(1+n.beta))]

  n = length(X)
  p = length(Theta)

  #Pr(Gmc|Gm)
  gt = matrix(0,2,3)
  gt[1,1] = 1
  gt[1,2] = 1/2
  gt[2,2] = 1/2
  gt[2,3] = 1

  #Pr(Gpc)
  gp = gpc.add(theta)
  gp.r = gpc.gr.add(theta)

  h.gr = matrix(0,n,p)

  r.fn = cbind(Im.dalog.add(0,X,theta,eta),Im.dalog.add(1,X,theta,eta),Im.dalog.add(2,X,theta,eta));
  r.gr = list();
  r.gr[[1]] = Im.dalog.gr.add(0,X,theta,eta);
  r.gr[[2]] = Im.dalog.gr.add(1,X,theta,eta);
  r.gr[[3]] = Im.dalog.gr.add(2,X,theta,eta);

  for(j in 1:3){
    for(k in 1:2){
      for (l in 1:2) {

        p.fn = Im.penetrance0(j-1,k-1,l-1,X,beta)
        p.gr = Im.penetrance0.gr(j-1,k-1,l-1,X,beta)

        h.gr.theta = p.fn*gt[k,j]*(gp.r[l]*r.fn[,j] + gp[l]*r.gr[[j]][,1])

        h.gr.eta = r.gr[[j]][,-1]*p.fn*gt[k,j]*gp[l]

        tmp = p.gr*r.fn[,j]*gt[k,j]*gp[l]

        h.gr[,1] = h.gr[,1] + h.gr.theta
        h.gr[,2:(1+n.beta)] = h.gr[,2:(1+n.beta)] + sweep(cbind(1,j-1,k+l-2,k-l,X),1,tmp,'*')
        h.gr[,-(1:(1+n.beta))] = h.gr[,-(1:(1+n.beta))] + h.gr.eta
      }
    }
  }

  return(h.gr)

}



#logistic P(Y|Gm,Gcm,Gcp,X;beta)

Im.penetrance0 = function(Gm,Gcm,Gcp,X,beta){
  E = cbind(1,Gm,Gcm+Gcp,Gcm-Gcp,X)
  tmp = exp(as.vector(E%*%beta))

  return(tmp/(1+tmp))
}

#d(pr(Y=1|Gm,Gcm,Gcp,X;beta))/d(beta)

Im.penetrance0.gr = function(Gm,Gcm,Gcp,X,beta){
  E = cbind(1,Gm,Gcm+Gcp,Gcm-Gcp,X)
  tmp = exp(as.vector(E%*%beta))

  return(tmp/(1+tmp)^2)
}

#logistic P(Y|Gm,Gc,Gcm,Gcp,X;beta)

Im.penetrance1 = function(E,Y,beta){

  tmp = exp(as.vector(E%*%beta))

  return(1-Y+(2*Y-1)*tmp/(1+tmp))
}

#d(pr(Y|Gm,Gc,Gcm,Gcp,X;beta))/d(beta)

Im.penetrance1.gr = function(E,Y,beta){

  tmp = exp(as.vector(E%*%beta))

  return((2*Y-1)*tmp/(1+tmp)^2)
}


#design matrix 1

Im.design.matrix1 = function(gm,gc,gcm,gcp,X){

  cbind(1,gm,gc,gcm-gcp,X)
}

#design matrix 2

Im.design.matrix2_1 = function(gm,gc,gcm,gcp,X,phi){

  gcm[-which(phi==1)] = 1
  gcp[-which(phi==1)] = 0
  cbind(1,gm,gc,gcm-gcp,X)
}

Im.design.matrix2_2 = function(gm,gc,gcm,gcp,X,phi){

  gcm[-which(phi==1)] = 0
  gcp[-which(phi==1)] = 1
  cbind(1,gm,gc,gcm-gcp,X)
}


#main function

imprinting = function(Y,gmm,gcc,X,loci,hap,f,ppi){

  F = distinguish(gmm, gcc, loci, hap)
  gm = F$gm
  gc = F$gc
  gcm = F$gcm
  gcp = F$gcp
  phi = F$phi

  gm.log = gm[phi==1]
  gc.log = gc[phi==1]
  gcm.log = gcm[phi==1]
  gcp.log = gcp[phi==1]
  Y.log = Y[phi==1]
  X.log = X[phi==1]

  Z = Im.design.matrix1(gm,gc,gcm,gcp,X)
  if(mean(phi==1)==1){
    Z1 = Z
    Z2 = Z
    Pm = rep(1,length(Y))
    Pp = rep(0,length(Y))
  } else{
    Z1 = Im.design.matrix2_1(gm,gc,gcm,gcp,X,phi)
    Z2 = Im.design.matrix2_2(gm,gc,gcm,gcp,X,phi)
    W = distinguish0(gmm[-which(phi==1),],gcc[-which(phi==1),],loci,hap,ppi)
    Pm = rep(1,length(Y))
    Pp = rep(0,length(Y))
    Pm[-which(phi==1)] = W$Pm
    Pp[-which(phi==1)] = W$Pp
  }
  P1 = Pm/(Pm+Pp)
  P2 = Pp/(Pm+Pp)

  Z.log = Im.design.matrix1(gm.log,gc.log,gcm.log,gcp.log,X.log)
  fit = glm(Y.log ~ 0 + Z.log,family = binomial)
  res = summary(fit)$coef
  est.log = as.vector(res[,1])
  cov.log = vcov(fit)
  pval.log = as.vector(res[,4])

  theta = 1-sqrt(mean(gm==0))
  n.beta = ncol(Z) # number of all regression parameters in the penetrance model
  n.eta = max(1,ncol(X)) # number of covariates possibly associated with Gm
  theta0 = log(f/(1-f))
  beta = rep(theta0,n.beta)
  beta[-1] = est.log[-1]
  Theta = c(theta,beta,rep(0,n.eta))
  n = length(Y)
  ns = c(sum(Y==1),sum(Y==0))
  lambda = ns[1]/n/f - ns[2]/n/(1-f)
  p = length(Theta)

  II = list()
  II[[1]] = which(Y==0)
  II[[2]] = which(Y==1)

  upper = lower = Theta
  upper[1] = 0.99
  lower[1] = 0.01
  upper[-1] = Theta[-1] + 2
  lower[-1] = Theta[-1] - 2

  #number of (gc,gm)
  m11 = sum(gc==0&gm==0)
  m12 = sum(gc==0&gm==1)
  m21 = sum(gc==1&gm==0)
  m22 = sum(gc==1&gm==1)
  m23 = sum(gc==1&gm==2)
  m32 = sum(gc==2&gm==1)
  m33 = sum(gc==2&gm==2)

  #likelihood
  fn.add = function(Theta){
    theta=Theta[1]
    beta = Theta[2:(1+n.beta)]
    eta = Theta[-(1:(1+n.beta))]

    tr = Im.tr.add(theta)
    h.fn = Im.H.add(n.beta,Theta,X)

    res = m11*log(tr[1,1])+m12*log(tr[1,2])+m21*log(tr[2,1])+m22*log(tr[2,2])+m23*log(tr[2,3])+m32*log(tr[3,2])+m33*log(tr[3,3])
    res = res + sum(log(Im.dalog.add(gm,X,theta,eta)))
    res = res + sum(log(P1*Im.penetrance1(Z1,Y,beta)+P2*Im.penetrance1(Z2,Y,beta)))
    res = res - sum(log(1+lambda*(h.fn-f)))
    return(-res);
  }

  #gradient
  gr.add = function(Theta){
    theta=Theta[1]
    beta = Theta[2:(1+n.beta)]
    eta = Theta[-(1:(1+n.beta))]

    logq.gr = Im.tr.gr.add(theta)/Im.tr.add(theta)

    logq.theta = rep(NA,n)
    logq.theta[gc==0&gm==0] = logq.gr[1,1]
    logq.theta[gc==1&gm==0] = logq.gr[2,1]
    logq.theta[gc==0&gm==1] = logq.gr[1,2]
    logq.theta[gc==1&gm==1] = logq.gr[2,2]
    logq.theta[gc==2&gm==1] = logq.gr[3,2]
    logq.theta[gc==1&gm==2] = logq.gr[2,3]
    logq.theta[gc==2&gm==2] = logq.gr[3,3]


    logr.gr = sweep(Im.dalog.gr.add(gm,X,theta,eta),1,Im.dalog.add(gm,X,theta,eta),'/')

    logr.theta = logr.gr[,1]
    logr.eta = logr.gr[,-1]

    logp.gr1 = P1*Im.penetrance1.gr(Z1,Y,beta)/(P1*Im.penetrance1(Z1,Y,beta)+P2*Im.penetrance1(Z2,Y,beta))
    logp.gr2 = P2*Im.penetrance1.gr(Z2,Y,beta)/(P1*Im.penetrance1(Z1,Y,beta)+P2*Im.penetrance1(Z2,Y,beta))
    logp.beta = sweep(Z1,1,logp.gr1,'*') + sweep(Z2,1,logp.gr2,'*')

    h.fn = Im.H.add(n.beta,Theta,X)
    h.gr = Im.H.gr.add(n.beta,Theta,X)
    others = sweep(lambda*h.gr,1,1+lambda*(h.fn-f),'/')

    all = cbind(logq.theta+logr.theta,logp.beta,logr.eta) - others
    -colSums(all)
  }


  fit = nloptr(x0=Theta,eval_f=fn.add,eval_grad_f=gr.add,lb=lower,ub=upper,opts=list("algorithm"="NLOPT_LD_LBFGS","xtol_rel"=1.0e-4) )
  hess = hessian(func=fn.add,x=fit$solution)

  est = fit$solution

  #covariance
  score.cov.add = function(Theta){

    theta=Theta[1]
    beta = Theta[2:(1+n.beta)]
    eta = Theta[-(1:(1+n.beta))]

    logq.gr = Im.tr.gr.add(theta)/Im.tr.add(theta)

    logq.theta = rep(NA,n)
    logq.theta[gc==0&gm==0] = logq.gr[1,1]
    logq.theta[gc==1&gm==0] = logq.gr[2,1]
    logq.theta[gc==0&gm==1] = logq.gr[1,2]
    logq.theta[gc==1&gm==1] = logq.gr[2,2]
    logq.theta[gc==2&gm==1] = logq.gr[3,2]
    logq.theta[gc==1&gm==2] = logq.gr[2,3]
    logq.theta[gc==2&gm==2] = logq.gr[3,3]


    logr.gr = sweep(Im.dalog.gr.add(gm,X,theta,eta),1,Im.dalog.add(gm,X,theta,eta),'/')

    logr.theta = logr.gr[,1]
    logr.eta = logr.gr[,-1]

    logp.gr1 = P1*Im.penetrance1.gr(Z1,Y,beta)/(P1*Im.penetrance1(Z1,Y,beta)+P2*Im.penetrance1(Z2,Y,beta))
    logp.gr2 = P2*Im.penetrance1.gr(Z2,Y,beta)/(P1*Im.penetrance1(Z1,Y,beta)+P2*Im.penetrance1(Z2,Y,beta))
    logp.beta = sweep(Z1,1,logp.gr1,'*') + sweep(Z2,1,logp.gr2,'*')

    h.fn = Im.H.add(n.beta,Theta,X)
    h.gr = Im.H.gr.add(n.beta,Theta,X)
    others = sweep(lambda*h.gr,1,1+lambda*(h.fn-f),'/')

    all = cbind(logq.theta+logr.theta,logp.beta,logr.eta) - others

    res = matrix(0,p,p)

    for(k in 1:2){
      tmp = all[II[[k]],]
      tt = as.matrix(scale(tmp,center=TRUE,scale=FALSE))*ns[k]/(ns[k]-1)
      res = res + t(tt)%*%tt;
    }

    return(res);
  }

  V = hess
  U = score.cov.add(est);

  V.inv = solve(V);
  cov = (V.inv%*%U)%*%V.inv;

  return(list(est =est,cov =cov,est.log=est.log,cov.log=cov.log,pval.log=pval.log))

}



