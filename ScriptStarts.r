##########################################
## Script for Quadvariate Starts Analysis
##  in paper 
## "Model Identification via Total Frobenius
##   Norm of Multivariate Spectra"
##########################################

## wipe
rm(list=ls())

library(devtools)
library(xtable)

# suppose directory is set to where sigex is located, e.g.
#setwd("H:\\SigEx\\sigex-master")
load_all(".")

setwd("H:\\FrobNorm")

##########################
### Part I: Load Data 

# automatic with the sigex package 


####################################
### Part II: Metadata Specifications  

start.date = c(1964,1)
period <- 12

## create ts object and plot
dataALL.ts <- sigex.load(starts,start.date,period,
                         c("South","West","NE","MW"),TRUE)

## all data with no transform
transform <- "none"
aggregate <- FALSE
subseries <- c(1,2,3,4)
begin.date <- c(2004,1)
end.date <- end(dataALL.ts)
range <- list(begin.date,end.date)
data.ts <- sigex.prep(dataALL.ts,transform,aggregate,subseries,range,TRUE)
 
 
###############################
### Part III: Model Declaration

N <- dim(data.ts)[2]
T <- dim(data.ts)[1]

## model construction
mdl <- NULL
mdl <- sigex.add(mdl,seq(1,N),"arma",c(0,0),NULL,"trend",c(1,-2,1))
mdl <- sigex.add(mdl,seq(1,N),"arma",c(0,0),NULL,"first seasonal",c(1,-sqrt(3),1))
mdl <- sigex.add(mdl,seq(1,N),"arma",c(0,0),NULL,"second seasonal",c(1,-1,1))
mdl <- sigex.add(mdl,seq(1,N),"arma",c(0,0),NULL,"third seasonal",c(1,0,1))
mdl <- sigex.add(mdl,seq(1,N),"arma",c(0,0),NULL,"fourth seasonal",c(1,1,1))
mdl <- sigex.add(mdl,seq(1,N),"arma",c(0,0),NULL,"fifth seasonal",c(1,sqrt(3),1))
mdl <- sigex.add(mdl,seq(1,N),"arma",c(0,0),NULL,"sixth seasonal",c(1,1))
mdl <- sigex.add(mdl,seq(1,N),"arma",c(0,0),NULL,"irregular",1)
mdl <- sigex.meaninit(mdl,data.ts,0)


#############################
### Part IV: Model Estimation

##############
## MOM fitting

mdl.mom <- mdl
constraint <- NULL
par.default <- sigex.default(mdl.mom,data.ts,constraint)
par.mom <- sigex.momfit(data.ts,par.default,mdl.mom)
psi.mom <- sigex.par2psi(par.mom,mdl.mom)
sigex.lik(psi.mom,mdl.mom,data.ts,FALSE)
# divergence  959.806

## residual analysis
resid.mom <- sigex.resid(psi.mom,mdl.mom,data.ts)[[1]]
resid.mom <- sigex.load(t(resid.mom),start(data.ts),
                        frequency(data.ts),colnames(data.ts),TRUE)
resid.acf <- acf(resid.mom,lag.max=4*period,plot=TRUE)$acf

## examine condition numbers
log(sigex.conditions(data.ts,psi.mom,mdl.mom))

## model checking
sigex.portmanteau(resid.mom,4*period,length(psi.mom))

## bundle for reduced span
analysis.mom <- sigex.bundle(data.ts,transform,mdl.mom,psi.mom)


##############
## MLE Fitting 

## initialize with default values
constraint <- NULL
par.mle <- sigex.default(mdl,data.ts,constraint)
psi.mle <- sigex.par2psi(par.mle,mdl)

## run fitting:   this takes a long time!
fit.mle <- sigex.mlefit(data.ts,par.mle,constraint,mdl,"bfgs",debug=TRUE)
# divergence 948.6092

## manage output
psi.mle <- sigex.eta2psi(fit.mle[[1]]$par,constraint)
par.mle <- fit.mle[[2]]

## another round!
fit.mle <- sigex.mlefit(data.ts,par.mle,constraint,mdl,"bfgs",debug=TRUE)
# divergence 923.4877

## manage output
psi.mle <- sigex.eta2psi(fit.mle[[1]]$par,constraint)
par.mle <- fit.mle[[2]]

## residual analysis
resid.mle <- sigex.resid(psi.mle,mdl,data.ts)[[1]]
resid.mle <- sigex.load(t(resid.mle),start(data.ts),
                        frequency(data.ts),colnames(data.ts),TRUE)
resid.acf <- acf(resid.mle,lag.max=4*period,plot=TRUE)$acf

## examine condition numbers
log(sigex.conditions(data.ts,psi.mle,mdl))

## model checking
sigex.portmanteau(resid.mle,4*period,length(psi.mle))

## bundle for reduced span
analysis.mle <- sigex.bundle(data.ts,transform,mdl,psi.mle)


##################################
### Part V: MOM fit on entire span

## all data with no transform
transform <- "none"
aggregate <- FALSE
subseries <- c(1,2,3,4)
begin.date <- start.date
end.date <- end(dataALL.ts)
range <- list(begin.date,end.date)
data.ts <- sigex.prep(dataALL.ts,transform,aggregate,subseries,range,TRUE)

N <- dim(data.ts)[2]
T <- dim(data.ts)[1]

## model construction
mdl <- NULL
mdl <- sigex.add(mdl,seq(1,N),"arma",c(0,0),NULL,"trend",c(1,-2,1))
mdl <- sigex.add(mdl,seq(1,N),"arma",c(0,0),NULL,"first seasonal",c(1,-sqrt(3),1))
mdl <- sigex.add(mdl,seq(1,N),"arma",c(0,0),NULL,"second seasonal",c(1,-1,1))
mdl <- sigex.add(mdl,seq(1,N),"arma",c(0,0),NULL,"third seasonal",c(1,0,1))
mdl <- sigex.add(mdl,seq(1,N),"arma",c(0,0),NULL,"fourth seasonal",c(1,1,1))
mdl <- sigex.add(mdl,seq(1,N),"arma",c(0,0),NULL,"fifth seasonal",c(1,sqrt(3),1))
mdl <- sigex.add(mdl,seq(1,N),"arma",c(0,0),NULL,"sixth seasonal",c(1,1))
mdl <- sigex.add(mdl,seq(1,N),"arma",c(0,0),NULL,"irregular",1)
mdl <- sigex.meaninit(mdl,data.ts,0)

mdl.mom <- mdl
constraint <- NULL
par.default <- sigex.default(mdl.mom,data.ts,constraint)
par.mom2 <- sigex.momfit(data.ts,par.default,mdl.mom)
psi.mom2 <- sigex.par2psi(par.mom2,mdl.mom)
sigex.lik(psi.mom2,mdl.mom,data.ts,FALSE)
# divergence   6329.107

## residual analysis
resid.mom2 <- sigex.resid(psi.mom2,mdl.mom,data.ts)[[1]]
resid.mom2 <- sigex.load(t(resid.mom2),start(data.ts),
                        frequency(data.ts),colnames(data.ts),TRUE)
resid.acf2 <- acf(resid.mom2,lag.max=4*period,plot=TRUE)$acf

## examine condition numbers
log(sigex.conditions(data.ts,psi.mom2,mdl.mom))

## model checking
sigex.portmanteau(resid.mom2,4*period,length(psi.mom2))

## bundle for full span
analysis.mom2 <- sigex.bundle(data.ts,transform,mdl.mom,psi.mom2)


######################################
### Part VI: Print Covariance Matrices

## print covariance matrices
for(i in 1:8)
{
  dat <- par.mle[[1]][[i]] %*% diag(exp(par.mle[[2]][[i]])) %*% t(par.mle[[1]][[i]])                                                                  
  dat <- cbind(dat,par.mom[[1]][[i]] %*% diag(exp(par.mom[[2]][[i]])) %*% t(par.mom[[1]][[i]]))
  dat <- cbind(dat,par.mom2[[1]][[i]] %*% diag(exp(par.mom2[[2]][[i]])) %*% t(par.mom2[[1]][[i]]))               
  rownames(dat) <- c("South","West","NE","MW")
  colnames(dat) <- c("South","West","NE","MW","South","West","NE","MW","South","West","NE","MW")
  tab <- xtable(dat,digits=3)
  print(tab)
}


