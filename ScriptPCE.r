##########################################
## Script for Bivariate Inflation Analysis
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


##########################
### Part I: Load Data 

setwd("H:\\FrobNorm")
pce.core <- read.table("PCECore86.txt")
pce.total <- read.table("PCETotal86.txt")
pce <- cbind(pce.core,pce.total) 
 

####################################
### Part II: Metadata Specifications  

start.date <- c(1986,1)
end.date <- c(2010,4)
period <- 4

## create ts object and plot
dataALL.ts <- sigex.load(pce,start.date,period,
                         c("Core","Total"),TRUE)

## all data with no transform
transform <- "none"
aggregate <- FALSE
subseries <- c(1,2)
begin.date <- start.date
end.date <- end(dataALL.ts)
range <- list(begin.date,end.date)
data.ts <- sigex.prep(dataALL.ts,transform,aggregate,subseries,range,TRUE)


###############################
### Part III: Model Declaration

N <- dim(data.ts)[2]
T <- dim(data.ts)[1]

## model construction: unrestricted
mdl1 <- NULL
mdl1 <- sigex.add(mdl1,seq(1,N),"arma",c(0,0),NULL,"trend",c(1,-1))
mdl1 <- sigex.add(mdl1,seq(1,N),"arma",c(0,0),NULL,"irregular",1)
mdl1 <- sigex.meaninit(mdl1,data.ts,0)
 
## model construction: common trends
mdl2 <- NULL
mdl2 <- sigex.add(mdl2,1,"arma",c(0,0),NULL,"trend",c(1,-1))
mdl2 <- sigex.add(mdl2,seq(1,N),"arma",c(0,0),NULL,"irregular",1)
mdl2 <- sigex.meaninit(mdl2,data.ts,0)


#############################
### Part IV: Model Estimation

##############
## MOM fitting

mdl.mom <- mdl1
constraint <- NULL
par.default <- sigex.default(mdl.mom,data.ts,constraint)
par.mom <- sigex.momfit(data.ts,par.default,mdl.mom)
psi.mom <- sigex.par2psi(par.mom,mdl.mom)
sigex.lik(psi.mom,mdl.mom,data.ts,FALSE)
# divergence  -1680.292

## residual analysis
resid.mom <- sigex.resid(psi.mom,mdl.mom,data.ts)[[1]]
resid.mom <- sigex.load(t(resid.mom),start(data.ts),
                        frequency(data.ts),colnames(data.ts),TRUE)
resid.acf <- acf(resid.mom,lag.max=4*period,plot=TRUE)$acf

## examine condition numbers
log(sigex.conditions(data.ts,psi.mom,mdl.mom))

## model checking
sigex.portmanteau(resid.mom,4*period,length(psi.mom))

## bundle for default span
analysis.mom <- sigex.bundle(data.ts,transform,mdl.mom,psi.mom)

# correlations
corr.trend <- par.mom[[1]][[1]][2,1]/sqrt(par.mom[[1]][[1]][2,1]^2+exp(par.mom[[2]][[1]][2] - par.mom[[2]][[1]][1]))
corr.irr <- par.mom[[1]][[2]][2,1]/sqrt(par.mom[[1]][[2]][2,1]^2+exp(par.mom[[2]][[2]][2] - par.mom[[2]][[2]][1]))

## print covariance matrices
xtable(10^6*(par.mom[[1]][[1]] %*% diag(exp(par.mom[[2]][[1]])) %*% t(par.mom[[1]][[1]])),digits=3)	# trend
xtable(10^6*(par.mom[[1]][[2]] %*% diag(exp(par.mom[[2]][[2]])) %*% t(par.mom[[1]][[2]])),digits=3)	# irreg

############################
## MLE Fitting: unrestricted

## initialize with MOM results
psi.mom <- analysis.mom[[4]]
par.mom <- sigex.psi2par(psi.mom,mdl1,data.ts)

#  Initialize with MOM estimates
constraint <- NULL
psi.mle <- sigex.par2psi(par.mom,mdl1)

## run fitting:   this takes a few minutes
fit.mle <- sigex.mlefit(data.ts,par.mom,constraint,mdl1,"bfgs",debug=TRUE)
# divergence -1695.502

## manage output
psi.mle <- sigex.eta2psi(fit.mle[[1]]$par,constraint)
par.mle <- fit.mle[[2]]

## another round!
fit.mle <- sigex.mlefit(data.ts,par.mle,constraint,mdl1,"bfgs",debug=TRUE)
# divergence -1702.763

## manage output
psi.mle <- sigex.eta2psi(fit.mle[[1]]$par,constraint)
par.mle <- fit.mle[[2]]

## residual analysis
resid.mle <- sigex.resid(psi.mle,mdl1,data.ts)[[1]]
resid.mle <- sigex.load(t(resid.mle),start(data.ts),
                        frequency(data.ts),colnames(data.ts),TRUE)
resid.acf <- acf(resid.mle,lag.max=4*period,plot=TRUE)$acf

## examine condition numbers
log(sigex.conditions(data.ts,psi.mle,mdl1))

## model checking
sigex.portmanteau(resid.mle,4*period,length(psi.mle))

# bundle for default span
analysis.mle <- sigex.bundle(data.ts,transform,mdl1,psi.mle)

# correlations
corr.trend <- par.mle[[1]][[1]][2,1]/sqrt(par.mle[[1]][[1]][2,1]^2+exp(par.mle[[2]][[1]][2] - par.mle[[2]][[1]][1]))
corr.irr <- par.mle[[1]][[2]][2,1]/sqrt(par.mle[[1]][[2]][2,1]^2+exp(par.mle[[2]][[2]][2] - par.mle[[2]][[2]][1]))

# print covariance matrices
xtable(10^6*(par.mle[[1]][[1]] %*% diag(exp(par.mle[[2]][[1]])) %*% t(par.mle[[1]][[1]])),digits=3)	# trend
xtable(10^6*(par.mle[[1]][[2]] %*% diag(exp(par.mle[[2]][[2]])) %*% t(par.mle[[1]][[2]])),digits=3)	# irreg


#############################
## MLE Fitting: common trends

## compute reduced rank trend covariance matrix
cov.mat <- par.mle[[1]][[1]] %*% diag(exp(par.mle[[2]][[1]])) %*%
  t(par.mle[[1]][[1]])
l.trend <- sqrt(cov.mat[2,2]/cov.mat[1,1])
d.trend <- cov.mat[1,1]

## parameter initialization based on prior info
constraint <- NULL
psi.mle2 <- c(l.trend,log(d.trend),psi.mle[4:8])
par.mle2 <- sigex.psi2par(psi.mle2,mdl2,data.ts)

## run fitting:   this takes a few minutes
fit.mle2 <- sigex.mlefit(data.ts,par.mle2,constraint,mdl2,"bfgs",debug=TRUE)
# divergence -1708.075

## manage output
psi.mle2 <- sigex.eta2psi(fit.mle2[[1]]$par,constraint)
par.mle2 <- fit.mle2[[2]]

## residual analysis
resid.mle2 <- sigex.resid(psi.mle2,mdl2,data.ts)[[1]]
resid.mle2 <- sigex.load(t(resid.mle2),start(data.ts),
                        frequency(data.ts),colnames(data.ts),TRUE)
resid.acf2 <- acf(resid.mle2,lag.max=4*period,plot=TRUE)$acf

## examine condition numbers
log(sigex.conditions(data.ts,psi.mle2,mdl2))

## model checking
sigex.portmanteau(resid.mle2,4*period,length(psi.mle2))

# bundle for default span
analysis.mle2 <- sigex.bundle(data.ts,transform,mdl2,psi.mle2)

# correlations
corr.irr <- par.mle2[[1]][[2]][2,1]/sqrt(par.mle2[[1]][[2]][2,1]^2+exp(par.mle2[[2]][[2]][2] - par.mle2[[2]][[2]][1]))

# print covariance matrices
xtable(10^6*(par.mle2[[1]][[1]] %*% exp(par.mle2[[2]][[1]]) %*% t(par.mle2[[1]][[1]])),digits=3)	# trend
xtable(10^6*(par.mle2[[1]][[2]] %*% diag(exp(par.mle2[[2]][[2]])) %*% t(par.mle2[[1]][[2]])),digits=3)	# irreg


####################################
### Part V: Tests for co-integration

compute_variance <- function(sigmu, sigit){
  
  G <- matrix(c(1,2,2,6),2,2)
  G11 <- G
  G12 <- matrix(c(2,6,6,20),2,2)
  G21 <- G12
  G22 <- matrix(c(6,20,20,70),nrow=2,ncol=2)
  G11 <- solve(G) %*% G11 %*% solve(G)
  G12 <- solve(G) %*% G12 %*% solve(G)
  G21 <- solve(G) %*% G21 %*% solve(G)
  G22 <- solve(G) %*% G22 %*% solve(G)
  var <- G11 %x% sigmu %x% sigmu + 
    G12 %x% sigmu %x% sigit + 
    G21 %x% sigit %x% sigmu +  
    G22 %x% sigit %x% sigit 
  return(2*var)
}

## compute test statistic for rank 1

mult <- diag(c(1,.5,.5,1))
mult[2,3] <- .5
mult[3,2] <- .5
mult <- diag(2) %x% mult

sigmuhat <- 10^6*(par.mom[[1]][[1]] %*% diag(exp(par.mom[[2]][[1]])) %*% t(par.mom[[1]][[1]]))
sigithat <- 10^6*(par.mom[[1]][[2]] %*% diag(exp(par.mom[[2]][[2]])) %*% t(par.mom[[1]][[2]]))
  
vhat <- compute_variance(sigmuhat,sigithat)
vhat <- mult %*% vhat

Det <- det(sigmuhat)
grad <- matrix(c(sigmuhat[2,2], -2*sigmuhat[2,1], sigmuhat[1,1]),1,3)
vdet <- grad %*% vhat[c(1,2,4), c(1,2,4)] %*% t(grad)/T
tdet <- Det/sqrt(vdet)



