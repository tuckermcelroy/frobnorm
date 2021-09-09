### POWER AND SIZE SIMULATION FOR COINTEGRATION RANK TEST
### SETTINGS ARE SIMILAR TO THAT OF TABLE 1


##### LOAD ECCE SIGNUM ######
setwd("R:\\DATA\\SHARE\\TimeSeries\\Frobnorm\\sigex-master\\sigex-master\\")
load_all(".")

##### GENERATE UC MODEL DATA ######

gen_UCdata<- function(N,T,D,sigit12, sigmu12, sigsa12, dist){
  mu = matrix(0,N,T)
  sa = mu
  it = mu
  sa1= matrix(0,N,(T+1000))

  if(dist == 1){
    for(t in 1: T){
      it[,t] = sigit12%*%rnorm(N,0,1)
    }
    mu[,1] = sigmu12%*%rnorm(N,0,1) ## Normal
    for(t in 2: T){
      mu[,t] = mu[,(t-1)] + sigmu12%*%rnorm(N,0,1) # Normal
    }

    for(t in 1:D){
      sa1[,t] =  sigsa12%*%rnorm(N,0,1) # normal
    }

    for(t in (D+1):(T+1000)){
      if(N ==1){sa1[,t] = -sum(sa1[,(t-1):(t-D+1)])+ sigsa12%*%rnorm(N,0,1) # normal
      } else {
        sa1[,t] = -apply(sa1[,(t-1):(t-D+1)],1,sum) + sigsa12%*%rnorm(N,0,1) # normal
      }
    }
  }else{
    tdf = 4
    tsd = sqrt(tdf/(tdf-2))
    for(t in 1: T){
      #it[,t] = sigit12%*%rt(N,df=tdf)/tsd
       it[,t] = sigit12%*%rnorm(N,0,1)
    }
    mu[,1] = sigmu12%*%rt(N,df=tdf)/tsd ## Student's t
    for(t in 2: T){
      mu[,t] = mu[,(t-1)] + sigmu12%*%rt(N,df=tdf)/tsd  # Student's t
    }

    for(t in 1:D){
      sa1[,t] =  sigsa12%*%rt(N,df=tdf)/tsd #students t
    }

    for(t in (D+1):(T+1000)){
      if(N ==1){sa1[,t] = -sum(sa1[,(t-1):(t-D+1)])+ sigsa12%*%rt(N,df=tdf)/tsd # student's t
      } else {
        sa1[,t] = -apply(sa1[,(t-1):(t-D+1)],1,sum) + sigsa12%*%rt(N,df=tdf)/tsd # student's t
      }
    }

  }
  sa = sa1[,1001:(T+1000)]

  data = mu + sa + it
  data = t(data)
  return(data)
}


########  RANK TEST SIMULATION #######

  
rep = 5000  # Monte Carlo Replications
nsamp = c(200,500)  # Sample Sizes
D = 4  ####  seasonal difference
options(warn=-1)

# Parameters for Variance Components
smu1 = 1
smu2 = .8
rho1 = c(0,.6,.9,.95,1)


ssa1 = 1
ssa2 = .6
rho2 = c(0,.4,.8,.9,1)

powermu<-powersa<-array(0,dim=c(length(rho1),length(rho2),length(nsamp)))

mult = diag(c(1,.5,.5,1))
mult[2,3] = .5
mult[3,2] = .5
mult = diag(3)%x%mult

#######################

f1 <- function(x) 16*cos(x)^2*cos(x/2)^2
f2 <- function(x) 2*(1 - cos(x))
f3 <- function(x) 2*(1 - cos(4*x))
f = c(f1,f2,f3)
G = matrix(0, 3,3)
G1 = array(0,dim=c(3,3,3,3))
for(i in 1:3){for(j in 1:3){
  G[i,j] = (1/(2*pi))*integrate(function(x) f[[i]](x)*f[[j]](x), -pi, pi)$value
}}
Ginv = solve(G)
for(i in 1:3){for(j in 1:3){
  for(l in 1:3){ for(k in 1:3){
    G1[i,j,l,k] <- (1/(2*pi))*integrate(function(x) f[[i]](x)*f[[j]](x)*f[[l]](x)*f[[k]](x), -pi, pi)$value
  }}
}}
for(l in 1:3){ for(k in 1:3){
  G1[,,l,k] = Ginv%*%G1[,,l,k]%*%Ginv
}}




for(nn in 1 :length(nsamp)){
  
  T = nsamp[nn]
  Detmu<-vdetmu<-tdetmu<-Detsa<-vdetsa<-tdetsa<-matrix(0,rep,1)

  ## model construction
  mdl <- NULL
  mdl <- sigex.add(mdl,seq(1,N),"arma",c(0,0),0,"trend",c(1,-1))   
  mdl <- sigex.add(mdl,seq(1,N),"arma",c(0,0),0,"seasonal",rep(1,D))    
  mdl <- sigex.add(mdl,seq(1,N),"arma",c(0,0),0,"irregular",1)
  constraint <- NULL
  
  for(tt in 1:length(rho1)){
    for(ss in 1:length(rho2)){
        rhomu = rho1[tt]
        rhosa = rho1[ss]
    
        sigmu = matrix(c(smu1^2, rhomu*smu1*smu2, rhomu*smu1*smu2, smu2^2),N,N)
        sigsa = matrix(c(ssa1^2, rhosa*ssa1*ssa2, rhosa*ssa1*ssa2, ssa2^2),N,N)
        sigit = diag(N)
        sigmu12 = sqrtm(sigmu)
        sigsa12 = sqrtm(sigsa)
        sigit12 = sqrtm(sigit)

    
          for (rr in 1:rep){
            
            dist = 1   ## for normal  
            #dist = 2   ## for students t with tdf=4
            data = gen_UCdata(N,T,D,sigit12, sigmu12, sigsa12, dist)
          
            mdl <- sigex.meaninit(mdl,data,0)  ### Needs data to be TxN ###
            par.default <- sigex.default(mdl,data,constraint)
            psi.default <- sigex.par2psi(par.default,mdl)
            mdl.mom <- mdl
            
            par.mom <- sigex.momfit(data,par.default,mdl.mom)
            
            l = par.mom[[1]][[1]]
            d = diag(exp(par.mom[[2]][[1]]))
            sigmuhat = l%*%d%*%t(l)
            l = par.mom[[1]][[2]]
            d = diag(exp(par.mom[[2]][[2]]))
            sigsahat = l%*%d%*%t(l)
            l = par.mom[[1]][[3]]
            d = diag(exp(par.mom[[2]][[3]]))
            sigithat = l%*%d%*%t(l)
          
            sig = array(0,dim = c(N,N,3))
            sig[,,1] = sigmuhat
            sig[,,2] = sigsahat
            sig[,,3] = sigithat
      
            vhat = array(0, dim=c(12,12,3,3))
            for(l in 1:3){for(k in 1:3){ vhat[,,l,k] = 
                            G1[,,l,k]%x%sig[,,l]%x%sig[,,k] }}
            vhat = 2*mult%*%apply(vhat,c(1,2),sum)
      
            Detmu[rr] = det(sigmuhat)
            gradmu = matrix(c(sigmuhat[2,2], -2*sigmuhat[2,1], sigmuhat[1,1]),1,3)
            vdetmu[rr] = gradmu%*%vhat[c(1,2,4), c(1,2,4)]%*%t(gradmu)/T
            tdetmu[rr] = Detmu[rr]/sqrt(vdetmu[rr])
      
            Detsa[rr] = det(sigsahat)
            gradsa = matrix(c(sigsahat[2,2], -2*sigsahat[2,1], sigsahat[1,1]),1,3)
            vdetsa[rr] = gradsa%*%vhat[c(5,6,8), c(5,6,8)]%*%t(gradsa)/T
            tdetsa[rr] = Detsa[rr]/sqrt(vdetsa[rr])
          }
        powermu[tt,ss,nn] = mean((tdetmu > qnorm(.95)),na.rm=T)
        powersa[tt,ss,nn] = mean((tdetsa > qnorm(.95)),na.rm=T)
    }
  }
}


powermu = matrix(powermu,length(rho1),length(rho2)*length(nsamp))
powersa = matrix(powersa,length(rho1),length(rho2)*length(nsamp))

