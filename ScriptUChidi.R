### EFFICIENCY OF MOM ESTIMATOR
### SETTINGS ARE SIMILAR TO THAT OF FIGURE 1


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
      if(N ==1){sa1[,t] = -sum(sa1[,(t-1):(t-D+1)]) + 
        sigsa12%*%rnorm(N,0,1) # normal
      } else {
        sa1[,t] = -apply(sa1[,(t-1):(t-D+1)],1,sum) + 
          sigsa12%*%rnorm(N,0,1) # normal
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
      if(N ==1){sa1[,t] = -sum(sa1[,(t-1):(t-D+1)]) + 
        sigsa12%*%rt(N,df=tdf)/tsd # student's t
      } else {
        sa1[,t] = -apply(sa1[,(t-1):(t-D+1)],1,sum) + 
          sigsa12%*%rt(N,df=tdf)/tsd # student's t
      }
    }

  }
  sa = sa1[,1001:(T+1000)]

  data = mu + sa + it
  data = t(data)
  return(data)
}


##### MAIN SIMULATION ######



dm = 25  # max limit of dimension
rep1 = 4 # number of cases for each dimension
rep = 5000 # Monte Carlo replications
nsamp = c(100,1000) # Sample sizes
D = 4  # seasonal difference
options(warn=-1)
out = array(0, dim = c(rep1*(dm-1),9,2))

#########################################
#  Generating different parameter configuration
#  For each dimension random variance components 
#  with random rank reduction
#########################################

for(aa in 2:dm){
  for(bb in 1:rep1){
    case = rep1*(aa-2) + bb
    N = aa
    rnkmu = sample(N,1)-1  ## Rank reduction in sigmu
    rnksa = sample(N,1)-1  ## Rank reduction in Sigsa
    effparmmu  = (N-rnkmu)*(N - ((N-rnkmu)-1)/2)
    effparmsa  = (N - rnksa)*(N - ((N- rnksa)-1)/2)
  
if(rnkmu ==0){ sigmu12 = matrix(rnorm(N^2,sd = .1),N,N)} else {
sigmu12 = cbind(matrix(rnorm(N*(N - rnkmu), 
                             sd = .1),N,(N - rnkmu)), matrix(0,N,rnkmu))}
if(rnksa ==0){ sigsa12 = matrix(rnorm(N^2, sd = .1),N,N)} else {
sigsa12 = cbind(matrix(rnorm(N*(N - rnksa), 
                             sd = .1),N,(N - rnksa)), matrix(0,N,rnksa))}
sigit12 = matrix(diag(N),N,N)
  

      for(nn in 1 :length(nsamp)){
  
        T = nsamp[nn]
        est1 = matrix(0, rep, N*N)
        est2 = matrix(0, rep, N*N)
        est3 = matrix(0, rep, N*N)
        
        ## model construction
        mdl <- NULL
        mdl <- sigex.add(mdl,seq(1,N),"arma",c(0,0),0,"trend",c(1,-1))   
        mdl <- sigex.add(mdl,seq(1,N),"arma",c(0,0),0,"seasonal",rep(1,D))    
        mdl <- sigex.add(mdl,seq(1,N),"arma",c(0,0),0,"irregular",1)
        constraint <- NULL
  
          for (rr in 1:rep){
            
            dist = 1   ## for normal  
            #dist = 2   ## for students t with tdf=4
            data = gen_UCdata(N,T,D,sigit12, sigmu12, sigsa12, dist)
            
            
            mdl <- sigex.meaninit(mdl,data,0)  ### Needs data to be TxN ###
            par.default <- sigex.default(mdl,data,constraint)
            psi.default <- sigex.par2psi(par.default,mdl)
            
            mdl.mom <- mdl
            
            ######################################
            # MOM estimation and reduced specification
            
            
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
            
            est1[rr,] = c(sigmuhat) - c(sigmu12%*%t(sigmu12))
            est2[rr,] = c(sigsahat) - c(sigsa12%*%t(sigsa12))
            est3[rr,] = c(sigithat) - c(sigit12%*%t(sigit12))
          }

    out[case, ,nn] = c(N, rnkmu, rnksa, effparmmu, effparmsa, N*(N+1)/2, 
                      c(sqrt(mean(est1^2)), sqrt(mean(est2^2)), 
                        sqrt(mean(est3^2)))*sqrt(T))    
    }
  }
}

