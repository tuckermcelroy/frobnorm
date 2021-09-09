## Generate data from VAR(2) model in Section 5 with either Gaussian or t errors
## Check size and power of Model Evaluation Statistic using Corollary 1



library(Matrix)
library(mvtnorm)
library(expm)
library(xtable)

sqrtm <- function(A){
return(eigen(A)$vectors%*%diag(sqrt(eigen(A)$value))%*%t(eigen(A)$vectors))
}

## Simulate VAR(p)
varp_sim<-function(T,N,p,matphi,sigma2,gp2){
e = matrix(rnorm(N*T),N,T)  ## Chnage to  rt for students-t
y = matrix(0,N,T)
y[,1:p] = matrix(gp2%*%c(e[,1:p]),N,p)
for(i in (p+1):T){
y[,i] = matphi%*%matrix(y[,(i-1):(i-p)],N*p,1)  + sigma2%*%e[,i]
}
return(y)
}


N = 2  ## DIMENSION
p = 2  ## TRUE VAR ORDER
sigma = diag(N)
sigma2 = sqrtm(sigma)
phi =  array(c(.3,0,-.3,.4,.01,-.1,-.1,.25),dim=c(2,2,2))
matphi = matrix(as.vector(phi),N,N*p)
bphi= rbind(matphi,cbind(diag(N*(p-1)), matrix(0,N*(p-1),N)))
gp0 = matrix(solve(diag(N^2*p^2) - 
        bphi%x%bphi)%*%as.vector(bdiag(sigma,diag(0,N*(p-1)))),N*p,N*p)
gp2 = sqrtm(gp0)


rep = 5000
pf = seq(1,8)  ## VAR ORDERS TO FIT
nsamp = c(200,500,1000)  ## SAMPLE SIZES
siz = matrix(0,length(pf),length(nsamp))
t<-s<-J<-f<-matrix(0,rep,length(pf))


## Model evaluation statistic: size and power for different 
## orders of misspecification and sample sizes
for(nn in 1:length(nsamp)){
  T = nsamp[nn]
    for(pp in 1:length(pf)){
      pfit = pf[pp]
        for(k in 1:rep){
            print(c(nn,pp,k))
          y = varp_sim(T,N,p,matphi,sigma2,gp2)
          out = ar.yw(t(y),aic=FALSE,order.max=pfit,
            intercept = FALSE,demean = FALSE)
          r =  out$resid[(pfit+1):T,]
          sig = out$var.pred
          pdg = apply(mvfft(r),1,function(x) x%*%Conj(t(x)))/T
          vec = diag(matrix(seq(N^2),N,N))
          rm = rowMeans(pdg)
          J[k,pp] = sqrt(T-pfit)*Re(sum(pdg*Conj(pdg))/T - 
             sum(rm*Conj(rm)) - (sum(rowMeans(pdg[vec,])))^2)
          f[k,pp] = sqrt(4*(sum(diag(sig%^%4)) + sum(diag(sig%^%2))^2))
        }
    siz[pp,nn] = mean(abs(J[,pp]/f[,pp])>qnorm(.975))
    }
}




