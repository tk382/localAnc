library(MCMCpack)
library(MASS)
library(mvtnorm)
R.utils::sourceDirectory("R/", modifiedOnly=FALSE)
Rcpp::sourceCpp('src/all_functions.cpp')
set.seed(1991)
n=400
T = 5
nu = T+10
Phi = diag(T)
sigmabeta = 1
missingprop = 0
dat = createData(n=400, Phi = diag(T), nu=T+10, T = 5, sigmabeta = 1, missingprop = 0)
X = dat$X; Y = dat$Y; Sigma = dat$Sigma

niter = 2000
beta = matrix(0, niter, T)
gamma = matrix(0, niter, T)
Sigma = array(dim=c(T, T, niter))
sigmabeta = rep(0,niter)
#initialize
beta[1,] = colMeans(X*Y, na.rm = TRUE)
gamma[1,] = as.numeric(abs(beta[1,])>0.5)
beta[1,] = beta[1,] * gamma[1,]
#Sigma[,,1] = em_with_zero_mean_c(Y, 100)$Sigma
Sigma[,,1] = dat$Sigma
sigmabeta[1] = var(beta[1,gamma[1,]==1]/diag(Sigma[,,1])[gamma[1,]==1])
Vbeta = 0.2
for (i in 2:niter){
  gam1 = gamma[i-1,]; beta1 = beta[i-1,]; Sigma1 = Sigma[,,i-1]; sigmabeta1 = sigmabeta[i-1]
  betagamma = update_betagam(X, Y, gam1, beta1, Sigma1, sigmabeta1, Vbeta)
  change= betagamma$change;     beta2 = betagamma$beta;
  gam2 = betagamma$gamma;       changeind = betagamma$changeind
  Sigma1 = round(Sigma1, 5)
  A = sum(betagam_accept(X, Y, sigmabeta1, Sigma1, gam1, beta1, gam2, beta2, changeind, change))
  check = runif(1,0,1)
  if(exp(A)>check){
    print('updated')
    gamma[i,] = gam2; beta[i,] = beta2; Sigma[,,i] = Sigma1; sigmabeta[i] = sigmabeta1
  }else{
    print('not updated')
    gamma[i,] = gam1; beta[i,] = beta1; Sigma[,,i] = Sigma1; sigmabeta[i] = sigmabeta1
  }
  print(beta2)
  print(beta1)
  print(exp(A))
  i=i+1
}

