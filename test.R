library(MCMCpack)
library(MASS)
library(mvtnorm)
R.utils::sourceDirectory("R/", modifiedOnly=FALSE)
Rcpp::sourceCpp('src/all_functions.cpp')
set.seed(1991)
#set.seed(1987)
n=400
T = 5
nu = T+10
Phi = diag(T)
sigmabeta = 1
missingprop = 0.6
dat = createData(n, Phi, nu, T, sigmabeta, missingprop)
X = dat$X; Y = dat$Y; Sigma = dat$Sigma

niter = 100
beta = matrix(0, niter, T)
gamma = matrix(0, niter, T)
Sigma = array(dim=c(T, T, niter))
sigmabeta = rep(0,niter)

#initialize
beta[1,] = colMeans(X*Y, na.rm = TRUE)
gamma[1,]= sample(c(0,1), T, replace=TRUE)
beta[1,] = beta[1,] * gamma[1,]
Sigma[,,1] = em_with_zero_mean_c(Y, 100)$Sigma
Sigma[,,1] = dat$Sigma
sigmabeta[1] = 1
Vbeta = 0.1
tar = rep(0,niter)
bgiter = 50
tar = matrix(0,niter, 3)
h = rep(0,niter)
hiter=50
plot(1:niter, rep(0,niter), pch=NA, ylim = c(-2500,-1800))
for (i in 2:niter){
  gam1     = gamma[i-1,]; beta1 = beta[i-1,];
  Sigma1   = Sigma[,,i-1]; sigmabeta1 = sigmabeta[i-1]
  h[i-1]   = get_h_from_sigmabeta(X,sigmabeta[i-1],Sigma1,gam1,n,T)
  h1       = h[i-1]
  bg       = update_betagam(X, Y, gam1, beta1, Sigma1, sigmabeta1, Vbeta, bgiter)
  gam2     = bg$gam[bgiter, ]
  beta2    = bg$beta[bgiter, ]
  Sigma2   = update_Sigma(n,nu,X,beta2,Phi,Y)
  hsig     = update_h(h1, hiter, gam2, beta2, Sigma2, X, T)

  h[i]     = hsig$h[hiter];   gamma[i,] = gam2;
  beta[i,] = beta2;          Sigma[,,i] = Sigma2;
  sigmabeta[i] = hsig$sigbeta[hiter]
  if(sigmabeta[i]==Inf){sigmabeta[i] = 1000}
  tar[i,] = get_target(X,Y,sigmabeta[i], Sigma[,,i], gamma[i,],beta[i,])
  points(i, sum(tar[i,]))
}



library(reshape2); library(ggplot2)

colors <- c("lightblue", "darkblue")
datgamma = as.data.frame(gamma[1:100, ])
colnames(datgamma) = paste0("var",0:4)
rownames(datgamma) = paste0("iter",1:100)
datgamma$id = rownames(datgamma)
datgamma$id = factor(datgamma$id, levels = datgamma$id)
head(melt(datgamma))
ggplot(melt(datgamma),
       aes(x = variable, y = id, fill = factor(value))) +
  geom_tile(color = 'white') +
  scale_fill_manual(values=colors) +
  ylab("iterations")+
  xlab("gamma")+
  theme(axis.text.y=element_blank())

datbeta = as.data.frame(beta[1:100,])
colnames(datbeta) = paste0("var",0:4)
rownames(datbeta) = paste0("iter",1:100)
datbeta$id = rownames(datbeta)
datbeta$id = factor(datbeta$id, levels=datbeta$id)
datbeta2= melt(datbeta)
datbeta2$xx = rep(1:100, 5)
ggplot(datbeta2, aes(x=xx, y=value, col=variable))+
  geom_line()+
  ylab('beta')+
  xlab('iteration')


plot(rowSums(tar)[-1], type = 'l', ylab='joint likelihood', xlab = 'iterations')
