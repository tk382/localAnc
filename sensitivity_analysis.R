#sensitivity analysis

library(MCMCpack)
library(MASS)
library(mvtnorm)
library(dplyr)
library(ggplot2)
library(data.table)
R.utils::sourceDirectory("R/", modifiedOnly=FALSE)
Rcpp::sourceCpp('src/all_functions.cpp')
set.seed(1991)
#set.seed(1987)
n=500
T = 5
nu = T+10
Phi = diag(T)
missingprop = 0
niter = 100
bgiter = 100
smallworlditer=100


sigbetavec = c(0.001, 0.01, 0.1, 1, 10)
ratiovec = c(1/2,1/4,1/10,1/20,1/50,1/100)

save_convergence_iter = matrix(0,length(sigbetavec), length(ratiovec))
truegamma = rep(0,5)
truebeta = rep(0,5)

for (k in 1:length(sigbetavec)){
  truesigmabeta = sigbetavec[k]
  dat = createData(n, Phi, nu, T, truesigmabeta, missingprop)
  X = dat$X; Y = dat$Y; Sigma = dat$Sigma
  truebeta = rbind(truebeta, dat$beta)
  truegamma = rbind(truegamma, dat$gamma)
  for (j in 1:length(ratiovec)){
    Vbeta = ratiovec[j] * truesigmabeta
    #saving
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
    Vbeta = 0.05
    tar = rep(0,niter)
    tar = matrix(0,niter, 3)
    h = rep(0,niter)
    hiter=50
    marcor = abs(colMeans(X*Y, na.rm=TRUE))
    for (i in 2:niter){
      gam1     = gamma[i-1,]; beta1 = beta[i-1,];
      Sigma1   = Sigma[,,i-1]; sigmabeta1 = sigmabeta[i-1]
      h[i-1]   = get_h_from_sigmabeta_c(X,sigmabeta[i-1],Sigma1,gam1,n,T)
      h1       = h[i-1]
      bg       = update_betagam_sw_c(X,
                                     Y,
                                     gam1,
                                     beta1,
                                     Sigma1,
                                     marcor,
                                     sigmabeta1,
                                     Vbeta,
                                     bgiter,
                                     100)
      gam2     = bg$gam[bgiter, ]
      beta2    = bg$beta[bgiter, ]
      Sigma2   = update_Sigma_c(n,nu,X,beta2,Phi,Y)
      hsig     = update_h_c(h1, hiter, gam2, beta2, Sigma2, X, T)
      h[i]     = hsig$h[hiter];   gamma[i,] = gam2;
      beta[i,] = beta2;          Sigma[,,i] = Sigma2;
      sigmabeta[i] = hsig$sigbeta[hiter]
      if(sigmabeta[i]==Inf){sigmabeta[i] = 1000}
      tar[i,] = get_target_c(X,Y,sigmabeta[i], Sigma[,,i], gamma[i,],beta[i,])
      if(i>3){
        colmean = colMeans(gamma[2:i, ])
        if(identical(as.numeric(colmean>0.5), as.numeric(dat$gamma==1))){
          save_convergence_iter[k,j] = i
        }
      }
      if(i%%10==0){print(i)}
    }
    library(reshape2); library(ggplot2)

    niter=100
    colors <- c("lightblue", "darkblue")
    datgamma = as.data.frame(gamma[1:niter, ])
    colnames(datgamma) = paste0("var",1:T)
    rownames(datgamma) = paste0("iter",1:niter)
    datgamma$id = rownames(datgamma)
    datgamma$id = factor(datgamma$id, levels = datgamma$id)
    colnames(datgamma) = factor(colnames(datgamma), levels = colnames(datgamma))
    head(melt(datgamma))
    mdatgamma = melt(datgamma)

    colnames(mdatgamma)[3] = 'gamma'
    mdatgamma$gamma = as.factor(mdatgamma$gamma)
    ggplot(mdatgamma,
         aes(x = variable, y = id, fill = gamma)) +
      geom_tile(color = 'white') +
      scale_fill_manual(values=colors) +
      ylab("iterations")+
      xlab("gamma")+
      theme(axis.text.y=element_blank())+
      ggtitle(paste(sigbetavec[k], ratiovec[j]))
    ggsave(filename=paste0(sigbetavec[k], '_', ratiovec[j],'_gamma.png'))

    datbeta = as.data.frame(beta[1:niter,])
    colnames(datbeta) = paste0("var",1:T)
    rownames(datbeta) = paste0("iter",1:niter)
    datbeta$id = rownames(datbeta)
    datbeta$id = factor(datbeta$id, levels=datbeta$id)
    datbeta2= melt(datbeta)
    datbeta2$xx = rep(1:niter, T)
    ggplot(datbeta2, aes(x=xx, y=value, col=variable))+
      geom_line()+
      ylab('beta')+
      xlab('iteration')+
      ggtitle(paste(sigbetavec[k], ratiovec[j]))
    ggsave(filename=paste0(sigbetavec[k], '_', ratiovec[j],'_beta.png'))
  }
}





