library(MCMCpack)
library(MASS)
library(mvtnorm)
library(dplyr)
library(broom)
library(ggplot2)
library(data.table)
library(tidyr)
library(gridExtra)
R.utils::sourceDirectory("R/", modifiedOnly=FALSE)
Rcpp::sourceCpp('src/localAnc.cpp')
set.seed(1991)
L = read.table('/Volumes/im-lab/nas40t2/tae/differentialeQTLs/multitissue/prelimtests/relevant_L.txt', header = TRUE, stringsAsFactors = FALSE)
Phi = as.matrix(read.table('/Volumes/im-lab/nas40t2/tae/differentialeQTLs/data/covariance_matrices_mean.txt'))
tissues = colnames(rawdat)
T=44

get_tidy = function(rawdat, X){
  Y = as.matrix(rawdat)
  T = ncol(Y)
  n = nrow(Y)
  marcor = colMeans(X*Y, na.rm=TRUE)
  tissues = colnames(Y)

  centeredX = X - mean(X)
  XX = rep(centeredX,T)
  tt = rep(colnames(Y), each=n)
  mod = (lm(as.numeric(scale(Y)) ~ tt:XX-1))
  tidy = tidy(mod)
  #plot(-log10(tidy$p.value), main='(-log10(pvalue))', ylab='-log10(p)', xlab='tissues')
  return(tidy)
}
plot_gamma = function(result, title){
  gamma = result$chain1$gamma
  nniter = nrow(gamma)
  colors <- c("lightblue", "darkblue")
  datgamma = as.data.frame(gamma[1:nniter, ])
  colnames(datgamma) = paste0('v',1:T)
  rownames(datgamma) = paste0("iter",1:nniter)
  datgamma$id = rownames(datgamma)
  datgamma$id = factor(datgamma$id, levels = datgamma$id)
  colnames(datgamma) = factor(colnames(datgamma), levels = colnames(datgamma))
  mdatgamma = melt(datgamma)
  colnames(mdatgamma)[3] = 'gamma'
  mdatgamma$gamma = as.factor(mdatgamma$gamma)
  gam1graph = ggplot(mdatgamma,
                     aes(x = variable, y = id, fill = gamma)) +
    geom_tile(color = 'white') +
    scale_fill_manual(values=colors) +
    ylab("iterations")+
    xlab("gamma")+
    theme(axis.text.y=element_blank())+ggtitle(title)

  gamma = result$chain2$gamma
  nniter = nrow(gamma)
  colors <- c("lightblue", "darkblue")
  datgamma = as.data.frame(gamma[1:nniter, ])
  colnames(datgamma) = paste0('v',1:T)
  rownames(datgamma) = paste0("iter",1:nniter)
  datgamma$id = rownames(datgamma)
  datgamma$id = factor(datgamma$id, levels = datgamma$id)
  colnames(datgamma) = factor(colnames(datgamma), levels = colnames(datgamma))
  mdatgamma = melt(datgamma)
  colnames(mdatgamma)[3] = 'gamma'
  mdatgamma$gamma = as.factor(mdatgamma$gamma)
  gam2graph = ggplot(mdatgamma,
                     aes(x = variable, y = id, fill = gamma)) +
    geom_tile(color = 'white') +
    scale_fill_manual(values=colors) +
    ylab("iterations")+
    xlab("gamma")+
    theme(axis.text.y=element_blank())+ggtitle(title)
  grid.arrange(gam1graph, gam2graph, nrow=1)
}
plot_beta = function(result, title){
  beta = result$chain1$beta
  nniter = nrow(beta)
  datbeta = as.data.frame(beta[1:nniter,])
  colnames(datbeta) = paste0("v",1:T)
  rownames(datbeta) = paste0("iter",1:nniter)
  datbeta$id = rownames(datbeta)
  datbeta$id = factor(datbeta$id, levels=datbeta$id)
  datbeta2= melt(datbeta)
  datbeta2$xx = rep(1:nniter, T)
  beta1graph = ggplot(datbeta2, aes(x=xx, y=value, col=variable))+
    geom_line(alpha=0.7)+
    ylab('beta')+
    xlab('iteration')+ggtitle(title)

  beta = result$chain2$beta
  datbeta = as.data.frame(beta[1:nniter,])
  colnames(datbeta) = paste0("v",1:T)
  rownames(datbeta) = paste0("iter",1:nniter)
  datbeta$id = rownames(datbeta)
  datbeta$id = factor(datbeta$id, levels=datbeta$id)
  datbeta2= melt(datbeta)
  datbeta2$xx = rep(1:nniter, T)
  beta2graph = ggplot(datbeta2, aes(x=xx, y=value, col=variable))+
    geom_line(alpha=0.7)+
    ylab('beta')+
    xlab('iteration')+ggtitle(title)
  grid.arrange(beta1graph, beta2graph, nrow=2)

}

rawdat = read.table('/Volumes/im-lab/nas40t2/tae/differentialeQTLs/data/xrra1_agesexraceadjusted.txt', header=TRUE)
AAind = match(L$sub, rownames(rawdat))
X = rep(0, nrow(rawdat)); X[AAind] = L$xrra1
xrra1_tidy = get_tidy(rawdat, X)

rawdat = read.table('/Volumes/im-lab/nas40t2/tae/differentialeQTLs/data/hlac_agesexraceadjusted.txt', header=TRUE)
AAind = match(L$sub, rownames(rawdat))
X = rep(0, nrow(rawdat)); X[AAind] = L$hlaC
hlac_tidy = get_tidy(rawdat, X)

rawdat = read.table('/Volumes/im-lab/nas40t2/tae/differentialeQTLs/data/brca1_agesexraceadjusted.txt', header=TRUE)
AAind = match(L$sub, rownames(rawdat))
X = rep(0, nrow(rawdat)); X[AAind] = L$brcaL
brca1_tidy = get_tidy(rawdat, X)


load('hla_c_run1.Rdata')
plot_gamma(result, 'hlaC')
plot_beta(result, 'hlaC')

load('xrra1_run3.Rdata')
plot_gamma(result, 'xrra1')
plot_beta(result, 'xrra1')

load('brca1_run1.Rdata')
plot_gamma(result, 'brca1')
plot_beta(result, 'brca1')
