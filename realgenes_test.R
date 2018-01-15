library(MCMCpack)
library(MASS)
library(mvtnorm)
library(dplyr)
library(broom)
library(ggplot2)
library(data.table)
R.utils::sourceDirectory("R/", modifiedOnly=FALSE)
Rcpp::sourceCpp('src/localAnc.cpp')
set.seed(1991)


#rawdat = read.table('/Volumes/im-lab/nas40t2/tae/differentialeQTLs/multitissue/prelimtests/hla-c-table.txt', header = TRUE)
#rawdat = read.table('/Volumes/im-lab/nas40t2/tae/differentialeQTLs/multitissue/prelimtests/brca1-table.txt', header = TRUE)
rawdat = read.table('/Volumes/im-lab/nas40t2/tae/differentialeQTLs/multitissue/prelimtests/xrra1-table.txt', header=TRUE)
L = read.table('/Volumes/im-lab/nas40t2/tae/differentialeQTLs/multitissue/prelimtests/relevant_L.txt', header = TRUE, stringsAsFactors = FALSE)
AAind = match(L$sub, rownames(rawdat))
X = rep(0, nrow(rawdat))
#X[AAind] = L$hlaC
#X[AAind] = L$brcaL
X[AAind] = L$xrra1
Y = as.matrix(rawdat)
T = ncol(Y)
n = nrow(Y)
Phi = as.matrix(read.table('/Volumes/im-lab/nas40t2/tae/differentialeQTLs/data/covariance_matrices_mean.txt'))
marcor = colMeans(X*Y, na.rm=TRUE)
tissues = colnames(Y)

#linear regression result
centeredX = X - mean(X)
XX = rep(centeredX,44)
tt = rep(colnames(Y), each=444)
mod = (lm(as.numeric(scale(Y)) ~ tt:XX-1))
tidy = tidy(mod)

signal_ind = sort(log10(tidy$p.value), index.return=TRUE)$ix
tmpgam = rep(0, T)
tmpgam[signal_ind[1:8]] = 1
initial_chain1 = list(beta=tmpgam*tidy$estimate,gamma=tmpgam,Sigma=diag(T),sigmabeta=1)
tmpgam = rep(0,T)
tmpgam[signal_ind[1:12]] = 1
initial_chain2 = list(beta=tmpgam*tidy$estimate,gamma=tmpgam,Sigma=diag(T),sigmabeta=1)

#marcor = abs(tidy$statistic)
#ix = sort(marcor, index.return=TRUE)$ix
#Y = Y[,ix]

result = run2chains_c(centeredX, scale(Y),
                      initial_chain1, initial_chain2,
                      Phi,
                      #colMeans(X*Y, na.rm=TRUE),
                      tidy$statistic,
                      niter=50, bgiter=1000,switer=50,burnin=5)
# run2chains argument
# X,Y,
# initial_chain1,
# initial_chain2,
# Phi,
# marcor,
# niter=1000,
# bgiter=500,
# hiter=50,
# switer=50,
# burnin=5

nu=T+4
tmpgam = sample(c(0,1),size=T,replace=TRUE)
result2 = doMCMC_c(centeredX,Y,n,T,Phi,nu,tmpgam*marcor, tmpgam, diag(T), var(marcor), abs(tidy$statistic), mean(marcor^2)*0.001,100,1000,50,50)

#save(result, file='hla_c_run2.Rdata')
#save(result, file='brca1_run1.Rdata')
save(result, file = 'xrra1_run1.Rdata')


gamma = result$chain2$gamma
nniter = 50
colors <- c("lightblue", "darkblue")
datgamma = as.data.frame(gamma[1:nniter, ])
colnames(datgamma) = paste0('v',1:44)
rownames(datgamma) = paste0("iter",1:nniter)
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
  theme(axis.text.y=element_blank())

beta = result$chain2$beta
datbeta = as.data.frame(beta[1:nniter,])
colnames(datbeta) = paste0("v",1:T)
rownames(datbeta) = paste0("iter",1:nniter)
datbeta$id = rownames(datbeta)
datbeta$id = factor(datbeta$id, levels=datbeta$id)
datbeta2= melt(datbeta)
datbeta2$xx = rep(1:nniter, T)
ggplot(datbeta2, aes(x=xx, y=value, col=variable))+
  geom_line(alpha=0.7)+
  ylab('beta')+
  xlab('iteration')

ggplot(melt(outSigma1[,,i-1]), aes(x=Var1, y=Var2, fill=value))+
  geom_tile()+
  scale_fill_gradient(high='indianred', low='antiquewhite2',limits=c(-0.3,1.5))+
  ggtitle('(c) Ours')
