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


rawdat = read.table('/Volumes/im-lab/nas40t2/tae/differentialeQTLs/data/xrra1_agesexraceadjusted.txt', header=TRUE)
# rawdat = read.table('/Volumes/im-lab/nas40t2/tae/differentialeQTLs/data/hlac_agesexraceadjusted.txt', header=TRUE)
# rawdat = read.table('/Volumes/im-lab/nas40t2/tae/differentialeQTLs/data/brca1_agesexraceadjusted.txt', header=TRUE)

L = read.table('/Volumes/im-lab/nas40t2/tae/differentialeQTLs/multitissue/prelimtests/relevant_L.txt', header = TRUE, stringsAsFactors = FALSE)
AAind = match(L$sub, rownames(rawdat))
X = rep(0, nrow(rawdat))
X[AAind] = L$xrra1
#X[AAind] = L$brcaL
#X[AAind] = L$hlaC
Y = as.matrix(rawdat)
T = ncol(Y)
n = nrow(Y)
Phi = as.matrix(read.table('/Volumes/im-lab/nas40t2/tae/differentialeQTLs/data/covariance_matrices_mean.txt'))
marcor = colMeans(X*Y, na.rm=TRUE)
tissues = colnames(Y)

#linear regression result
centeredX = X - mean(X)
XX = rep(centeredX,T)
tt = rep(colnames(Y), each=n)
mod = (lm(as.numeric(scale(Y)) ~ tt:XX-1))
tidy = tidy(mod)
plot(-log10(tidy$p.value), main='brca1_(-log10(pvalue))', ylab='-log10(p)', xlab='tissues')

tmpgam = sample(c(0,1),size=T,replace=TRUE)
#tmpgam = rep(0, T)
initial_chain1 = list(beta=tmpgam*tidy$estimate,gamma=tmpgam,Sigma=diag(T),sigmabeta=1)
tmpgam = sample(c(0,1), size=T, replace=TRUE)
#tmpgam = rep(1, T)
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

#save(result, file='hla_c_run1.Rdata')
#save(result, file='xrra1_run2.Rdata')
#save(result, file = 'brca1_run1.Rdata')

load('brca1_run1.Rdata')

gamma = result$chain1$gamma
nniter = nrow(gamma)
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
gam1graph = ggplot(mdatgamma,
       aes(x = variable, y = id, fill = gamma)) +
  geom_tile(color = 'white') +
  scale_fill_manual(values=colors) +
  ylab("iterations")+
  xlab("gamma")+
  theme(axis.text.y=element_blank())+ggtitle('brca1_gam1')

gamma = result$chain2$gamma
nniter = nrow(gamma)
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
gam2graph = ggplot(mdatgamma,
                   aes(x = variable, y = id, fill = gamma)) +
  geom_tile(color = 'white') +
  scale_fill_manual(values=colors) +
  ylab("iterations")+
  xlab("gamma")+
  theme(axis.text.y=element_blank())+ggtitle('brca1_gam2')

grid.arrange(gam1graph, gam2graph, nrow=1)


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

ggplot(melt(result$chain1$Sigma[,,nniter]), aes(x=Var1, y=Var2, fill=value))+
  geom_tile()+
  scale_fill_gradient(high='indianred', low='antiquewhite2')
ggplot(melt(Phi), aes(x=Var1,y=Var2,fill=value)) + geom_tile() +
  scale_fill_gradient(high='indianred', low='antiquewhite2')
