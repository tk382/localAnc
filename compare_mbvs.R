library(MCMCpack)
library(MASS)
library(mvtnorm)
library(dplyr)
library(ggplot2)
library(data.table)
R.utils::sourceDirectory("R/", modifiedOnly=FALSE)
Rcpp::sourceCpp('src/all_functions.cpp')
set.seed(1991)


#create data with the following beta
beta = c(seq(0.05, 0.5, by=0.05), rep(0,10))
n = 400; T = 20; nu = T+5
Sigma = matrix(0.3, T, T); diag(Sigma) = 1
X = sample(c(0,1,2), size=n, replace = TRUE)
error = mvrnorm(n, rep(0,T), Sigma)
gamma = c(rep(1,10), rep(0,10))
Y = X %*% t(beta) + error
dat = list(Y=Y,X=X,beta=beta,gamma=gamma,Sigma=Sigma) #true data: save this


nu = T+5
initialbet = rep(0,T)
initialgam = rep(0,T)
initialSigma = diag(T)
initialsb = 1
marcor = colMeans(X*Y, na.rm = TRUE)
Vbeta = mean(marcor^2) * 0.01
niter=500; bgiter=500; hiter=50; switer=30
Y = as.matrix(Y)
print('starting the first algorithm..')
res = doMCMC(X=X,Y=Y,n=n,T=T,Phi=diag(T),nu=nu,
             initialbet,
             initialgam,
             initialSigma,
             initialsb,
             marcor,
             Vbeta,
             niter,
             bgiter,
             hiter,
             switer)
out = list(res, dat)
save(out, file = 'different_betas.Rdata')

library(mBvs)
data <- data.frame(L=X)
form1 <- as.formula( ~ L)
form2 <- as.formula( ~ 1)
lin.pred <- list(form1, form2)
p  <- dim(data)[2]
p_adj  <- 0
q  <- dim(Y)[2]
eta = 0.1
v = rep(10, q)
omega = rep(log(0.5/(1-0.5)), p-p_adj)
common.beta0 <- c(rep(0, q), 10^6)
rho0  <- q + 4
Psi0 <- diag(3, q)
US.Sigma <- c(rho0, Psi0)
hyperParams <- list(eta=eta, v=v, omega=omega, beta0=common.beta0,
                    US=list(US.Sigma=US.Sigma))
numReps    <- 1000
thin       <- 1
burninPerc <- 0.5
mhProp_beta_var  <- matrix(0.5, p+p_adj, q)
mhrho_prop <- 1000
mhPsi_prop <- diag(1, q)
mcmc.US <- list(run=list(numReps=numReps, thin=thin, burninPerc=burninPerc),
                tuning=list(mhProp_beta_var=mhProp_beta_var,
                            mhrho_prop=mhrho_prop, mhPsi_prop=mhPsi_prop))

beta0 <- rep(0, q)
B <- matrix(sample(x=c(0.3, 0), size=q, replace = TRUE), p+p_adj, q)
gamma <- B
gamma[gamma != 0] <- 1
Sigma <- diag(1, q)
lambda <- rep(0.5, q)
sigmaSq <- 1
startValues <- vector("list", 2)
startValues[[1]] <- as.vector(c(beta0, B, gamma, Sigma))
beta0 <- rep(0.2, q)
Sigma <- diag(0.5, q)
startValues[[2]] <- as.vector(c(beta0, B, gamma, Sigma))

mBvsresult = mvnBvs(Y,
                   lin.pred,
                   data, model = "unstructured",
                   hyperParams, startValues,
                   mcmc.US)

save(mBvsresult, file='mBvsresult.Rdata')


gamma = res$gam
gamma = t(apply(mBvsresult$chain1$gamma.p, 3, as.numeric))

nniter = 500
colors <- c("lightblue", "darkblue")
datgamma = as.data.frame(gamma[1:nniter, ])
colnames(datgamma) = paste0('var',1:20)
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
  theme(axis.text.y=element_blank())+
  ggtitle('(a) Ours' )








nniter = 500
beta = res$beta
beta = mBvsresult$chain1$beta0.p

datbeta = as.data.frame(beta[1:nniter,])
colnames(datbeta) = paste0("var",1:T)
rownames(datbeta) = paste0("iter",1:nniter)
datbeta$id = rownames(datbeta)
datbeta$id = factor(datbeta$id, levels=datbeta$id)
datbeta2= melt(datbeta)
datbeta2$xx = rep(1:nniter, T)
ggplot(datbeta2, aes(x=xx, y=value, col=variable))+
  geom_line(alpha=0.7)+
  ylab('beta')+
  xlab('iteration')+
  ylim(-0.2,0.75)+
  ggtitle('(b) mBvs')

unname(outSigma1)
outSigma1 = apply(out[[1]]$Sigma[,,-(1:100)], c(1,2), mean)
outSigma2 = apply(mBvsresult$chain1$Sigma.p[,,-(1:100)], c(1,2),mean)
ggplot(melt(outSigma1), aes(x=Var1, y=Var2, fill=value))+
  geom_tile()+
  scale_fill_gradient(high='indianred', low='antiquewhite2')+
  ggtitle('(c) Ours')
ggplot(melt(outSigma2), aes(x=Var1, y=Var2, fill=value))+
  geom_tile()+
  scale_fill_gradient(high='indianred', low='antiquewhite2')+
  ggtitle('(c) mBvs')


for (i in 1:20){
  ind = gamma[-(1:100),i]==1
  betatmp = beta[-(1:100),i]
  print(round(sd(betatmp[ind]),3))
}
