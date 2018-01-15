library(MCMCpack)
library(MASS)
library(mvtnorm)
library(dplyr)
library(ggplot2)
library(data.table)
R.utils::sourceDirectory("R/", modifiedOnly=FALSE)
Rcpp::sourceCpp('src/localAnc.cpp')
set.seed(1991)

##create data##
beta = c(seq(0.05, 0.5, by=0.05), rep(0,10))
n = 400; T = 20; nu = T+5
Sigma = matrix(0.3, T, T); diag(Sigma) = 1
X = sample(c(0,1,2), size=n, replace = TRUE)
error = mvrnorm(n, rep(0,T), Sigma)
gamma = c(rep(1,10), rep(0,10))
Y = X %*% t(beta) + error
dat = list(Y=Y,X=X,beta=beta,gamma=gamma,Sigma=Sigma) #true data: save this
colnames(Y) = paste0('t', 1:T)

centeredX = X - mean(X)
XX = rep(centeredX,T)

tt = rep(colnames(Y), each=n)
tt = factor(tt, levels=colnames(Y))
mod = (lm(as.numeric(scale(Y)) ~ tt:XX-1))
tidy = tidy(mod)

#############################################


##initialize our algorithm##
nu = T+5
initialbet = rep(0,T)
initialgam = rep(0,T)
initialSigma = diag(T)
initialsb = 1
marcor = colMeans(X*Y, na.rm = TRUE)
Vbeta = mean(marcor^2) * 0.01
niter=60; bgiter=500; hiter=50; switer=30
Phi = diag(T)
Y = as.matrix(Y)
initial_chain1 = list(beta=initialbet, gamma = initialgam, Sigma = diag(T), sigmabeta = initialsb)
#initial_chain1 = list(beta = rep(0,T), gamma = rep(0,T), Sigma = diag(T), sigmabeta = var(marcor))
tmpgam = sample(c(0,1), size=T, replace=TRUE)
initial_chain2 = list(beta = tmpgam*marcor, gamma = tmpgam, Sigma = em_with_zero_mean_c(Y, 100), sigmabeta = var(marcor))

res2 = run2chains_c(X,Y,initial_chain1,initial_chain2,diag(T), colMeans(X*Y, na.rm=TRUE),niter=60,bgiter=500,hiter=50,switer=30)

##initialize mbvs##
library(mBvs)
data <- data.frame(L=X)
form1 <- as.formula( ~ L)
form2 <- as.formula( ~ 1)
lin.pred <- list(form1, form2)
p  <- dim(data)[2]
p_adj  <- 0
q  <- dim(Y)[2]
eta = 0.1
v = rep(var(marcor)*5, q)
omega = rep(log(0.5/(1-0.5)), p-p_adj)
common.beta0 <- c(rep(0, q), 0.1)
rho0 	<- q + 4
Psi0	<- diag(1, q)
US.Sigma <- c(rho0, Psi0)
hyperParams <- list(eta=eta, v=v, omega=omega, beta0=common.beta0,
                    US=list(US.Sigma=US.Sigma))
numReps    <- 60
thin       <- 1
burninPerc <- 0.5
mhProp_beta_var  <- matrix(0.001, p+p_adj, q)
mhrho_prop <- 500
mhPsi_prop <- diag(1, q)
mcmc.US <- list(run=list(numReps=numReps, thin=thin, burninPerc=burninPerc),
                tuning=list(mhProp_beta_var=mhProp_beta_var,
                            mhrho_prop=mhrho_prop, mhPsi_prop=mhPsi_prop))
lambda = rep(0.5, q)
sigmaSq = 1
startValues = vector("list", 1)
beta0 = rep(0, q)
gamma = sample(c(0,1),size=q, replace=TRUE)
B = matrix(marcor*gamma, nrow=1)
Sigma = diag(1, q)
startValues[[1]] = as.vector(c(beta0, B, gamma, Sigma))


nsim = 100
gam_ours = matrix(0,nsim,T)
beta_ours = matrix(0,nsim,T)
Sigma_ours = array(0, dim=c(T,T,nsim))

gam_mbvs = matrix(0,nsim,T)
beta_mbvs = matrix(0,nsim,T)
Sigma_mbvs = array(0,dim=c(T,T,nsim))
for (sim in 1:nsim){
  res = doMCMC_c(X=X,Y=Y,n=n,T=T,Phi=diag(T),nu=nu,
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
  mBvsresult = mvnBvs(Y,
                    lin.pred,
                    data, model = "unstructured",
                    hyperParams, startValues,
                    mcmc.US)
  ourgam = res$gam[31:60,]
  theirgam = t(apply(mBvsresult$chain1$gamma.p, 3, as.numeric))
  ourbeta = res$beta[31:60,]
  theirbeta = t(apply(mBvsresult$chain1$B.p, 3, as.numeric))
  gam_ours[sim, ] = colMeans(ourgam)
  gam_mbvs[sim, ] = colMeans(theirgam)
  for (i in 1:20){
    ind1 = ourgam[,i]==1
    ind2 = theirgam[,i]==1
    beta_ours[sim, i] = mean(ourbeta[ind1, i])
    beta_mbvs[sim, i] = mean(theirbeta[ind2, i])
  }
  Sigma_ours[,,sim] = apply(res$Sigma[,,31:60], c(1,2), mean)
  Sigma_mbvs[,,sim] = apply(mBvsresult$chain1$Sigma.p, c(1,2), mean)
}
save(gam_ours, file = 'gam_ours.Rdata')
save(beta_ours, file = 'beta_ours.Rdata')
save(Sigma_ours, file = 'Sigma_ours.Rdata')
save(gam_mbvs, file = 'gam_mbvs.Rdata')
save(beta_mbvs, file = 'beta_mbvs.Rdata')
save(Sigma_mbvs, file = 'Sigma_mbvs.Rdata')


gam_ours = gam_ours[1:50, ]
beta_ours = beta_ours[1:50, ]
Sigma_ours = Sigma_ours[,,1:50]
gam_mbvs = gam_mbvs[1:50, ]
beta_mbvs = beta_mbvs[1:50, ]
Sigma_mbvs = Sigma_mbvs[,,1:50]






save(mBvsresult, file='mBvsresult_idealsetting.Rdata')
save(out, file = 'different_betas.Rdata')
load('mBvsresult.Rdata')
load('different_betas.Rdata')
res = out[[1]]

gamma = result2$gam
gamma = t(apply(mBvsresult$chain1$gamma.p, 3, as.numeric))
gamma = t(apply(mBvsresult$chain2$gamma.p, 3, as.numeric))

nniter = 100
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
  theme(axis.text.y=element_blank())+
  ggtitle('(a) mBvs' )

nniter = 30
beta = res$beta[31:60, ]
beta = t(apply(mBvsresult$chain1$B.p, 3, as.numeric))
beta2 = t(apply(mBvsresult$chain2$B.p, 3, as.numeric))

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
  xlab('iteration')+
  ylim(-0.2,1)+
  ggtitle('(b) mBvs')

outSigma1 = result2$Sigma[,,100]
outSigma1 = apply(res$Sigma[,,31:60], c(1,2), mean)
outSigma2 = apply(mBvsresult$chain1$Sigma.p[,,], c(1,2),mean)
ggplot(melt(outSigma1), aes(x=Var1, y=Var2, fill=value))+
  geom_tile()+
  scale_fill_gradient(high='indianred', low='antiquewhite2',limits=c(-0.3,1.5))+
  ggtitle('(c) Ours')
ggplot(melt(outSigma2), aes(x=Var1, y=Var2, fill=value))+
  geom_tile()+
  scale_fill_gradient(high='indianred', low='antiquewhite2', limits=c(-0.3,1.5))+
  ggtitle('(c) mBvs')



for (i in 1:20){
  #print(round(mean(beta_mbvs[,i]*gam_mbvs[,i], na.rm=TRUE),3))
  print(round(sd(beta_mbvs[,i] * gam_mbvs[,i], na.rm=TRUE)*1000,3))
}
