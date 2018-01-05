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

n=500
T = 10
nu = T+5
Phi = diag(T)
sigmabeta = 1
missingprop = 0.3
dat = createData(n, Phi, nu, T, sigmabeta, missingprop)
X = dat$X; Y = dat$Y; Sigma = dat$Sigma
niter = 50
bgiter=200
hiter=50
switer=30
marcor = colMeans(X*Y, na.rm=TRUE)
Vbeta = mean(marcor^2) * 0.002







res = doMCMC_c(X=X,Y=Y,n=n,T=T,Phi=Phi,nu=nu,
             initialbeta,
             initialgamma,
             initialSigma,
             initialsb,
             marcor,
             Vbeta,
             niter=30,
             bgiter,
             hiter,
             switer)


save(res, file = "start_zeros.Rdata")


res2 = doMCMC(X=X,Y=Y,n=n,T=T,Phi=Phi,nu=nu,
             rep(0,T),
             rep(1,T),
             initialSigma,
             initialsb,
             marcor,
             Vbeta,
             niter,
             bgiter,
             hiter,
             switer)
save(res2, file="start_ones.Rdata")

initialgam = sample(c(0,1),size=T,replace=TRUE)
initialbet = initialgam*marcor
res3 = doMCMC(X=X,Y=Y,n=n,T=T,Phi=Phi,nu=nu,
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
save(res3, file="start_random.Rdata")



library(reshape2); library(ggplot2)


##compare with R


niter=2; bgiter=100;hiter=100;switer=100;
res = doMCMC_c(X=X,Y=Y,n=n,T=T,Phi=Phi,nu=nu,
             initialbeta,
             initialgamma,
             initialSigma,
             initialsb,
             marcor,
             Vbeta,
             niter=2,
             bgiter,
             hiter,
             switer)

res2 = doMCMC_R(X=X,Y=Y,n=n,T=T,Phi=Phi,nu=nu,
                initialbeta,
                initialgamma,
                initialSigma,
                initialsb,
                marcor,
                Vbeta,
                niter=2,
                bgiter,
                hiter,
                switer)

microbenchmark(res = doMCMC_c(X=X,Y=Y,n=n,T=T,Phi=Phi,nu=nu,
                              initialbeta,
                              initialgamma,
                              initialSigma,
                              initialsb,
                              marcor,
                              Vbeta,
                              niter=2,
                              bgiter,
                              hiter,
                              switer),
               res2 = doMCMC_R(X=X,Y=Y,n=n,T=T,Phi=Phi,nu=nu,
                               initialbeta,
                               initialgamma,
                               initialSigma,
                               initialsb,
                               marcor,
                               Vbeta,
                               niter=2,
                               bgiter,
                               hiter,
                               switer),
               times = 20)


#Unit: milliseconds
#expr     min       lq      mean    median        uq     max    neval cld
#res  593.7042  598.722  625.7765  613.6157  634.5351  702.2796   20   a
#res2 1510.7107 1554.319 1618.2936 1636.3927 1667.9950 1734.1371  20   b


#after edit
#Unit: milliseconds
#expr       min        lq     mean   median       uq      max neval cld
#res  967.3445  989.4781 1018.904 1018.312 1040.978 1096.271    20  a
#res2 4670.5110 4740.5762 4828.640 4830.155 4874.727 5052.816    20   b


initial_chain1 = list(beta = rep(0,T), gamma = rep(0,T), Sigma = diag(T), sigmabeta = var(marcor))
tmpgam = sample(c(0,1), size=T, replace=TRUE)
initial_chain2 = list(beta = tmpgam*marcor, gamma = tmpgam, Sigma = em_with_zero_mean_c(Y, 100), sigmabeta = var(marcor))
result = run2chains_c(X,Y,
                      initial_chain1, initial_chain2,
                      Phi, colMeans(X*Y, na.rm=TRUE))


centeredX = X - mean(X)
XX = rep(centeredX,10)
tt = rep(colnames(Y), each=500)
mod = (lm(as.numeric(scale(Y)) ~ tt:XX-1))
tidy = tidy(mod)


round(colMeans(result$chain1$gamma[-(1:5),]))
round(colMeans(result$chain2$gamma[-(1:5),]))
dat$gamma
