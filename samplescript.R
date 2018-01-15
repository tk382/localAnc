#sample script
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
beta = c(rep(0.5, 22), rep(0,22))
n = 444; T = length(beta); nu = T+5
Sigma = matrix(0.8, T, T); diag(Sigma) = 1
X = sample(c(0,1,2), size=n, replace = TRUE)
error = mvrnorm(n, rep(0,T), Sigma)
gamma = c(rep(1,22), rep(0,22))
Y = X %*% t(beta) + error
dat = list(Y=Y,X=X,beta=beta,gamma=gamma,Sigma=Sigma) #true data: save this
# Y[sample(1:400, 10),9] = NA
# Y[sample(1:400, 300), 10] = NA
# for (i in c(1:8, 11:20)){
#   Y[sample(1:400, 200), i] = NA
# }
for (i in 1:44){
  Y[which(is.na(tempY[,i])), i] = NA
}

##set up the parameters##
nu = T+5
marcor = colMeans(X*Y, na.rm = TRUE)
Phi = diag(T)
Y = as.matrix(Y)


colnames(Y) = paste0('t',1:44)
XX = rep(X,T)
tt = rep(colnames(Y), each=n)
mod = (lm(as.numeric(scale(Y)) ~ tt:XX-1))
tidy = tidy(mod)
tidy$term = colnames(Y); tidy[,2:5] = round(tidy[,2:5],3)
plot(-log10(tidy$p.value))


#starting points -
#two chains start from different random gammas
tmpgam = sample(c(0,1), size=T, replace=TRUE)
initial_chain1 = list(beta = tmpgam*marcor, gamma = tmpgam, Sigma = Phi, sigmabeta = var(marcor))
tmpgam = sample(c(0,1), size=T, replace=TRUE)
initial_chain2 = list(beta = tmpgam*marcor, gamma = tmpgam, Sigma = Phi, sigmabeta = var(marcor))

result = run2chains_c(X,Y,initial_chain1,initial_chain2,
                      Phi,colMeans(X*Y, na.rm=TRUE),
                      niter=100,bgiter=500,hiter=50,switer=30)

# run2chains argument
# X,Y,initial_chain1,initial_chain2,
# Phi, marcor,
# niter=1000,bgiter=500,hiter=50,switer=50,burnin=5

#let's look at the first chain's result
get_target(X,Y,1,result$chain2$Sigma[,,50], dat$gam, dat$beta, 20)
get_target(X,Y,1,result$chain2$Sigma[,,50], rep(0,20), rep(0,20), 20)

sum(get_target(X,Y,1,result$chain2$Sigma[,,50], dat$gam, dat$beta, 20))
sum(get_target(X,Y,1,result$chain2$Sigma[,,50], rep(0,20), rep(0,20), 20))

get_target(X,Y,1,result$chain2$Sigma[,,50], c(rep(1,9),rep(0,11)), c(dat$beta[1:9], rep(0,11)), 20)
get_target(X,Y,1,result$chain2$Sigma[,,50], dat$gamma, dat$beta, 20)


ggplot(melt(result$chain1$Sigma[,,50]), aes(x=Var1, y=Var2, fill=value)) + geom_tile() +
  scale_color_gradient()

save(result, file = 'samplescript_result.Rdata')
save(dat, file='samplescript_trudata.Rdata')





##visualize result
#gamma
gamma = result$chain1$gamma
nniter = 35
colors <- c("lightblue", "darkblue")
datgamma = as.data.frame(gamma[1:nniter, ])
colnames(datgamma) = paste0('v',1:T)
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

#beta
beta = result$chain1$beta
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

#Sigma
outSigma1 = apply(result$chain2$Sigma[,,6:35], c(1,2), mean)
ggplot(melt(outSigma1), aes(x=Var1, y=Var2, fill=value))+
  geom_tile()+
  scale_fill_gradient(high='indianred', low='antiquewhite2')
