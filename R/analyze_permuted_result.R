library(ggplot2)
library(reshape2)
library(gridExtra)
#analyze result
plot_gamma = function(result, title, T=44){
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
plot_beta = function(result, title, T=44){
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

load('/Volumes/im-lab/nas40t2/tae/differentialeQTLs/multitissue/localAnc/permuted_result.Rdata')
dat = summarized_result1

#ggplot(melt(dat$gamma1, id='genename'),aes(x=variable,y=genename,fill=value)) +
#  geom_tile() + scale_fill_gradient(high='indianred',low='antiquewhite1')
nullpip = c(as.numeric(as.matrix(dat$gamma2[,-1])), as.numeric(as.matrix(dat$gamma1[,-1])))
hist(nullpip, breaks=100, freq=FALSE)
#plot(nullpip)
#length(nullpip)


load('/Volumes/im-lab/nas40t2/tae/differentialeQTLs/multitissue/localAnc/1000genes_summarized_result1.Rdata')
one = summarized_result1; rm(summarized_result1)
allpip = c(as.numeric(as.matrix(one$gamma1[,-1])), as.numeric(as.matrix(one$gamma2[,-1])))
hist(allpip, breaks=100, freq=FALSE)
#plot(allpip, pch='.')
