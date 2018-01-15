library(reshape2)
library(ggplot2)
library(gridExtra)
tmp = read.table(
  '/Volumes/im-lab/nas40t2/tae/differentialeQTLs/bytissues/Muscle_Skeletal/temp/tests/e~A/significant.result.txt',
  header = TRUE,
  stringsAsFactors = FALSE)
names = tmp$name

#'ENSG00000163466.11' : pretty cool result - about half selected in 200 iterations
genename = 'ENSG00000204666.3'
T = 44
setwd('/Volumes/im-lab/nas40t2/tae/differentialeQTLs/multitissue/localAnc')
load(paste0(genename,'.Rdata'))


plot_gamma = function(result, title){
  gamma = result$chain1$gamma
  nniter = nrow(gamma)
  colors <- c("lightblue", "darkblue")
  datgamma = as.data.frame(gamma[1:nniter, ])
  colnames(datgamma) = paste0('v',1:44)
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
  colnames(datgamma) = paste0('v',1:44)
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
    geom_line(alpha=0.3)+
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
    geom_line(alpha=0.3)+
    ylab('beta')+
    xlab('iteration')+ggtitle(title)
  grid.arrange(beta1graph, beta2graph, nrow=2)

}

plot_gamma(final$result, genename)
plot_beta(final$result, genename)
ind1 = which(colMeans(final$result$chain1$gamma[-(1:20), ])>0.5)
ind2 = which(colMeans(final$result$chain2$gamma[-(1:20), ])>0.5)
plot(-log10(final$tidy$p.value))
points(ind1, -log10(final$tidy$p.value[ind1]), col = 'red')
