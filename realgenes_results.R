load('/Volumes/im-lab/nas40t2/tae/differentialeQTLs/multitissue/localAnc/brca1_run1.Rdata')
brca1 = result
load('/Volumes/im-lab/nas40t2/tae/differentialeQTLs/multitissue/localAnc/hla_c_run2.Rdata')
hlac = result
load('/Volumes/im-lab/nas40t2/tae/differentialeQTLs/multitissue/localAnc/xrra1_run1.Rdata')
xrra = result
rm(result)

xrra_exp = read.table('/Volumes/im-lab/nas40t2/tae/differentialeQTLs/multitissue/prelimtests/xrra1-table.txt',header=TRUE)
tissues = colnames(tissuenames)


plotgamma = function(gene,chain,gamma){
  #gamma = gene$chain$gamma
  nniter = 204
  colors <- c("lightblue", "darkblue")
  datgamma = as.data.frame(gamma[1:nniter, ])
  colnames(datgamma) = tissues
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
  theme(axis.text.y=element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1))+
  ggtitle(paste(gene,chain))
}
plotgamma('brca1','chain1', brca1$chain1$gamma)
plotgamma('brca1','chain2', brca1$chain2$gamma)
plotgamma('hlac','chain1', hlac$chain1$gamma)
plotgamma('hlac','chain2', hlac$chain2$gamma)
plotgamma('xrra1','chain1', xrra$chain1$gamma)
plotgamma('xrra1','chain2', xrra$chain2$gamma)


plotbeta = function(gene, chain, beta){
#beta = brca1$chain1$beta
datbeta = as.data.frame(beta[1:nniter,])
colnames(datbeta) = paste0('v',1:44)
rownames(datbeta) = paste0("iter",1:nniter)
datbeta$id = rownames(datbeta)
datbeta$id = factor(datbeta$id, levels=datbeta$id)
datbeta2= melt(datbeta)
datbeta2$xx = rep(1:nniter, T)
ggplot(datbeta2, aes(x=xx, y=value, col=variable))+
  geom_line(alpha=0.7)+
  ylab('beta')+
  xlab('iteration')+
  ggtitle(paste(gene, chain))
}
plotbeta('brca1','chain1', brca1$chain1$beta)
plotbeta('brca1','chain2', brca1$chain2$beta)
plotbeta('hlac','chain1', hlac$chain1$beta)
plotbeta('hlac','chain2', hlac$chain2$beta)
plotbeta('xrra1','chain1', xrra$chain1$beta)
plotbeta('xrra1','chain2', xrra$chain2$beta)


#Sigma
plotsigma = function(gene,chain,Sigma){
#Sigma = apply(hlac$chain1$Sigma[,,101:200], c(1,2), mean)
ggplot(melt(Sigma), aes(x=Var1, y=Var2, fill=value))+
  geom_tile()+
  scale_fill_gradient(high='indianred', low='antiquewhite2')+
  ggtitle(paste(gene,chain))
}
plotsigma('brca1','chain1',apply(brca1$chain1$Sigma[,,101:200],c(1,2),mean))
plotsigma('hlac','chain1',apply(hlac$chain1$Sigma[,,101:200],c(1,2),mean))
plotsigma('xrra','chain1',apply(xrra$chain1$Sigma[,,101:200],c(1,2),mean))

xrra_exp$id = rownames(xrra_exp)
head(melt(xrra_exp))
ggplot(melt(xrra_exp), aes(x=variable, y=id, fill=value))+
  geom_tile()+
  scale_fill_gradient(high='indianred', low='antiquewhite2')
cov = em_with_zero_mean_c(as.matrix(xrra_exp[,-45]),100)
head(melt(cov))
ggplot(melt(cov), aes(x=Var1, y=Var2, fill=value))+
  geom_tile()+
  scale_fill_gradient(high='indianred', low='antiquewhite2')
