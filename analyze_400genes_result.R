library(ggplot2)
library(reshape2)
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

#### DATA SET ONE ####
load('/Volumes/im-lab/nas40t2/tae/differentialeQTLs/multitissue/localAnc/400genes_summarized_result1.Rdata')
one = summarized_result1
rm(summarized_result1)
ggplot(melt(one$finalgamma, id='genename'),aes(x=variable,y=genename,fill=value)) +
         geom_tile() + scale_fill_gradient(high='indianred',low='antiquewhite1')


converged_one = which(one$niter!=1000)
unconverged_one = which(one$niter == 1000)
convergedgamma = one$finalgamma[converged_one, ]

##distribution of 'how many tissues' are selected
newdf = data.frame(tissue_count = as.factor(rowSums(convergedgamma[,2:45])/2))
ggplot(newdf, aes(x=tissue_count)) + geom_bar()

#plot the per-tissue number of selected genes
newdf = data.frame(tissues = colnames(one$p)[-1], counts = colSums(convergedgamma[,2:45])/2)
ggplot(newdf, aes(x=tissues, y=counts)) + geom_point()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

gam1_uncon = one$gamma1[unconverged_one,]
gam2_uncon = one$gamma2[unconverged_one, ]

for (i in 1:length(unconverged_one)){
  if(max(abs(as.numeric(gam1_uncon[i, 2:45]) - as.numeric(gam2_uncon[i, 2:45])))>0.5){print(i)}
}
tmpname = one$finalgamma$genename[unconverged_one[28]]
load(paste0('/Volumes/im-lab/nas40t2/tae/differentialeQTLs/multitissue/localAnc/400genes_result1/',tmpname,'Rdata'))
plot_gamma(final$result, tmpname)
plot_beta(final$result, tmpname)

#### DATA SET 2 ####
load('/Volumes/im-lab/nas40t2/tae/differentialeQTLs/multitissue/localAnc/400genes_summarized_result2.Rdata')
two = summarized_result2
rm(summarized_result2)
ggplot(melt(two$finalgamma, id='genename'), aes(x=variable,y=genename,fill=value)) +
  geom_tile() + scale_fill_gradient(high='indianred',low='antiquewhite1')

converged_two = which(two$niters < 1000) #103 out of 125 converged
convergedgamma2 = two$finalgamma[converged_two, ]


##distribution of 'how many tissues' are selected
##since there are NA's, use the ratio
apply(convergedgamma2[,2:45], 1, function(x) mean(x, na.rm=TRUE)/2)

newdf = data.frame(tissue_selection_ratio = as.factor(apply(convergedgamma2[,2:45], 1, function(x) mean(x, na.rm=TRUE)/2)))
ggplot(newdf, aes(x=tissue_selection_ratio)) + geom_bar()

#plot the per-tissue number of selected genes
newdf = data.frame(tissues = colnames(one$p)[-1], counts = colSums(convergedgamma2[,2:45], na.rm=TRUE)/2)
ggplot(newdf, aes(x=tissues, y=counts)) + geom_point()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


