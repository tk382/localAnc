update_betagam = function(X, Y, gam1, beta1, sigma, sigmabeta, Vbeta){
  temp = update_gamma(X, Y, gam1, 1)
  gam2 = temp$newgamma
  changeind = temp$changeind
  beta2 = beta1
  beta2 = beta2*gam2
  ind = which(gam2==1)
  beta2[ind] = beta1[ind] + rnorm(length(ind), 0, Vbeta)
  change = gam2[changeind]
  return(list(beta = beta2, gamma = gam2, change = change, changeind = changeind))
}


update_gamma = function(X, Y, gamma, mag){
  newgamma = gamma
  T = length(gamma)
  p1 = p2 = 0.5;
  for (j in 1:mag){
    ind0 = which(gamma==0)
    ind1 = which(gamma!=0)
    if(length(ind0)==T){p1=1; p2=0}else if(length(ind1)==T){p1=0; p2=1}
    case = sample(c(1,2), size=1, prob = c(p1,p2))
    if (case==1){
      marcor = abs(colMeans(Y*X, na.rm=TRUE)[ind0])
      k = length(ind0)
      if(k>1){
        add = sample(ind0, size=1, prob=marcor)
      }else{add=ind0}
      newgamma[add] = 1
      changeind = add
    }else if (case==2){
      remove = sample(ind1, size=1)
      newgamma[remove] = 0
      changeind = remove
    }
  }
  return(list(newgamma=newgamma, changeind=changeind))
}
