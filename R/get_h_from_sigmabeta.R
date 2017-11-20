get_h_from_sigmabeta = function(X,sigmabeta,Sigma,gamma,n,T){
  Sigma2 = Sigma
  Sigma2[gamma==0, ] = Sigma2[, gamma==0] = 0
  num = sigmabeta * sum(diag(Sigma2)) * sum(X^2)/n
  denom = num + sum(diag(Sigma))
  return(num/denom)
}
