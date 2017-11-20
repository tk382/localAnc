get_sigmabeta_from_h = function(h, gamma, Sigma, X, T){
  n = length(X)
  num = h * sum(diag(Sigma))
  ds = diag(Sigma)
  denom = (1-h)*sum(ds[gamma!=0]) * sum(X^2)/n
  return(num/denom)
}

