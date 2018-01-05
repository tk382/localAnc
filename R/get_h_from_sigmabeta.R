get_h_from_sigmabeta = function(X,
                                sigmabeta,
                                Sigma,
                                gamma,
                                n,
                                T){
  ds = diag(Sigma)
  ind = which(gamma==1)
  num = sum(X^2)/n * sum(ds[ind]) * sigmabeta
  denom = num + sum(ds)
  return(num/denom)
}
