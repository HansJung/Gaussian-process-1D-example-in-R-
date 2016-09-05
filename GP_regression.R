GP_regression = function(Ky,K_star, y){
  Ky.inv = solve(Ky)
  if (det(Ky)>0.1){
    f_mean = t(K_star) %*% Ky.inv %*% y
  }
  else{
    Ky.inv = solve(Ky + 1e-6*diag(length(y)))
    f_mean = t(K_star) %*% Ky.inv %*% y
  }
  return (f_mean)
}