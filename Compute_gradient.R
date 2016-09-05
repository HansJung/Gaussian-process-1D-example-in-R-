Compute_gradient = function(param, Ky, Alpha, X_train, y, n_iter){
  source("Kernel_matrix.R")
  source("GP_regression.R")
  source("Compute_Alpha.R")
  
  # Construct the gradient matrix 
  ## initial val. 
  l = param[1]
  sigma_f = param[2] 
  sigma_y = param[3]
  
  Grad_l = matrix(0,dim(Ky)[1],dim(Ky)[2]) 
  Grad_f = matrix(0,dim(Ky)[1],dim(Ky)[2]) 
  Grad_y = matrix(0,dim(Ky)[1],dim(Ky)[2]) 
  
  len_x = length(X_train)
  Ky.inv = solve(Ky)
  
  for (iter_idx in 1:n_iter){
    Comp_1 = Alpha %*% t(Alpha) - Ky.inv
    for (i in 1:len_x){
      for (j in 1:len_x){
        part_1 = exp( -0.5*(1/(l^2))*((X_train[i] - X_train[j])^2))
        Grad_l[i,j] = (sigma_f^2)*part_1*((1/l)^3)*((X_train[i] - X_train[j])^2)
        Grad_f[i,j] = 2*sigma_f*part_1
        if (i==j){
          Grad_y[i,j] = 2*sigma_y
        }
      }
    }
    
    Comp_l = Comp_1 %*% Grad_l
    Comp_f = Comp_1 %*% Grad_f
    Comp_y = Comp_1 %*% Grad_y
  
    nll_l = -0.5*sum(diag(Comp_l))
    nll_f = -0.5*sum(diag(Comp_f))
    nll_y = -0.5*sum(diag(Comp_y))
    
    if (l - 0.01*nll_l > 0){
      l = l - 0.01*nll_l
    }
    if (sigma_f - 0.01*nll_f > 0){
      sigma_f = sigma_f - 0.01*nll_f 
    }
    if (sigma_y - 0.01*nll_y > 0){
      sigma_y = sigma_y - 0.01*nll_y
    }
    param = c(l,sigma_f,sigma_y)
    Ky = Kernel_matrix(X_train,X_train,param,T)
    if (det(Ky) > 0.1){
      Ky.inv = solve(Ky)
    }
    else{
      Ky.inv = solve(Ky + 1e-6*diag(length(y)))
    }
    K_star = Kernel_matrix(X_train,X_test, param,F)
    Alpha = Compute_Alpha(Ky,y)
  }
  return (list(param=param,Ky=Ky,K_star=K_star,Alpha=Alpha))
  

}