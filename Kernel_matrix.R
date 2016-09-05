Kernel_matrix = function(X1,X2,param,same){
  # case of 1D 
  len_1 = length(X1)
  len_2 = length(X2)
  KM = matrix(0,length(X1),length(X2))
  
  l = param[1]
  sigma_f = param[2] 
  sigma_y = param[3]
  
  for (i in 1:len_1){
    xp = X1[i]
    for (j in 1:len_2){
      xq = X2[j]
      if ( (i ==j) & (same==T) ){
        KM[i,j] = sigma_f^2* exp(-1/(2*(l^2)) * (xp-xq)^2 ) + sigma_y^2
      }
      else{
        KM[i,j] = sigma_f^2* exp(-1/(2*(l^2)) * (xp-xq)^2 )
      }
    }
  }
  
  return (KM)
  
}