source("Kernel_matrix.R")
source("GP_regression.R")
source("Compute_Alpha.R")
source("Compute_gradient.R")


# toy data generating 
set.seed(1234)
n_obs = 8
sampled = (runif(n_obs,min=-10,max=10))
X_train = seq(min(sampled),max(sampled),length=n_obs)
y = runif(n_obs,min=-5,max=5)
X_test = seq(min(X_train)-1,max(X_train)+1,length=300)

# Initial parameter setting 
param = c(1,1,0.1)
Ky = Kernel_matrix(X_train,X_train,param,T)
K_star = Kernel_matrix(X_train,X_test, param,F)
K_starstar = Kernel_matrix(X_test,X_test,param,F)
Alpha = Compute_Alpha(Ky,y)

n_iter = 10
result_gradient = Compute_gradient(param, Ky, Alpha, X_train, y, n_iter)
f_mean = GP_regression(result_gradient$Ky,result_gradient$K_star,y)

plot(X_train,y,col="blue")
lines(X_test,f_mean,col="red")