Compute_Alpha = function(Ky,y){
  Alpha = solve(Ky) %*% y
  return (Alpha)
}