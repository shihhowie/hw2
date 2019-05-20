#' elnet coord
#'
#'@description Use coordinate descent on elastic net
#' @param x design matrix
#' @param y target vector
#' @param lambda penalty weight
#' @param alpha distribution of penalty
#' @param iter number of iterations
#'
#' @return
#' @export
#'
#' @examples
#' library(MASS)
#' n=50
#' p=10
#' beta = c(2,0,-2,0,1,0,-1,0,0,0)
#' sigma = diag(1,p,p)
#' sigma[1,2]=0.8
#' sigma[2,1]=0.8
#' sigma[5,6]=0.8
#' sigma[6,5]=0.8

#' x = mvrnorm(n, mu=rep(0,p),Sigma=sigma)
#' y = x%*%beta + rnorm(n,0,1)

#' elnet_coord(x,y, 4, 0.5, 20)
elnet_coord = function(x,y,lambda, alpha, iter){
  n = length(y)
  p = dim(x)[2]
  b = rnorm(p, 0, 1)
  y = y-mean(y)
  for(l in 1:iter){
    for(k in 1:p){
      U = 0
      V = 2*sum(x[,k]**2)/n+2*lambda*(1-alpha)
      for(i in 1:n){
        temp = 0
        for(j in 1:p){
          if(j!=k) temp = temp+x[i,j]*b[j]
        }
        U = U + (y[i]-temp)*x[i,k]
      }
      U = 2*U/n
      if(U > lambda*alpha){
        b[k]=(U-lambda*alpha)/V
      }else if(U < -lambda*alpha){
        b[k]=(U+lambda*alpha)/V
      }else{
        b[k]=0
      }
    }
  }
  return(b)
}
