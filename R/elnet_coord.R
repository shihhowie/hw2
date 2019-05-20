#' elnet coord
#'
#' @param x
#' @param y
#' @param lambda
#' @param alpha
#' @param iter
#'
#' @return
#' @export
#'
#' @examples elnet_coord(x,y,0.1,0.5,20)
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
