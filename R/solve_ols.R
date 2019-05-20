
library(doParallel)
library(magrittr)

#' Solve OLS
#'
#' @param A
#' @param b
#' @param iter
#'
#' @return
#' @export
#'
#' @examples solve_ols(A,b,20)
solve_ols = function(A,b,iter){
  n = dim(A)[1]
  x = rep(0, n)
  errs = NULL
  times = NULL
  ptm <- proc.time()
  clusterExport(cluster, varlist = c("A", "b", "n"), environment())
  for(k in 1:iter){
    x = parSapply(cluster, 1:n, function(i){
      temp = 0
      for(j in 1:n){
        if(j!=i){
          temp = temp + A[i,j]*x[j]
        }
      }
      res = 1/A[i,i]*(b[i]-temp)
      res
    })
    errs = c(errs, error(x))
    t = proc.time()-ptm
    times = c(times, t[3])
  }
  return(list(x, errs, times))
}
