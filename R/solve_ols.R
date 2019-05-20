library(doParallel)
library(magrittr)

error = function(x){
  err = norm(x-v, "2")/norm(v, "2")
  return(err)
}

#' Solve OLS
#' @description Solve ordinary least square using iterative methods
#'
#'
#' @param A design matrix
#' @param b target vector
#' @param iter number of iteration for the iterative method
#' @param method Jacobi(Jacobi), or Gauss-Seidal(GS)
#' @param n_cores Number of parallel processes
#'
#' @return list(x, run_time)
#' @export
#'
#' @examples
#' A = diag(2, 10, 10)
#' b = 1:10

#' Jacobi_res = solve_ols(A,b,method="Jacobi",n_cores=2)
#' GS_res = solve_ols(A,b,method="GS")
#' print(Jacobi_res[1])
#' print(GS_res[1])
solve_ols = function(A,b,method="Jacobi",n_cores=1,iter=10){
  n = dim(A)[1]
  x = rep(0, n)
  times = NULL
  ptm <- proc.time()
  if(method=="Jacobi"){
    no_cores = n_cores
    if(n_cores > detectCores())  no_cores = detectCores()
    cluster = makeCluster(no_cores)
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
      t = proc.time()-ptm
      times = c(times, t[3])
    }

  }

  if(method=="GS"){
    for(k in c(1:iter)){
      for(i in 1:n){
        temp = 0
        for(j in 1:n){
          if(j!=i){
            temp = temp + A[i,j]*x[j]
          }
          x[i] = 1/A[i,i]*(b[i]-temp)
        }
      }
      t = proc.time()-ptm
      times = c(times, t[3])
    }
  }

  return(list(x, times))
}

