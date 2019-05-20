#' algo_leverage
#'
#' @description Use leverage sampling to solve least square problem.  Returns the intercept and slope
#'
#'
#' @param A Design matrix
#' @param b Target vector
#' @param n_samples Number of samples to use for leverage sampling.
#'
#' @return
#' @export
#'
#' @examples
#' A = rt(100,6)
#' e = rnorm(100, 0,1)
#' b = -A+e
#' res = algo_leverage(A,b,100)
#' print(res)

algo_leverage = function(A,b,n_samples){
  A_hat = (A%*%solve(t(A)%*%A)%*%t(A))
  leverage = diag(A_hat)
  prob = leverage/norm(leverage, "2")
  n = length(b)

  lev_sample = sample(1:n, n_samples, prob=prob)
  A_l = A[lev_sample]
  b_l = b[lev_sample]
  fit_lev = lm(b_l ~ A_l)
  b_lev = fit_lev$coefficients[2]
  a_lev = fit_lev$coefficients[1]

  return(list(a_lev, b_lev))
}
