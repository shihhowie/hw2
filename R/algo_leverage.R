#' algo_leverage
#'
#' @param A
#' @param b
#' @param n_samples
#'
#' @return
#' @export
#'
#' @examples algo_leverage(A,b,100)
algo_leverage = function(A,b,n_samples){
  A_hat = (A%*%solve(t(A)%*%A)%*%t(A))
  leverage = diag(A_hat)
  prob = leverage/norm(leverage, "2")

  lev_sample = sample(1:n_samples, i, prob=prob)
  x_l = x[lev_sample]
  y_l = y[lev_sample]
  fit_lev = lm(y_l ~ x_l)
  b_lev = fit_lev$coefficients[2]
  a_lev = fit_unif$coefficients[1]

  return(list(a_lev, b_lev))
}
