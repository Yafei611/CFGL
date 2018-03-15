#' Condition-adaptive fused graphical lasso
#'
#' The function jointly construct gene co-expression network for multiple class using Condition-adaptive Fused Graphical Lasso. Pairwise screening matrics are required to adjust between-condition lasso penalty.
#'
#' @param Y A list expression data which are n*p matrices. all matrices should have a same n and p.
#' 
#' @param lambda1 The tuning parameter for the graphical lasso penalty.
#' 
#' @param lambda2 The tuning parameter for the between condition group lasso penalty.
#' 
#' @param btc.screening A list of screening matrices (p*p) for between condition penalty. Can be obtained using the function \code{get_scr_mat}. When setting as NULL, the function will perform a standard fused graphical lasso.
#' 
#' @param penalize.diag Binary variables that determine whether lambda1 and lambda2 are applied to the diagonal of inverse matrices. 
#' 
#' @param weight Experimental features that assigning weights to each class. Leaving it as default (NULL) is suggested.
#' 
#' @param rho Step size parameter for ADMM algorithm. Large values decrease the step size.
#' 
#' @param rho.increment Adjustment for rho. In each ADMM iteration, rho will be updated as rho=rho*rho.increment.
#' 
#' @param maxiter The maximum number of ADMM interactions.
#' 
#' @param tol The criterion for ADMM convergence.
#' 
#' @param truncate All value in the estimated inverse convenience below this number will be set to 0.
#' 
#' @param loglik.trace Store trace of the likelihood of estimation in each iteration.
#' 
#' @return \code{CFGL} produces a list that contains estimated inverse matrices and other necessary components.
#' \itemize{
#'  \item{$theta} {The estimation of inverse matrices}
#'  \item{$iters} {The numebr of ADMM iterations}
#'  \item{$loglik.trace} {Trace of log-likelihood}
#' }
#' @details Please refer \bold{An adaptive procedure for inferring condition-specific gene co-expression network }
#' @export
#' @examples
#' x = expr
#' plot(x[,1],x[,2],xlim=c(-8,8),ylim=c(-8,8),cex = .4);
#' idr.out = IDR.3component(x = x)

CFGL <- function(Y, lambda1, lambda2, btc.screening=NULL, penalize.diag = c(TRUE, TRUE),
                 weights=NULL, rho=1, rho.increment=1, maxiter=500, tol=1e-4, truncate=1e-05, loglik.trace=FALSE){
  K <- length(Y)
  p <- c()
  n <- c()
  S <- list()
  for (k in 1:K) {
    p[k] <- dim(Y[[k]])[2]
    n[k] <- dim(Y[[k]])[1]
    S[[k]] <- cov(Y[[k]]) * (n[k] - 1)/n[k]
  }
  
  if (!(K==2||K==3)) stop("K must be 2 or 3 in this version")
  if (p[1]!=p[2]) stop("p must be the same for each k")
  
  p <- p[1]
  if (is.null(weights)) weights <- rep(1, K)
  
  lam1.m <- get_lam_mat(lambda1, p, penalize.diagonal = penalize.diag[1])
  lam2.m0 <- get_lam_mat(lambda2, p, penalize.diagonal = penalize.diag[2])
  
  if ((K==2)&&(!is.null(btc.screening))) lam2.m = lam2.m0 * btc.screening else lam2.m = lam2.m0 
  
  if (K==3){
    if (!is.null(btc.screening)){
      lam2.m <- list()
      lam2.m[[1]] <- btc.screening[[1]] * lam2.m0 #t12
      lam2.m[[2]] <- btc.screening[[2]] * lam2.m0 #t13
      lam2.m[[3]] <- btc.screening[[3]] * lam2.m0 #t23
    } 
    else{
      lam2.m <- list()
      lam2.m[[1]] <- lam2.m0 #t12
      lam2.m[[2]] <- lam2.m0 #t13
      lam2.m[[3]] <- lam2.m0 #t23
    }
  }
  
  # admm start
  out.admm = admm.iter(S = S,n = n,lam1.m = lam1.m,lam2.m = lam2.m, weights = weights,
                       rho = rho,rho.increment = rho.increment,
                       maxiter = maxiter,tol = tol,loglik.trace = loglik.trace)
  
  # assumed theta == Z, see the difference
  diff_theta_z = 0
  for (k in 1:K) {
    diff_theta_z = diff_theta_z + sum(abs(out.admm$theta[[k]] - out.admm$Z[[k]]))
  }
  
  # round down
  theta = list()
  for (k in 1:K) {
    rounddown = abs(out.admm$Z[[k]]) < truncate
    diag(rounddown) = FALSE
    theta[[k]] = out.admm$Z[[k]] * (1 - rounddown)
  }
  
  out.CFGL= list(theta = theta ,iters = out.admm$iter)
  if (loglik.trace) {
    out.V$loglik.trace = out.admm$loglik.trace
  }
  
  return(out.CFGL)
}






