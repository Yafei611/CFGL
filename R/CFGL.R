#' Condition-adaptive fused graphical lasso
#'
#' The function estimates differences between 2 precision matrixes. Then a multiple testing procedure with false discovery rate control will be applied to determine different entries between 2 matrixes. The testing result will be turned into a binary matrix which can be used as screening matrix for CFGL.
#'
#' @param expr1 A n*p matrix or data frame of normalized gene expression data. The rows correspond samples (n), the columns correspond genes (p). 
#' 
#' @param expr2 The second gene expression data, should be in the same format, size as expr1.
#' 
#' @param s The tuning parameter for matrixes differences estimation, leave it as NULL to automatically select. 
#' 
#' @param s.seq The candidates for s selection.
#' 
#' @param alpha Prespecified level of false discovery rate. A relatively loose criterion is suggested for determines screening matrix.
#' 
#' @param verbose Set verbose to TURE to show details of s selection.
#' 
#' @details Please refer \bold{Yin et.al (2016). Testing differential networks with applications to the detection of gene-gene interactions. Biometrika(2015),pp. 1-20}
#' @export



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
  
  out.JGLS = list(theta = theta , diff_theta_z = diff_theta_z, iters = out.admm$iter)
  if (loglik.trace) {
    out.JGLS$loglik.trace = out.admm$loglik.trace
    out.JGLS$diffval.trace = out.admm$DiffVal.tr
  }
  
  return(out.JGLS)
}






