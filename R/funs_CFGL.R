
centering_data <- function(dat){
  dat.ct = list()
  for (k in 1:length(dat)) {
    dat.ct[[k]] = matrix(0,nr=dim(dat[[k]])[1],nc=dim(dat[[k]])[2])
    dat.ct[[k]] = apply(X = dat[[k]], MARGIN = 2, function(x) x-mean(x))
  }
  return(dat.ct)
}

get_lam_mat <- function (lambda, p, penalize.diagonal) {
  if (is.matrix(lambda)) {
    if (sum(lambda != t(lambda)) > 0) stop("error: penalty matrix is not symmetric")
    if (sum(abs(dim(lambda) - p)) != 0) stop("error: penalty matrix has wrong dimension")
  }
  if (length(lambda) == 1) lambda = matrix(lambda, p, p)
  if (!penalize.diagonal) diag(lambda) = 0
  return(lambda)
}

cal_loglik <- function(theta,S,n,lam1.m,lam2.m){
  loglik = 0
  K = length(theta)
  for (k in 1:K){
    loglik = loglik + n[k] * log(det(theta[[k]])) - n[k] * sum(diag(S[[k]] %*% theta[[k]])) - sum(lam1.m * abs(theta[[k]]))
  }
  if (K==2) loglik = loglik - sum(lam2.m*abs(theta[[1]]-theta[[2]]))
  if (K==3){
    loglik = loglik - sum(lam2.m[[1]]*abs(theta[[1]]-theta[[2]]))
    - sum(lam2.m[[2]]*abs(theta[[1]]-theta[[3]]))
    - sum(lam2.m[[3]]*abs(theta[[2]]-theta[[3]]))
  }
  return(loglik)
}

soft_thresholding <- function (a, lam.m) {
  out <- sign(a) * pmax(0, abs(a) - lam.m)
  return(out)
}

cal_z_L1_t2 <- function (A, rho, lam1.m, lam2.m){ 
  S1 <- abs(A[[1]] - A[[2]]) <= 2 * lam2.m/rho
  Z1_1 <- (A[[1]] + A[[2]])/2
  Z2_1 <- Z1_1
  S2 <- (A[[1]] > A[[2]] + 2 * lam2.m/rho)
  Z1_2 <- A[[1]] - lam2.m/rho
  Z2_2 <- A[[2]] + lam2.m/rho
  S3 <- (A[[2]] > A[[1]] + 2 * lam2.m/rho)
  Z1_3 <- A[[1]] + lam2.m/rho
  Z2_3 <- A[[2]] - lam2.m/rho
  Z1 <- soft_thresholding(a = S1 * Z1_1 + S2 * Z1_2 + S3 * Z1_3, lam.m = lam1.m/rho)
  Z2 <- soft_thresholding(a = S1 * Z2_1 + S2 * Z2_2 + S3 * Z2_3, lam.m = lam1.m/rho)
  return(list(Z1, Z2))
}

cal_z_L1_t3 <- function(A, rho, lam1.m, lam2.m){
  p <- dim(A[[1]])[1]
  vl <- dim(A[[1]])[1]^2
  
  # turn matrix list to vector
  A_vec <- cbind(as.vector(A[[1]]),as.vector(A[[2]]),as.vector(A[[3]])) #1,2,3
  P_vec <- cbind(as.vector(lam2.m[[3]]),as.vector(lam2.m[[2]]),as.vector(lam2.m[[1]])) / rho  # 23,13,12
  
  condi <- rep(0,vl)
  Z_vec <- matrix(0,nc=3,nr=vl)
  
  # cond 1 z1=z2=z3
  fcd.temp <- find_condi_1(av = A_vec,pv = P_vec )
  condi[fcd.temp$id] <- 1
  Z_vec[fcd.temp$id,] <- fcd.temp$zv
  
  # cond 2 z1>z2>z3
  index1 <- matrix(c(1,2,3,1,3,2,
                     2,1,3,2,3,1,
                     3,1,2,3,2,1),nc=3,byrow = T)
  for (indexi in 1:6){
    inputi <- which(condi==0)
    fcd.temp <- find_condi_2(av = A_vec[inputi,],pv = P_vec[inputi,], ind = index1[indexi,] )
    if (sum(condi[inputi][fcd.temp$id]>0)>0) stop("error conflect condition")  # for test
    condi[inputi][fcd.temp$id] <- 2
    Z_vec[inputi[fcd.temp$id],] <- fcd.temp$zv
  }
  
  # cond 3 z1=z2>z3
  index2 <- matrix(c(1,2,3,2,3,1,3,1,2),nc=3,byrow = T)
  for (indexi in 1:3){
    inputi <- which(condi==0)
    fcd.temp <- find_condi_3(av = A_vec[inputi,],pv = P_vec[inputi,], ind = index2[indexi,] )
    if (sum(condi[inputi][fcd.temp$id]>0)>0) stop("error conflect condition")  # for test
    condi[inputi][fcd.temp$id] <- 3
    Z_vec[inputi[fcd.temp$id],] <- fcd.temp$zv
  }
  
  # cond 4 z1=z2<z3
  for (indexi in 1:3){
    inputi <- which(condi==0)
    fcd.temp <- find_condi_4(av = A_vec[inputi,],pv = P_vec[inputi,], ind = index2[indexi,] )
    if (sum(condi[inputi][fcd.temp$id]>0)>0) stop("error conflect condition")  # for test
    condi[inputi][fcd.temp$id] <- 4
    Z_vec[inputi[fcd.temp$id],] <- fcd.temp$zv
  }
  
  if (sum(condi==0)) stop("error condi==0")
  
  Z <- list(matrix(Z_vec[,1],nc=p,byrow = F),
            matrix(Z_vec[,2],nc=p,byrow = F),
            matrix(Z_vec[,3],nc=p,byrow = F))
  
  Z[[1]] <- soft_thresholding(Z[[1]], lam.m = lam1.m/rho)
  Z[[2]] <- soft_thresholding(Z[[2]], lam.m = lam1.m/rho)
  Z[[3]] <- soft_thresholding(Z[[3]], lam.m = lam1.m/rho)
  return(Z)
}

find_condi_1 <- function(av,pv){
  av <- matrix(av,ncol=3)
  pv <- matrix(pv,ncol=3)
  z <- (av[,1]+av[,2]+av[,3])/3
  c_temp <- which( (z<=av[,1]+pv[,2]+pv[,3])&(z>=av[,1]-pv[,2]-pv[,3])&
                     (z<=av[,2]+pv[,1]+pv[,3])&(z>=av[,2]-pv[,1]-pv[,3])&
                     (z<=av[,3]+pv[,1]+pv[,2])&(z>=av[,3]-pv[,1]-pv[,2]))
  zv <- cbind(z[c_temp],z[c_temp],z[c_temp])
  return(list(zv=zv,id=c_temp))
}

find_condi_2 <- function(av,pv,ind){
  av <- matrix(av,ncol=3)
  pv <- matrix(pv,ncol=3)
  c_temp <- which(((av[,ind[1]]-pv[,ind[2]]-pv[,ind[3]])>(av[,ind[2]]+pv[,ind[3]]-pv[,ind[1]])) & 
                    ((av[,ind[2]]+pv[,ind[3]]-pv[,ind[1]])>(av[,ind[3]]+pv[,ind[1]]+pv[,ind[2]]))==T)
  zv <- cbind(av[c_temp,ind[1]]-pv[c_temp,ind[2]]-pv[c_temp,ind[3]],
              av[c_temp,ind[2]]+pv[c_temp,ind[3]]-pv[c_temp,ind[1]],
              av[c_temp,ind[3]]+pv[c_temp,ind[1]]+pv[c_temp,ind[2]])
  zv[,c(ind[1],ind[2],ind[3])] <- zv
  return(list(zv=zv,id=c_temp))
}

find_condi_3 <- function(av,pv,ind){
  av <- matrix(av,ncol=3)
  pv <- matrix(pv,ncol=3)
  c_temp <- which( (abs(av[,ind[1]]-av[,ind[2]]-pv[,ind[2]]+pv[,ind[1]])<=2*pv[,ind[3]]) &
                     ((av[,ind[1]]+av[,ind[2]]-pv[,ind[2]]-pv[,ind[1]])/2>=av[,ind[3]]+pv[,ind[1]]+pv[,ind[2]]) )
  temp <- (av[c_temp,ind[1]]+av[c_temp,ind[2]]-pv[c_temp,ind[2]]-pv[c_temp,ind[1]])/2
  
  zv <- cbind(temp,
              temp,
              av[c_temp,ind[3]]+pv[c_temp,ind[1]]+pv[c_temp,ind[2]])
  zv[,c(ind[1],ind[2],ind[3])] <- zv
  return(list(zv=zv,id=c_temp))
}

find_condi_4 <- function(av,pv,ind){
  av <- matrix(av,ncol=3)
  pv <- matrix(pv,ncol=3)
  c_temp <- which( (abs(av[,ind[1]]-av[,ind[2]]+pv[,ind[2]]-pv[,ind[1]])<=2*pv[,ind[3]]) &
                     ((av[,ind[1]]+av[,ind[2]]+pv[,ind[2]]+pv[,ind[1]])/2<=av[,ind[3]]-pv[,ind[1]]-pv[,ind[2]]) )
  temp <- (av[c_temp,ind[1]]+av[c_temp,ind[2]]+pv[c_temp,ind[2]]+pv[c_temp,ind[1]])/2
  
  zv <- cbind(temp,
              temp,
              av[c_temp,ind[3]]-pv[c_temp,ind[1]]-pv[c_temp,ind[2]])
  zv[,c(ind[1],ind[2],ind[3])] <- zv
  return(list(zv=zv,id=c_temp))
}

admm.iter <- function(S,n,lam1.m,lam2.m,weights=NULL,rho=1,rho.increment=1,maxiter=500,tol=1e-05,loglik.trace=FALSE){
  K <- length(S)
  p <- dim(S[[1]])[1]
  
  theta <- list()
  Z <- list()
  U <- list()
  for (k in 1:K) {
    theta[[k]] <- diag(1/diag(S[[k]]))
    Z[[k]] <- matrix(0, p, p)
    U[[k]] <- matrix(0, p, p)
  }
  
  iter <- 0
  diff_value <- 1
  loglik.tr <- rep(0,maxiter)
  DiffVal.tr <- rep(0,maxiter)
  
  while (iter < maxiter && diff_value > tol) {
    
    theta.prev <- theta
    iter <- iter + 1
    if (loglik.trace) loglik.tr[iter] <- cal_loglik(theta = Z,S = S,n = n,lam1.m = lam1.m,lam2.m = lam2.m)
    
    # step a update theta
    for (k in 1:K) {
      edecomp <- eigen(S[[k]] - rho * Z[[k]]/weights[k] + rho * U[[k]]/weights[k])
      D <- edecomp$values
      V <- edecomp$vectors
      D2 <- weights[k]/(2 * rho) * (-D + sqrt(D^2 + 4 * rho/weights[k]))
      theta[[k]] <- V %*% diag(D2) %*% t(V)
    }
    
    # step b update Z
    A <- list()
    for (k in 1:K) {
      A[[k]] <- theta[[k]] + U[[k]]
    }
    if (K == 2) Z <- cal_z_L1_t2(A, rho, lam1.m, lam2.m)
    if (K == 3) Z <- cal_z_L1_t3(A, rho, lam1.m, lam2.m)
    
    # step c update U
    for (k in 1:K) {
      U[[k]] <- U[[k]] + (theta[[k]] - Z[[k]])
    }
    
    # difference between est
    diff_value <- 0
    for (k in 1:K) {
      diff_value <- diff_value + sum(abs(theta[[k]] - theta.prev[[k]]))/sum(abs(theta.prev[[k]]))
    }
    
    if (loglik.trace) DiffVal.tr[iter] <- diff_value
    
    rho <- rho * rho.increment
  }
  
  out.admm = list(theta = theta, Z = Z, iters = iter)
  if (loglik.trace) {
    out.admm$loglik.trace = loglik.tr[1:iter]
    out.admm$DiffVal.trace = DiffVal.tr[1:iter]
  }
  return(out.admm)
}





