
# selecting tuning parameter s automatically
s_selection <- function(expr1,expr2,ss=seq(0.1,2,0.1),verbose=F){
  p <- dim(expr1)[2]
  temp1 <- pnorm(sqrt(log(p)))
  d <- matrix(0,nr=length(ss),nc=10)
  
  for (i in 1:length(ss)){
    W <- get_diff_W(expr1,expr2,s = ss[i], tri = T)
    for (j in 1:10){
      temp2 <- (1-temp1)*j/10
      nomi <- sum(abs(W)>qnorm(1-temp2))
      deno <- temp2*p*(p-1)
      d[i,j] <- (nomi/deno-1)^2
    }
    if(verbose) print(paste("ss =",ss[i],"/",tail(ss,1),"   d =",sum(d[i,])))
  }
  s.slected <- ss[which.min(rowSums(d))]
  if(verbose) print(paste("selected s =",s.slected))
  return(s.slected)
}

# estimate the different between precision matrixes
get_diff_W <- function(expr1,expr2,s=2,tri=F){
  
  n1 <- dim(expr1)[1]
  n2 <- dim(expr2)[1]
  p <- dim(expr1)[2]
  
  covm1 <- cov(expr1)
  covm2 <- cov(expr2)
  sigma1 <- diag(covm1)
  sigma2 <- diag(covm2)
  expr1.t <- t(expr1)
  expr2.t <- t(expr2)
  
  b1 <- matrix(0,nr=p-1,nc=p)
  b2 <- matrix(0,nr=p-1,nc=p)
  c1 <- matrix(0,nr=n1,nc=p)
  c2 <- matrix(0,nr=n2,nc=p)
  T1 <- matrix(0,nr=p,nc=p)
  T2 <- matrix(0,nr=p,nc=p)
  W <- matrix(0,nr=p,nc=p)
  
  for (i in 1:p){
    y1 <- expr1[,i]
    x1 <- expr1[,-i]
    lam1 <- s*sqrt(sigma1[i]*log(p)/n1)
    temp1 <- glmnet(x1,y1,family = "gaussian",alpha = 1, lambda = lam1)
    b1[,i] <- as.matrix(coefficients(temp1)[-1])
    
    y2 <- expr2[,i]
    x2 <- expr2[,-i]
    lam2 <- s*sqrt(sigma2[i]*log(p)/n2)
    temp2 <- glmnet(x2,y2,family = "gaussian",alpha = 1, lambda = lam2)
    b2[,i] <- as.matrix(coefficients(temp2)[-1])
    
    c1[,i] <- y1-mean(y1) - (x1–rep(1,nrow(x1))%*%t(colMeans(x1))) %*%b1[,i]
    c2[,i] <- y2-mean(y2) - (x2-rep(1,nrow(x2))%*%t(colMeans(x2))) %*%b2[,i]
  }
  
  r1 <- t(c1)%*%c1/n1
  r2 <- t(c2)%*%c2/n2
  s1 <- colMeans(c1^2)
  s2 <- colMeans(c2^2)
  
  for (i in 1:(p-1)){
    for (j in (i+1):p){
      T1[i,j] <- (r1[i,j]+s1[i]*b1[i,j]+s1[j]*b1[j-1,i])/(r1[i,i]*r1[j,j])
      T2[i,j] <- (r2[i,j]+s2[i]*b2[i,j]+s2[j]*b2[j-1,i])/(r2[i,i]*r2[j,j])
      W[i,j]  <- (T1[i,j]-T2[i,j])/
        sqrt(  (1+b1[i,j]^2*r1[i,i]/r1[j,j]) / (r1[i,i]*r1[j,j]*n1) + (1+b2[i,j]^2*r2[i,i]/r2[j,j]) / (r2[i,i]*r2[j,j]*n2) )
    }
  }
  if (!tri) W <- W+t(W)
  return(W)
}

# multiple testing procedure by controlling FPR
get_W_theshold <- function(W,alpha){
  p <- dim(W)[1]
  t.upper <- 2*sqrt(log(p))
  t0 <- abs(W[upper.tri(W)])
  t1 <- t0[which(t0<t.upper)]
  t2 <- sort(t1,decreasing = T)
  temp <- (p^2-p)/2
  
  thes <- NULL
  use.t.upper <- F
  x=NULL
  for(i in 1:length(t2)){
    x[i] <- 2*(1-pnorm(t2[i]))*temp/i
    if (x[i]>=alpha) { 
      if (i>1) {thes<-t2[i-1]}
      if (i==1) {thes<-t.upper; use.t.upper<-T}
      break 
    }
  }
  W_thes <- abs(W)>=thes
  return(list(W_thes=W_thes, thes=thes, use.t.upper=use.t.upper))
}

