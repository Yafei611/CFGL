
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
#' @return \code{CFGL} produces a list that contains estimated inverse matrices and other necessary components.
#' @export

theta2rmat <- function(theta,top_edge=NULL,min_edge=0,keep.diag=F,verbose=F){
  temp0 <- NULL
  temp1 <- NULL
  temp2 <- NULL
  
  rmat <- list()
  rmat2 <- list()
  
  for (i in 1:length(theta)) {
    temp <- diag(theta[[i]])
    rmat[[i]] <- -theta[[i]]/sqrt(temp%*%t(temp))
  }
  
  if (!is.null(top_edge)){
    for (i in 1:length(rmat)) {
      temp0 <- rmat[[i]]
      diag(temp0) <- 0
      temp1 <- c(temp1,abs(as.vector(temp0)))
    }
    temp1 <- temp1[which(temp1!=0)]
    min_edge0 <- temp1[which(rank(-temp1,ties.method = "random")==top_edge)]
    if (min_edge0>min_edge) print("Min_edge was overrided")
    min_edge <- max(min_edge0,min_edge)
  }
  
  if (verbose) print(paste("min_edge is",min_edge))
  
  for (i in 1:length(rmat)) {
    temp2 <- rmat[[i]]
    diag(temp2) <- 0
    temp2[abs(temp2)<min_edge] <- 0
    rmat2[[i]] <- temp2
    if (keep.diag) diag(rmat2[[i]]) <- diag(rmat[[i]])
  }
  return(rmat2)
}


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
#' @return \code{CFGL} produces a list that contains estimated inverse matrices and other necessary components.
#' @export


get_sp_net_3t <- function(rmat){
  nm <- matrix(0,dim(rmat[[1]])[1],dim(rmat[[1]])[1])
  nt <- list()
  nt$t1 <- rmat[[1]]
  nt$t2 <- rmat[[2]]
  nt$t3 <- rmat[[3]]
  
  nt$t1s <- nm;lb <- which((nt$t1!=0)&(nt$t2==0)&(nt$t3==0))
  nt$t1s[lb] <- rmat[[1]][lb]
  
  nt$t2s <- nm;lb <- which((nt$t1==0)&(nt$t2!=0)&(nt$t3==0))
  nt$t2s[lb] <- rmat[[2]][lb]
  
  nt$t3s <- nm;lb <- which((nt$t1==0)&(nt$t2==0)&(nt$t3!=0))
  nt$t3s[lb] <- rmat[[3]][lb]
  
  nt$t12s <- nm; lb <- which((nt$t1!=0)&(nt$t2!=0)&(nt$t3==0))
  nt$t12s[lb] <- (abs(rmat[[1]][lb])+abs(rmat[[2]][lb]))/2
  
  nt$t13s <- nm; lb <- which((nt$t1!=0)&(nt$t2==0)&(nt$t3!=0))
  nt$t13s[lb] <- (abs(rmat[[1]][lb])+abs(rmat[[3]][lb]))/2
  
  nt$t23s <- nm; lb <- which((nt$t1==0)&(nt$t2!=0)&(nt$t3!=0))
  nt$t23s[lb] <- (abs(rmat[[2]][lb])+abs(rmat[[3]][lb]))/2
  
  nt$t123s <- nm; lb <- which((nt$t1!=0)&(nt$t2!=0)&(nt$t3!=0))
  nt$t123s[lb] <- (abs(rmat[[1]][lb])+abs(rmat[[2]][lb])+abs(rmat[[3]][lb]))/3
  
  return(nt)
}


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
#' @return \code{CFGL} produces a list that contains estimated inverse matrices and other necessary components.
#' @export


get_sp_net_2t <- function(rmat){
  nm <- matrix(0,dim(rmat[[1]])[1],dim(rmat[[1]])[1])
  nt <- list()
  nt$t1 <- rmat[[1]]
  nt$t2 <- rmat[[2]]
  
  nt$t1s <- nm;lb <- which((nt$t1!=0)&(nt$t2==0))
  nt$t1s[lb] <- rmat[[1]][lb]
  
  nt$t2s <- nm;lb <- which((nt$t1==0)&(nt$t2!=0))
  nt$t2s[lb] <- rmat[[2]][lb]
  
  nt$t12 <- nm; lb <- which((nt$t1!=0)&(nt$t2!=0))
  nt$t12[lb] <- (abs(rmat[[1]][lb])+abs(rmat[[2]][lb]))/2
  
  return(nt)
}


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
#' @return \code{CFGL} produces a list that contains estimated inverse matrices and other necessary components.
#' @export


show_net <- function(mat,gname=c(1:dim(mat)[1])){
  nt <- mat
  colnames(nt) <- gname
  lb <- rowSums(nt|nt)!=0
  nt <- nt[which(lb),which(lb)]
  ntp <- graph_from_adjacency_matrix(nt,weighted = T)
  
  V(ntp)$color = "firebrick1"
  V(ntp)$size = 4
  V(ntp)$label=V(ntp)$name
  V(ntp)$label.cex=0.7
  V(ntp)$label.color="black"
  
  E(ntp)$arrow.size =0
  E(ntp)$color = "gray60"
  E(ntp)$width=1.5
  
  plot(ntp,layout=layout.fruchterman.reingold)
}


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
#' @return \code{CFGL} produces a list that contains estimated inverse matrices and other necessary components.
#' @export


get_top_node <- function(mat,topn,gname){
  temp <- rowSums(abs(mat)>0)
  names(temp) <- gname
  return(sort(temp,decreasing = T)[1:topn])
}

