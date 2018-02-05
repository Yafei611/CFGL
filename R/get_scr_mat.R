#' Determine screening matrix for CFGL by testing differential entries between 2 precision matrixes
#'
#' The function estimates differences between 2 precision matrixes. Then a multiple testing procedure with false discovery rate control will be applied to determine different entries between 2 matrixes. The testing result will be turned into a binary matrix which can be used as screening matrix for CFGL.
#'
#' @import glmnet
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

get_scr_mat <- function(expr1,expr2,s=NULL,s.seq=seq(0.2,2,0.2),alpha=0.4,verbose=F){
  if (is.null(s)){
    s.selected <- s_selection(expr1,expr2,ss = s.seq,verbose = verbose)
    if (verbose) print("S was not assigned, will be auto-selected")
  }
  W <- get_diff_W(expr1 = expr1,expr2 = expr2,s = s.selected)
  temp <- get_W_theshold(W,alpha = alpha)
  W.diff.m <- !temp$W_thes
  return(list(scr.mat=W.diff.m,s=s.selected,t.threshold=temp$thes))
}



