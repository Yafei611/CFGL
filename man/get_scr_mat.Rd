% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_scr_mat.R
\name{get_scr_mat}
\alias{get_scr_mat}
\title{Determine screening matrix for CFGL by testing differential entries between 2 precision matrixes}
\usage{
get_scr_mat(expr1, expr2, s = NULL, s.seq = seq(0.2, 2, 0.2), alpha = 0.4,
  verbose = F)
}
\arguments{
\item{expr1}{A n*p matrix or data frame of normalized gene expression data. The rows correspond samples (n), the columns correspond genes (p).}

\item{expr2}{The second gene expression data, should be in the same format, size as expr1.}

\item{s}{The tuning parameter for matrixes differences estimation, leave it as NULL to automatically select.}

\item{s.seq}{The candidates for s selection.}

\item{alpha}{Prespecified level of false discovery rate. A relatively loose criterion is suggested for determines screening matrix.}

\item{verbose}{Set verbose to TURE to show details of s selection.}
}
\description{
The function estimates differences between 2 precision matrixes. Then a multiple testing procedure with false discovery rate control will be applied to determine different entries between 2 matrixes. The testing result will be turned into a binary matrix which can be used as screening matrix for CFGL.
}
\details{
Please refer \bold{Yin et.al (2016). Testing differential networks with applications to the detection of gene-gene interactions. Biometrika(2015),pp. 1-20}
}
