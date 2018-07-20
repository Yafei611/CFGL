# Condition-adaptive Fused Graphical Lasso

Condition-adaptive fused graphical lasso (CFGL) is a data-driven approach to incorporate condition specificity in the estimation of co-expression networks. More details can be found in  https://www.biorxiv.org/content/early/2018/03/28/290346


The CFGL package is an implementation of condition-adaptive fused graphical lasso method. It allows you:

- Build co-expression network for multiple biological conditions simultaneously
- Build tissue-specific or shared networks
- Visualize co-expression network and get hub genes

This document introduces you an example of CFGL application. It shows you how to build a co-expression network with CFGL using multi-condition expression data. Also, several functions for co-expression network exploratory analysis and visualization are presented.

## Install the package

The package can be installed from Github. R package **devtools** is required. By specifying `build_vignettes = T` to create a vignettes for the CFGL. Note that the CFGL depends two other packages: **glmnet** and **igraph**.

```{r, message = FALSE}
library(devtools)
#install_github("Yafei611/CFGL", build_vignettes = T)
library(CFGL)
```

## Input data

The input data for the CFGL should be a list of multiple matrices. Each matrix corresponding an expression data set for one biological condition.  The row of the matrix corresponds to the samples while the column corresponds to the features (genes/transcripts). Also, the data should satisfy the following requirements:

- The matrix columns must be a same set of features across conditions and in the same sequence.
- The expression data should be approximately normally distributed. It can be either array intensity or log scaled read counts.
- The expression data should be properly cleaned and normalized before the analysis.

In the example, we use rat expresion data for brain and heart. The data contain 19 rat strain for each tissue and 100 gene are randomly selected for the analysis.

```{r,include = T}
x <- expr                       # Get the preloaded data
gname <- colnames(x[[1]])       # Get the name of genes
str(x)                          # Check the data
```

## Determine the screening matrix

CFGL allows users to provide external information of condition-specificity which named as screening matrix. The screening matrix is a binary matrix controls whether similarity should or should not be encouraged between each pair of condition for each edge. 

The screening matrix can be obtained by using the external information from public databases, like KEGG, COXPRESdb or MSigDB depend on the research question and expression dataset. When such information is not available, we also provided a function `get_scr_mat` to calculate a screening matrix.

In this example, we use the `get_scr_mat` to determine a screening matrix.

```{r}
temp <- get_scr_mat(expr1 = x[[1]],expr2 = x[[2]])
scr.mat <- temp$scr.mat
```
CFGL can work without screening matrix. When `scr.mat` is absent, CFGL is equavalent to Fused Graphical Lasso (FGL). 

## Perform CFGL

Now, we perform network analysis using function `CFGL`. CFGL estimates the conditional dependency (precision matrix) among features. To run 'CFGL', one need to specify two tuning parameter `lam1` and `lam2` and a screening matrix (optional).

```{r}
lam1=0.0006
lam2=0.0004
temp = CFGL(expr, lambda1 = lam1, lambda2 = lam2, btc.screening = scr.mat)
theta <- temp$theta
```

`CFGL` returns a list of objects, and the `theta` is the estimated precision matrices which describe dependency among features. Since the example data contains two set of expression data for heart and brain, the `theta` is a list of two matrices correspondingly.

Then, we calculated the conditional correlation based on the `theta`. The conditional correlation quantifies the conditional dependence between the features.

```{r}
rmat <- theta2rmat(theta) 
```

## Tunning parameters selection

Runing CFGL analysis involves tuning parameter seceltion. We suggest to perfrom CFGL with a series tunning paraters and select result according to minimal BIC certieron.



