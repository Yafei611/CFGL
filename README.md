# Condition-adaptive Fused Graphical Lasso

Condition-adaptive fused graphical lasso (CFGL) is a data-driven approach to incorporate 
condition specificity in the estimation of co-expression networks. More details can be found in 
 https://www.biorxiv.org/content/early/2018/03/28/290346


The CFGL package is an implementation of condition-adaptive fused graphical lasso method. It 
allows you:

- Build co-expression network for multiple biological conditions simultaneously
- Build tissue-specific or shared networks
- Visualize co-expression network and get hub genes

This document introduces you an example of CFGL application. It shows you how to build a co-
expression network with CFGL using multi-condition expression data. Also, several functions for 
co-expression network exploratory analysis and visualization are presented.

## Install the package

The package can be installed from Github. R package **devtools** is required. By specifying 
`build_vignettes = T` to create a vignettes for the CFGL. Note that the CFGL depends two other 
packages: **glmnet** and **igraph**.

```{r, message = FALSE,eval = FALSE}
library(devtools)
install_github("Yafei611/CFGL", build_vignettes = T)
```

Before performing analysis, we also need to load the package.

```{r}
library(CFGL)
```

## Input data

The input data for the CFGL should be a list of multiple matrices. Each matrix corresponding an 
expression data set for one biological condition.  The row of the matrix corresponds to the 
samples while the column corresponds to the features (genes/transcripts). Also, the data should 
satisfy the following requirements:

- The matrix columns must be a same set of features across conditions and in the same sequence.
- The expression data should be approximately normally distributed. It can be either array 
intensity or log scaled read counts.
- The expression data should be properly cleaned and normalized before the analysis.

In the example, we use rat expresion data for brain and heart. The data contain 19 rat strain 
for each tissue and 100 gene are randomly selected for the analysis.

```{r,include = T}
x <- expr                       # Get the preloaded data
gname <- colnames(x[[1]])       # Get the name of genes
str(x)                          # Check the data
```

## Determine the screening matrix

CFGL allows users to provide external information of condition-specificity which named as 
screening matrix. The screening matrix is a binary matrix controls whether similarity should or 
should not be encouraged between each pair of condition for each edge. 

The screening matrix can be obtained by using the external information from public databases, 
like KEGG, COXPRESdb or MSigDB depend on the research question and expression dataset. When 
such information is not available, we also provided a function `get_scr_mat` to calculate a 
screening matrix.

In this example, we use the `get_scr_mat` to determine a screening matrix.

```{r}
temp <- get_scr_mat(expr1 = x[[1]],expr2 = x[[2]])
scr.mat <- temp$scr.mat
```
CFGL can work without screening matrix. When `scr.mat` is absent, CFGL is equavalent to Fused 
Graphical Lasso (FGL). 

## Perform CFGL

Now, we perform network analysis using function `CFGL`. CFGL estimates the conditional 
dependency (precision matrix) among features. To run 'CFGL', one need to specify two tuning 
parameter `lam1` and `lam2` and a screening matrix (optional).

```{r}
lam1=0.0006
lam2=0.0004
temp = CFGL(expr, lambda1 = lam1, lambda2 = lam2, btc.screening = scr.mat)
theta <- temp$theta
```

`CFGL` returns a list of objects, and the `theta` is the estimated precision matrices which 
describe dependency among features. Since the example data contains two set of expression data 
for heart and brain, the `theta` is a list of two matrices correspondingly.

Then, we calculated the conditional correlation based on the `theta`. The conditional 
correlation quantifies the conditional dependence between the features.

```{r}
rmat <- theta2rmat(theta) 
```

## Tunning parameters selection

Runing CFGL analysis involves tuning parameter seceltion. We suggest to perfrom CFGL with a 
series tunning paraters and select result according to minimal BIC certieron.


## Network visualization

One can visualize the constructed network using R package **igraph**. Here, we show the codes 
that visualize networking using **igraph**. Now, let check the network for heart tissue.


```{r, message = FALSE, fig.width=7, fig.height=4}
library("igraph")

nt <- rmat[[1]]
colnames(nt) <- gname

lb <- rowSums(nt|nt)!=0              
nt <- nt[which(lb),which(lb)]                           # removing unlinked nodes

ntp <- graph_from_adjacency_matrix(nt,weighted = T)     # turning the matrix into network

V(ntp)$color = "firebrick1"                             # defining attributions of nodes
V(ntp)$size = 4
V(ntp)$label=V(ntp)$name
V(ntp)$label.cex=0.7
V(ntp)$label.color="black"

E(ntp)$arrow.size =0                                    # defining attributions of edges
E(ntp)$color = "gray60"
E(ntp)$width=1.5

par(mar=c(0,0,0,0))
plot(ntp,layout=layout.graphopt)
```

For a quick check, one can also use the wrapper function in **CFGL** package. Now we check the 
network for the brain tissue.

```{r, message = FALSE, fig.width=7, fig.height=4}
show_net(rmat[[2]],gname = gname)
```

Check http://igraph.org/redirect.html to find more information of **igraph**.

## Tissue-specific/-shared network

Once can obtain tissue-specific or shared network using the function `get_sp_net_2t`. This 
function takes a list of two matrices (network) as input and output tissue-specific or shared 
network.

```{r}
rmat.sp <- get_sp_net_2t(rmat) 
names(rmat.sp)
```

- 't1'  : network for tissue 1
- 't2'  : network for tissue 2
- 't1s' :  network contains edges that only appeared in the tissue 1
- 't2s' :  network contains edges that only appeared in the tissue 2
- 't12' :  network contains edges that shared by 2 tissues

Let us check the brain specific network

```{r, message = FALSE, fig.width=7, fig.height=4}
show_net(rmat.sp$t1s,gname)
```

For a list of three networks. One can get tissue-specific network using `get_sp_net_3t`. For a 
list of four or more tissues, the functions for listing all tissue-specific/shared network are 
not provided, since there are too many combinations. Writing your own code to get network 
shared by particular tissues are suggested.


## Get network hubs

Nodes that have many connections are called hubs. Network hubs can be obtained using function 
`get_top_node`. The function will count the connections for each node and lists nodes have the 
most connections.

Here, we list top hubs for brain and heart-specific networks.

```{r}
get_top_node(rmat.sp$t1s,topn = 5,gname)
get_top_node(rmat.sp$t2s,topn = 5,gname)
```

## Contributing
We are continuing to add new features. Any kind of contribution, like bug reports or feature 
requests,is welcome.


