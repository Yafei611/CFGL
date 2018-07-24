# Condition-adaptive Fused Graphical Lasso

Condition-adaptive fused graphical lasso (CFGL) is a data-driven approach to incorporate 
condition specificity in the estimation of co-expression networks. For thorough details, see the preprint:
 https://www.biorxiv.org/content/early/2018/03/28/290346


The CFGL package is an implementation of condition-adaptive fused graphical lasso method. With CFGL, you can:

- Build co-expression networks for multiple biological conditions simultaneously
- Build tissue-specific or shared networks
- Visualize co-expression network and get hub genes

This document introduces an example application of CFGL. It demonstrates how to build a co-expression network with CFGL using multi-condition expression data. Also, several functions for co-expression network exploratory analysis and visualization are presented.

## Install the package

The package can be installed from Github. R package `devtools` is required. CFGL dependencies include two other 
packages: `glmnet` and `igraph`.

```{r, message = FALSE, eval = FALSE}
library(devtools)
install_github("Yafei611/CFGL", build_vignettes = TRUE)
library(CFGL)
```


## Input data

The input data for CFGL should be a list of multiple matrices. Each matrix corresponds to an expression data set for one biological condition.  The rows of the matrix correspond to the samples while the columns correspond to the features (genes/transcripts). Also, the data should satisfy the following requirements:

- The matrix columns must be the same set of features across conditions and in the same sequence.
- The expression data should be approximately normally distributed. It can be either array intensity or log scaled read counts.
- The expression data should be properly cleaned and normalized before the analysis.

In the example, we use rat expresion data for brain and heart. The data contain 19 rat strains for each tissue and 100 genes are randomly selected for the analysis.

```{r, include = TRUE, echo = 2:10}
library(CFGL)
# Get the preloaded data
data('expr')
x <- expr            

# Get the name of genes
gname <- colnames(x[[1]])       

# Check the data
str(x)                          
```

## Determine the screening matrix

CFGL allows users to provide external information concerning condition-specificity. This information provides the condition adaptivity of the method, and it contained in something termed the *screening matrix*. The screening matrix is a binary matrix controls whether or not similarity should be encouraged between each pairs of conditions.

The screening matrix can be obtained by using the external information from public databases, like KEGG, COXPRESdb or MSigDB, depending on the user's research question and expression dataset. When external information is not available, we also provide a function `get_scr_mat` to calculate a screening matrix for the specific problem.

In this example, we use `get_scr_mat` to determine a screening matrix.

```{r}
temp <- get_scr_mat(expr1 = x[[1]],expr2 = x[[2]])
scr.mat <- temp$scr.mat
```

CFGL can work without a screening matrix. When `scr.mat` is absent, CFGL is equavalent to Fused Graphical Lasso (FGL). 

## Perform CFGL

Now, we perform network analysis using the function `CFGL`. `CFGL` estimates the conditional dependency (precision matrix) among features. To run `CFGL`, one needs to specify two tuning parameters, `lam1` and `lam2`, and optionally, a screening matrix.

```{r}
lam1 <- 0.0006
lam2 <- 0.0004
temp <- CFGL(expr, lambda1 = lam1, lambda2 = lam2, btc.screening = scr.mat)
theta <- temp$theta
```

`CFGL` returns a list of objects. `theta` is the estimated precision matrices which describe dependency among features. Since the example data contain two sets of expression data 
for heart and brain, the `theta` is a list of two corresponding precision matrices.

Next, we calculate the conditional correlation based on `theta`. The conditional correlation quantifies the conditional dependence between the features.

```{r}
rmat <- theta2rmat(theta) 
```

## Tuning parameter selection

Running the CFGL analysis involves tuning parameter seceltion. We suggest to perform CFGL with a series of tuning parameters and select optimal parameters to minimize Bayesian Information Criterion (BIC).


## Network visualization

Visualization of the constructed network can be done using the R package **igraph**. Here, we show code that accomplishes this.

Let's first check the network constructed for heart tissue.


```{r, message = FALSE, fig.width=7, fig.height=4, warning=FALSE}
library(igraph)

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

For a quick check, one can also use the wrapper function, `show_net`, provided in the `CFGL` package. Now we check the network for the brain tissue.

```{r, message = FALSE, fig.width=7, fig.height=4}
show_net(rmat[[2]], gname = gname)
```

The functionality of `igraph` is extensive. More information is available at http://igraph.org/redirect.html.

## Tissue-specific/-shared network

Once can obtain a tissue-specific or shared network using the function `get_sp_net_2t`. This function takes a list of two matrices (network) as input and outputs a tissue-specific or shared network.

```{r}
rmat.sp <- get_sp_net_2t(rmat) 
names(rmat.sp)
```

- 't1'  : network for tissue 1
- 't2'  : network for tissue 2
- 't1s' :  network contains edges that only appeared in the tissue 1
- 't2s' :  network contains edges that only appeared in the tissue 2
- 't12' :  network contains edges that shared by 2 tissues

Let's check the brain specific network.

```{r, message = FALSE, fig.width=7, fig.height=4}
show_net(rmat.sp$t1s,gname)
```

For a list of three networks, one can get tissue-specific network using `get_sp_net_3t`. For a list of four or more tissues, the functions for listing all tissue-specific/shared network are 
not provided, due to combinatorial challenges.


## Get network hubs

Nodes that have many connections are called hubs. Network hubs can be obtained using function `get_top_node`. The function will count the connections for each node and list nodes which have the most connections.

Here, we list top hubs for brain and heart-specific networks.

```{r}
get_top_node(rmat.sp$t1s,topn = 5, gname)
get_top_node(rmat.sp$t2s,topn = 5, gname)
```

## Contributing
We are continuing to add new features. Any kind of contribution, like bug reports or feature requests, are welcome.

