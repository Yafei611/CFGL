---
title: "Put the title of your vignette here"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Put the title of your vignette here}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---

# Condition-adaptive Fused Graphical Lasso

CFGL is a tool to jointly construct co-expression networks for gene expression profiles from multiple conditions. By using a data-driven approach to capturing condition-specific co-expression patterns, this method is effective in identifying both co-expression patterns that are specific to a condition and that are common across conditions. 

## Installing

The package can be installed from github. R package **devtools** is required.

```
library(devtools)
install_github("Yafei611/CFGL")
library(CFGL)
```

## Example

### Determine a screening matrix
The screening matrix can be determine by using external information. If external information is absent, the screening matrix can be calculated as following:

```
x <- expr
temp <- get_scr_mat(expr1 = x[[1]],expr2 = x[[2]])
scr.mat <- temp$scr.mat
```

### Perform CFGL

```
temp = CFGL(expr, lambda1 = 0.0008, lambda2 = 0.0008, btc.screening = scr.mat)
theta <- temp$theta
```





