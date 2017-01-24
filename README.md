---
title: "Eigen vector selection"
author: "Sudeep Sahadevan"
date: "09/17/2015"
output:
 BiocStyle::html_document:
  toc: yes
  theme: united
  highlight: tango
---

Functions to perform informative eigen vector selection based on the algorithm proposed by:  
Tao Xiang and Shaogang Gong (2008). Spectral clustering with eigenvector selection. Pattern Recogn. 41(3), 1012-1029.   
[Article DOI: 10.1016/j.patcog.2007.07.023](http://dx.doi.org/10.1016/j.patcog.2007.07.023)  
[Pdf link](http://www.eecs.qmul.ac.uk/~sgg/papers/XiangGong-PR08.pdf)  

##**compute.params**

##### Description
Compute posterior probability for the gaussian mixture model. Given an eigen vector, compute the posterior probabilities that the given vector is a gaussian mixture under the given parameters. The variable names in this function follows the pattern adopted by Xiang and Gong (2008) in their manuscript.   
This function solves the following equations from the paper (Expectation part):

$h^{1}_{kn}=\frac{(1-R_{ek})\, \mathscr{N}(e_{kn}|\mu_{k1},\, \sigma_{k1})}{(1-R_{ek})\, \mathscr{N}(e_{kn}|\mu_{k1},\, \sigma_{k1})+\, \mathit{w}_{k}\, R_{ek}\, \mathscr{N}(e_{kn}|\mu_{k2},\, \sigma_{k2})+\, (1-\mathit{w}_{k})\, R_{ek}\, \mathscr{N}(e_{kn}|\mu_{k1},\, \sigma_{k1})}$

$h^{2}_{kn}=\frac{\mathit{w}_{k}\, R_{ek}\, \mathscr{N}(e_{kn}|\mu_{k2},\, \sigma_{k2})}{(1-R_{ek})\, \mathscr{N}(e_{kn}|\mu_{k1},\, \sigma_{k1})+\, \mathit{w}_{k}\, R_{ek}\, \mathscr{N}(e_{kn}|\mu_{k2},\, \sigma_{k2})+\, (1-\mathit{w}_{k})\, R_{ek}\, \mathscr{N}(e_{kn}|\mu_{k1},\, \sigma_{k1})}$

$h^{3}_{kn}=\frac{(1-\mathit{w}_{k})\, R_{ek}\, \mathscr{N}(e_{kn}|\mu_{k3},\, \sigma_{k3})}{(1-R_{ek})\, \mathscr{N}(e_{kn}|\mu_{k1},\, \sigma_{k1})+\, \mathit{w}_{k}\, R_{ek}\, \mathscr{N}(e_{kn}|\mu_{k2},\, \sigma_{k2})+\, (1-\mathit{w}_{k})\, R_{ek}\, \mathscr{N}(e_{kn}|\mu_{k1},\, \sigma_{k1})}$

The maximization part updates the parameters as:

$R_{ek}^{new}=1-\frac{1}{N}\sum_{n=1}^{N}h_{kn}^1$

$\mathit{w}_{ek}^{new}=\frac{1}{R_{ek}^{new}N}\sum_{n=1}^{N}h_{kn}^{2}$

$\mu_{k2}^{new}=\frac{\sum_{n=1}^{N}h_{kn}^{2}e_{kn}}{\sum_{n=1}^{N}h_{kn}^{2}}$

$\sigma_{k2}^{new}=\frac{\sum_{n=1}^{N}h_{kn}^{2}(e_{kn}-\mu_{kn}^{2})^{2}}{\sum_{n=1}^{N}h_{kn}^{2}}$

$\mu_{k3}^{new}=\frac{\sum_{n=1}^{N}h_{kn}^{3}e_{kn}}{\sum_{n=1}^{N}h_{kn}^{3}}$

$\sigma_{k3}^{new}=\frac{\sum_{n=1}^{N}h_{kn}^{3}(e_{kn}-\mu_{kn}^{3})^{2}}{\sum_{n=1}^{N}h_{kn}^{3}}$

##### Parameters
* vec: input eigen vector $e_{kn}$
* rel: numeric variable, $0\, {\leq}\,R_{ek}\,\leq\,1$, relevance of the vector. (default value: 0.50)
* mean2: mean of the first gaussian mixture, $\mu_{k2}$, if NULL, this parameter is  estimated based on `init` option 
* mean3: mean of the second gaussian mixture, $\mu_{k3}$, if NULL, this parameter is randomly estimated based on `init` option
* var2: variance of the first gaussian mixture, $\sigma_{k2}$, if NULL, this parameter is randomly estimated based on `init` option
* var3: variance of the second gaussian mixture, $\sigma_{k3}$, if NULL, this parameter is randomly estimated based on `init` option
* w: weight of the gaussian mixture, $\mathit{w}_{k}$, if NULL, weight is randomly estimated as `w <- runif(1,min = 0, max = 1)`
* <a name="init1"></a> init: initialization options "random" or "cluster", "random" random estimation of parameters and "cluster" use cluster mean from k-means clustering with centers = 2. For details on kmeans clustering see `kmeans` R function

##### Usage
This function is not expected to be used as such, but rather as a part of [compute.relevance](#comprel) function

##### Return
A list of many things:

* rnk: estimated $R_{ek}^{new}$
* wnk: estimated $\mathit{w}_{k}^{new}$
* m2nk: estimated $\mu_{k2}^{new}$
* v2nk: estimated $\sigma_{k2}^{new}$
* m3nk: estimated $\mu_{k3}^{new}$
* v2nk: estimated $\sigma_{k3}^{new}$


##**compute.relevance**

#####<a name="comprel"></a>Description

Given an eigenvector, compute the relevance of the vector according to the expectation maximization algorithm proposed by Xiang and Gong (2008).

##### Parameters
* vec: input eigen vector $e_{k}$ in the equations
* tol: tolerance level for convergence (default: $1e^{-6}$)
* maxit: maximum number of iterations for the expectation maximization step before convergence (default: 2500)
* maxtrials: maximum number of multiple runs (default: 25)
* init: see [init](#init1) description

##### Usage
Example usage: create random dataset 
```{r}
testdata <- matrix(runif(36,0,1),6,6)
testdata
```
Make it symmetric, and assign 0 to diagonal elements
```{r}
testdata <- testdata %*% t(testdata)
diag(testdata) <- 0
testdata
```
Compute Laplacian as $L\, =\, D-W$ and compute eigen decomposition of the laplacian
```{r}
testlap <- diag(rowSums(testdata))-testdata
testlap
testeig <- eigen(testlap,isSymmetric(testlap))
testeig
```

The eigenvectors can be used for relevance estimation like:
```{r comprel,warning=FALSE,message=FALSE}
testrel <- compute.relevance(testeig$vectors[,2],tol=1e-6,maxit=2500,maxtrials=2 )
testrel
```
The relevant eigenvectors with rnk > 0.50 can be used for further downstream processing

##### Return
A list of values:

* rnk: a vector of n $R_{ek}^{new}$ values, where `n==maxtrials`
* problist: a list of mixture probabilities, calculated as $p(e_{kn}|\theta_{ekn}^{2})\,=\,\mathit{w}_{k}\, \mathscr{N} (e_{kn}|\mu_{k2},\, \sigma_{k2})\, +\, (1-\mathit{w}_{k})\, \mathscr{N}(e_{kn}|\mu_{k3},\, \sigma_{k3})$

## **wrapper.compute.relevance**

##### Description

Wrapper for the function [compute.relevance](#comprel), instead of using a single eigenvector as input, use the eigenvector matrix. 

##### Parameters

* mat: input eigen vector matrix
* tol: tolerance level for convergence (default: $1e^{-6}$)
* maxit: maximum number of iterations for the expectation maximization step (default: 2500)
* maxtrials: maximum number of multiple runs (default: 25)
* init: see [init](#init1) description
* ncpus: number of cores to use, requires [doMC](https://cran.r-project.org/web/packages/doMC/index.html) and [foreach](https://cran.r-project.org/web/packages/foreach/index.html) packages for ncpus>1

##### Usage

```{r wcomprel,warning=FALSE,message=FALSE}
testrel <- wrapper.compute.relevance(testeig$vectors[,c(2:4)])
testrel
```

##### Return

A vector, with same length as `ncol(mat)` and each value is the relevance value that maximizes the log likelihood function $ln\, p(e_{kn}|\theta_{ekn}^{2})\,=\, ln\, \sum_{i=1}^{N}\mathit{w}_{k}\, \mathscr{N} (e_{ki}|\mu_{k2},\, \sigma_{k2})\, +\, (1-\mathit{w}_{k})\, \mathscr{N}(e_{ki}|\mu_{k3},\, \sigma_{k3})$ after `maxtrial` repeats

