---
title: "Eigen vector selection"
author: "Sudeep Sahadevan"
date: "09/17/2015"
output: html_document
---

Functions to perform informative eigen vector selection based on the algorithm proposed by:

Tao Xiang and Shaogang Gong. 2008. Spectral clustering with eigenvector selection. Pattern Recogn. 41, 3 (March 2008), 1012-1029. DOI=10.1016/j.patcog.2007.07.023 

[Article DOI](http://dx.doi.org/10.1016/j.patcog.2007.07.023)

[Pdf link](http://www.eecs.qmul.ac.uk/~sgg/papers/XiangGong-PR08.pdf)

#### <span style="color:blue">**mixture.prob**</span>


##### Description
Compute posterior probability for the mixture. Given an eigen vector, compute the posterior probabilities that the given vector is a gaussian mixture under the given parameters. The variable names in this function follows the pattern adopted by Xiang and Gong (2008) in their manuscript. This function solves the following equations from the paper:

$h^{1}_{kn}=\frac{(1-R_{ek})\, \mathscr{N}(e_{kn}|\mu_{k1},\, \sigma_{k1})}{(1-R_{ek})\, \mathscr{N}(e_{kn}|\mu_{k1},\, \sigma_{k1})+\, \mathit{w}_{k}\, R_{ek}\, \mathscr{N}(e_{kn}|\mu_{k2},\, \sigma_{k2})+\, (1-\mathit{w}_{k})\, R_{ek}\, \mathscr{N}(e_{kn}|\mu_{k1},\, \sigma_{k1})}$

$h^{2}_{kn}=\frac{\mathit{w}_{k}\, R_{ek}\, \mathscr{N}(e_{kn}|\mu_{k2},\, \sigma_{k2})}{(1-R_{ek})\, \mathscr{N}(e_{kn}|\mu_{k1},\, \sigma_{k1})+\, \mathit{w}_{k}\, R_{ek}\, \mathscr{N}(e_{kn}|\mu_{k2},\, \sigma_{k2})+\, (1-\mathit{w}_{k})\, R_{ek}\, \mathscr{N}(e_{kn}|\mu_{k1},\, \sigma_{k1})}$

$h^{3}_{kn}=\frac{(1-\mathit{w}_{k})\, R_{ek}\, \mathscr{N}(e_{kn}|\mu_{k3},\, \sigma_{k3})}{(1-R_{ek})\, \mathscr{N}(e_{kn}|\mu_{k1},\, \sigma_{k1})+\, \mathit{w}_{k}\, R_{ek}\, \mathscr{N}(e_{kn}|\mu_{k2},\, \sigma_{k2})+\, (1-\mathit{w}_{k})\, R_{ek}\, \mathscr{N}(e_{kn}|\mu_{k1},\, \sigma_{k1})}$

##### Parameters
* vec: input eigen vector $e_{kn}$
* rel: numeric variable, $0\, {\leq}\,R_{ek}\,\leq\,1$, relevance of the vector. default value: 0.50
* mean2: mean of the first gaussian mixture, $\mu_{k2}$, if NULL, estimated as `jitter(mean(vec),jitter(runif(1,max=25)))`
* mean3: mean of the second gaussian mixture, $\mu_{k3}$,if NULL, estimated as `jitter(mean(vec),jitter(runif(1,max=25)))`
* var2: variance of the first gaussian mixture, $\sigma_{k2}$,if NULL, estimated as `jitter(var(vec),jitter(runif(1,max=25)))`
* var3: variance of the second gaussian mixture, $\sigma_{k3}$,if NULL, estimated as `jitter(var(vec),jitter(runif(1,max=25)))`
* w: weight of the first gaussian, $\mathit{w}_{k}$,if NULL, weight is randomly drawn from a uniform distribution with min = 0, max = 1

##### Usage
This function is supposed to be the expectation part of the algorithm proposed by Xiang and Gong (2008) and it is not expected to run this independently.

##### Return
A list of many things:

* mean2: randomly estimated $\mu_{k2}$
* var2: randomly estimated $\sigma_{k2}$
* mean3: randomly estimated $\mu_{k3}$
* var3: randomly estimated $\sigma_{k3}$
* w: randomly estimated $\mathit{w}_{k}$
* h1: $h^{1}_{kn}$ calculated using the first formula
* h2: $h^{2}_{kn}$ calculated using the second formula
* h3: $h^{3}_{kn}$ calculated using the third formula 
* mixtpdf: $\mathit{w}_{k}\, \mathscr{N}(e_{kn}|\mu_{k2},\, \sigma_{k2})+\, (1-\mathit{w}_{k}) \mathscr{N}(e_{kn}|\mu_{k3},\, \sigma_{k3})$ the combined mixture probability

#### <span style="color:blue">**compute.vars**</span>

##### Description

Compute parameters 
This function performs the maximization step (part) of  the algorithm proposed by Xiang and Gong (2008) in their manuscript. The variable names in this function also follows the pattern adopted by Xiang and Gong (2008).

The parameters are updated according to:
<!---
#' rn <- 1-mean(h1kn)
#' wn <- mean(h2kn)*1/rn
#' m2n <- sum(h2kn*vec)/sum(h2kn)
#' v2n <- sum(h2kn*(vec-m2n)^2)/sum(h2kn)
#' m3n <- sum(h3kn*vec)/sum(h3kn)
#' v3n <- sum(h3kn*(vec-m3n)^2)/sum(h3kn) --->

$R_{ek}^{new}=1-\frac{1}{N}\sum_{n=1}^{N}h_{kn}^1$

$\mathit{w}_{ek}^{new}=\frac{1}{R_{ek}^{new}N}\sum_{n=1}^{N}h_{kn}^{2}$

$\mu_{k2}^{new}=\frac{\sum_{n=1}^{N}h_{kn}^{2}e_{kn}}{\sum_{n=1}^{N}h_{kn}^{2}}$

$\sigma_{k2}^{new}=\frac{\sum_{n=1}^{N}h_{kn}^{2}(e_{kn}-\mu_{kn}^{2})^{2}}{\sum_{n=1}^{N}h_{kn}^{2}}$

$\mu_{k3}^{new}=\frac{\sum_{n=1}^{N}h_{kn}^{3}e_{kn}}{\sum_{n=1}^{N}h_{kn}^{3}}$

$\sigma_{k3}^{new}=\frac{\sum_{n=1}^{N}h_{kn}^{3}(e_{kn}-\mu_{kn}^{3})^{2}}{\sum_{n=1}^{N}h_{kn}^{3}}$