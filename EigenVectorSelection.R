# TODO: Add comment
# 
# Author: Sudeep Sahadevan
# 
# Functions to perform informative eigen vector selection 
# Based on the algorithm proposed by:
# Tao Xiang and Shaogang Gong. 2008. Spectral clustering with eigenvector selection. 
# Pattern Recogn. 41, 3 (March 2008), 1012-1029. DOI=10.1016/j.patcog.2007.07.023 
# http://dx.doi.org/10.1016/j.patcog.2007.07.023 
# Pdf link: http://www.eecs.qmul.ac.uk/~sgg/papers/XiangGong-PR08.pdf
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
###############################################################################

#' Compute posterior probability for the mixture
#' Given an eigen vector, compute the posterior probabilities that the given vector is a gaussian mixture under the given parameters.
#' The variable names in this function follows the pattern adopted by Xiang and Gong (2008) in their manuscript.
#' 
#' Expectation part
#' h1kn <- ((1-rel)*pnorm(vec,mean1,sqrt(var1)))/(((1-rel)*pnorm(vec,mean1,sqrt(var1)))+(w*rel*pnorm(vec,mean2,sqrt(var2)))+((1-w)*rel*pnorm(vec,mean3,sqrt(var3)))) # eqn1
#' h2kn <- (w*rel*pnorm(vec,mean2,sqrt(var2)))/(((1-rel)*pnorm(vec,mean1,sqrt(var1)))+(w*rel*pnorm(vec,mean2,sqrt(var2)))+((1-w)*rel*pnorm(vec,mean3,sqrt(var3)))) # eqn2
#' h3kn <- ((1-w)*rel*pnorm(vec,mean3,sqrt(var3))))/(((1-rel)*pnorm(vec,mean1,sqrt(var1)))+(w*rel*pnorm(vec,mean2,sqrt(var2)))+((1-w)*rel*pnorm(vec,mean3,sqrt(var3)))) # eqn3
#'
#' Maximization part
#' rn <- 1-mean(h1kn)
#' wn <- mean(h2kn)*1/rn
#' m2n <- sum(h2kn*vec)/sum(h2kn)
#' v2n <- sum(h2kn*(vec-m2n)^2)/sum(h2kn)
#' m3n <- sum(h3kn*vec)/sum(h3kn)
#' v3n <- sum(h3kn*(vec-m3n)^2)/sum(h3kn)
#' 
#' @param vec: input eigen vector
#' @param rel: numeric variable, 0 <= rek <=1, relevance of the vector. default value: 0.50
#' @param mean2: mean of the first gaussian mixture, if NULL, estimated as jitter(mean(v),jitter(runif(1,max=10)))
#' @param mean3: mean of the second gaussian mixture, if NULL, estimated as jitter(mean(v),jitter(runif(1,max=10)))
#' @param var2: variance of the first gaussian mixture, if NULL, estimated as jitter(var(v),jitter(runif(1,max=10)))
#' @param var3: variance of the second gaussian mixture, if NULL, estimated as jitter(var(v),jitter(runif(1,max=10)))
#' @param w: weight of the first gaussian, if NULL, weight is randomly drawn from a uniform distribution with min = 0, max = 1
#' @param init: initialization options "random" or "cluster", "random" random estimation of parameters and "cluster" use cluster mean from k means clustering with c = 2
#' @return a list of many things
#' 
compute.params <- function(vec,rel=0.5,mean2=NULL,mean3=NULL,var2=NULL,var3=NULL,w=NULL,init="cluster"){
	check.na <- function(invec){
		invec[is.na(invec)] <- 0
		return(invec)
	}
	# 	Expectation step begin
	if(is.null(vec)){stop("ERROR! input eigen vector cannot be empty\n")}
	else if(length(vec)<2){"ERROR! input eigen vector must have more than one value\n"}
	init <- match.arg(arg=init,choices = c("random","cluster"),several.ok = FALSE)
	vec <- na.omit(vec)
	mean1 <- mean(vec)
	var1 <- mean((vec-mean1)^2) # WARNING! dnorm|rnorm uses standard devation NOT variance as input!, since var() gives the unbiased variance, calculate the biased (population version)
	if((is.null(mean2)||is.null(mean3)) && (init=="random")){ # if mean is NULL, use random estimate for variance as well
		nsamp <- length(vec)*runif(1,0.1,0.9)
		vec2 <- sample(vec,max(3,nsamp),FALSE) # avoid sampling <3 values 
		mean2 <- mean(vec2)
		var2 <- mean((vec2-mean2)^2)
		vec3 <- setdiff(vec,vec2)
		mean3 <- mean(vec3)
		var3 <- mean((vec3-mean3)^2)
	}
	else if((is.null(mean2)||is.null(mean3)) && (init=="cluster")){ # use k-means clustering based parameter initialization
		k2 <- kmeans(vec,centers=2,iter.max=2500,nstart=10)
		if(k2$ifault==4){
			k2 <- kmeans(x=vec,centers=2,iter.max=2500,nstart=10,algorithm="MacQueen")
		}
#		k2$centers <- jitter(k2$centers) # add small amount of noise to the cluster means
		k21 <- vec[which(k2$cluster==1)]
		k22 <- vec[which(k2$cluster==2)]
		mean2 <- k2$centers[1]
		var2 <- mean((k21-mean2)^2)
		mean3 <- k2$centers[2]
		var3 <- mean((k22-mean3)^2)
	}
	if(is.null(w)){
		w <- runif(1,min = 0, max = 1)
	}
	pdf1 <- pdf2 <- pdf3 <- rep(0,length(vec))
	if(round(var1,15)>0){
		pdf1 <- (1-rel)*dnorm(vec,mean1,sqrt(var1))
	}
	if(round(var2,15)>0){
		pdf2 <- w*rel*dnorm(vec,mean2,sqrt(var2))
	} 
	if(round(var3,15)>0){
		pdf3 <- (1-w)*rel*dnorm(vec,mean3,sqrt(var3))
	}
	cdn <- pdf1 + pdf2 + pdf3
	h1kn <- check.na(pdf1/cdn)
	h2kn <- check.na(pdf2/cdn)
	h3kn <- check.na(pdf3/cdn)
	# 	Maximization step
	rn <- 1 
	wn <- 0 # updated relevance and weight
	m2n <- v2n <- 0 # updated mean and variance for the first mixture
	m3n <- v3n <- 0 # updated mean and variance for the second mixture
	# 	update params
	if(round(sum(h1kn),15)>0){ # sanity check to avoid crashes and bizzare results and round of sum to avoid rounding off errors
		rn <- 1-mean(h1kn)
	}
	if(round(sum(h2kn),15)>0){
		wn <- mean(h2kn)
		if(rn>0){ # avoid wn/0
			wn <- wn/rn
		} else{ wn <- 0}
		m2n <- sum(h2kn*vec)/sum(h2kn)
		v2n <- sum(h2kn*((vec-m2n)^2))/sum(h2kn)
	}
	if(round(sum(h3kn),15)>0){
		m3n <- sum(h3kn*vec)/sum(h3kn)
		v3n <- sum(h3kn*((vec-m3n)^2))/sum(h3kn)
	}
#	if(is.na(var1)||is.nan(var1)){ 
#		return(pdf1=pdf1,cdn=cdn)
#	}
#	if(is.na(v2n)||is.nan(v2n)){ 
#		return(list(pdf2=pdf2,cdn=cdn)) }
#	if(is.na(v3n)||is.nan(v3n)){ 
#		return(list(pdf3=pdf3,cdn=cdn)) 
#	}
	return(list(rnk=rn,wnk=wn,m2nk=m2n,v2nk=v2n,m3nk=m3n,v3nk=v3n))
}

#' Compute relevance of an eigen vector using expectation maximization algorithm
#' Algorithm as proposed by Xiang and Gong (2008)
#' @param vec: input eigen vector
#' @param tol: tolerance level for convergence
#' @param maxit: maximum number of iterations for the expectation maximization step
#' @param init: initialization options "random" or "cluster", "random" random estimation of parameters and "cluster" use cluster mean from k means clustering with c = 2
#' @param maxtrials: maximum number of multiple runs
compute.relevance <- function(vec,tol=1e-6,maxit=2500,maxtrials=25,init="cluster"){
	vec <- na.omit(vec)
#	function for convergence check
#	p1: a vector of densities
#	p2: another vector of densities
#	ln: Boolean, whether to convert values to natural logarithm
	conv <- function(p1,p2,ln=TRUE){
		if(ln){
			p1 <- sum(log(p1))
			p2 <- sum(log(p2))
		}
		if((sum(abs(p1-p2)))<tol){ return(TRUE)}
		else{FALSE}
	}
#	some variables to store stuff
	rnvec <- rep(0,maxtrials)
	problist <- list()
#	data init
	orig.mean <- mean(vec) # need it here
	orig.var <- mean((vec-orig.mean)^2)
	orig.prob <- dnorm(vec,orig.mean,sqrt(orig.var))
	for(r in 1:maxtrials){
#		init vars 
		vec <- na.omit(vec)
		vecvars <- vector()
		rn <- 0.5
		wn <- NULL
		m2 <- v2 <- NULL 
		m3 <- v3 <- NULL
		prev.prob <- cur.prob <- rep(1,length(vec))
		it <- 1
		convergence <- FALSE
# 		do some work		
		while((!convergence)){
			em.vars <- compute.params(vec = vec,rel = rn,mean2 = m2,mean3 = m3,var2 = v2,var3 = v3,w = wn,init=init)
			rn <- em.vars$rnk
			wn <- em.vars$wnk
			m2 <- em.vars$m2nk
			v2 <- em.vars$v2nk
			m3 <- em.vars$m3nk
			v3 <- em.vars$v3nk			
			if((abs(rn)==0)){
				convergence <- TRUE
				cur.prob <- orig.prob
			}
			if((v2==0) && (v3==0)){
				convergence <- TRUE
				warning("The second and third variance values are zero\n")
				cur.prob <- orig.prob
			}else{
				cur.prob <-  ((1-rn)*orig.prob)+.helper.mixture.pdf(vec,rn,wn,m2,v2,m3,v3)
				cur.prob[which(cur.prob==0)] <- 1e-200 # if the probability is 0, add an extremely small value so that log likelihood is not -Inf
				convergence <- conv(cur.prob,prev.prob,TRUE)
			}
			prev.prob <- cur.prob
			if(is.nan(cur.prob[1])){ return (pre.vars)}
			pre.vars <- em.vars
			if(it==maxit){
				cnv <- sum(abs(sum(log(cur.prob)) - sum(log(prev.prob))))
				warning("The algorithm did not converge, but maximum number of iterations reached\n Difference: ",cnv,"\n")
				convergence <- TRUE
			}
			it <- it + 1
		}
		rnvec[r] <- rn
		tmpprob <- .helper.mixture.pdf(vec,NULL,wn,m2,v2,m3,v3)
		tmpprob[which(tmpprob==0)]<- 1e-200
		problist[[r]] <- tmpprob
	}
	return(list(rnk=rnvec,problist=problist))
}

#' helper function calculate mixture pdf according to Xiang and Gong (2008)
#' 
.helper.mixture.pdf <- function(vec,rn=NULL,wn,m2n,v2n,m3n,v3n){
	if(is.null(rn)){ rn <- 1 }
	p2 <- (rn*wn*dnorm(vec,m2n,sqrt(v2n)))
	p3 <- (rn*(1-wn)*dnorm(vec,m3n,sqrt(v3n)))
	return(p2+p3)
}

#' Wrapper for compute.relevance function, given a matrix of eigen vectors, compute the relevance
#' and return the rnk value that maximizes the log likelyhood function
#' @param mat: input eigen vector matrix
#' @param tol: tolerance level for convergence
#' @param maxit: maximum number of iterations for the expectation maximization step
#' @param maxtrials: maximum number of multiple runs
#' @param init: initialization options "random" or "cluster", "random" random estimation of parameters and "cluster" use cluster mean from k means clustering with c = 2
#' @param ncpus: number of cores to use, requires doMC and foreach packages 
#' @return vector, with same length as ncol(mat) and each value is the relevance value that maximizes the log likelyhood after maxtrial repeats
wrapper.compute.relevance <- function(mat,tol=1e-6,maxit=2500,maxtrials=25,init="cluster",ncpus=1){
	mat <- na.omit(mat)
	ncm <- ncol(mat)
	if(ncpus>1){
		library(doMC)
		library(foreach)
		registerDoMC(cores=ncpus)
		system(sprintf(paste("taskset -p 0x",paste(rep("f",ncpus),collapse="")," %d",sep=""), Sys.getpid())) # makes sure that doMC can use all the cores in Openblas environment
#		do work
		rellist <- foreach(r=1:ncm)%dopar%{
			if((r%%50==0)||r==ncm){ message(r,"/",ncm,"\n") }
			tmp.list <- compute.relevance(mat[,r],tol,maxit,maxtrials,init)
			maxprob <- which.max(unlist(lapply(tmp.list$problist,function(x) sum(log(x)))))
			return(tmp.list$rnk[maxprob])
			
		}
		return(unlist(rellist))
	}else{
		relvec <- rep(0,ncm)
		for(i in 1:ncm){
			rel.list <- compute.relevance(mat[,i],tol,maxit,maxtrials,init)
#			find the relevance variable that maximizes the log likelyhood function
			maxprob <- which.max(unlist(lapply(rel.list$problist,function(x) sum(log(x)))))  
			relvec[i] <- rel.list$rnk[maxprob]
			message(i,"/",ncm)
		}
		return(relvec)
	}
}

#' helper function, given a vector of values, mean and variance (optional)
#' return probability density function for the normal distribution
.helper.multivariate.pdf <- function(vec,mux=NULL,varx=NULL){
	vec <- na.omit(vec)
	if(is.null(mux)){
		mux <- mean(vec)
	}
	if(is.null(varx)){
		varx <- mean((vec-mux)^2)
	}
	if(round(varx,15)>0){
#		return((1/(sqrt(((pi*2)^length(vec))*varx)))*(exp(-(vec-mux)^2/(2*varx))))
		return((1/(sqrt(varx)*((pi*2)^(length(vec)/2))))*(exp(-(vec-mux)^2/(2*varx))))
	}else{
		return(rep(0,length(vec)))
	}
}

#' order eigen vectors accroding to relevance
#' given an eigen vector matrix, a corresponding relevance vector and threshold cut off,  prune and order the eigen vectors according to relevance
#' @param evec: eigen vector matrix
#' @param rel: relevance vector
#' @param thresh: relevance threshold
eigenvec.order <- function(evec,rel,thresh=0.50){
	if(ncol(evec)<length(rel)){ stop("ERROR! ncol(evec) must be >= length(rel)\n") }
	prel <- rel[order(rel,decreasing = TRUE)]
	pvec <- evec[,order(rel,decreasing = TRUE)]
	pvec <- pvec[,which(prel>thresh)]
	return(pvec)
}