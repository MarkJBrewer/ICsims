#' Run simulations allowing for between-sample heterogeneity
#' 
#' \code{simulation.study} implements a simulation framework sampling repeatedly
#' from linear regression models and GLMs, allowing for between-sample heterogeneity.
#' The purpose is to allow the study of AIC and related statistics in the
#' context of model selection, with prediction quality as target.
#'
#' @param type Character string determining what type of model to fit. At present,
#' available model types are "lm" and "glm", with the former the default.
#' @param nsims Number of simulated data sets to analyse for each sample size
#' @param nsamples Vector of integers containing the sample sizes
#' @param nX Number of "real" covariates
#' @param nZ Number of "spurious" covariates
#' @param alpha Intercept for simulation model
#' @param beta.x Either: vector of slopes for the X covariates; or a single numeric values
#' for a constant slope for all X's
#' @param meanX Either: vector of means for the X covariates; or a single numeric value
#' for a constant mean across all X's
#' @param meanZ As for meanX but for the Z covariates
#' @param XZCov Covariance matrix of the X's and Z's. Must be of dimension (nX+nZ) by
#' (nX+nZ). Ignored if \code{simulate.from.data==TRUE} or if \code{is.null(rho)==FALSE}
#' @param varmeanX Either: vector of variances for the means of the X covariates; or a
#' single numeric value for a constant mean across all X's. Non-zero
#' values will produce a different set of covariate means for each
#' individual simulated data set
#' @param varmeanZ As for varmeanX but for the Z covariates
#' @param simulate.from.data Logical. If \code{TRUE}, function takes actual covariate data to
#' use as the basis of simulations; if \code{FALSE} (the default) the function uses the
#' distributions defined by the model parameters given as input to the function
#' @param X Matrix of "real" covariates; only used if \code{simulate.from.data==TRUE}
#' @param Y Vector of "real" response variables; only used if \code{simulate.from.data==TRUE},
#' and if given the values of alpha and beta above will be ignored,
#' but instead derived from a regression model of Y against X
#' @param var.res Residual variance of the simulation model
#' @param var.RE.Intercept Random effect variance for the intercept
#' @param var.RE.X Either: vector of random effect variances for the X covariate slopes;
#' or a single numeric value for no random slopes in the X's
#' @param rho A numeric constant specifying the mean correlation between the X's and the
#' Z's
#' @param epsilon A numeric constant specifying the level of variability around the mean
#' correlation rho; note that a necessary condition is \code{max(abs(rho+epsilon))<=1}, so
#' combinations of rho and epsilon which break this constraint will cause and error (as
#' the simcor function would produce correlations outside the range [-1,1]). If not
#' supplied but \code{is.null(rho)==TRUE}, then epsilon is set to zero
#' @param corsim.var If generating the covariance matrices using rho and epsilon, we
#' to specify the variances (which otherwise are in the leading
#' diagonal of XZCov)
#' @param noise.epsilon A numeric constant used to specify whether XZCov is to vary from
#' sample to sample. Higher values indicate more variability; note that this cannot be
#' greater than 1 minus the largest absolute value of (off-diagonal) correlations in the
#' corresponding correlation matrix
#' @param step.k Numeric value of the AIC criterion in the stepwise analysis; defaults
#' to about 3.84, corresponding to a p-value of 0.05 for both adding and removing variables
#' @param keep.dredge Logical constant on whether to keep the dredge outputs
#' (\code{TRUE==yes}); required if \code{simulate.from.data==TRUE}
#' @param Xin.or.out Vector of length nX (or \code{nrow(X)}) of logicals, specifying whether
#' an X is made available as data (\code{TRUE} for yes; \code{FALSE} for no)
#' @param glm.family If a GLM is to be fitted, the error distribution must be supplied
#' (to the standard family argument to glm).
#' @param glm.offset An (optional) offset can be supplied if fitting a GLM. (Not
#' currently implemented.)
#' @param filename Character string providing the root for the output files. Intermediate files are
#' saved as "filenameX.RData" where X is an incremental count from 1 to length(nsamples). The final
#' output is in "filename.RData".
#' @param binomial.n If fitting a binomial GLM, the number of trials per sample. Must be
#' either a scalar (in which case the same number of trials are used for each sample) or a
#' vector of length nsamples. (Default is 1)
#' @return If \code{keep.dredge==FALSE} (the default), the output is a list of length equal to the length
#' of nsamples, each containing two matrices, \code{reg.bias} (prediction bias for each
#' sample) and \code{reg.rmse} (root mean square error of prediction for each sample). Each
#' of these two matrices has length \code{nsims} and four columns, corresponding to model 
#' selection by AICc, AIC, BIC and stepwise regression. If \code{keep.dredge==TRUE}, then the output
#' is a list of lists, with a top level list with length equal to the length of nsamples as before, and
#' with the next level having length equal to \code{nsims}; this inner list contains the full model set
#' output from \code{dredge}, converted to a matrix for storage efficiency.
#' @export
simulation.study <- function(type="lm",nsims=1000,
    nsamples=c(20,50,100,200,500,1000,2000,5000,10000),
    alpha=0,beta.x=1,nX=10,nZ=5,meanX=0,meanZ=0,
    XZCov=diag(nX+nZ),varmeanX=0,varmeanZ=0,simulate.from.data=FALSE,
    X=NULL,Y=NULL,var.res=1,var.RE.Intercept=0,var.RE.X=0,
    rho=NULL,epsilon=NULL,corsim.var=NULL,noise.epsilon=NULL,
    step.k=qchisq(0.05,1,lower.tail=FALSE),keep.dredge=FALSE,
    Xin.or.out=rep(TRUE,nX),glm.family=NULL,glm.offset=NULL,
    binomial.n=1,filename="results"){

    if(!(type%in%c("lm","glm"))){
        stop("Error: type not recognised; should be either \"lm\" or \"glm\".")
    }
    if(type=="glm"){
        if(is.null(glm.family)){
            stop("Error: please supply \"glm.family\" to specify error distribution for GLM.")
        }
        if(!(substr(glm.family,1,3)%in%c("poi","bin"))){
            stop("Error: GLM errors currently restricted to \"binomial\" and \"poisson\"")
        }
    }

    if( type=="glm" && substr(glm.family,1,3)=="bin" ){
        rmse.calc <- function(model,newdata,y){
	       sqrt(mean((predict(model,newdata=newdata,type="response")-y)^2))
        }
        if(length(binomial.n)>1.5){
            if(length(nsamples)>1.5){
                stop("Error: can only run simulations for a single value of nsamples if a vector binomial.n is supplied.")
            }
            if(length(binomial.n)!=nsamples){
                stop("Error: if a vector binomial.n is supplied, it must have length nsamples")
            }
        }
    }else{
        rmse.calc <- function(model,newdata,y){
	       sqrt(mean((predict(model,newdata=newdata,type="response")-y)^2))
        }
    }

    if(!is.null(glm.offset)){
        if(length(glm.offset)!=n){
            glm.offset <- rep(glm.offset[1],n)
        }
    }
    if(simulate.from.data){
        if(is.null(X)){
            stop("Error: X required if simulating from data.\n")
        }
        X <- as.matrix(X)
        dimnames(X) <- NULL
        nX <- ncol(X)
        meanX <- colMeans(X)
        XZCov <- cov(X)
        corsim.var <- diag(XZCov)
        nZ <- 0
        if(!is.null(Y)){
            if(type=="lm"){
                reg.temp <- lm(Y~X)
            }else{
                if( is.null(glm.offset) ){
                    reg.temp <- glm(Y~X,family=glm.family)
                }else{
                    reg.temp <- glm(Y~X,family=glm.family,offset=glm.offset)
                }
            }
            alpha <- coef(reg.temp)[1]
            beta.x <- coef(reg.temp)[-1]
        }
        keep.dredge <- TRUE
    }

    if(keep.dredge){
        dredge.out <- list()
    }

    XZCor <- cov2cor(XZCov)
    #print(XZCor)

    if(!is.null(rho)){
        if(is.null(epsilon)){
            epsilon <- 0
        }else{
            if(abs(epsilon+rho)>1){
                stop("Error: abs(epsilon+rho)>1\n")
            }
        }
        if(is.null(corsim.var)){
            corsim.var <- 1
        }
        if(length(corsim.var)==1){
            corsim.var <- rep(corsim.var,nX+nZ)
        }
        XZCor <- simcor(k=1,size=nX+nZ,rho=rho,epsilon=epsilon)
        XZCov <- sqrt(corsim.var)*XZCor*rep(sqrt(corsim.var),each=nX+nZ)
    }
    if(!is.null(noise.epsilon)){
        if(is.null(corsim.var)){
            corsim.var <- 1
        }
        if(length(corsim.var)==1){
            corsim.var <- rep(corsim.var,nX+nZ)
        }
        maxeps <- 0.9999-max(abs(XZCor*(upper.tri(XZCor)+lower.tri(XZCor))))
        noise.epsilon <- min(noise.epsilon,maxeps)
    }

    n.sample.sizes <- length(nsamples)
    criteria <- c("AICc","AIC","BIC","stepwise")
    ncriteria <- length(criteria)

    sd.res <- sqrt(var.res)
    sdmeanX <- sqrt(varmeanX)
    sdmeanZ <- sqrt(varmeanZ)

    beta.z <- rep(0,nZ)
    Sigma.RE.X <- var.RE.X*diag(nX)
    sd.RE.Intercept <- sqrt(var.RE.Intercept) # should perhaps be correlated with RE.X?

    if(length(beta.x)==1){
        beta.x <- rep(beta.x,nX)
    }

    results <- list()
    for(i in 1:n.sample.sizes){
        print(date())
        if(keep.dredge){
            dredge.out[[i]] <- list()
        }
        n <- nsamples[i]
        all.x.names <- paste("x.",1:nX,sep="")
        z.names <- paste("z.",1:nZ,sep="")
        x.names <- all.x.names[Xin.or.out]
        reg.terms <- vector("list",nsims)
        reg.bias <- array(NA,dim=c(nsims,ncriteria))
        colnames(reg.bias) <- criteria
        reg.rmse <- array(NA,dim=c(nsims,ncriteria))
        colnames(reg.rmse) <- criteria
        if(nZ > 0.5){
            reg.eqn <- as.formula(paste("y",paste(paste("",x.names,sep="",collapse="+"),
                paste("",z.names,sep="",collapse="+"),sep="+",collapse="+"),sep="~"))
        }else{
            reg.eqn <- as.formula(paste("y",paste(paste("",x.names,sep="",collapse="+"),
                sep="+",collapse="+"),sep="~"))
        }
        options(na.action = "na.fail") # Necessary
        reg.preds <- array(NA,dim=c(n,ncriteria))
        for(j in 1:nsims){
            cat(paste(i,":",j,":"))
            # Generate "true" response
            x.sim <- rnorm(nX,meanX,sdmeanX)
            z.sim <- rnorm(nZ,meanZ,sdmeanZ)
            if(!is.null(noise.epsilon)){
                XZCor.temp <- noisecor(XZCor,epsilon=noise.epsilon)
                XZCov <- sqrt(corsim.var)*XZCor.temp*rep(sqrt(corsim.var),each=nX+nZ)
            }
            x.and.z <- mvrnorm(n=n,mu=c(x.sim,z.sim),Sigma=XZCov)
            x <- x.and.z[,1:nX]
            if(nZ > 0.5){
                z <- x.and.z[,(1+nX):(nZ+nX)]
            }
            alpha.sim <- rnorm(1,alpha,sd.RE.Intercept)
            beta.sim <- mvrnorm(n=1,mu=beta.x,Sigma=Sigma.RE.X)
            if(type=="lm"){
                y <- rnorm(n,0,sd.res)
                y <- as.numeric(y+alpha.sim+x%*%beta.sim)
            }else{
                if(substr(glm.family,1,3)=="bin"){
                    y1.prob <- as.numeric(alpha.sim+x%*%beta.sim)
                    y1.prob <- exp(y1.prob)/(1+exp(y1.prob))
                    y1 <- rbinom(n,binomial.n,y1.prob)
                    y2 <- binomial.n-y1
                    y <- cbind(y1,y2)
                }else{ # hence Poisson
                    y.mean <- exp(as.numeric(alpha.sim+x%*%beta.sim))
                    if(!is.null(glm.offset)){
                        y.mean <- exp(glm.offset)*y.mean
                    }
                    y <- rpois(n,y.mean)
                }
            }
            if(nZ > 0.5){
                reg.data <- data.frame(y=y,x=x,z=z)
            }else{
                reg.data <- data.frame(y=y,x=x)
            }
            if(type=="lm"){
                reg.model <- lm(reg.eqn,data=reg.data)
            }else{
                if(is.null(glm.offset)){
                    reg.model <- glm(reg.eqn,family=glm.family,data=reg.data)
                }else{
                    reg.model <- glm(reg.eqn,family=glm.family,offset=glm.offset,data=reg.data)
                }
            }
            reg.dredge <- dredge(reg.model,extra=c("AIC","BIC"))
            capture.output(reg.stepwise <- step(reg.model,k=step.k))
            # Get second, test set
            x.sim <- rnorm(nX,meanX,sdmeanX)
            z.sim <- rnorm(nZ,meanZ,sdmeanZ)
            if(!is.null(noise.epsilon)){
                XZCor.temp <- noisecor(XZCor,epsilon=noise.epsilon)
                XZCov <- sqrt(corsim.var)*XZCor.temp*rep(sqrt(corsim.var),each=nX+nZ)
            }
            x.and.z <- mvrnorm(n=n,mu=c(x.sim,z.sim),Sigma=XZCov)
            x.new <- x.and.z[,1:nX]
            if(nZ > 0.5){
                z.new <- x.and.z[,(1+nX):(nZ+nX)]
            }
            alpha.sim <- rnorm(1,alpha,sd.RE.Intercept)
            beta.sim <- mvrnorm(n=1,mu=beta.x,Sigma=Sigma.RE.X)
            if(type=="lm"){
                y.new <- rnorm(n,0,sd.res)
                y.new <- as.numeric(y.new+alpha.sim+x.new%*%beta.sim)
            }else{
                if(substr(glm.family,1,3)=="bin"){
                    y1.prob <- as.numeric(alpha.sim+x.new%*%beta.sim)
                    y1.prob <- exp(y1.prob)/(1+exp(y1.prob))
                    y1 <- rbinom(n,binomial.n,y1.prob)
                    y2 <- binomial.n-y1
                    y.new <- cbind(y1,y2)
                }else{ # hence Poisson
                    y.mean <- exp(as.numeric(alpha.sim+x.new%*%beta.sim))
                    if(!is.null(glm.offset)){
                        y.mean <- exp(glm.offset)*y.mean
                    }
                    y.new <- rpois(n,y.mean)
                }
            }
            if(nZ > 0.5){
                newdata <- data.frame(x=x.new,z=z.new)
            }else{
                newdata <- data.frame(x=x.new)
            }
            if( type=="glm" && substr(glm.family,1,3)=="poi" ){
                newdata$glm.offset <- glm.offset
            }
            if( type=="glm" && substr(glm.family,1,3)=="bin" ){
                y.calc <- y.new[,1]/rowSums(y.new)
            }else{
                y.calc <- y.new
            }
            reg.best <- get.models(reg.dredge,subset=which.min(AICc))[[1]]
            reg.preds[,1] <- predict(reg.best,newdata=newdata,type="response")
            reg.best <- get.models(reg.dredge,subset=which.min(AIC))[[1]]
            reg.preds[,2] <- predict(reg.best,newdata=newdata,type="response")
            reg.best <- get.models(reg.dredge,subset=which.min(BIC))[[1]]
            reg.preds[,3] <- predict(reg.best,newdata=newdata,type="response")
            reg.preds[,4] <- predict(reg.stepwise,newdata=newdata,type="response")
            reg.bias[j,] <- colMeans(reg.preds-y.calc)
            reg.rmse[j,] <- sqrt(colMeans((reg.preds-y.calc)^2))
            if(keep.dredge){
                capture.output(dredge.out[[i]][[j]] <- as.matrix(model.sel(reg.dredge,rank=rmse.calc,rank.args=list(newdata=newdata,y=y.calc),extra=alist(AICc,AIC,BIC))))
            }
        }
        options(na.action = "na.omit")
        print(date())

        results[[i]] <- list(reg.bias,reg.rmse)
        names(results[[i]]) <- c("reg.bias","reg.rmse")
        save(results[[i]],file=paste(filename,i,".RData",sep=""))
    }
    save(results,file=paste(filename,".RData",sep=""))
    if(keep.dredge){
        return(dredge.out)
    }else{
        return(results)
    }
}

