# File:     pslm.R
# Author:   Olli Saarela (olli.saarela@utoronto.ca)
# Date:     2016-02-06
# Summary:  Produce the simulation example of Section 7 in
#  Saarela, O., Belzile, L. R. and D. A. Stephens. A Bayesian view of doubly robust causal inference,
#  Biometrika (2016), 103 (3): 667-681, doi:10.1093/biomet/asw025

# Set output directory
outpath <- ""

##Load necessary libraries and cast warnings to errors
rm(list=ls())
library(mvtnorm)
library(splines)
library(MCMCpack)
library(parallel)
# options(warn=2)

##Canonical link function, inverse link function and error function - used in DGM
expit <- function(x) {1/(1+exp(-x))}
logit <- function(p) {log(p)-log(1-p)}
erf <- function(x) {2 * pnorm(x * sqrt(2)) - 1}

# from <- 1; to <- 1; numobs <- 1000; ncov <- 4; coeffa <- 0.4; coeffb <- -1.0
# a0 <- 0.0; b0 <- 0.0; effect <- 1.0; nrep <- 1000; i=1; j=1; no=1; msor=TRUE; msps=FALSE

pssim <- function(no, outpath=outpath, numobs=1000, ncov=4, coeffa=0.4, coeffb=-1.0, a0=0.0, b0=0.0, effect=1.0, nrep=1000,
                  msor=TRUE, msps=FALSE) {
    #Number of replications for the simulation study
    numsim <- 150
    ncores <- 7
    tot <- ncores * numsim
    from <- seq(1,tot,by=numsim)[no]
    to <- seq(numsim,tot,by=numsim)[no]

    ##Weighted log-likelihood function used in call to optim - includes an extra parameter
    ##log-likelihood proportional to that of gaussian model wrt param[1:(idx-1)], aka beta
    ##additional param [idx] is precision, parametrized as exp(-param[idx])
    ##this restriction ensures positive variance
	loglik <- function(param, resp, mmat, weight) {
        idx <- length(param)
        lp <- mmat %*% param[1:(idx-1)]
        return(-param[idx] * sum(weight)/2.0 - 1.0/(2.0 * exp(param[idx])) * sum(weight * (resp - lp)^2))
    }
    #Gradient function of loglik
    gradient <- function(param, resp, mmat, weight) {
        idx <- length(param)
        lp <- mmat %*% param[1:(idx-1)]
        return(c(1.0/exp(param[idx]) * colSums(matrix(weight * (resp - lp), numobs, ncol(mmat)) * mmat),
                 -sum(weight)/2.0 + 1.0/(2.0 * exp(param[idx])) * sum(weight * (resp - lp)^2)))
    }
    ##Covariance matrix  - Fisher information matrix (as quadratic form of the score function, i.e. observed)
    scorecovariance <- function(param, resp, mmat, weight) {
        idx <- length(param)
        lp <- mmat %*% param[1:(idx-1)]
        scorem <- cbind(1.0/exp(param[idx]) * matrix(weight * (resp - lp), numobs, ncol(mmat)) * mmat,
                        -weight/2.0 + 1.0/(2.0 * exp(param[idx])) * weight * (resp - lp)^2)
        return(t(scorem) %*% scorem)
    }

    scorecrossprod <- function(param1, param2, resp1, resp2, mmat1, mmat2, weight1, weight2) {
        idx <- length(param1)
        lp1 <- mmat1 %*% param1[1:(idx-1)]
        scorem1 <- cbind(1.0/exp(param1[idx]) * matrix(weight1 * (resp1 - lp1), numobs, ncol(mmat1)) * mmat1,
                         -weight1/2.0 + 1.0/(2.0 * exp(param1[idx])) * weight1 * (resp1 - lp1)^2)
        lp2 <- mmat2 %*% param2
        scorem2 <- mmat2 * matrix(weight2 * (resp2 - 1.0/(1.0 + exp(-lp2))), nrow(mmat2), ncol(mmat2), byrow=FALSE)
        return(t(scorem1) %*% scorem2)
    }
    ##Numerical derivative of the gradient wrt gamma, based on finite differences
    gd <- function(param1, param2, adjusted=FALSE) {
        delta <- 0.0001
        dv <- matrix(NA, length(param1), length(param2))
        for (j in 1:length(param2)) {
            e <- rep(0, length(param2))
            e[j] <- 1
            if (adjusted) {
                dv[,j] <- (gradient(param1, y, cbind(1.0, z, x, predict(psbasis, expit(cbind(1, lfx) %*% (param2 + delta * e)))), rep(1.0, numobs)) -
                           gradient(param1, y, cbind(1.0, z, x, predict(psbasis, expit(cbind(1, lfx) %*% (param2 - delta * e)))), rep(1.0, numobs)))/(2 * delta)
            }   else {
                dv[,j] <- (gradient(param1, y, cbind(1.0, z, predict(psbasis, expit(cbind(1, lfx) %*% (param2 + delta * e)))), rep(1.0, numobs)) -
                           gradient(param1, y, cbind(1.0, z, predict(psbasis, expit(cbind(1, lfx) %*% (param2 - delta * e)))), rep(1.0, numobs)))/(2 * delta)
            }
        }
        return(dv)
    }
    ##Joint log-likelihood function: binary exposure, normal response
    ##param consists of parameters for outcome, variance of the latter model and ps model
    ##ps model parametrized viz lfx (transformed covariates)
	jointloglik <- function(param, adjusted=TRUE) {
        idx <- 3 + adjusted * ncol(x) + ncol(psbasis)
        mmat2 <- cbind(1.0, lfx)
        lp2 <- mmat2 %*% param[(idx+1):length(param)]
        if (adjusted) {
            mmat1 <- cbind(1.0, z, x, predict(psbasis, expit(lp2)))
        }         else {
            mmat1 <- cbind(1.0, z, predict(psbasis, expit(lp2)))
        }
        lp1 <- mmat1 %*% param[1:(idx-1)]

        return(-param[idx] * numobs/2.0 - 1.0/(2.0 * exp(param[idx])) * sum((y - lp1)^2)
               + sum(z * lp2 - log(1.0 + exp(lp2))))
        #previous line includes contribution of binomial exposure model with canonical link function
    }

    ##Joint log-likelihood function, version 1
    ##with fitted values derived from fixed gamma (parameter of the fitted ps via glm)
    jointloglik1 <- function(param, adjusted=TRUE) {
        idx <- 3 + adjusted * ncol(x) + ncol(psbasis)
        mmat2 <- cbind(1.0, lfx)
        lp2 <- mmat2 %*% gamma
        if (adjusted) {
            mmat1 <- cbind(1.0, z, x, predict(psbasis, expit(lp2)))
        }         else {
            mmat1 <- cbind(1.0, z, predict(psbasis, expit(lp2)))
        }
        lp1 <- mmat1 %*% param[1:(idx-1)]
        return(-param[idx] * numobs/2.0 - 1.0/(2.0 * exp(param[idx])) * sum((y - lp1)^2)
               + sum(z * lp2 - log(1.0 + exp(lp2))))
    }
    ##Joint log-likelihood function, version 2
    ##with fixed variance for the normal model
    jointloglik2 <- function(param, adjusted=TRUE) {
        idx <- 3 + adjusted * ncol(x) + ncol(psbasis)
        mmat2 <- cbind(1.0, lfx)
        lp2 <- mmat2 %*% param
        if (adjusted) {
            mmat1 <- cbind(1.0, z, x, predict(psbasis, expit(lp2)))
        } else {
            mmat1 <- cbind(1.0, z, predict(psbasis, expit(lp2)))
        }
        lp1 <- mmat1 %*% phi[1:(idx-1)]
        return(-phi[idx] * numobs/2.0 - 1.0/(2.0 * exp(phi[idx])) * sum((y - lp1)^2)
               + sum(z * lp2 - log(1.0 + exp(lp2))))
    }
    ##Simulation study
    nsim <- to - (from - 1)
    ##Correlation matrix
    s <- matrix(0.0, ncov, ncov)
    diag(s) <- 1.0
    ##Matrix used to store results - fitted parameters and variance of the estimates
    results <- matrix(NA, nsim, 89)

    for (i in 1:nsim) {
        iter <- from + (i - 1)
        set.seed(iter)

        ##Covariates are multinormal with unit variance and covariance 0, mean zero
        x <- rmvnorm(numobs, rep(0.0, ncov), s)
        lfx <- x
        ##Nonlinear transformation of the covariates to create confounding
        for(j in 1:ncov) {
            assign(paste("x", j, sep=""), x[,j])
            # assign(paste("fx", j, sep=""), as.vector(erf(abs(x[,j])/sqrt(2))))
            # assign(paste("lfx", j, sep=""), qnorm(get(paste("fx", j, sep=""))))
            assign(paste("lfx", j, sep=""), abs(x[,j])/sqrt(1-2/pi))
            # x[,j] <- get(paste("x", j, sep=""))
            # lfx[,j] <- get(paste("lfx", j, sep=""))
        }
        # Predictors for the treatment:
        {
        if (msps)
            lfx <- cbind(x1, x2, x3)
        else
            lfx <- cbind(lfx1, x2, x3)
        }
        # Predictors for the outcome:
        {
        if (msor)
            x <- cbind(x1, x2, x4)
        else
            x <- cbind(lfx1, x2, x4)
        }
        #apply(x, 2, sd)
        #apply(lfx, 2, sd)
        ##Exposure model and generation of assignment to treatment variable z
        pz <- expit(a0 + coeffa * lfx1 + coeffa * x2 + 2.0 * coeffa * x3)
        z <- (runif(numobs) < pz)
        ##Outcome model and generation of response variable y, via a linear model
        y <- b0 + effect * z + coeffb * lfx1 + coeffb * x2 + coeffb * x4 + rnorm(numobs)
        # py <- expit(b0 + effect * z + coeffb * lfx1 + coeffb * x2 + coeffb * x4)
        # y <- (runif(numobs) < py)

        ##[1] Average treatment effect (ATE)
        results[i,1] <- mean(b0 + effect * 1 + coeffb * lfx1 + coeffb * x2 + coeffb * x4) - mean(b0 + coeffb * lfx1 + coeffb * x2 + coeffb * x4)
        ##Propensity score model, using true DGM
        psmodel <- glm(z ~ lfx, family=binomial(link=logit))
        ##Naive propensity score using mean level
        pssmodel <- glm(z ~ 1, family=binomial(link=logit))
        # pssmodel <- glm(z ~ x, family=binomial(link=logit))
        ps <- as.numeric(expit(cbind(1, lfx) %*% coef(psmodel)))
        pss <- as.numeric(expit(rep(1, numobs) * coef(pssmodel)))
        # pss <- as.numeric(expit(cbind(1, x) %*% coef(pssmodel)))
        ##[2] IPTW estimator
        results[i,2] <- mean(z * (y/ps) - (1 - z) * (y/(1 - ps)))
 	    ##Default spline basis for regression:
        psbasis <- bs(ps, Boundary.knots=c(0,1))
        ##Clever-covariate
        cc <- z/ps - (1 - z)/(1.0 - ps)
        iptw <- z/ps + (1 - z)/(1.0 - ps)
        ##Unadjusted (naive) estimation
        model1 <- lm(y ~ z)
        results[i,3] <- coef(model1)[2]
        results[i,4] <- vcov(model1)[2,2]

        # plot(iptw)
        # coef(model1)[2]
	    ##Adjusted - with covariates
        model2 <- lm(y ~ z + x)
        results[i,5] <- coef(model2)[2]
        results[i,6] <- vcov(model2)[2,2]

        # model1 <- glm(y ~ z, family=binomial(link=logit))
        # summary(model1)
        # model1w <- glm(y ~ z, family=binomial(link=logit), weights=iptw)
        # summary(model1w)

        # model2 <- glm(y ~ z + x, family=binomial(link=logit))
        # summary(model2)
        # model2w <- glm(y ~ z + x, family=binomial(link=logit), weights=iptw)
        # summary(model2w)

	    ##Adjusted - with ps spline
        model3 <- lm(y ~ z + psbasis)
        results[i,7] <- coef(model3)[2]
        results[i,8] <- vcov(model3)[2,2]
	    ##Adjusted - doubly robust
        model4 <- lm(y ~ z + x + psbasis)
        results[i,9] <- coef(model4)[2]
        results[i,10] <- vcov(model4)[2,2]
	    ##Adjusted - clever covariates
        model5 <- lm(y ~ z + cc)
        y1 <- cbind(1, 1, 1/ps) %*% coef(model5)
        y0 <- cbind(1, 0, -1/(1-ps)) %*% coef(model5)
        results[i,11] <- mean(y1 - y0)
	    ##Adjusted - covariates plus clever covariates
        model6 <- lm(y ~ z + x + cc)
        y1 <- cbind(1, 1, x, 1/ps) %*% coef(model6)
        y0 <- cbind(1, 0, x, -1/(1-ps)) %*% coef(model6)
        results[i,12] <- mean(y1 - y0)

        # Adjusted variances:
        ##Optimize with constant weights,using MLE estimates as starting value for maximization
        ##Model: y~z+psbasis
        maxim <- optim(c(coef(model3), log(summary(model3)$sigma)), fn=loglik, gr=gradient, resp=y, mmat=cbind(1.0, z, psbasis), weight=rep(1.0, numobs),
                       control=list(fnscale=-1), method='BFGS', hessian=TRUE)
        ##Estimated variance-covariance matrix of the outcome model
        ##E[-U_i^{\phi \phi}(phi_0; \gamma_0)]^{-1}
        vcm <- solve(-maxim$hessian)
        ##E[-U_i^{\gamma \gamma}(\gamma_0)]^{-1}
        psvcm <- vcov(psmodel)
        ##E[U_i^{\phi \gamma}(\phi_0; \gamma_0)]
        gdm <- gd(maxim$par, coef(psmodel), adjusted=FALSE)
        scov <- scorecovariance(maxim$par, y, cbind(1.0, z, psbasis), weight=rep(1.0, numobs))
        ##E[U_i^{\phi}(\phi_0; \gamma_0)U_i^{\gamma}(\gamma_0)^\top]
        scp <- scorecrossprod(maxim$par, coef(psmodel), y, z, cbind(1.0, z, psbasis), cbind(1.0, lfx), rep(1.0, numobs), rep(1.0, numobs))
        ##estimated variance based on Hessian
        vars <- diag(vcm)
        ##robust variance
        rvars <- diag(vcm %*% scov %*% vcm)
        ##adjusted robust variance (see Appendix)
        # arvars <- diag(vcm) +
        arvars <- diag(vcm %*% scov %*% vcm)
                  diag(vcm %*% scp %*% psvcm %*% t(gdm) %*% vcm) +
                  diag(vcm %*% gdm %*% psvcm %*% t(scp) %*% vcm) +
                  diag(vcm %*% gdm %*% psvcm %*% t(gdm) %*% vcm)
        results[i,13] <- maxim$par[2]
        results[i,14] <- vars[2]
        results[i,15] <- rvars[2]
        results[i,16] <- arvars[2]

        ##Similar as above
        ##Model: y~z+x+psbasis
        maxim <- optim(c(coef(model4), log(summary(model4)$sigma)), fn=loglik, gr=gradient, resp=y, mmat=cbind(1.0, z, x, psbasis), weight=rep(1.0, numobs),
                       control=list(fnscale=-1), method='BFGS', hessian=TRUE)
        vcm <- solve(-maxim$hessian)
        psvcm <- vcov(psmodel)
        gdm <- gd(maxim$par, coef(psmodel), adjusted=TRUE)
        scov <- scorecovariance(maxim$par, y, cbind(1.0, z, x, psbasis), weight=rep(1.0, numobs))
        scp <- scorecrossprod(maxim$par, coef(psmodel), y, z, cbind(1.0, z, x, psbasis), cbind(1.0, lfx), rep(1.0, numobs), rep(1.0, numobs))
        vars <- diag(vcm)
        rvars <- diag(vcm %*% scov %*% vcm)
        ##asymptotic variance of Bayesian estimator discussed in appendix
        # arvars <- diag(vcm) +
        arvars <- diag(vcm %*% scov %*% vcm) +
                  diag(vcm %*% scp %*% psvcm %*% t(gdm) %*% vcm) +
                  diag(vcm %*% gdm %*% psvcm %*% t(scp) %*% vcm) +
                  diag(vcm %*% gdm %*% psvcm %*% t(gdm) %*% vcm)

        results[i,17] <- maxim$par[2]
        results[i,18] <- vars[2]
        results[i,19] <- rvars[2]
        results[i,20] <- arvars[2]
        # Joint estimation:
        ##Iterative estimation of the exposure and outcome models
        ##phi is a placeholder for const+z+psbasis+var
        ##gamma is vector of ps model parameters
        phi <- rep(0.0, 3 + ncol(x) + ncol(psbasis))
        gamma <- rep(0.0, 1 + ncol(lfx))
        # phi <- coef(model5)
        # gamma <- coef(psmodel)
        # jointloglik(c(phi, gamma))

        maxopiter <- 100
        eps <- 1e-05
        for (opiter in 1:maxopiter) {
            maxim <- optim(phi, fn=jointloglik1, control=list(fnscale=-1), method='BFGS')
            maxfound <- (maxim$convergence == 0)
            if (!maxfound)
                stop('No convergence.')
            phi <- maxim$par

            maxim <- optim(gamma, fn=jointloglik2, control=list(fnscale=-1), method='BFGS')
            maxfound <- (maxim$convergence == 0)
            if (!maxfound)
                stop('No convergence.')
            gamma <- maxim$par

            lik <- jointloglik(c(phi, gamma))
            # cat('opiter=', opiter, '\n', sep='')
            # cat('phi:', phi, '\n', sep=' ')
            # cat('gamma:', gamma, '\n', sep=' ')
            # cat('lik:', lik, '\n', sep=' ')
            if (opiter > 1) {
                if (abs(likold - lik) < eps)
                    break
            }
            likold <- lik
        }
        hessian <- optim(par=c(phi, gamma), fn=jointloglik, gr=NULL,
                        control=list(fnscale=-1, maxit=0), method='BFGS', hessian=TRUE)$hessian
        # cbind(c(pi, alpha, beta, theta), sqrt(diag(solve(-hessian))))
        results[i,21] <- phi[2]
        results[i,22] <- diag(solve(-hessian))[2]

        # Bayesian estimators and bootstrap:

        set.seed(1)
        pe <- matrix(NA, nrep, 2)
        ve <- matrix(NA, nrep, 2)
        pp <- matrix(NA, nrep, 2)
        bt <- matrix(NA, nrep, 2)
        btiptw <- matrix(NA, nrep, 2)
        pd <- matrix(NA, nrep, 2)
        btcc <- matrix(NA, nrep, 2)
        btdr <- matrix(NA, nrep, 2)
        dr <- rep(NA, nrep)
        for (j in 1:nrep) {
            xi <- as.numeric(rdirichlet(1, rep(1.0, numobs)))

            psgmodel <- glm(z ~ lfx, family=binomial(link=logit), weight=numobs * xi)
            ##Naive propensity score using mean level
            # pssgmodel <- glm(z ~ 1, family=binomial(link=logit), weight=xi)
            pssgmodel <- glm(z ~ 1, family=binomial(link=logit))

            psg <- as.numeric(expit(cbind(1, lfx) %*% coef(psgmodel)))
            pssg <- as.numeric(expit(rep(1, numobs) * coef(pssgmodel)))

            # alpha <- as.numeric(rmvnorm(1, coef(pssmodel), vcov(pssmodel)))
            # gamma <- as.numeric(rmvnorm(1, coef(psmodel), vcov(psmodel)))
            # pssg <- expit(cbind(1, x) %*% alpha)
            ##mean proportion of treated
            # pssg <- expit(rep(1.0, numobs) * alpha)
        	##individual predicted propensity scores
            # psg <- expit(cbind(1, lfx) %*% gamma)
        	##numerator and denominator of likelihood weights
            jw <- (z * pssg + (1 - z) * (1.0 - pssg))/(z * psg + (1 - z) * (1.0 - psg))

			##spline basis on population proportion of treated in sample
            psbasisg <- predict(psbasis, psg)

            model3g <- lm(y ~ z + psbasisg)
            #point estimate and variance estimate of the Normal approximation
            pe[j,1] <- coef(model3g)[2]
            ve[j,1] <- vcov(model3g)[2,2]
            ##draw from posterior predictive distribution
            pp[j,1] <- as.numeric(rmvnorm(1, coef(model3g), vcov(model3g)))[2]

            model4g <- lm(y ~ z + x + psbasisg)
            pe[j,2] <- coef(model4g)[2]
            ve[j,2] <- vcov(model4g)[2,2]
            pp[j,2] <- as.numeric(rmvnorm(1, coef(model4g), vcov(model4g)))[2]

            # Importance sampling estimators (Dirichlet):

            jw <- xi * jw
            model11 <- lm(y ~ z, weights=jw)
            y1 <- cbind(rep(1.0, numobs), 1) %*% coef(model11)
            y0 <- cbind(rep(1.0, numobs), 0) %*% coef(model11)
            pd[j,1] <- mean(y1 - y0)

            model11 <- lm(y ~ z + x, weights=jw)
            y1 <- cbind(1, 1, x) %*% coef(model11)
            y0 <- cbind(1, 0, x) %*% coef(model11)
            pd[j,2] <- mean(y1 - y0)

            model11 <- lm(y ~ z + x, weights=numobs * xi)
            m <- cbind(1, z, x) %*% coef(model11)
            y1 <- cbind(1, 1, x) %*% coef(model11)
            y0 <- cbind(1, 0, x) %*% coef(model11)
            cc <- z/psg - (1 - z)/(1.0 - psg)
            dr[j] <- sum(xi * (y - m) * cc) + sum(xi * (y1 - y0))

	        ##nonparametric bootstrap: drawn indices and subsamples
            bootidx <- sample(1:numobs, numobs, replace=TRUE)
            lfxb <- lfx[bootidx,]
            xb <- x[bootidx,]
            zb <- z[bootidx]
            yb <- y[bootidx]

            psbmodel <- glm(zb ~ lfxb, family=binomial(link=logit))
            pssbmodel <- glm(zb ~ 1, family=binomial(link=logit))
            psb <- expit(cbind(1, lfxb) %*% coef(psbmodel))
            pssb <- as.numeric(expit(rep(1, numobs) * coef(pssbmodel)))
            bw <- (zb * pssb + (1 - zb) * (1.0 - pssb))/(zb * psb + (1 - zb) * (1.0 - psb))

            model2b <- lm(yb ~ zb + xb)
            psbasisb <- predict(psbasis, psb)
            model3b <- lm(yb ~ zb + psbasisb)
            bt[j,1] <- coef(model3b)[2]
            model4b <- lm(yb ~ zb + xb + psbasisb)
            bt[j,2] <- coef(model4b)[2]

            psbasisb <- predict(psbasis, psb)
            model3b <- lm(yb ~ zb, weights=bw)
            btiptw[j,1] <- coef(model3b)[2]
            model4b <- lm(yb ~ zb + xb, weights=bw)
            btiptw[j,2] <- coef(model4b)[2]

            # Clever covariate estimators:

            ccb <- zb/psb - (1 - zb)/(1.0 - psb)

            model5b <- lm(yb ~ zb + ccb)
            y1 <- cbind(1, 1, 1/psb) %*% coef(model5b)
            y0 <- cbind(1, 0, -1/(1-psb)) %*% coef(model5b)
            btcc[j,1] <- mean(y1 - y0)

            model6b <- lm(yb ~ zb + xb + ccb)
            y1 <- cbind(1, 1, xb, 1/psb) %*% coef(model6b)
            y0 <- cbind(1, 0, xb, -1/(1-psb)) %*% coef(model6b)
            btcc[j,2] <- mean(y1 - y0)

            # IPT-weighted estimators:

            btdr[j,1] <- mean(zb * (yb/psb) - (1 - zb) * (yb/(1 - psb)))
            y1 <- cbind(1, 1, xb) %*% coef(model2b)
            y0 <- cbind(1, 0, xb) %*% coef(model2b)
            btdr[j,2] <- mean((yb * zb - (zb - psb) * y1)/psb - (yb * (1 - zb) - ((1 - zb) - (1 - psb)) * y0)/(1 - psb))

            # print(j)
        }
        ##model 3 (spline basis) and model 4 (doubly robust results
        ##Normal approximation
        results[i,23:24] <- colMeans(pe)
        ##variance decomposition for point estimate
        results[i,25:26] <- colMeans(ve) + apply(pe, 2, var)
        ##posterior predictive
        results[i,27:28] <- colMeans(pp)
        results[i,29:30] <- apply(pp, 2, var)
        results[i,31:34] <- as.numeric(apply(pp, 2, quantile, probs=c(0.025,0.975)))
        ##nonparametric bootstrap, propensity score
        results[i,35:36] <- colMeans(bt)
        results[i,37:38] <- apply(bt, 2, var)
        results[i,39:42] <- as.numeric(apply(bt, 2, quantile, probs=c(0.025,0.975)))
        ##nonparametric bootstrap, IPTW regression
        results[i,57:58] <- colMeans(btiptw)
        results[i,59:60] <- apply(btiptw, 2, var)
        results[i,61:64] <- as.numeric(apply(btiptw, 2, quantile, probs=c(0.025,0.975)))
        ##importance sampling estimator
        results[i,49:50] <- colMeans(pd)
        results[i,51:52] <- apply(pd, 2, var)
        results[i,53:56] <- as.numeric(apply(pd, 2, quantile, probs=c(0.025,0.975)))
        ##nonparametric bootstrap, clever covariate:
        results[i,70:71] <- colMeans(btcc)
        results[i,72:73] <- apply(btcc, 2, var)
        results[i,74:77] <- as.numeric(apply(btcc, 2, quantile, probs=c(0.025,0.975)))
        ##nonparametric bootstrap, IPTW estimators:
        results[i,78:79] <- colMeans(btdr)
        results[i,80:81] <- apply(btdr, 2, var)
        results[i,82:85] <- as.numeric(apply(btdr, 2, quantile, probs=c(0.025,0.975)))
        ##importance sampling estimator/DR:
        results[i,86] <- mean(dr)
        results[i,87] <- var(dr)
        results[i,88:89] <- as.numeric(quantile(dr, probs=c(0.025,0.975)))

        # Models with true propensity scores known:

        pzbasis <- predict(psbasis, pz)
        cc0 <- z/pz - (1 - z)/(1.0 - pz)

        model7 <- lm(y ~ z + pzbasis)
        results[i,43] <- coef(model7)[2]
        results[i,44] <- vcov(model7)[2,2]

        model8 <- lm(y ~ z + x + pzbasis)
        results[i,45] <- coef(model8)[2]
        results[i,46] <- vcov(model8)[2,2]

        model9 <- lm(y ~ z + cc0)
        y1 <- cbind(1, 1, 1/pz) %*% coef(model9)
        y0 <- cbind(1, 0, -1/(1-pz)) %*% coef(model9)
        results[i,47] <- mean(y1 - y0)

        model10 <- lm(y ~ z + x + cc0)
        y1 <- cbind(1, 1, x, 1/pz) %*% coef(model10)
        y0 <- cbind(1, 0, x, -1/(1-pz)) %*% coef(model10)
        results[i,48] <- mean(y1 - y0)

        # Oracle outcome model:

        model12 <- lm(y ~ z + lfx1 + x2 + x4)
        results[i,65] <- coef(model12)[2]
        results[i,66] <- vcov(model12)[2,2]

        # IPT weighted doubly robust estimator:

        y1 <- cbind(1, 1, x) %*% coef(model2)
        y0 <- cbind(1, 0, x) %*% coef(model2)
        results[i,67] <- mean((y * z - (z - ps) * y1)/ps - (y * (1 - z) - ((1 - z) - (1 - ps)) * y0)/(1 - ps))
        #Empirical sandwich method - variance estimate (formula (18) in Lunceford and Davidian (2004)
        Einv<-apply(sapply(1:numobs, function(x){ps[x]*(1-ps[x])*cbind(1,lfx)[x,]%*%t(cbind(1,lfx)[x,])},simplify="array"),c(1,2),mean)
        H<-colMeans((z*y*(1-ps)/ps+(1-z)*y*ps/(1-ps))*cbind(1,lfx))
        HEinv<-H%*%Einv
        I.IPW<-y*z/ps-(1-z)*y/(1-ps)-results[i,2]-(z-ps)*sapply(1:numobs,function(x){HEinv%*%cbind(1,lfx)[x,]})
        #Variance estimate for IPTW
        results[i,68] <- (1/numobs^2)*sum(I.IPW^2)
        #Empirical sandwich method - variance estimate (formula (21) in Lunceford and Davidian (2004)
        I.DR<-(z*y-y1*(z-ps))/ps-((1-z)*y+y0*(z-ps))/(1-ps)-results[i,67]
        #Variance estimate for IPTW-DR
        results[i,69] <- (1/numobs^2)*sum(I.DR^2)
        if (i %% 1 == 0) {
            print(iter)
        }
    }
    write.table(cbind(from:to, results), paste(outpath, 'results', numobs, no, '.txt', sep=''), row.names=FALSE, col.names=FALSE)
	return(NULL)
}
##End of function
##Number of core (must be changed above also)
ncores <- 7

##Parallel computations for the function
# mclapply(1:ncores,pssim,mc.cores=ncores)
# mclapply(1:ncores, pssim, outpath="~/Dropbox/work/ps_olli/data/model1/", msor=TRUE, msps=FALSE, mc.cores=ncores)
mclapply(1:ncores, pssim, outpath=outpath, msor=FALSE, msps=TRUE, mc.cores=ncores)



