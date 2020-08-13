# Assume we have n observations from N(0, \Sigma_n)
# Aim: find MLE:
# \Sigma_n is a symmetric Toeplitz matrix of the form
# S = (s_{ij}) and 
# s_{ij} = a^{\lvert j - i \rvert} + b \mathbbm{1}_{i = j}

# have observations x_1, ..., x_n
# obtain covariance matrix \Sigma_n
# optimise to find gH and s2H

.sigmaMatrix.baOU <- function(g, s2, n, delta.t) {
        # form: c a^{abs(i - j)} + b 1_{i = j}
        c <- (s2 / (2 * g^3)) * (exp(- g * delta.t) + exp(g * delta.t) - 2)
        a <- exp(-g * delta.t)
        b <- (s2 / (2 * g^3)) * (exp(- g * delta.t) - 
                               exp(g * delta.t) + 2 * g * delta.t)
        rowOne.1 <- (s2 / g^3) * (exp(-g * delta.t) + g * delta.t - 1)
        # powers
        j <- 1:(n - 1)
        rowOne <- c * a^j
        rowOne <- c(rowOne.1, rowOne)
        # old: numerical errors this way,
        # as c and b and v. large and c + b \approx 0
        #rowOne <- rowOne + c(b, rep(0, times = n - 1))
        res <- toeplitz(rowOne)
        res
}

.logLikelihoodGaussian <- function(x, S.int) {
        # x is single observation from multivariate Gaussian
        # S.int is .sigmaMatrix.baOU
        U.chol <- chol(S.int)
        log_detS <- 2 * sum(log(diag(U.chol)))
        r2 <- sum(backsolve(U.chol, x, transpose = TRUE) ^ 2)
        #r1 <- t(x) %*% solve(A, x)
        #print(abs(r1 - r2)) # should be zero
        ll <- -0.5 * log_detS - 0.5 * r2
        names(ll) <- "ll"
        ll
}

mle.baOU <- function(x, delta.t, gamma0, sigma20) {
        # x is matrix, with columns time series
        if (is.null(nrow(x))) x <- as.matrix(x)
        n <- nrow(x)
        res <- optim(par = c(gamma0, sigma20), function(gs = c(g, s2)) {
                             g <- abs(gs[1])
                             if (g == 0) g <- 0.01
                             s2 <- gs[2]
                             #message(sprintf("g = %.2f; s2 = %.2f", g, s2))
                             S.int <- .sigmaMatrix.baOU(g = g, s2 = s2,
                                                       n = n, delta.t = delta.t)
                             ll <- .logLikelihoodGaussian(x = x, S.int = S.int)
                             #S.mle <- .sigmaMLEMatrix(x = x)
                             #opt.min <- .integratedOptimFunction(S = S.int,
                                                                 #S.mle = S.mle)
                             -ll
                                },
                        #method = "CG", control = list(maxit = 500))
                        method = "Nelder-Mead")
                        #method = "L-BFGS-B", lower = c(1e-4, 1e-4), upper = c(NA, NA))
        res$par[1] <- abs(res$par[1])
        res
}

baML.forecast <- function(x, n, g, s2, delta.t, alpha = 0.05) {
        # find point forecast and PI for time (n \dt + \dt)
        # assuming:
        # x is realisation from block-average OU
        # g, s2 are parameter estimates using x from 'mle.baOUFunctions.R'
        # delta.t is time interval of x realisations
        # alpha is 100 (1 - alpha)% PI
        x <- x[1:n]
        S <- .sigmaMatrix.baOU(g = g, s2 = s2, n = (n + 1), delta.t = delta.t)
        S11 <- S[1:n, 1:n]
        S21 <- S[n + 1, 1:n]
        S12 <- S21 # as Toeplitz matrix
        #S12 <- S[1:n, n + 1]
        S22 <- S[n + 1, n + 1]
        U.chol <- chol(S11)
        # point forecast
        y_1 <- backsolve(U.chol, S21, transpose = TRUE)
        y_2 <- backsolve(U.chol, x, transpose = TRUE)
        x.nh <- sum(y_1 * y_2)
        #x.nh_test <- sum(S21 * solve(S11, x)
        #print(abs(x.nh_test - x.nh)) # should be zero
        # prediction interval
        # y_3 <- backsolve(U.chol, S12, transpose = TRUE)
        # y_3 = y_1 as S12 = S21, as S Toeplitz
        pi.nh_var <- S22 - sum(y_1^2)
        z <- qnorm(alpha / 2, lower.tail = FALSE)
        pi.nh <- c(x.nh - z * sqrt(pi.nh_var), x.nh + z * sqrt(pi.nh_var))
        # out
        res <- list(pF = x.nh, PI = pi.nh)
        res
}

TEST <- exists("TEST") && isTRUE(TEST)
if (isTRUE(TEST)) { # baML.forecast test {{{
        g <- 1
        s2 <- 100
        delta_t <- 1
        n <- 10
        S <- .sigmaMatrix.baOU(g = g, s2 = s2, n = n + 1, delta.t = delta_t)
        S11 <- S[1:n, 1:n]
        S21 <- S[n + 1, 1:n]
        S12 <- S[1:n, n + 1]
        S22 <- S[n + 1, n + 1]
        message("S12 equal to S21...")
        switch(isTRUE(all.equal(S12, S21)) + 1, message("FALSE"), message("TRUE: S12 == S21"))
        # realisation
        x <- MASS::mvrnorm(n = 1, mu = rep(0, n), Sigma = S11)

        ### test
        # test pF inverse correct
        U.chol <- chol(S11)
        # point forecast
        y_1 <- backsolve(U.chol, S21, transpose = TRUE)
        y_2 <- backsolve(U.chol, x, transpose = TRUE)
        x.nh <- sum(y_1 * y_2)
        x.nh_test <- sum(S21 * solve(S11, x))
        message("Inverse using Choleski and backsolve, compared to 'solve' directly")
        message("\tequal to zero if true...")
        print(abs(x.nh_test - x.nh)) # should be zero

        # prediction interval
        str.tmp <- "backsolve(U.chol, S12, transpose = TRUE)"
        message(sprintf("checking whether %s required", str.tmp))
        y_3 <- backsolve(U.chol, S12, transpose = TRUE)
        message("\tif 'TRUE' then not required")
        switch(isTRUE(all.equal(y_1, y_3)) + 1, message("FALSE"), message("TRUE: y_1 == y_3"))
}#}}}

