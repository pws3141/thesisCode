# forecasting block-average Ornstein-Uhlenbeck rocess,
# given some drift and diffusion values


# require .sigmaMatrix.baOU function to obtain covariance matrix
# given some drift and diffusion
source("../ch2/parameterFitting/mle_ornsteinU_functions.R")

# one-step-ahead point forecast and prediction interval
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
