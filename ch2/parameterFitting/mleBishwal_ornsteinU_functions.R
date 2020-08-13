.stratonovichIntegral <- function(X, fnX) {
        n <- length(X)
        midpointf <- 0.5 * (fnX[1:(n-1)] + fnX[2:n])
        differenceX <- X[2:n] - X[1:(n-1)]
        S1n <- sum(midpointf * differenceX)
        S1n
}

.mleBishwalGammaOU <- function(X, t) {
        # estimate MLE of gamme in 
        # dX = \gamma (\mu - X)dt + dW
        # with discrete obs from t0 = 0 to tn = T
        stopifnot(length(t) == length(X), 
                    t[1] == 0)
        n <- length(X)
        tn <- t[n]
        XT <- X[n]
        x0 <- X[1]
        stratInt <- .stratonovichIntegral(X = X, fnX = X)
        topOne <- (XT - x0) * sum(X[1:(n-1)] * (t[2:n] - t[1:(n-1)]))
        topTwo <- tn * (stratInt - 0.5 * tn) 
                        
        bottomOne <- tn * sum(X[1:(n-1)]^2 * (t[2:n] - t[1:(n-1)]))
        bottomTwo <- sum(X[1:(n-1)] * (t[2:n] - t[1:(n-1)]))^2
        gammaHat <- (topOne - topTwo) / (bottomOne - bottomTwo)
        gammaHat
}

.mleBishwalMuOU <- function(X, t) {
        # estimate MLE of gamme in 
        # dX = \gamma (\mu - X)dt + dW
        # with discrete obs from t0 = 0 to tn = T
        stopifnot(length(t) == length(X), 
                    t[1] == 0)
        n <- length(X)
        tn <- t[n]
        XT <- X[n]
        x0 <- X[1]
        stratInt <- .stratonovichIntegral(X = X, fnX = X)
        topOne <- (XT - x0) * sum(X[1:(n-1)]^2 * (t[2:n] - t[1:(n-1)]))
        topTwoOne <- stratInt - 0.5 * tn
        topTwoTwo <- sum(X[1:(n-1)] * (t[2:n] - t[1:(n-1)]))
        topTwo <- topTwoOne * topTwoTwo
        bottomOne <- (XT - x0) * sum(X[1:(n-1)] * (t[2:n] - t[1:(n-1)]))
        bottomTwo <- tn * (stratInt - 0.5 * tn)
        muHat <- (topOne - topTwo) / (bottomOne - bottomTwo)
        muHat
}


mleBishwalOUParameterEstimation <- function(X, t) {
        stopifnot(length(t) == length(X), 
                    t[1] == 0)
        n <- length(X)
        tn <- t[n]
        x0 <- X[1]
        XT <- X[n]
        gammaHat <- .mleBishwalGammaOU(X = X, t = t)
        muHat <- .mleBishwalMuOU(X = X, t = t)
        res <- list(gammaHat = gammaHat, 
                    muHat = muHat)
        res
}

