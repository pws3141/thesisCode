# Aim:
# Simulate SDE using Euler approximation.

# Input:
# - breaks
# either: 
        # - spline (a, b, c) values for each break
        # - values of m(t) for each time points t_i
# - gamma
# - sigma^2

# Output:
# - simulated SDE path

.eulerSimulationOUSpline <- function(x0, t, g, s2, 
                                            spline = NULL, breaks = NULL,
                                            muT = NULL, muTDer = NULL) {
        # approximate Orstein-Uhlenbeck process
        # dX = -\gamma X dt + \sigma dW
        # Y = X + m(t)
        # dY = -gamma (Y - m(t)) dt + m'(t) dt + sigma dW
        # by Euler-Maruyama approximation
        # t is time discretisation in [t0, T]
        # t0 = t_0 < t_1 < ... < t_N = T
        N <- length(t) - 1 
        y <- numeric(length = (N + 1))
        y[1] <- x0
        maxT <- max(t)
        deltaT <- t[2:(N+1)] - t[1:N]
        deltaW <- rnorm(N, mean = 0, sd = sqrt(deltaT))
        # split t up into relevant breaks
        if (is.null(muT) | is.null(muTDer)) {
                n <- length(breaks) - 1
                numBreaks <- nrow(spline)
                if (n != numBreaks) stop("breaks is not compatible with spline")
                muTList <- .muT(t = t, spline = spline, breaks = breaks, deriv = TRUE)
                muT <- muTList$muT
                muTDer <- muTList$muTDer
        } 
        ###
        s <- sqrt(s2)
        for(j in 1:N) {
                termOne <- -g * (y[j] - muT[j] - muTDer[j] / g) * deltaT[j]
                termTwo <- s * deltaW[j]
                y[j + 1] <- y[j] + termOne + termTwo
        }
        list(t = t, y = y, mu = muT)
}

.eulerSimulationOU <- function(x0, t, g, s2, mu) {
        # approximate Orstein-Uhlenbeck process
        # dX = \gamma (\mu - X) dt + \sigma dW
        # by Euler-Maruyama approximation
        # tau is time discretisation in [t0, T]
        # t0 = \tau_0 < \tau_1 < ... < \tau_N = T
        N <- length(t) - 1 
        y <- numeric(length = (N + 1))
        y[1] <- x0
        deltaT <- t[2:(N+1)] - t[1:N]
        deltaW <- rnorm(N, mean = 0, sd = sqrt(deltaT))
        s <- sqrt(s2)
        for(n in 1:N) {
                termOne <- g * (mu - y[n]) * deltaT[n]
                if (abs(termOne) > abs(y[n])) {
                        message(sprintf("y = %.3f, g*y*dt = %.3f", y[n], termOne))
                }
                termTwo <- s * deltaW[n]
                y[n+1] <- y[n] + termOne + termTwo
        }
        list(t = t, y = y)
}

eulerSimulationOU <- function(x0, t, g, s2, mu = 0, spline = NULL, breaks = NULL,
                                    muT = NULL, muTDer = NULL) {
        if (!(is.null(spline) & is.null(breaks)) | !(is.null(muT) & is.null(muTDer))) {
                res <- .eulerSimulationOUSpline(x0 = x0, t = t, g = g, s2 = s2,
                                                spline = spline, breaks = breaks,
                                                muT = muT, muTDer = muTDer)
        } else {
                res <- .eulerSimulationOU(x0 = x0, t = t, g = g, s2 = s2, mu = mu)
        }
        res
}
