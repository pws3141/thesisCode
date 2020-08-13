# Aim: calulate likelihood and estimate MLE using optim
# of OU process
# dX = - gamma X dt + sigma dB

# - 1D time series
# - breaks
# - spline (a, b, c) values for each break
# - gamma0
# - sigma0^2
# Output: 
# - gamma and sigma^2 estimates
likelihoodOU <- function(X, t, g, s2) {
        n <- length(X)
        deltaT <- t[2] - t[1] # assume equal spacing
        termOne <- (g / (pi * s2 * (1 - exp(-2 * g * deltaT))))^(n / 2)
        termTwo1 <- (-g / (s2 * (1 - exp(-2 * g * deltaT))))
        termTwo2 <- sum((X[2:n] - X[1:(n-1)] * exp(-g * deltaT))^2)
        termTwo <- exp(termTwo1 * termTwo2)
        res <- termOne * termTwo
        res
}

logLikelihoodOU <- function(X, t, g, s2) {
        n <- length(X)
        deltaT <- t[2] - t[1] # assume equal spacing
        termOne <- (n / 2) * (log(g) - 
                              log(pi * s2 * (1 - exp(-2 * g * deltaT))))
        termTwo <- g / (s2 * (1 - exp(-2 * g * deltaT))) *
                    sum((X[2:n] - X[1:(n-1)] * exp(-g * deltaT))^2)
        logL <- termOne - termTwo
        logL
}

if (FALSE) {
.dldgamma <- function(X, t, g, s2) {
        n <- length(X)
        deltaT <- t[2] - t[1] # assume equal spacing
        termOne <- (n / 2) * (1 / g * 
                              (2 * deltaT * exp(-2 * g * deltaT)) /
                              (1 - exp(-2 * g * deltaT)))
        termTwo1 <- (1 - (1 + 2 * deltaT * g) * exp(-2 * deltaT * g)) * 
                        sum((X[2:n] - X[1:(n - 1)] * exp(- deltaT * g))^2)
        termTwo2 <- 2 * g * deltaT * exp(-g * deltaT) *
                        sum((X[2:n] - X[1:(n - 1)] * exp(-g * deltaT)) * X[1:(n - 1)])
        termTwo <- 1 / (s2 * (1 - exp(-2 * g * deltaT))) * 
                       (termTwo1 + termTwo2)
        res <- termOne - termTwo
        res
}
}

.dldgamma <- function(X, t, g, s2) {
        n <- length(X)
        deltaT <- t[2] - t[1] # assume equal spacing
        termOne <- (n / 2) * 
                        (2 * g * deltaT + 1 - exp(2 * g * deltaT)) /
                        (g * (1 - exp(2 * g * deltaT)))
        termTwo1 <- exp(2 * g * deltaT) *
                        (exp(2 * g * deltaT) - 1 - 2 * g * deltaT) /
                        (s2 * (exp(2 * g * deltaT) - 1)^2) *
                        sum((X[2:n] - X[1:(n - 1)] * exp(-g * deltaT))^2)
        termTwo2 <- (g / (s2 * (1 - exp(-2 * g * deltaT)))) *
                        2 * deltaT * exp(-2 * g * deltaT) *
                        sum(X[1:(n - 1)] * (X[2:n] * exp(g * deltaT) - X[1:(n - 1)]))
        res <- termOne - (termTwo1 + termTwo2)
        res
}

.dldsigma2 <- function(X, t, g, s2) {
        n <- length(X)
        deltaT <- t[2] - t[1] # assume equal spacing
        termOne <- - n / (2 * s2)
        termTwo <- g / (s2^2 * (1 - exp(-2 * g * deltaT))) *
                        sum((X[2:n] - X[1:(n - 1)] * exp(-g * deltaT))^2)
        res <- termOne + termTwo
        res
}

.sigmaHat2 <- function(X, t, g) {
        n <- length(X)
        deltaT <- t[2] - t[1] # assume equal spacing
        res <- 2 * g / (n * (1 - exp(-2 * g * deltaT))) * 
                sum((X[2:n] - X[1:(n - 1)] * exp(-g * deltaT))^2)
        names(res) <- "s2H"
        res
}

# TODO: make mleOptim work for non constant Delta t
.mleOptim <- function(X, t, gamma0) {
        n <- length(X)
        deltaT <- t[2] - t[1] # assume equal spacing
        if(length(gamma0) == 1) {
                opt <- optim(par = gamma0,
                             fn = function(g) {
                                     s2 <- .sigmaHat2(X = X, t = t, g = g)
                                     resTmp <- logLikelihoodOU(X = X, t = t, 
                                                               g = g, s = s2)
                                     #cat(g, " : ", s2, " : ")
                                     #cat(resTmp, "\n")
                                     resTmp
                              }, 
                             #gr = function(g) {
                                     #s2 <- .sigmaHat2(X = X, t = t, g = g)
                                     #dldg <- .dldgamma(X = X, t = t,
                                                       #g = g, s = s2)
                                     #cat(s2, " : ")
                                     #cat(dldg, "\n")
                                     #dldg
                             #}, 
                             method = "L-BFGS-B", lower = 1e-10, upper = Inf, 
                             control = list(fnscale = -1, maxit = 500))
        } else {
                optAll <- lapply(gamma0, function(par) {
                                opt <- optim(par = par, 
                                     fn = function(g) {
                                             s2 <- .sigmaHat2(X = X, t = t, g = g)
                                             resTmp <- logLikelihoodOU(X = X, t = t, 
                                                                       g = g, s = s2)
                                             resTmp
                                      }, 
                                     #gr = function(g) {
                                             #s2 <- .sigmaHat2(X = X, t = t, g = g)
                                             #dldg <- .dldgamma(X = X, t = t,
                                                               #g = g, s = s2)
                                             #dldg
                                     #}, 
                                     method = "L-BFGS-B", lower = 1e-10, upper = Inf, 
                                     control = list(fnscale = -1, maxit = 500))
                              })
                optValue <- sapply(optAll, function(o) o$value)
                whichMinValue <- which.min(optValue)
                opt <- optAll[[whichMinValue]]
        }
        sHat2 <- .sigmaHat2(X = X, t = t, g = opt$par)
        opt$par <- c(opt$par, sHat2) 
        names(opt$par) <- c("gammaHat", "sigmaHat2")
        opt
}

.mleSplineOptim <- function(Y, t, gamma0, spline, breaks) {
        n <- length(Y)
        resMuT <- .muT(t = t, spline = spline, breaks = breaks, deriv = FALSE)
        X <- Y - resMuT$muT
        opt <- .mleOptim(X = X, t = t, gamma0 = gamma0)
        opt
}

mleOUParameterEstimation <- function(X, t, spline = NULL, breaks = NULL, gamma0) {
        # if spline and breaks NULL, then observations from process
        # dX = -gamma X dt + sigma dB.
        # otherside, from Y = X + m(t), where m(t) quadratic spline
        if (is.null(spline) & is.null(breaks)) {
                res <- .mleOptim(X = X, t = t, gamma0 = gamma0)
        } else
        {
                if (!is.null(spline) & !is.null(breaks)) {
                res <- .mleSplineOptim(Y = X, t = t, gamma0 = gamma0,
                                       spline = spline, breaks = breaks)
                } else {
                        stop("Need to specify 'spline' and 'breaks', or neither")
                }
        }
        res
}

#### additional info and test code

TEST <- exists("TEST") && isTRUE(TEST)

if (TEST) {
message("
Sourced \'mleOUFunctions.R\', to find MLE of Ornstein-Uhlenbeck process.
")

#source("mleOUFunctions.R")
source("../eulerApproximation_ornsteinU_functions.R")
lenT <- 1e5 + 1
thinBy <- 100
t <- seq(from = 0, to = 100, length = lenT)
seqThin <- seq(from = 1, to = lenT, by = thinBy)
tThin <- t[seqThin]

g <- 1
s2 <- 1
sim <- .eulerSimulationOU(x0 = 0, t = t, g = g, s2 = s2, mu = 0)
yThin <- sim$y[seqThin]

message(sprintf("checking derivative dldg function and dlds2 for g = %.1f and s2 = %.1f...", g, s2))
# plot ll against gamma
gammaSeq <- seq(from = 0.1, to = 10, length = 1000)
llgSim <- sapply(gammaSeq, function(gam) {
                        logLikelihoodOU(X = yThin, t = tThin, g = gam, s2 = s2)
                })
dldgSim <- sapply(gammaSeq, function(gam) {
                .dldgamma(X = yThin, t = tThin, g = gam, s2 = s2)
                })

#message(sprintf("true gamma equals %d", g))
whichMaxLLg <- which.max(llgSim)
whichdldgNeg <- which(dldgSim < 0)
message(sprintf("LL maximum for gammaSeq equal %.3f", gammaSeq[whichMaxLLg]))
message(sprintf(".dldg crosses x-axis at gammaSeq equal %.3f", 
                gammaSeq[min(whichdldgNeg)]))

# plot ll against s2
s2Seq <- seq(from = 0.1, to = 10, length = 1000)
lls2Sim <- sapply(s2Seq, function(ss) {
                        logLikelihoodOU(X = yThin, t = tThin, g = g, s2 = ss)
                })
dlds2Sim <- sapply(s2Seq, function(ss) {
                .dldsigma2(X = yThin, t = tThin, g = g, s2 = ss)
                })

whichMaxLLs2 <- which.max(lls2Sim)
whichdlds2Neg <- which(dlds2Sim < 0)
#message(sprintf("true sigma2 equals %d", s2))
message(sprintf("LL maximum for s2Seq equal %.3f", s2Seq[whichMaxLLs2]))
message(sprintf(".dlds2 crosses x-axis at s2Seq equal %.3f", 
                s2Seq[min(whichdlds2Neg)]))

PLOT <- exists("PLOT") && isTRUE(PLOT)
message(sprintf("set PLOT <- TRUE for plot: currently PLOT = %s", PLOT))
if (PLOT) {
        par(mfrow = c(2, 2))
        plot(gammaSeq, llgSim, xlim = c(0, 5),
             main = "log-likelihood vs. gamma", t = "l")
        abline(v = gammaSeq[min(whichMaxLLg)], lty = 2)
        plot(gammaSeq, dldgSim, main = "dldg vs. gamma", 
             xlim = c(0, 5), ylim = c(-150, 150), t = "l")
        abline(h = dldgSim[min(whichdldgNeg)], lty = 2)
        abline(v = gammaSeq[min(whichdldgNeg)], lty = 2)

        plot(s2Seq, lls2Sim, xlim = c(0, 5), 
             main = "log-likelihood vs. sigma2", t = "l")
        abline(v = s2Seq[min(whichMaxLLs2)], lty = 2)
        plot(s2Seq, dlds2Sim, main = "dlds2 vs. sigma2",
             xlim = c(0, 5), ylim = c(-150,150), t = "l")
        abline(h = dlds2Sim[min(whichdlds2Neg)], lty = 2)
        abline(v = s2Seq[min(whichdlds2Neg)], lty = 2)
} #end if (PLOT)

#############
#########
#####
##
message("
checking maximum of log likelihood and zero of derivative dldg function match up... 
        (this may take some time)")
s2 <- 1
gg <- seq(from = 1, to = 101, by = 10)
ggLen <- length(gg)
ggInt <- gg[2] - gg[1]
message(sprintf("set sigma2 = %d and vary gamma from %.1f to %.1f by %.1f",
                s2, min(gg), max(gg), ggInt))
llMax <- numeric(ggLen)
dldgZero <- numeric(ggLen)
iter <- 100
message(sprintf("for each gamma between %.1f and %.1f, perform %d iterations:",
                min(gg), max(gg), iter))
message(sprintf("simulate %d OU processes using EM approximation,
thin chain, then calculate dldg for varying gamma.", iter))
message(sprintf("   for LL: choose gHat to be the value that maximises log-likelihood"))
message(sprintf("   for dldg: choose gHat to be the last value before dldg crosses x-axis
                "))

for(i in 1:ggLen) {
        gTmp <- gg[i]
        llMaxGG <- numeric(iter)
        dldgZeroGG <- numeric(iter)
        for (j in 1:iter) {
                sim <- .eulerSimulationOU(x0 = 0, t = t, g = gTmp, s2 = s2, mu = 0)
                yThin <- sim$y[seqThin]
                gammaSeq <- seq(from = 0.1, to = (max(gg) + 20), length = 1000)
                llgSim <- sapply(gammaSeq, function(gam) {
                                logLikelihoodOU(X = yThin, t = tThin, g = gam, s2 = s2)
                                })
                dldgSim <- sapply(gammaSeq, function(gam) {
                                .dldgamma(X = yThin, t = tThin, g = gam, s2 = s2)
                                })
                whichMaxLLg <- which.max(llgSim)
                whichdldgNeg <- which(dldgSim < 0)
                llMaxGG[j] <- gammaSeq[whichMaxLLg]
                dldgZeroGG[j] <- gammaSeq[min(whichdldgNeg)]
        }
        llMax[i] <- mean(llMaxGG)
        dldgZero[i] <- mean(dldgZeroGG)
}
resG <- cbind(gg, llMax, dldgZero)
colnames(resG) <- c("g", "ll: gHat", "dldg: gHat")
message("approximations for gamma are given by:")
print(resG)

#############
#########
#####
##
message("
checking maximum of log likelihood and zero of derivative dlds2 function match up... 
        (this may take some time)")
g <- 1
ss2 <- seq(from = 1, to = 101, by = 10)
ss2Len <- length(ss2)
ss2Int <- ss2[2] - ss2[1]
message(sprintf("set gamma = %d, and vary sigma2 from %.1f to %.1f by %.1f",
                g, min(ss2), max(ss2), ss2Int))
llMax <- numeric(ss2Len)
dlds2Zero <- numeric(ss2Len)
iter <- 100
message(sprintf("for each sigma2 between %.1f and %.1f, perform %d iterations:",
                min(ss2), max(ss2), iter))
message(sprintf("simulate %d OU processes using EM approximation,
thin chain, then calculate dlds2 for varying sigma2", iter))
message(sprintf("   for LL: choose s2Hat to be the value that maximises log-likelihood"))
message(sprintf("   for dlds2: choose s2Hat to be the last value before dlds2 crosses x-axis
                "))

for(i in 1:ss2Len) {
        s2Tmp <- ss2[i]
        llMaxSS <- numeric(iter)
        dlds2ZeroSS <- numeric(iter)
        for (j in 1:iter) {
                sim <- .eulerSimulationOU(x0 = 0, t = t, g = g, s2 = s2Tmp, mu = 0)
                yThin <- sim$y[seqThin]
                s2Seq <- seq(from = 0.1, to = (max(ss2) + 20), length = 1000)
                lls2Sim <- sapply(s2Seq, function(ss) {
                                logLikelihoodOU(X = yThin, t = tThin, g = g, s2 = ss)
                })
                dlds2Sim <- sapply(s2Seq, function(ss) {
                                .dldsigma2(X = yThin, t = tThin, g = g, s2 = ss)
                })
                whichMaxLLs2 <- which.max(lls2Sim)
                whichdlds2Neg <- which(dlds2Sim < 0)
                llMaxSS[j] <- s2Seq[whichMaxLLs2]
                dlds2ZeroSS[j] <- s2Seq[min(whichdlds2Neg)]
        }
        llMax[i] <- mean(llMaxSS[j])
        dlds2Zero[i] <- mean(dlds2ZeroSS)
}
resS <- cbind(ss2, llMax, dlds2Zero)
colnames(resS) <- c("s2", "ll: s2Hat", "dlds2: s2Hat")
message("approximations for sigma2 are given by:
        ")
print(resS)

##########
######
###
lenT <- 1e6 + 1
thinBy <- 100
t <- seq(from = 0, to = 100, length = lenT)
seqThin <- seq(from = 1, to = lenT, by = thinBy)
tThin <- t[seqThin]

s2 <- 1
gg <- seq(from = 1, to = 101, by = 10)
message(sprintf("
Checking convergence of mleOptim: gamma from %d to %d, by %d
        ", min(gg), max(gg), gg[2] - gg[1]))
lenGG <- length(gg)
cvg <- numeric(lenGG)
gHat <- numeric(lenGG)
for (i in 1:lenGG) {
        gTmp <- gg[i]
        sim <- .eulerSimulationOU(x0 = 0, t = t, g = gTmp, s2 = s2, mu = 0)
        yThin <- sim$y[seqThin]
        res <- .mleOptim(X = yThin, t = tThin, gamma0 = 1)
        cvg[i] <- res$convergence
        gHat[i] <- res$par[1]
}
message("optim convergence (0 == converged):
        ")
res <- rbind(gg, gHat, cvg)
rownames(res) <- c("gamma", "gammaHat", "convergence")
print(res)

} ## end if (TEST)
