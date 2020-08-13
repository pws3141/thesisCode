# Aim: estimate parameters of OU process
# dX = - gamma X dt + sigma dB
# using 'method of moments'
# i.e. given variance and lag 1 correlation

# Output: 
# - gamma and sigma^2 estimates

varOUProcess <- function(g, s2, t1, t2) {
        res <- (s2 / (2 * g)) * (1 - exp(-2 * g * (t2 - t1)))
        res
}

covarOUProcess <- function(g, s2, t1, t2) {
        # OU proce from time t1 to t2.
        res <- (s2 /  (2 * g)) * 
                (exp(-g * abs(t2 - t1)) - exp(-g * (t2 - t1)))
        res
}

corrOUProcess <- function(g, s2, t1, t2) {
        # OU proce from time t1 to t2.
        covarOU <- covarOUProcess(g = g, s2 = s2, t1 = t1, t2 = t2)
        varOU <- varOUProcess(g = g, s2 = s2, t1 = t1, t2 = t2)
        res <- covarOU / varOU
        res
}

moments <- function(X, t, verbose = FALSE) {
        #
        t0 <- t[1]
        stopifnot(t0 == 0)
        #
        tInt <- findInterval(t, 0:1)
        stopifnot(2 %in% tInt)
        tSize <- min(which(tInt == 2))
        t1 <- t[tSize]
        stopifnot(t1 == 1)
        if (verbose == TRUE) {
                message(sprintf("'t' in increments of %g, from %d to %.1f", 
                                t[2] - t0, t0, max(t)))
        }
        #
        len <- length(X)
        # correlation between points at distance t = 1
        X1Seq <- seq(from = 1, to = (len - tSize + 1), by = 1)
        X1 <- X[X1Seq]
        X2Seq <- seq(from = tSize, to = len, by = 1)
        X2 <- X[X2Seq]
        # variance and correlation
        v <- var(X)
        r <- cor(X1, X2)
        if (verbose == TRUE) {
                message(sprintf("corr = %.4f; var = %.4f", r, v))
        }
        resVR <- c(v, r)
        names(resVR) <- c("var", "cor")
        resX1X2 <- rbind(X1, X2)
        rownames(resX1X2) <- c("X1", "X2")
        list(vr = resVR, X = resX1X2)
}

momOUParameterEstimation <- function(X, t, verbose = FALSE, upper.range = FALSE) {
        #
        vr_X <- moments(X = X, t = t, verbose = verbose)
        vr <- vr_X$vr
        # variance and correlation
        v <- vr[1]
        r <- vr[2]
        if (isTRUE(upper.range)) {
                if (r <= 0) {
                        if (verbose == TRUE) {
                                message(sprintf("Correlation of chain is %.3f:", r),
                                sprintf(" must be positive for method of moment estimation..."),
                                sprintf("\n\tcalculating 95%% CI for correlation..."))
                        }
                        X1 <- vr_X$X[1, ]
                        X2 <- vr_X$X[2, ]
                        rr <- cor.test(X1, X2)
                        if (verbose == TRUE) {
                                message(sprintf(
                                "\tcorrelation = %.4f with CI (%.4f, %.4f)", 
                                rr$estimate, rr$conf.int[1], rr$conf.int[2]))
                                message(sprintf("\tchoosing correlation as upper value of CI, r = %.4f",
                                rr$conf.int[2]))
                                message(sprintf("\t\twe expect gammaHat and sigma2Hat to be less than the true values"))
                        }
                        #r <- rr$conf.int[2] 
                }
        } # end if upper range TRUE
        if (r <= 0) {
                mss <- paste0(sprintf("Upper CI value is %.4f:", r),
                                " must be positive for method of moment estimation...",
                                sprintf("\n\treturning NA for gammaHat and sigma2Hat") )
                warning(mss)
                res <- c(NA, NA)
        } 
        if (r > 0) {
                # parameter estimation 
                gHat <- - log(r)
                sHat2 <- -2 * v * log(r)
                res <- c(gHat, sHat2)
        } 
        names(res) <- c("gHat", "sHat2")
        res <- list(params = res, var = v, corr = r)
        res
}

### Testing

TEST <- exists("TEST") && isTRUE(TEST)

if (TEST) {
message("
Sourced \'methodOfMomentsParameterEstimationFunctions\', to estimate parameters
of an OU model: dX = -gamma X dt + sigma dB
")
library(testthat) 
context("momOUParameterEstimationFunction tests")
X <- rnorm(100)

message("Testing exected errors in 't'")
testthat::test_that("momOUParameterEstimation gives error if t0 != 0",{
      invalid_t = seq(from = 1, to = 5, length = 100)
      expect_error(momOUParameterEstimation(X = X, t = invalid_t), "t0 == 0 is not TRUE")
})
testthat::test_that("momOUParameterEstimation gives error if all t < 1",{
      invalid_t = seq(from = 0, to = 0.5, length = 100)
      expect_error(momOUParameterEstimation(X = X, t = invalid_t), "2 %in% tInt is not TRUE")
})
testthat::test_that("momOUParameterEstimation gives error if t1 != 1",{
      invalid_t = seq(from = 0, to = 5, length = 100)
      expect_error(momOUParameterEstimation(X = X, t = invalid_t), "t1 == 1 is not TRUE")
})

message("Testing output for specific simulation
        ")
v <- var(X)
r <- cor(X[2:100], X[1:99])
while (r < 0) {
        X <- rnorm(100)
        v <- var(X)
        r <- cor(X[2:100], X[1:99])
}
gHat <- - log(r)
sHat2 <- -2 * v * log(r)
gsHat <- suppressMessages(momOUParameterEstimation(X = X, t = 0:99))
testthat::expect_equivalent(
    object = c("first" = gsHat$params[1], "second" = gsHat$params[2]),
    expected = c("third" = gHat, "fourth" = sHat2))

###################
#############
#######
# test for varying sigma2, constant gamma
source("../eulerApproximation_ornsteinU_functions.R")
lenT <- 1e5 + 1
thinBy <- 100
t <- seq(from = 0, to = 100, length = lenT)
seqThin <- seq(from = 1, to = lenT, by = thinBy)
tThin <- t[seqThin]

g <- 1
s2 <- seq(from = 1, to = 101, by = 10)
lenS <- length(s2)
iter <- 100

message(sprintf("Looking at MoM for gamma = %d and sigma2 from %.1f to %.1f by %.1f:",
        g, min(s2), max(s2), s2[2] - s2[1]))
message(sprintf("   performing %d iterations for each combination...", iter))

gs <- cbind(rep(g, lenS), s2)
gsHat <- matrix(ncol = 2, nrow = lenS)
for (i in 1:lenS) {
        ss2 <- s2[i]
        gsHatSS <- matrix(ncol = 2, nrow = 100)
        for (j in 1:iter) {
                sim <- .eulerSimulationOU(x0 = 0, t = t, g = g, s2 = ss2, mu = 0)
                yThin <- sim$y[seqThin]
                gsHatTmp <- momOUParameterEstimation(X = sim$y, t = t)
                gsHatSS[j, ] <- gsHatTmp$params
        }
        gsHat[i, ] <- apply(gsHatSS, 2, mean)
}
gsHat <- cbind(gs[,1], gsHat[,1], gs[,2], gsHat[,2])
colnames(gsHat) <- c("g", "gHat", "s2", "s2Hat")
print(gsHat)

#######
# test for varying sigma2, constant gamma
#source("../eulerApproximation_ornsteinU_functions.R")
lenT <- 1e5 + 1
thinBy <- 100
t <- seq(from = 0, to = 100, length = lenT)
seqThin <- seq(from = 1, to = lenT, by = thinBy)
tThin <- t[seqThin]

g <- seq(from = 1, to = 10, by = 1)
s2 <- 1
lenG <- length(g)
iter <- 100

message(sprintf("Looking at MoM for sigma2 = %d and gamma from %.1f to %.1f by %.1f:",
        s2, min(g), max(g), g[2] - g[1]))
message(sprintf("   performing %d iterations for each combination...", iter))

gs <- cbind(g, rep(s2, lenG))
colnames(gs) <- c("g", "s2")
gsHat <- matrix(ncol = 2, nrow = lenG)
colnames(gsHat) <- c("gHat", "s2Hat")
for (i in 1:lenG) {
        gg <- g[i]
        gsHatGG <- matrix(ncol = 2, nrow = 100)
        for (j in 1:iter) {
                sim <- .eulerSimulationOU(x0 = 0, t = t, g = gg, s2 = s2, mu = 0)
                #yThin <- sim$y[seqThin]
                gsHatTmp <- momOUParameterEstimation(X = sim$y, t = t)
                gsHatGG[j, ] <- gsHatTmp$params
        }
        gsHat[i, ] <- apply(gsHatGG, 2, mean)
}
gsHat <- cbind(gs[,1], gsHat[,1], gs[,2], gsHat[,2])
colnames(gsHat) <- c("g", "gHat", "s2", "s2Hat")
print(gsHat)

# create error: r < 0 -> log(r) = NA
lenT <- 1e5 + 1
t <- seq(from = 0, to = 100, length = lenT)
g <- 20
s2 <- 1
set.seed(1)
sim <- .eulerSimulationOU(x0 = 0, t = t, g = g, s2 = s2, mu = 0)
gsHatTmp <- momOUParameterEstimation(X = sim$y, t = t)

} ## end if (TEST)
