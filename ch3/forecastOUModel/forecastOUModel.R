#! /usr/bin/env Rscript

#### top stuff

PROD <- !exists("PROD")

if (PROD) {

        ## stuff that happens the first time
        ## etc etc
        suppressPackageStartupMessages(library(ggplot2))
        suppressPackageStartupMessages(library(reshape2))
        suppressPackageStartupMessages(library(latex2exp))

        ## set the seed
        set.seed(1001)
        message("Using specified random seed (1001)")

} else {

        ## stuff that happens NOT the first time
        ## -- I don't usually have anything in here

        ## set the seed
        set.seed(NULL)
        message("Using random random seed")
}

my.Random.seed <- .Random.seed


# simulate Ornstein-Uhlenbeck  process
# d X = - \gamma dt + \sigma dW
# with \gamma = 1 and \sigma^2 = 100
# by Euler-Marayama approximation

# use 3/4 of trajectory to find parameter estimates
# find prediction interval for remaining 1/4 of trajectory
# by simulating N trajectories from x_n using OU model with estimated
# parameters

# load functions
dname <- "../ch2/"
source(paste0(dname, "eulerApproximation_ornsteinU_functions.R"))
source(paste0(dname, "parameterFitting/mle_ornsteinU_functions.R"))

# initialise
tBy = 1e-3
tMax <- 1000
t <- seq(from = 0, to = tMax, by = tBy)
lenT <- length(t)

tThinBy <- 1e-2
tThinBySeq <- tThinBy / tBy
tThinSeq <- seq(from = 1, to = lenT, by = tThinBySeq)
tThin <- t[tThinSeq]

gamma <- 1
sigma2 <- 100

x.0 <- rnorm(1, mean = 0, sd = sqrt(sigma2 / (2 * gamma)))
sim <- eulerSimulationOU(x0 = x.0, t = t, g = gamma, s2 = sigma2, mu = 0)
x <- sim$y
xThin <- x[tThinSeq]

tTrainingTrue <- t <= 750
tTraining <- t[tTrainingTrue]
tTest <- t[!tTrainingTrue]
xTraining <- x[tTrainingTrue]
xTest <- x[!tTrainingTrue]

lenTraining <- length(xTraining)
lenTest <- length(xTest)

tThinTrainingTrue <- tThin <= 750
tThinTraining <- tThin[tThinTrainingTrue]
tThinTest <- tThin[!tThinTrainingTrue]
xThinTraining <- xThin[tThinTrainingTrue]
xThinTest <- xThin[!tThinTrainingTrue]

lenThinTraining <- length(xThinTraining)
lenThinTest <- length(xThinTest)


# estimate gamma, sigma2
x.mle <- mleOUParameterEstimation(X = xThinTraining, t = tThinTraining, gamma0 = gamma)
gHat <- x.mle$par[1]
s2Hat <- x.mle$par[2]

message(sprintf("gH = %.4f; s2Hat = %.4f.",
                gHat, s2Hat))

xn <- xTraining[lenTraining] # same as last obs val of thinned training obs

# horizon
h <- seq(from = 0, to = tMax - max(tTraining), by = tThinBy)
# calculate point forecasts
x.nh <- xn * exp(- gHat * h)
# calculate variance of OU over horizon
var.h <- (s2Hat / (2 * gHat)) * (1 - exp(-2 * gHat * h))
# calculate 95% prediction interval
pi.matrix <- matrix(nrow = 2, ncol = length(h))
a <- 0.05
z <- qnorm(a / 2, lower.tail = FALSE)
pi.matrix <- sapply(seq_len(length(h)), function(i) {
                        x.nh_i <- x.nh[i]
                        var.h_i <- var.h[i]
                        res <- c(x.nh_i - z * sqrt(var.h_i), x.nh_i + z * sqrt(var.h_i))
                        res
                })


# look at which values of test set lie outside PI
whichAbove <- which(xThinTest > pi.matrix[2, ])
whichBelow <- which(xThinTest < pi.matrix[1, ])
propOutside <- (length(whichAbove) + length(whichBelow)) / length(xThinTest)
message(sprintf("%.1f%% of xTest values outside of %.1f%% P.I.",
                propOutside * 100, (1 - a) * 100))

# plotting
# sparse results for plotting
message("plotting...")

h.max <- 50
# sparse everything to make plotting nicer
sparsePlot.trainingSeq <- seq(from = 1, to = lenThinTraining, by = 100)
sparsePlot.testSeq <- seq(from = 1, to = lenThinTest, by = 100)

tTraining.plot <- tThinTraining[sparsePlot.trainingSeq]
xTraining.plot <- xThinTraining[sparsePlot.trainingSeq]

tTest.plot <- tThinTest[sparsePlot.testSeq]
xTest.plot <- xThinTest[sparsePlot.testSeq]

x.nh.plot <- x.nh[sparsePlot.testSeq]
pi.matrix.plot <- pi.matrix[, sparsePlot.testSeq]

par(mai = c(.6, .3, .1, .1), mgp = c(2, 0.7, 0), cex.lab = 0.8,
  cex.axis = 0.8, cex.main = 1, font.main = 1, las = 1)
op <- par(no.readonly = TRUE)

pdf("forecastOUModel.pdf", 6, 3, family = "serif", pointsize = 10,
                bg = "transparent")
par(op, bty = "n") # reset par values

plot(x = tTraining.plot, y = xTraining.plot, 
     ylim = c(-max(abs(xThin)), max(abs(xThin))), 
     xlim = c(700, (max(tThinTraining) + h.max)), 
     t = 'l', xlab = "t", ylab = "y")
lines(x = tTest.plot, y = x.nh.plot, col = "grey", lwd = 1)
lines(x = tTest.plot, y = xTest.plot, lty = 3)
lines(x = tTest.plot, y = pi.matrix.plot[1, ], lty = 2)
lines(x = tTest.plot, y = pi.matrix.plot[2, ], lty = 2)

invisible(dev.off())


