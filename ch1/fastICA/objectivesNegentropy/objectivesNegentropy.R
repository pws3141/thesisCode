#! /usr/bin/env Rscript

cat("start ...\n")
library(fastICA)
library(jvcoords)
library(sfsmisc) # integrate.xy function

source("./Jhatf0Functions.R")
source("../entropy.R")

set.seed(1)

x <- as.matrix(read.csv("../points/points.csv", header = FALSE))
n <- nrow(x)

##################################
# helper functions

cat("setting up fastICA objective function ...\n")
G1 <- function(x) log(cosh(x))
G1.base <- mean(G1(rnorm(1e7)))

# systematically try a grid of directions to project onto
M <- 1000
theta <- seq(0, pi, length.out = M + 1)[1:M]
W <- rbind(cos(theta), sin(theta))

######################################################################
# plot the m-spacing and fastICA objective functions

cat("obtaining m-spacing and fastICA directions...\n")

# find the best direction for fastICA
fi <- fastICA(x, n.comp = 2, alg.typ = "deflation")
stopifnot(isTRUE(all.equal(cov(fi$S) * (n - 1) / n, diag(2))))
stopifnot(isTRUE(all.equal(x %*% fi$K %*% fi$W, fi$S)))
stopifnot(isTRUE(all.equal(t(fi$W) %*% fi$W, diag(2))))
w <- fi$K %*% fi$W %*% c(1, 0) * sqrt((n - 1) / n)
stopifnot(abs(sum(w^2) - 1) < 1e-6)
xf <- as.vector(x %*% w)
stopifnot(isTRUE(all.equal(xf * sqrt(n / (n - 1)), fi$S[, 1])))
if (w[1] > 0) {
    w <- -w
}

# find the best direction using m-spacing
# take M projections in half circle, choose best one
x_rot <- x %*% W
x_rot.entropy <- mSpacingEntropy(t(x_rot), m = sqrt(n))
v <- W[, which.min(x_rot.entropy)]
stopifnot(abs(sum(v^2) - 1) < 1e-6)
if (v[1] < 0) {
    v <- -v
}

# ploting negentropy with 0 marked on y-axis
cat("plotting objective functions ...\n")

entropyGaussian <- 0.5 * (1 + log(2 * pi))

# function to scale negentropy for plotting
ficaScale <- function(negFica, negEnt) {
        negFica <- negFica - min(negFica)
        minFica <- min(negFica)
		maxFica <- max(negFica)
		maxNegentropy <- max(negEnt)
		scaleConst <- maxNegentropy / maxFica
		res <- scaleConst * negFica
		res
}

#dev.new(width = 5.5, height = 3) # in inches
#acrossthetop <- dev.cur()
## plus other sizes ...
## this figure goes across the top
#dev.set(acrossthetop)
#
par(mar = c(4, 4, 3, 1), mgp = c(2, 0.7, 0), cex.lab = 0.8,
  cex.axis = 0.8, cex.main = 1, font.main = 1, las = 1)
op.1 <- par(no.readonly = TRUE)
par(mai = c(.6, .3, .1, .1), mgp = c(2, 0.7, 0), cex.lab = 0.8,
  cex.axis = 0.8, cex.main = 1, font.main = 1, las = 1)
op.2 <- par(no.readonly = TRUE)


# Jf0 objective function
intLower <- -10
intUpper <- 10
message(sprintf("calculating Jf0 over %d directions from theta = 0 to pi...", M))

# find vals for K
Ga <- 1
Kopt <- optim(par = c(0, 0, 0, 1), 
        function(init) .gToKOptimFunction(init, a = Ga),
        control = list(maxit = 1000))
Kvals <- Kopt$par

A <- kappa <- eta <- a <- numeric(length = M)
cRot <- numeric(length = M)
Jf0 <- numeric(length = M)
for(i in 1:M) {
    sRot <- x %*% W[, i]
    denn <- density(sRot, bw = "sj")
    cRot[i] <- integrate.xy(denn$x, .K(denn$x, vals = Kvals, a = Ga) * denn$y)
    f0Opt <- optim(par = c(0, -0.5, 1), 
                function(init) {
                        .f0Optim(init, c = cRot[i], K.vals = Kvals, G.a = Ga)
                }, control = list(maxit = 1000))
    f0vals <- f0Opt$par
    kappa[i] <- f0vals[1]
    eta[i] <- f0vals[2]
    a[i] <- f0vals[3]
    Ainv <- integrate(function(x) {
              exp(kappa[i] * x + eta[i] * x^2 + a[i] * .K(x, vals = Kvals, a = Ga))
                }, lower = intLower, upper = intUpper)
    A[i] <- 1 / Ainv$value

    Jf0Tmp <- Jf0Function(A = A[i], f0.vals = f0vals, K.vals = Kvals, G.a = Ga,
                  int.lower = intLower, int.upper = intUpper)
    Jf0[i] <- Jf0Tmp
}

# find negentropy of all x_rot projections
m <- floor(sqrt(n))
negEnt <- entropyGaussian - x_rot.entropy
negEnt <- negEnt - min(negEnt) # make min(negEnt) = 0

# scale Jf0 w.r.t. negEnt s.t. range is the same
Jf0Scaled <- ficaScale(negFica = Jf0, negEnt = negEnt)
# smooth a little
Jf0Smoothed <- smooth.spline(theta, Jf0Scaled, spar = 0.4)$y

# fICA objective function (here sorting is irrelavent)
fica1 <- -(colMeans(G1(x_rot)) - G1.base)^2
ficaScaled <- ficaScale(negFica=-fica1, negEnt=negEnt)

pdf("objectivesNegentropy.pdf", 5.5, 3.7, family = "serif", pointsize = 10,
                bg = "transparent")
par(op.2, bty = "n") # reset par values

scaleFactor <- max(negEnt) + 1e-1 # put functions on top of e.o.
plot(theta, negEnt + 2*scaleFactor, type = "l", ylim = c(0, 1.3),
     xaxt = "n",
    yaxt = "n",
     xlab = expression(theta ~~ "(projection angle)"), ylab = NA)
title(ylab = "non-Gaussianity", mgp = c(0.7, 0, 0))
axis(1, at = (0:4)*pi/4, padj = c(-.8, -.4, -.4, -.4, -1.2), cex.axis = .8,
     labels = c(expression(0),
		expression(1/4 ~ pi),
		expression(1/2 ~ pi),
		expression(3/4 ~ pi),
		expression(pi)))
#axis(2, at = 0, padj = 0.6, cex.axis = 0.8, labels = expression(0))
abline(v = atan2(v[2], v[1]) %% pi)

lines(x = theta, y = Jf0Smoothed + scaleFactor, lty = 2)

lines(x = theta, y = ficaScaled, lty = 3)
abline(v = atan2(w[2], w[1]) %% pi, lty = 3)

invisible(dev.off())

pdf("objectivesNegentropy_Jf.pdf", 5, 2.4, family = "serif", pointsize = 10,
                bg = "transparent")
par(op.2) # reset par values
plot(theta, negEnt, type = "l", ylim = c(0, 0.45),
     xaxt = "n",
    yaxt = "n",
     xlab = expression(theta ~~ "(projection angle)"), ylab = NA)
title(ylab = "non-Gaussianity", mgp = c(0.7, 0, 0))
axis(1, at = (0:4)*pi/4, padj = c(-.8, -.4, -.4, -.4, -1.2), cex.axis = .8,
     labels = c(expression(0),
		expression(1/4 ~ pi),
		expression(1/2 ~ pi),
		expression(3/4 ~ pi),
		expression(pi)))
axis(2, at = 0, padj = 0.6, cex.axis = 0.8, labels = expression(0))
abline(v = atan2(v[2], v[1]) %% pi)

invisible(dev.off())

pdf("objectivesNegentropy_Jf0.pdf", 5, 2.4, family = "serif", pointsize = 10,
                bg = "transparent")
par(op.2) # reset par values

plot(theta, Jf0Smoothed, type = "l", ylim = c(0, 0.45),
     xaxt = "n",
    yaxt = "n",
     xlab = expression(theta ~~ "(projection angle)"), ylab = NA)
title(ylab = "non-Gaussianity", mgp = c(0.7, 0, 0))
axis(1, at = (0:4)*pi/4, padj = c(-.8, -.4, -.4, -.4, -1.2), cex.axis = .8,
     labels = c(expression(0),
		expression(1/4 ~ pi),
		expression(1/2 ~ pi),
		expression(3/4 ~ pi),
		expression(pi)))
axis(2, at = 0, padj = 0.6, cex.axis = 0.8, labels = expression(0))

maxJf0 <- which.max(Jf0Smoothed)
abline(v = theta[maxJf0], lty = 3)

invisible(dev.off())

pdf("objectivesNegentropy_hatJf0.pdf", 5, 2.4, family = "serif", pointsize = 10,
                bg = "transparent")
par(op.2) # reset par values

plot(theta, ficaScaled, type = "l", ylim = c(0, 0.45),
     xaxt = "n",
    yaxt = "n",
     xlab = expression(theta ~~ "(projection angle)"), ylab = NA)
title(ylab = "non-Gaussianity", mgp = c(0.7, 0, 0))
axis(1, at = (0:4)*pi/4, padj = c(-.8, -.4, -.4, -.4, -1.2), cex.axis = .8,
     labels = c(expression(0),
		expression(1/4 ~ pi),
		expression(1/2 ~ pi),
		expression(3/4 ~ pi),
		expression(pi)))
axis(2, at = 0, padj = 0.6, cex.axis = 0.8, labels = expression(0))
abline(v = atan2(w[2], w[1]) %% pi, lty = 3)

invisible(dev.off())

FOURTH <- exists("FOURTH") && isTRUE(FOURTH)
if (FOURTH) {
        # fICA fourth moments objective function
        # see Mieetinen (2015)
        fm.ica <- abs(colMeans(x_rot^4) - 3)
        fm.ica.scaled <- ficaScale(negFica = fm.ica, negEnt = negEnt)
        pdf("objectivesNegentropy_fourth.pdf", 5, 2.4, family = "serif", pointsize = 10,
                        bg = "transparent")

        plot(theta, fm.ica.scaled, type = "l", ylim = c(0, 0.45),
             xaxt = "n",
            yaxt = "n",
             xlab = expression(theta ~~ "(projection angle)"), ylab = NA)
        title(ylab = "non-Gaussianity", mgp = c(0.7, 0, 0))
        axis(1, at = (0:4)*pi/4, padj = c(-.8, -.4, -.4, -.4, -1.2), cex.axis = .8,
             labels = c(expression(0),
                expression(1/4 ~ pi),
                expression(1/2 ~ pi),
                expression(3/4 ~ pi),
                expression(pi)))
        axis(2, at = 0, padj = 0.6, cex.axis = 0.8, labels = expression(0))
        abline(v = atan2(w[2], w[1]) %% pi, lty = 3)

        invisible(dev.off())
}

invisible(dev.off())

cat("done\n")
