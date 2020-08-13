#! /usr/bin/env Rscript

cat("start ...\n")
library(fastICA)
library(jvcoords)
library(sfsmisc) # integrate.xy function

source("../entropy.R")

set.seed(1)

x <- as.matrix(read.csv("../points/points.csv", header = FALSE))

# find f0 density for the two directions
# add plot that show density of f and f_0 on
#   m-spacing and fastICA projections

cat("Computing densities of m-spacing and fastICA projections...\n")
# Hyvin denn
# find vals for KX
gx.to.kx <- function(initial) {
	# function to be optimised to find parameters of K given G
	if(missing(initial)) {
		initial <- c(1, 1, 1, 1)
	}
	alpha <- initial[1]
	beta <- initial[2]
	gamma <- initial[3]
	delta <- initial[4]
	p <- integrate(function(x) (1 / sqrt(2 * pi)) * exp(-x^2 / 2) *
		(log(cosh(x)) + alpha * x^2 + beta * x + gamma) / delta,
		lower = -10, upper = 10)
	q <- integrate(function(x) x * (1 / sqrt(2 * pi)) * exp(-x^2 / 2) *
		(log(cosh(x)) + alpha * x^2 + beta * x + gamma) / delta,
		lower = -10, upper = 10)
	r <- integrate(function(x) x^2 * (1 / sqrt(2 * pi)) * exp(-x^2 / 2) *
		(log(cosh(x)) + alpha * x^2 + beta * x + gamma) / delta,
		lower = -10, upper = 10)
	s <- integrate(function(x) (1 / sqrt(2 * pi)) * exp(-x^2 / 2) *
		((log(cosh(x)) + alpha * x^2 + beta * x + gamma) / delta)^2,
		lower = -10, upper = 10)
	return(p$value^2 + q$value^2 + r$value^2 + (1 - s$value)^2)
}


Kx.fn <- function(x, opt.vals = Kx.opt.vals) {
	# function outputs K given K params
	# with G(x) = log(cosh(x))
		alpha <- opt.vals[1]
		beta <- opt.vals[2]
		gamma <- opt.vals[3]
		delta <- opt.vals[4]
		Gx <- log(cosh(x))
		Kx <- (Gx + alpha * x^2 + beta * x + gamma) / delta
		return(Kx)
}

f0hat.fn <- function(x, c, Kx.opt.vals) {
	# returns \hat f_0(x) fastICA density approximation
	# given K parameters
	f0hatx <- dnorm(x) * (1 + c * Kx.fn(x, opt.vals = Kx.opt.vals))
	return(f0hatx)
}

xLim <- 3
optKx <- optim(par = c(0, 0, 0, 1),
		function(init) gx.to.kx(init),
		control = list(maxit = 1000))

Kx.opt.vals <- optKx$par

x_plot <- seq(-xLim, xLim, length = 1000)

# f0 density in m-spacing direction
# find the best direction using m-spacing
# take M projections in half circle, choose best one
M <- 5000
theta <- seq(0, pi, length.out = M + 1)[1:M]
W <- rbind(cos(theta), sin(theta))
n <- nrow(x)
x_rot <- x %*% W
x_rot.entropy <- mSpacingEntropy(t(x_rot), m = sqrt(n))
v <- W[, which.min(x_rot.entropy)]

xc <- as.vector(x %*% v)
dc <- density(xc, bw = "SJ")
c_dc <- integrate.xy(dc$x, Kx.fn(dc$x) * dc$y)
f0hat_xc <- f0hat.fn(x_plot, c = c_dc, Kx.opt.vals = Kx.opt.vals)

# f0 density in fastICA directions
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
    xf <- -xf
}
df <- density(xf, bw = "SJ")

c_df <- integrate.xy(df$x, Kx.fn(df$x) * df$y)
f0hat_xf <- f0hat.fn(x_plot, c = c_df, Kx.opt.vals = Kx.opt.vals)

# add denns of f and f0hat together
#dc <- density(resClust$S[,1], bw="SJ")
#df <- density(resFast$S[,1], bw="SJ")

cat("Plotting densities...\n")

par(mai = c(.4, .7, .1, .1), mgp = c(2, 0.7, 0), cex.lab = 0.8,
  cex.axis = 0.8, cex.main = 1, font.main = 1, las = 1)
op <- par(no.readonly = TRUE)

maxPlot <- max(dc$y, f0hat_xc, 
               df$y, f0hat_xf)
pdf("densityClusterICADirection.pdf", width = 2.75, height = 2.3,
    family = "serif", pointsize = 10)
par(op)
plot(dc, main = "",
     ylim = c(0, maxPlot), xlim = c(-xLim, xLim), lty = 1, xlab = "")
lines(x_plot, f0hat_xc, lty = 3)
invisible(dev.off())

pdf("densityFastICADirection.pdf", width=2.75, height = 2.3,
    family = "serif", pointsize = 10)
par(op)
plot(df, main = "",
     ylim = c(0, maxPlot), xlim = c(-xLim, xLim), lty = 1, xlab = "")
lines(x_plot, f0hat_xf, lty = 3)
invisible(dev.off())

cat("done\n")
