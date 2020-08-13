#! /usr/bin/env Rscript

cat("start ...\n")
library(fastICA)
library(jvcoords)

source("../entropy.R")

# good values 6, 10, 11, 14, ...
set.seed(14)

# rotate the points anti-clockwise, so that the scatter plot looks nice
rotate.deg = 10

##################################
# helper functions

cat("setting up fastICA objective functions ...\n")
G1 <- function(x) log(cosh(x))
G2 <- function(x) -exp(-x^2/2)
G1.base <- mean(G1(rnorm(1e7)))
G2.base <- -1/sqrt(2)

###################################
# create a stripy, whitened cloud of points

# systematically try a grid of directions to project onto
M <- 100
theta <- seq(0, pi, length.out = M + 1)[1:M]
W <- rbind(cos(theta), sin(theta))

# direction of the gaps is theta[mid1]
mid1 <- 14

gap <- function(x) {
    s <- sin(6*x)
    .4 * s^6 * (s > 0)
}

random.vec <- function(n) {
    m <- floor(sqrt(n))

    x <- matrix(0, 0, 2)
    while (TRUE) {
        x <- rbind(x, matrix(rnorm(2 * (n - nrow(x))), ncol =  2))

        x <- whiten(x)$y
        all.proj <- apply(x %*% W, 2, sort)
        spacing <- all.proj[(m + 1):n, ] - all.proj[1:(n - m), ]
        ent <- colSums(log(spacing*(n/m)))
        dt <- theta[which.min(ent)] - theta[mid1]
        U <- matrix(c(cos(dt), sin(dt), -sin(dt), cos(dt)), 2, 2)
        x <- x %*% U
        if (mean(gap(x %*% W[, mid1]) > 0) > .5) {
            x <- -x
        }

        gg <- gap(x %*% W[, mid1])
        rows <- which(gg > 1e-3)
        if (length(rows) == 0) {
            break
        }
        x <- x[-rows, , drop = FALSE]
    }
    x
}

cat("creating random points ...\n")
n <- 2000
x <- random.vec(n)

# rotate the points for the scatter plot
phi <- rotate.deg / 180 * pi
U <- matrix(c(cos(phi), -sin(phi), sin(phi), cos(phi)), 2, 2)
x <- x %*% U

fname <- "points.csv"
cat("writing", fname, "...\n")
write.table(x, fname, sep = ",", row.names = FALSE, col.names = FALSE)

######################################################################
# create a scatter plot together with the marginals

cat("creating scatter plot ...\n")

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

# find the best direction using m-spacing
# take M projections in half circle, choose best one
M <- 5000
theta <- seq(0, pi, length.out = M + 1)[1:M]
W <- rbind(cos(theta), sin(theta))

n <- nrow(x)
x_rot <- x %*% W
x_rot.entropy <- mSpacingEntropy(t(x_rot), m = sqrt(n))
v <- W[, which.min(x_rot.entropy)]
stopifnot(abs(sum(v^2) - 1) < 1e-6)
xc <- as.vector(x %*% v)
if (v[1] < 0) {
    v <- -v
    xc <- -xc
}
dc <- density(xc, bw = "SJ")

pdf("points.pdf", 4.5, 5.5, family = "serif", pointsize = 10)
par(mai = c(0, 0, 0, 0))

plot(x, xlim = c(-6, 4), ylim = c(-6, 6), asp = 1, pch = 20, cex = .6,
     col = rgb(0, 0, 0, .5),
     axes = FALSE, xlab = "", ylab = "")

d.dist <- 4

Q <- matrix(c(w[1], -w[2], w[2], w[1]), 2, 2)
lines(c(-4 * w[1], 4 * w[1]), c(-4 * w[2], 4 * w[2]), lty = 3)
for (pos in seq(-pi, pi, l = 7)) {
    a <- c(pos, 0) %*% Q
    b <- c(pos, d.dist - .3) %*% Q
    c <- c(pos, d.dist) %*% Q
    lines(rbind(a, b), lty = 3)
    arrows(b[1], b[2], c[1], c[2], length = .08, lwd = .6)
}
z <- cbind(df$x, 4 * df$y + d.dist)
A <- z %*% Q
polygon(A, border = NA, col = "grey")
z <- rbind(c(-4, d.dist), c(4, d.dist))
lines(z %*% Q)
lines(A)

Q <- matrix(c(v[1], -v[2], v[2], v[1]), 2, 2)
lines(c(-4 * v[1], 4 * v[1]), c(-4 * v[2], 4 * v[2]))
for (pos in seq(-pi, pi, l = 7)) {
    a <- c(pos, 0) %*% Q
    b <- c(pos, d.dist) %*% Q
    arrows(a[1], a[2], b[1], b[2], length = .08, lwd = .6)
}
z <- cbind(dc$x, 4 * dc$y + d.dist)
A <- z %*% Q
polygon(A, border = NA, col = "grey")
z <- rbind(c(-4, d.dist), c(4, d.dist))
lines(z %*% Q)
lines(A)

invisible(dev.off())

cat("done\n")

