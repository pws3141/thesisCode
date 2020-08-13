.G <- function(x, a = 1) {
        stopifnot(a >= 1 & a <= 2) # following Hyv
        (1 / a) * log(cosh(a * x))
}

.gToKOptimFunction <- function(initial = c(1,1,1,1), a = 1,
                                int.lower = -10, int.upper = 10) {
    alpha <- initial[1]
    beta <- initial[2]
    gamma <- initial[3]
    delta <- initial[4]
    p <- integrate(function(x) (1 / sqrt(2 * pi)) * exp(-x^2 / 2) *
        (.G(x, a=a) + alpha * x^2 + beta * x + gamma) / delta,
        lower = int.lower, upper = int.upper) 
    q <- integrate(function(x) x * (1 / sqrt(2 * pi)) * exp(-x^2 / 2) *
        (.G(x, a=a) + alpha * x^2 + beta * x + gamma) / delta,
        lower = int.lower, upper = int.upper)
    r <- integrate(function(x) x^2 * (1 / sqrt(2 * pi)) * exp(-x^2 / 2) *
        (.G(x, a=a) + alpha * x^2 + beta * x + gamma) / delta,
        lower = int.lower, upper = int.upper)
    s <- integrate(function(x) (1 / sqrt(2 * pi)) * exp(-x^2 / 2) * 
        ((.G(x, a=a) + alpha * x^2 + beta * x + gamma) / delta)^2,
        lower = int.lower, upper = int.upper) 
    return(p$value^2 + q$value^2 + r$value^2 + (1 - s$value)^2)
}

.K <- function(x, vals, a = 1) {
        stopifnot(length(vals) == 4)
        alpha <- vals[1]
        beta <- vals[2]
        gamma <- vals[3]
        delta <- vals[4]
        Kx <- (.G(x, a=a) + alpha * x^2 + beta * x + gamma) / delta
        Kx
}

.f0Optim <- function(initial = c(0, -0.5, 1), c, K.vals, G.a = 1,
                        int.lower = -10, int.upper = 10) {
    stopifnot(!missing(c))
    kappa <- initial[1]
    eta <- initial[2]
    a <- initial[3]
    Ainv <- integrate(function(x) {
                      exp(kappa * x + eta * x^2 + a * .K(x, vals = K.vals, a = G.a))
                }, lower = int.lower, upper = int.upper)
    A <- 1 / Ainv$value
    p <- integrate(function(x) {
                   x * A * exp(kappa * x + eta * x^2 + a * .K(x, vals = K.vals, a = G.a))
                }, lower = int.lower, upper = int.upper)
    q <- integrate(function(x) {
                   x^2 * A * exp(kappa * x + eta * x^2 + a * .K(x, vals = K.vals, a = G.a))
                }, lower = int.lower, upper = int.upper)
    r <- integrate(function(x) {
                   .K(x, vals = K.vals) * A *
                           exp(kappa * x + eta * x^2 + a * .K(x, vals = K.vals, a = G.a))
                }, lower = int.lower, upper = int.upper)
    return(p$value^2 + (1 - q$value)^2 + (r$value - c)^2)
}

.f0 <- function(x, A, vals, K.vals, G.a = 1) {
        stopifnot(A > 0)
        stopifnot(length(vals) == 3)
        kappa <- vals[1]
        eta <- vals[2]
        a <- vals[3]
        y <- A * exp(kappa * x + eta * x^2 +
                a * .K(x, vals = K.vals, a = G.a))
        y
}

.logf0 <- function(x, A, vals, K.vals, G.a = 1) {
        stopifnot(A > 0)
        stopifnot(length(vals) == 3)
        kappa <- vals[1]
        eta <- vals[2]
        a <- vals[3]
        log_y <- log(A) + (kappa * x + eta * x^2 +
                        a * .K(x, vals = K.vals, a = G.a))
        log_y
}

Hf0Function <- function(A, f0.vals, K.vals, G.a = 1,
                int.lower = -10, int.upper = 10) {
        res <- integrate(function(x) {
                         .f0(x = x, A = A, vals = f0.vals, K.vals = K.vals, G.a = G.a) * 
                            .logf0(x = x, A = A, vals = f0.vals, K.vals = K.vals, G.a = G.a)
                 }, lower = int.lower, upper = int.upper)
        res <- -res$value
}

Jf0Function <- function(A, f0.vals, K.vals, G.a = 1,
                int.lower = -10, int.upper = 10) {
        Hnu <- (1 / 2) * (log(2 * pi) + 1)
        H <- Hf0Function(A = A, f0.vals = f0.vals, K.vals = K.vals, G.a = 1,
                        int.lower = int.lower, int.upper = int.upper)
        res <- Hnu - H
}
