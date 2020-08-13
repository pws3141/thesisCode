# basic pointwise observations
pointwiseMeans <- function(X, t, breaks) {
        # mean for each partition of all
        # X_t values where t \in partition
        #stopifnot(t[1] == 0)
        n <- length(t)
        numBreaks <- length(breaks) - 1
        maxBreak <- breaks[numBreaks + 1]
        # split 't' up into relevant 'breaks'
        tModBreaks <- t %% maxBreak
        # find which lower bound of 'breaks' each element of 't' belongs to
        tWhichBreak <- findInterval(tModBreaks, breaks)
        means <- numeric(length = numBreaks)
        for (i in 1:numBreaks) {
                # split t to just elements that belong to break i
                tSplitSeq <- which(tWhichBreak == i)
                # choose just these X values
                Xtmp <- X[tSplitSeq]
                means[i] <- mean(Xtmp)
        }
        means
}

# input 1: partition of time interval; means of intervals
# out 1: values of a, b, c, for a t^2 + b^t + c for each spline

# input 1: partition of time interval; values of a, b, c; time points
# out 2: the values of the quadratic function at each time points

splineOptimFunction <- function(par, means, breaks) {
        # par = list(a, b, c)
        # partition is the right hand x value
        len <- length(breaks)
        a <- par$a
        b <- par$b
        c <- par$c
        forwardSeq <- c(2:len, 1)
        aF <- a[forwardSeq]
        bF <- b[forwardSeq]
        cF <- c[forwardSeq]
        x1 <- breaks
        x0 <- breaks[c(len, 1:(len - 1))]
        a <- (3 / 2) * 
                ((bF - 2 * aF) * (x1^2 - x0^2) + (means + aF - cF) * (x1 - x0)) /
                ((x1^3 - x0^3) - 3 * (x1^2 - x0^2) + 3 * (x1 - x0))
        b <- bF - 2 * (a - aF)
        c <- a - aF + cF
        res <- list(a = a, b = b, c = c)
}

splineOptimCriteria <- function(par, means, breaks) {
        # par = list(a, b, c)
        # partition is the right hand x value
        len <- length(breaks)
        a <- par$a
        b <- par$b
        c <- par$c
        forwardSeq <- c(2:len, 1)
        aF <- a[forwardSeq]
        bF <- b[forwardSeq]
        cF <- c[forwardSeq]
        x1 <- breaks
        x0 <- breaks[c(len, 1:(len - 1))]
        minOne <- (1 / (x1 - x0)) * ((1 / 3) * a * (x1^3 - x0^3) +
                                     (1 / 2) * b * (x1^2 - x0^2) +
                                     c * (x1 - x0)) - means
        minTwo <- (a - aF) * x1^2 + (b - bF) * x1 + (c - cF)
        minThree <- 2 * (a - aF) * x1 + (b - bF)
        min <- minOne^2 + minTwo^2 + minThree^2
        sum(min)
}

systemMatrixA <- function(breaks, i) {
        n <- length(breaks) - 1
        pBackTmp <- breaks[i]
        pTmp <- breaks[i + 1]
        deltaP <- pTmp - pBackTmp
        if (i < n) {
                rowOne <- c(pTmp^2, pTmp, 1, -pTmp^2, -pTmp, -1)
                rowTwo <- c(2 * pTmp, 1, 0, -2 * pTmp, -1, 0)
                rowOne <- c(rep(0, length = 3 * (i - 1)), rowOne, 
                            rep(0, length = 3 * (n - (i + 1))))
                rowTwo <- c(rep(0, length = 3 * (i - 1)), rowTwo, 
                            rep(0, length = 3 * (n - (i + 1))))
        }
        if (i == n) {
                pOne <- breaks[1]
                rowOne1 <- c(-pOne^2, -pOne, -1)
                rowOne2 <- c(pTmp^2, pTmp, 1)
                rowTwo1 <- c(-2 * pOne, -1, 0)
                rowTwo2 <- c(2 * pTmp, 1, 0)
                rowOne <- c(rowOne1, rep(0, length = 3 * (n - 2)),
                            rowOne2)
                rowTwo <- c(rowTwo1, rep(0, length = 3 * (n - 2)),
                            rowTwo2)
        }
        rowThree <- (1 / (pTmp - pBackTmp)) * 
                        c((1 / 3) * (pTmp^3 - pBackTmp^3),
                          (1 / 2) * (pTmp^2 - pBackTmp^2), 1)
        rowThree <- c(rep(0, length = 3 * (i - 1)), rowThree, 
                    rep(0, length = 3 * (n - i)))
        matTmp <- cbind(rowOne, rowTwo, rowThree)
        matTmp
}

# out 1: values of a, b, c, for a t^2 + b^t + c for each spline
quadraticSplines <- function(means, breaks, method = c("system", "optim")) {
        # partition is the right hand x value
        # approximate n quadratic functions 
        # using mean of each function as constraint
        if (all(method == c("system", "optim"))) method = "system"
        n <- length(means)
        if (missing(breaks)) breaks <- (0:n) 
        if (method == "optim") {
                init <- numeric(3 * n)
                opt <- optim(par = init, function(abc) {
                                     abc <- list(a = abc[1:n], b = abc[(n+1):(2*n)],
                                                 c = abc[(2*n + 1):(3*n)])
                                     abcNew <- splineOptimFunction(par = abc, means = means, 
                                                                 breaks = breaks[-1])
                                     min <- splineOptimCriteria(par = abcNew, means = means,
                                                               breaks = breaks[-1])
                                     min
                                             }, method = "BFGS", control = list(maxit = 5000))
                optPar <- opt$par
                optPar <- list(a = optPar[1:n], b = optPar[(n+1):(2*n)],
                         c = optPar[(2*n + 1):(3*n)])
                res <- splineOptimFunction(par = optPar, means = means, 
                                   breaks = breaks[-1])
                res <- sapply(seq_len(n), function(i) {
                                      tmpSeq <- seq(from = i, by = n, length = 3)
                                      res[tmpSeq]
                                             })
                res <- do.call(cbind, res)
        }
        if (method == "system") {
                A <- lapply(seq_len(n), function(i) {
                                    systemMatrixA(breaks = breaks, 
                                                  i = i)
                                             })
                A <- t(do.call(cbind, A))
                rownames(A) <- NULL
                b <- unlist(lapply(seq_len(n), function(i) {
                                           mm <- means[i]
                                           c(rep(0, 2), mm)
                                             }))
                #b <- c(b, 0, 0, 0)
                res <- solve(a = A, b = b)
                res <- matrix(res, ncol = 3, byrow=TRUE)
                colnames(res) <- c("a", "b", "c")
                #res <- res[-(n+1),]
                opt <- list()
        }
        rownames(res) <- paste0("spline", seq_len(n))
        opt$par <- res
        opt$breaks <- breaks
        opt
}

# splinePlot: plot round the circle
splinePlot <- function(spline, r = 5, circle = TRUE, savePlot = TRUE) {
        # r = diameter of circle in plot
        degree <- c(2, 1, 0)
        par <- spline$par
        breaks <- spline$breaks
        # assume equal partition spacing
        deltaPartition <- breaks[2] - breaks[1]
        n <- length(breaks) - 1
        lenSeq <- 100
        xy <- lapply(seq_len(n), function(i) {
                             pTmp <- breaks[i+1]
                             pBackTmp <- breaks[i]
                        tmpSeq <- seq(from = pBackTmp, to = pTmp, length = lenSeq)
                        tmpSeq <- tmpSeq[-lenSeq]
                        xQuadratic <- outer(tmpSeq, degree, "^")
                        splineTmp <- par[i,]
                        yTmp <- as.vector(tcrossprod(splineTmp, xQuadratic))
                        #yTmp
                        list(x = tmpSeq, y = yTmp)
                                             })
        x <- unlist(lapply(seq_len(n), function(i) xy[[i]]$x))
        y <- unlist(lapply(seq_len(n), function(i) xy[[i]]$y))
        if (circle == TRUE) {
        theta <- 2 * pi * x / max(breaks)
        xCircle <- cos(theta) * r
        yCircle <- sin(theta) * r
        plotScale <- 1/15
        R <- r + plotScale * y
        pname <- "splineCirclePlot.pdf"
        if (savePlot == TRUE) pdf(pname)
                par(mai = c(0, 0, 0, 0))
                plot(xCircle, yCircle, t = "l", asp = 1, #col = "light gray",
                                xlim = c(-10, 10), ylim = c(-10, 10),
                                col = rgb(0, 0, 0, .5),
                                axes = FALSE, xlab = "", ylab = "")
                lines(cos(theta) * R, sin(theta) * R, lwd = 2)
                fanLinesTheta <- 2 * pi * breaks / max(breaks)
                xFanLines <- cos(fanLinesTheta) * r
                yFanLines <- sin(fanLinesTheta) * r
                lenFan <- 2
                for (i in 1:n) {
                        lines(x = c(0, lenFan * xFanLines[i]),
                              y = c(0, lenFan * yFanLines[i]), 
                                col = rgb(0, 0, 0, .5))
        }
        }
        if (circle == FALSE) {
                pname <- "splinePlot.pdf"
                if (savePlot == TRUE) pdf(pname)
                #par(mai = c(1, 1, 1, 1))
                plot(x, y, t = "l", lwd = 2,
                                col = rgb(0, 0, 0, 1),
                                axes = TRUE, xlab = "", ylab = "")
                abline(v = breaks, col = rgb(0, 0, 0, .5))
        }
        if (savePlot == TRUE) invisible(dev.off())
}

# input 1: partition of time interval; values of a, b, c; time points
# out 2: the values of the quadratic function at each time points
.muT <- function(t, spline, breaks, deriv=TRUE) {
        n <- length(t)
        numBreaks <- nrow(spline)
        maxBreak <- breaks[numBreaks + 1]
        deltaT <- t[2:n] - t[1:(n - 1)]
        maxT <- max(t)
        # split 't' up into relevant 'breaks'
        tModBreaks <- t %% maxBreak
        # find which lower bound of 'breaks' each element of 't' belongs to
        tWhichBreak <- findInterval(tModBreaks, breaks)
        degree <- c(2, 1, 0)
        derivativeDegree <- c(1, 0)
        muT <- muTDer <- numeric(length = n)
        for (i in 1:numBreaks) {
                # split t to just elements that belong to break i
                tSplitSeq <- which(tWhichBreak == i)
                tSplit <- t[tSplitSeq]
                tModSplit <- tModBreaks[tSplitSeq]
                # choose spline associated with break i
                splineTmp <- spline[i, ]
                # a t^2 + b t + c
                # first, find t^2
                xQuadratic <- outer(tModSplit, degree, "^")
                # then, do a t^2 + b t + c
                muTmp <- as.vector(tcrossprod(splineTmp, xQuadratic))
                # put result into muT vector in correct place
                muT[tSplitSeq] <- muTmp
                if (deriv == TRUE) {
                        xDerivative <- outer(tModSplit, derivativeDegree, "^")
                        xDerivative[, 1] <- 2 * xDerivative[, 1]
                        muDerivativeTmp <- as.vector(tcrossprod(splineTmp[-3], 
                                                          xDerivative))
                        muTDer[tSplitSeq] <- muDerivativeTmp
                }
        }
        if (deriv == TRUE) {
                res <- list(t = t, muT = muT, muTDer = muTDer)
        } else {
                res <- list(t = t, muT = muT)
        }
        res
}

