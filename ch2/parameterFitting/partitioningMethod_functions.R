# use box method to estimate drift
# in 1D, so box is partition of real line
# box = n, then [min(X), max(X)] split into n equal intervals
.sdeBoxDriftEstimation1D <- function(X, t, box) {
        n <- length(X)
        deltaT <- t[2] - t[1] # assuming t equal partition
        maxX <- max(X)
        minX <- min(X)
        if (length(box) == 1) {
                partition <- seq(from = minX - 1e-2, to = maxX + 1e-2, length = box)
        } else {
                partition <- box
                if (minX < min(box)) {
                        partition <- c(minX - 10e-3, partition)
                }
                if (maxX > max(box)) {
                        partition <- c(partition, maxX + 10e-3)
                }
                box <- length(partition)
        }
        deltaPartition <- partition[-1] - partition[-box]
        bHat <- sapply(seq_len(box-1), function(i) {
                               p <- partition[i]
                               deltaP <- deltaPartition[i]
                               whichXBox <- which(X[1:(n-1)] >= p & 
                                                  X[1:(n-1)] < (p + deltaP))
                               if (length(whichXBox) > 0) {
                                       XTmp <- X[whichXBox]
                                       XMoveTmp <- X[whichXBox + 1]
                                       bTmp <- sum(XMoveTmp - XTmp) / 
                                               (deltaT * length(whichXBox))
                                       bTmp
                               } else {
                                       NA
                               }
                        })
        bHat
}

.sdeBoxDiffusionEstimation1D <- function(X, t, box) {
        n <- length(X)
        deltaT <- t[2] - t[1] # assuming t equal partition
        maxX <- max(X)
        minX <- min(X)
        if (length(box) == 1) {
                partition <- seq(from = minX - 1e-2, to = maxX + 1e-2, length = box)
        } else {
                partition <- box
                if (minX < min(box)) {
                        partition <- c(minX - 10e-3, partition)
                         }
                if (maxX > max(box)) {
                        partition <- c(partition, maxX + 10e-3)
                         }
                box <- length(partition)
        }
        deltaPartition <- partition[-1] - partition[-box]
        aHat <- sapply(seq_len(box-1), function(i) {
                               p <- partition[i]
                               deltaP <- deltaPartition[i]
                               whichXBox <- which(X[1:(n-1)] >= p & 
                                                  X[1:(n-1)] < (p + deltaP))
                               if (length(whichXBox) > 0) {
                                       XTmp <- X[whichXBox]
                                       XMoveTmp <- X[whichXBox + 1]
                                       aTmp <- sum((XTmp - XMoveTmp)^2) / 
                                               (deltaT * length(whichXBox))
                                       aTmp
                               } else {
                                       NA
                               }
                        })
        aHat
}

# use box method to estimate drift in 2D, 
# box = n, then [min(X_1), max(X_1)] \times [min(X_2), max(X_2)] split into n equal intervals
.sdeBoxDriftEstimation2D <- function(X, t, box) {
        n <- nrow(X)
        deltaT <- t[2] - t[1] # assuming t equal partition
        maxX1 <- max(X[,1])
        minX1 <- min(X[,1])
        maxX2 <- max(X[,2])
        minX2 <- min(X[,2])
        partition1 <- seq(from = minX1 - 1e-3, to = maxX1 + 1e-3, length = box)
        partition2 <- seq(from = minX2 - 1e-3, to = maxX2 + 1e-3, length = box)
        deltaPartition1 <- partition1[2] - partition1[1]
        deltaPartition2 <- partition2[2] - partition2[1]
        partition <- cbind(rep(partition1[1:(box-1)], each = (box - 1)),
                           rep(partition2[1:(box-1)], times = (box - 1)))
        bHat <- apply(partition, 1, function(p) {
                               whichXBox <- which((X[1:(n-1), 1] >= p[1] & 
                                                  X[1:(n-1), 1] < (p[1] + deltaPartition1))
                                               & (X[1:(n-1), 2] >= p[2] & 
                                                  X[1:(n-1), 2] < (p[2] + deltaPartition2)))
                               len_XBox <- length(whichXBox)
                               if (len_XBox > 0) {
                                       XTmp <- X[whichXBox, ,drop=FALSE]
                                       XMoveTmp <- X[(whichXBox + 1), ,drop=FALSE]
                                       apply(XMoveTmp - XTmp, 2, function(y) sum(y)) /
                                               (deltaT * len_XBox)
                               } else {
                                       warning("No points in a partition...")
                                       c(0, 0)
                               }
                           })
        # see Notes: 19.05.28 - 2D partition matrix convention 
        # for matrix layout wrt partition
        bHat <- t(bHat)
        partitionMatX <- matrix(rep(partition1, each = box), ncol = box, 
                                nrow = box, byrow = TRUE)
        partitionMatY <- matrix(rep(partition2, times = box), ncol = box, 
                                nrow = box, byrow = TRUE)
        bMatX <- matrix(bHat[,1], ncol = (box - 1), nrow = (box - 1), byrow = TRUE)
        bMatY <- matrix(bHat[,2], ncol = (box - 1), nrow = (box - 1), byrow = TRUE)
        list(partX = partitionMatX, partY = partitionMatY,
             bX = bMatX, bY = bMatY)
}

sdeBoxDriftEstimation <- function(X, t, box) {
        p <- ncol(as.matrix(X))
        stopifnot(p == 1 | p == 2)
        if (p == 1) {
                res <- .sdeBoxDriftEstimation1D(X = X, t = t, box = box)
        }
        if (p == 2) {
                res <- .sdeBoxDriftEstimation2D(X = X, t = t, box = box)
        }
        res
}

sdeBoxDiffusionEstimation <- function(X, t, box) {
        p <- ncol(as.matrix(X))
        stopifnot(p == 1 | p == 2)
        if (p == 1) {
                res <- .sdeBoxDiffusionEstimation1D(X = X, t = t, box = box)
        }
        if (p == 2) {
                res <- NA
                warning("Diffusion only aproximated for p = 1")
        }
        res
}

boxParameterEstimation <- function(X, t, box) {
                bHat <- sdeBoxDriftEstimation(X = X, t = t, box = box)
                aHat <- sdeBoxDiffusionEstimation(X = X, t = t, box = box)
                if (is.vector(X)) {
                        n <- length(X)
                        deltaT <- t[2] - t[1] # assuming t equal partition
                        maxX <- max(X)
                        minX <- min(X)
                        if (length(box) == 1) {
                                partition <- seq(from = minX - 1e-2, to = maxX + 1e-2, 
                                                 length = box)
                                deltaPartition <- partition[2] - partition[1]
                                partitionMid <- partition[-box] + 0.5 * deltaPartition
                        } else {
                                partition <- box
                                if (minX < min(box)) {
                                        partition <- c(minX - 10e-3, partition)
                                                 }
                                if (maxX > max(box)) {
                                        partition <- c(partition, maxX + 10e-3)
                                                 }
                                box <- length(partition)
                                deltaPartition <- partition[-1] - partition[-box]
                                partitionMid <- partition[-box] + 0.5 * deltaPartition
                        }
                        res <- cbind(partition[-box], partition[-1], 
                                     partitionMid, bHat, aHat)
                        colnames(res) <- c("lowerP", "upperP", "midP", "bHat", "aHat")
                } #end if vector X
                else {
                        if (ncol(X) == 2) {
                                res <- cbind(bHat, aHat)
                        }
                }
                res
}

### Testing

TEST <- exists("TEST") && isTRUE(TEST)

if (TEST) {
message("
Sourced \'boxParameterEstimationFunctions\', to estimate parameters
of an OU model: dX = -gamma X dt + sigma dB
")
#library(testthat)
source("../eulerApproximation_ornsteinU_functions.R")

############
########
####
##
# test how error behaves for varying box size
g <- 1
s2 <- 1
iter <- 100
boxes <- seq(from = 5, to = 105, by = 10)
lenBox <- length(boxes)
message(sprintf("testing box method for varying number of boxes:
    from %d to %d by %d
    take mean of all iterations", min(boxes), max(boxes), boxes[2] - boxes[1])) 

deltaT <- 1e-3
toT <- 100
t <- seq(from = 0, to = toT, by = deltaT)

message(sprintf("for each number of boxes: 
    simulate %d processes with t from %.1f to %.1f and deltaT = %.4f",
    iter, min(t), max(t), deltaT))

listAllMean <- vector(mode = "list", length = lenBox)
listAllMedian <- vector(mode = "list", length = lenBox)
matOutMean <- matrix(ncol = 2, nrow = lenBox)
matOutMedian <- matrix(ncol = 2, nrow = lenBox)
for (i in 1:lenBox) {
        boxTmp <- boxes[i]
        resBoxMedian <- matrix(ncol = 2, nrow = iter)
        resBoxMean <- matrix(ncol = 2, nrow = iter)
        for (j in 1:iter) {
                sim <- .eulerSimulationOU(x0 = 0, t = t, g = g, s2 = s2, mu = 0)
                res <- boxParameterEstimation(X = sim$y, t = sim$t, box = boxTmp)
                # dX = -g X dt + s dW
                resOU <- t(apply(res, 1, function(x) {
                                       x1 <- - x[2] / x[1]
                                       x2 <- x[3]
                       c(x1, x2)
                }))
                resOUMedian <- apply(resOU, 2, median)
                resBoxMedian[j ,] <- resOUMedian
                resOUMean <- apply(resOU, 2, mean)
                resBoxMean[j ,] <- resOUMean
        }
        listAllMean[[i]] <- resBoxMean
        matOutMean[i ,] <- apply(resBoxMean, 2, mean)
        listAllMedian[[i]] <- resBoxMedian
        matOutMedian[i ,] <- apply(resBoxMedian, 2, mean)
}
matOut <- cbind(boxes, matOutMean, matOutMedian)
matOutMean <- cbind(boxes, matOutMean)
matOutMedian <- cbind(boxesMedian, matOutMedian)
colnames(matOutMean) <- colnames(matOutMedian) <- c("boxes", "gHat", "s2Hat")
colnames(matOut) <- c("boxes", "gHatMean", "s2HatMean",
                      "gHatMedian", "s2HatMedian")
print(matOut)

####
##
# test with user defined boxes
g <- 1
s2 <- 1
iter <- 100
boxes <- seq(from = 5, to = 105, by = 10)
lenBox <- length(boxes)
deltaT <- 1e-3
toT <- 100
t <- seq(from = 0, to = toT, by = deltaT)

sim <- .eulerSimulationOU(x0 = 0, t = t, g = g, s2 = s2, mu = 0)
partition <- c(-2, 1, 2)
res <- boxParameterEstimation(X = sim$y, t = sim$t, box = partition)


############
########
####
##
# test how error behaves as \Delta t \to \infty, t \to \infty
g <- 1
s2 <- 1
lenS <- length(s2)
iter <- 100

N <- 5000 # keep length of chain the same as \Delta t increases
thinBy <- seq(from = 10, to = 1010, by = 100)
len <- length(thinBy)

for (i in 1:len) {
        thinn <- thinBy[i]
        for (j in 1:iter) {
                seqThin <- seq(from = 1, to = N * thinBy, by = thinBy)
                deltaT <- 1e-3
                toT <- thinBy * deltaT * N 
                t <- seq(from = 0, to = toT, by = deltaT)
                tThin <- t[seqThin]
                sim <- .eulerSimulationOU(x0 = 0, t = t, g = g, s2 = s2, mu = 0)
                yThin <- sim$y[seqThin]
                res <- boxParameterEstimation(X = yThin, t = tThin, box = 15)
                # dX = -g X dt + s dW
                resOU <- t(apply(res, 1, function(x) {
                                       x1 <- - x[2] / x[1]
                                       x2 <- x[3]
                                       c(x1, x2)
                                }))
                resOU <- apply(resOU, 2, mean)
        }
}

}
