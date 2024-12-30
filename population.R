set.seed(123)

beta <- c(-1,-.75,1,0,0,.5)
beta.alt <- c(-0.5, -0.5, 0.75, 0.25, 0.25, 0.25)

## validation sample DGM1
npop <- 100000
n <- npop
p <- 2
trt <- rep(c(0,1), each=n/2)
val1.ind.A <- which(trt == 0)
val1.ind.B <- which(trt == 1)
X <- cbind(trt, matrix(rnorm(p*n), n, p))
colnames(X) <- c("trt", "x1", "x2")
X <- model.matrix( ~ trt * (x1+x2), data=data.frame(X))
ctrlX <- trtX <- X
ctrlX[ ,c("trt", "trt:x1", "trt:x2")] <- 0
trtX[ ,c("trt", "trt:x1", "trt:x2")] <- cbind(1, X[ ,"x1"], X[ ,"x2"])

p0 <- plogis(ctrlX %*% beta)
p1 <- plogis(trtX %*% beta)

y0 <- rbinom(n,1,p0)
y1 <- rbinom(n,1,p1)

y <- ifelse(trt==0, y0, y1)

val1.out <- data.frame(y, y0, y1, p0, p1, trt)
val1.X <- X
val1.ctrlX <- ctrlX
val1.trtX <- trtX

## validation sample DGM2
p <- 2
trt <- rep(c(0,1), each=n/2)
val2.ind.A <- which(trt == 0)
val2.ind.B <- which(trt == 1)
X <- cbind(trt, matrix(rnorm(p*n), n, p))
colnames(X) <- c("trt", "x1", "x2")
X <- model.matrix( ~ trt * (x1+x2), data=data.frame(X))
ctrlX <- trtX <- X
ctrlX[ ,c("trt", "trt:x1", "trt:x2")] <- 0
trtX[ ,c("trt", "trt:x1", "trt:x2")] <- cbind(1, X[ ,"x1"], X[ ,"x2"])

p0 <- plogis(ctrlX %*% beta.alt)
p1 <- plogis(trtX %*% beta.alt)

y0 <- rbinom(n,1,p0)
y1 <- rbinom(n,1,p1)

y <- ifelse(trt==0, y0, y1)

val2.out <- data.frame(y, y0, y1, p0, p1, trt)
val2.X <- X
val2.ctrlX <- ctrlX
val2.trtX <- trtX

# marginal event rate controls and treated
# DGM1
events <- list()
events[[1]] <- mean(val1.out$y[val1.out$trt == 0])
events[[2]] <- mean(val1.out$y[val1.out$trt == 1])

# DGM2
events[[3]] <- mean(val2.out$y[val2.out$trt == 0])
events[[4]] <- mean(val2.out$y[val2.out$trt == 1])
round(unlist(events), 3)

save(val1.out, val1.X, val1.ctrlX, val1.trtX, val1.ind.A, val1.ind.B,
     val2.out, val2.X, val2.ctrlX, val2.trtX, val2.ind.A, val2.ind.B,
     beta, beta.alt, file="population.RData")

