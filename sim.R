library(Hmisc)
library(rms)
library(iteval)

source("support_functions.R")

load("population.RData")

results <- list()

app.model.coefs <- list()

list.of.ite.discr.stats <- c("cforbenefit", "cforbenefit.y0hat",
                         "mbcforbenefit",
                         "true.diff.y0y1.mbcben",
                         "true.diff.y0y1.pop")

list.of.ite.cal.stats <- c("oe", "empirical", "true.oe", "true.cal",
                       "pop.oe", "pop.empirical", "pop.cal")

list.of.risk.cal.stats <- c("ctrl.app", "trt.app", "overall.app",
                            "ctrl.val1", "trt.val1", "overall.val1",
                            "ctrl.val2", "trt.val2", "overall.val2")

list.of.true.rmspe.stats <- c("rmspe.p01.app", "rmspe.p01.val1", "rmspe.p01.val2",
                              "rmspe.delta.app", "rmspe.delta.val1", "rmspe.delta.val2")

sample.sizes <- c(500, 750, 1000)

g=1
for(g in 1:3){

  ptm <- proc.time()

  apparent.discr <- vector("list", length = length(list.of.ite.discr.stats))
  names(apparent.discr) <- list.of.ite.discr.stats
  ext1.app.discr <- ext1.total.discr <- ext2.app.discr <- ext2.total.discr <-
    boot0.632.discr <- boot.opt.discr <- apparent.discr

  apparent.cal <- vector("list", length = length(list.of.ite.cal.stats))
  names(apparent.cal) <- list.of.ite.cal.stats
  ext1.app.cal <- ext1.total.cal <- ext2.app.cal <- ext2.total.cal <-
    boot0.632.cal <- boot.opt.cal <- apparent.cal

  risk.cal <- vector("list", length = length(list.of.risk.cal.stats))
  names(risk.cal) <- list.of.risk.cal.stats

  rmspe <- vector("list", length = length(list.of.true.rmspe.stats))
  names(rmspe) <- list.of.true.rmspe.stats

  n.dev <- sample.sizes[g]
  n.val <- 1000

  Sample <- c(sample(which(val1.out$trt==0), n.dev/2, replace=FALSE),
              sample(which(val1.out$trt==1), n.dev/2, replace=FALSE))

  ctrlX <- val1.ctrlX[Sample, ]
  trtX <- val1.trtX[Sample, ]
  X <- val1.X[Sample, ]
  trt <- val1.out$trt[Sample]
  ind.A <- which(trt == 0)
  ind.B <- which(trt == 1)

  p0 <- plogis(ctrlX %*% beta)
  p1 <- plogis(trtX %*% beta)
  p01 <- ifelse(trt==0, p0, p1)
  delta <- p1 - p0

  y0 <- rbinom(length(p0),1,p0)
  y1 <- rbinom(length(p1),1,p1)
  ben <- y1-y0 < 0

  y <- ifelse(trt==0, y0, y1)
  mod <- glm(y ~ -1 + X, family = "binomial")
  app.model.coefs <- coef(mod)

  y0hat <- plogis(ctrlX %*% coef(mod))
  y1hat <- plogis(trtX %*% coef(mod))
  deltahat <- y1hat - y0hat
  deltahatlp <- qlogis(y1hat) - qlogis(y0hat)
  yhat <- ifelse(trt==0, y0hat, y1hat)

  #######################
  ## Apparent          ##
  #######################

  ## Discrimination ##
  apparent.discr$cforbenefit <- cbendelta(deltahat, y, ind.A, ind.B, nresample=1000)
  apparent.discr$cforbenefit.y0hat <- cbeny0(y0hat, y1hat, y, ind.A, ind.B, match="y0hat", nresample=1000)
  apparent.discr$mbcforbenefit <- mbcb(y0hat, y1hat)
  apparent.discr$true.diff.y0y1.mbcben <- mbcb(y0hat, y1hat, p0, p1)

  # population
  y0hat.ref <- plogis(val1.ctrlX %*% coef(mod))
  y1hat.ref <- plogis(val1.trtX %*% coef(mod))
  deltahat.ref <- y1hat.ref-y0hat.ref
  apparent.discr$true.diff.y0y1.pop <- rcorr.cens(deltahat.ref, I(val1.out$y1-val1.out$y0))["C Index"]

  ## Calibration ##
  apparent.cal$oe <- oe(y0hat, y1hat, y, ind.A, ind.B)
  apparent.cal$empirical <- coef(cal(y0hat,y1hat,y,ind.B))
  apparent.cal$true.oe <- oe(y0hat, y1hat, p01, ind.A, ind.B) # so for Inf replications of Y for this sample of X
  y1hat.updated <- plogis(qlogis(p0) + deltahatlp)
  apparent.cal$true.cal <- coef(cal(p0,y1hat.updated,p1,ind.B)) # pure deltahat assessment
  # that is, if you were to succeed in getting y0hat right
  # (the only sensible interpretations seems the be when y0hat is accurate)

  # population
  y0hat.ref <- plogis(val1.ctrlX %*% coef(mod))
  y1hat.ref <- plogis(val1.trtX %*% coef(mod))
  deltahat.ref <- y1hat.ref-y0hat.ref
  deltahatlp.ref <- (val1.trtX - val1.ctrlX) %*% coef(mod)
  apparent.cal$pop.oe <- oe(y0hat.ref, y1hat.ref, val1.out$y, val1.ind.A, val1.ind.B)
  apparent.cal$pop.empirical <- coef(cal(y0hat.ref,y1hat.ref,val1.out$y,val1.ind.B))
  y1.pop.updated <- plogis(qlogis(val1.out$p0) + deltahatlp.ref)
  apparent.cal$pop.cal  <- coef(cal(val1.out$p0,y1.pop.updated,val1.out$p1))

  # risk.cal app
  val.prob.cols <- c("C (ROC)", "R2", "Brier", "Intercept", "Slope", "Emax", "E90", "Eavg")
  risk.cal$ctrl.app <- val.prob(y0hat[ind.A], y[ind.A], pl=FALSE)[val.prob.cols]
  risk.cal$trt.app <- val.prob(y1hat[ind.B], y[ind.B], pl=FALSE)[val.prob.cols]
  risk.cal$overall.app <- val.prob(yhat, y, pl=FALSE)[val.prob.cols]

  # rmspe app
  rmspe$rmspe.p01.app <- sqrt(mean((p01 - yhat)^2))
  rmspe$rmspe.delta.app <- sqrt(mean((delta - deltahat)^2))

  #######################
  ## bootstrap correct ##
  #######################
  B <- 250

  ## dicsr
  # for 0.632
  cben.boot <- cbeny0.boot <- ben.boot <- list()
  # for opt
  cben.boot.opt.discr <- cbeny0.boot.opt.discr <- ben.boot.opt.discr <- list()
  cben.orig.opt.discr <- cbeny0.orig.opt.discr <- ben.orig.opt.discr <- list()

  ## cal
  # for 0.632
  oe.boot.0.632 <- rep(NA, B)
  cal.mod.boot <- matrix(NA, B, 2)
  # for opt
  oe.boot <- oe.boot.orig <- rep(NA, B)
  coef.cal.mod.boot <- coef.cal.mod.orig <- matrix(NA, B, 2)

  b=1
  for(b in 1:B){
    # sample
    ind <- sample(1:n.dev, n.dev, replace=TRUE)

    # get boot matrices
    X.in <- X[ind, ]
    y.in <- y[ind]
    ctrlX.in <- trtX.in <- X.in
    ctrlX.in[ ,c("trt", "trt:x1", "trt:x2")] <- 0
    trtX.in[ ,c("trt", "trt:x1", "trt:x2")] <- cbind(1, X.in[ ,"x1"], X.in[ ,"x2"])
    ind.A.boot <- which(X.in[ ,"trt"] == 0)
    ind.B.boot <- which(X.in[ ,"trt"] == 1)

    mod.boot <- glm(y.in ~ -1 + X.in, family = "binomial")

    # get oos matrices
    X.out <- X[-ind, ]
    y.out <- y[-ind]
    ctrlX.out <- trtX.out <- X.out
    ctrlX.out[ ,c("trt", "trt:x1", "trt:x2")] <- 0
    trtX.out[ ,c("trt", "trt:x1", "trt:x2")] <- cbind(1, X.out[ ,"x1"], X.out[ ,"x2"])
    ind.A.out <- which(X.out[ ,"trt"] == 0)
    ind.B.out <- which(X.out[ ,"trt"] == 1)

    # 0.632+ discr
    # predict out-of-sample
    y0hat.boot.oos <- plogis(ctrlX.out %*% coef(mod.boot))
    y1hat.boot.oos <- plogis(trtX.out %*% coef(mod.boot))
    deltahat.boot.oos <- y1hat.boot.oos-y0hat.boot.oos

    cben.boot[[b]] <- cbendelta(deltahat.boot.oos, y.out,
                           ind.A.out, ind.B.out, nresample=1000)
    cbeny0.boot[[b]] <- cbeny0(y0hat.boot.oos, y1hat.boot.oos, y.out,
                                   ind.A.out, ind.B.out, match="y0hat", nresample=1000)
    #!#! ben.boot[[b]] <-  mbcb(y0hat.boot.oos, y1hat.boot.oos, y0hat[-ind], y1hat[-ind])
    mod.oos <- glm(y.out ~ -1 + X.out, family = "binomial")
    y0hat.true.oos <- plogis(ctrlX.out %*% coef(mod.oos))
    y1hat.true.oos <- plogis(trtX.out %*% coef(mod.oos))
    ben.boot[[b]] <-  mbcb(y0hat.boot.oos, y1hat.boot.oos, y0hat.true.oos, y1hat.true.oos)

    # opt corrected discr
    y0hat.boot <- plogis(ctrlX.in %*% coef(mod.boot))
    y1hat.boot <- plogis(trtX.in %*% coef(mod.boot))
    deltahat.boot <- y1hat.boot - y0hat.boot

    y0hat.orig <- plogis(ctrlX %*% coef(mod.boot))
    y1hat.orig <- plogis(trtX %*% coef(mod.boot))
    deltahat.boot.orig <- y1hat.orig - y0hat.orig

    cben.boot.opt.discr[[b]] <- cbendelta(deltahat.boot, y.in,
                                     ind.A.boot, ind.B.boot, nresample=1000)
    cben.orig.opt.discr[[b]] <- cbendelta(deltahat.boot.orig, y,
                                     ind.A, ind.B, nresample=1000)
    cbeny0.boot.opt.discr[[b]] <- cbeny0(y0hat.boot, y1hat.boot, y.in,
                                             ind.A.boot, ind.B.boot, match="y0hat", nresample=1000)
    cbeny0.orig.opt.discr[[b]] <- cbeny0(y0hat.orig, y1hat.orig, y,
                                             ind.A, ind.B, match="y0hat", nresample=1000)
    ben.boot.opt.discr[[b]] <- mbcb(y0hat.boot, y1hat.boot)
    ben.orig.opt.discr[[b]] <- mbcb(y0hat.orig, y1hat.orig, y0hat, y1hat)

    # 0.632 cal
    # predict out-of-sample
    y0hat.boot.oos <- plogis(ctrlX.out %*% coef(mod.boot))
    y1hat.boot.oos <- plogis(trtX.out %*% coef(mod.boot))
    deltahat.boot.oos <- y1hat.boot.oos-y0hat.boot.oos
    oe.boot.0.632[b] <- oe(y0hat.boot.oos, y1hat.boot.oos, y.out, ind.A.out, ind.B.out)
    # cal.mod.boot[b, ] <- coef(cal(y0hat.boot.oos, y1hat.boot.oos, y.out, ind.B.out))

    deltahat.boot.oos.lp <- qlogis(y1hat.boot.oos) - qlogis(y0hat.boot.oos)
    mod.oos.ctrl <- glm(y.out[ind.A.out] ~ 1 + X.out[ind.A.out, c("x1", "x2")], family = "binomial")
    y0hat.boot.oos <- plogis(cbind(1, ctrlX.out[ , c("x1", "x2")]) %*% coef(mod.oos.ctrl))
    y1hat.boot.oos <- plogis(qlogis(y0hat.boot.oos) +  deltahat.boot.oos.lp)
    cal.mod.boot[b, ] <- coef(cal(y0hat.boot.oos, y1hat.boot.oos, y.out, ind.B.out))

    # optimism corrected cal
    y0hat.boot <- plogis(ctrlX.in %*% coef(mod.boot))
    y1hat.boot <- plogis(trtX.in %*% coef(mod.boot))
    deltahat.boot <- y1hat.boot-y0hat.boot
    oe.boot[b] <- oe(y0hat.boot,y1hat.boot, y.in, ind.A.boot,ind.B.boot)
    coef.cal.mod.boot[b, ] <- coef(cal(y0hat.boot, y1hat.boot, y.in, ind.B.boot))

    y0hat.boot.orig <- plogis(ctrlX %*% coef(mod.boot))
    y1hat.boot.orig <- plogis(trtX %*% coef(mod.boot))
    deltahat.boot.orig <- y1hat.boot.orig-y0hat.boot.orig
    oe.boot.orig[b] <- oe(y0hat.boot.orig,y1hat.boot.orig,y,ind.A,ind.B)
    # coef.cal.mod.orig[b, ] <- coef(cal(y0hat.boot.orig, y1hat.boot.orig, y, ind.B))

    deltahat.boot.orig.lp <- qlogis(y1hat.boot.orig) - qlogis(y0hat.boot.orig)
    y1hat.boot.orig <- plogis(qlogis(y0hat) + deltahat.boot.orig.lp)
    coef.cal.mod.orig[b, ] <- coef(cal(y0hat, y1hat.boot.orig, y, ind.B))
  }

  # 0.632+ discr
  app <- apparent.discr$cforbenefit
  boot.oos <- mean(unlist(cben.boot))
  boot0.632.discr$cforbenefit <- boot0.632plus_cstat(app, boot.oos)

  app <- apparent.discr$cforbenefit.y0hat
  boot.oos <- mean(unlist(cbeny0.boot))
  boot0.632.discr$cforbenefit.y0hat <- boot0.632plus_cstat(app, boot.oos)

  app <- apparent.discr$mbcforbenefit
  boot.oos <- mean(unlist(ben.boot))
  boot0.632.discr$mbcforbenefit <- boot0.632plus_cstat(app, boot.oos)

  boot0.632.discr$true.diff.y0y1.mbcben <- apparent.discr$true.diff.y0y1.mbcben
  boot0.632.discr$true.diff.y0y1.pop <- apparent.discr$true.diff.y0y1.pop

  # opt discr
  boot.opt.discr$cforbenefit <- apparent.discr$cforbenefit - mean(unlist(cben.boot.opt.discr) - unlist(cben.orig.opt.discr))
  boot.opt.discr$cforbenefit.y0hat <- apparent.discr$cforbenefit.y0hat - mean(unlist(cbeny0.boot.opt.discr) - unlist(cbeny0.orig.opt.discr))
  boot.opt.discr$mbcforbenefit <- apparent.discr$mbcforbenefit - mean(unlist(ben.boot.opt.discr) - unlist(ben.orig.opt.discr))
  boot.opt.discr$true.diff.y0y1.mbcben <- apparent.discr$true.diff.y0y1.mbcben
  boot.opt.discr$true.diff.y0y1.pop <- apparent.discr$true.diff.y0y1.pop

  # 0.632+ cal
  app <- apparent.cal$empirical[2]
  boot.oos <- mean(cal.mod.boot[ ,2])
  boot0.632.cal$empirical <- c(NA, boot0.632plus_cal_slope(app, boot.oos))

  # opt
  boot.opt.cal$oe <- apparent.cal$oe + mean(oe.boot.orig - oe.boot)
  boot.opt.cal$empirical <- apparent.cal$empirical + apply(coef.cal.mod.orig - coef.cal.mod.boot, 2, function(x) mean(x))

  # boot cal references
  boot0.632.cal$true.oe <- boot.opt.cal$true.oe <- apparent.cal$true.oe
  boot0.632.cal$true.cal <- boot.opt.cal$true.cal  <- apparent.cal$true.cal
  boot0.632.cal$pop.oe <- boot.opt.cal$pop.oe  <- apparent.cal$pop.oe
  boot0.632.cal$pop.empirical <- boot.opt.cal$pop.empirical  <- apparent.cal$pop.empirical
  boot0.632.cal$pop.cal <- boot.opt.cal$pop.cal  <- apparent.cal$pop.cal

  #######################
  ## new data same DGM ##
  #######################
  Sample <- c(sample(which(val1.out$trt==0), n.val/2, replace=FALSE),
              sample(which(val1.out$trt==1), n.val/2, replace=FALSE))

  ctrlX <- val1.ctrlX[Sample, ]
  trtX <- val1.trtX[Sample, ]
  X <- val1.X[Sample, ]
  trt <- val1.out$trt[Sample]
  ind.A <- which(trt == 0)
  ind.B <- which(trt == 1)

  p0 <- plogis(ctrlX %*% beta)
  p1 <- plogis(trtX %*% beta)
  p01 <- ifelse(trt==0, p0, p1)
  delta <- p1 - p0

  y0 <- rbinom(length(p0),1,p0)
  y1 <- rbinom(length(p1),1,p1)
  ben <- y1-y0 < 0

  y <- ifelse(trt==0, y0, y1)

  y0hat <- plogis(ctrlX %*% coef(mod))
  y1hat <- plogis(trtX %*% coef(mod))
  deltahat <- y1hat - y0hat
  deltahatlp <- qlogis(y1hat) - qlogis(y0hat)
  yhat <- ifelse(trt==0, y0hat, y1hat)

  # apparent discr
  ext1.app.discr$cforbenefit <- cbendelta(deltahat, y, ind.A, ind.B, nresample=1000)
  ext1.app.discr$cforbenefit.y0hat <- cbeny0(y0hat, y1hat, y, ind.A, ind.B, match="y0hat", nresample=1000)
  ext1.app.discr$mbcforbenefit <- mbcb(y0hat, y1hat)
  ext1.app.discr$true.diff.y0y1.mbcben <- mbcb(y0hat, y1hat, p0, p1)
  ext1.app.discr$true.diff.y0y1.pop <- apparent.discr$true.diff.y0y1.pop

  # total.update discr
  mod.updated <- glm(formula = y ~ -1 + X, family = "binomial")
  y0hat.updated <- plogis(ctrlX %*% coef(mod.updated))
  y1hat.updated <- plogis(trtX %*% coef(mod.updated))

  ext1.total.discr$cforbenefit <- NA
  ext1.total.discr$cforbenefit.y0hat <- cbeny0(y0hat.updated, I(y0hat.updated + deltahat), y, ind.A, ind.B, match="y0hat", nresample=1000)
  ext1.total.discr$mbcforbenefit <- mbcb(y0hat, y1hat, y0hat.updated, y1hat.updated)
  ext1.total.discr$true.diff.y0y1.mbcben <- ext1.app.discr$true.diff.y0y1.mbcben
  ext1.total.discr$true.diff.y0y1.pop <- ext1.app.discr$true.diff.y0y1.pop

  # apparent cal
  ext1.app.cal$oe <- oe(y0hat,y1hat,y,ind.A,ind.B)
  ext1.app.cal$empirical <- coef(cal(y0hat,y1hat,y,ind.B))
  ext1.app.cal$true.oe <- oe(y0hat,y1hat,p01,ind.A,ind.B)
  y1hat.updated <- plogis(qlogis(p0) + deltahatlp)
  ext1.app.cal$true.cal <- coef(cal(p0,y1hat.updated,p1,ind.B))

  # population (just a copy from apparent)
  ext1.app.cal$pop.oe <- apparent.cal$pop.oe
  ext1.app.cal$pop.empirical <- apparent.cal$pop.empirical
  ext1.app.cal$pop.cal <- apparent.cal$pop.cal

  # total.update
  mod.updated <- glm(formula = y[ind.A] ~ 1 + X[ind.A, c("x1", "x2") ], family = "binomial")
  y0hat.updated <- plogis(cbind(1, ctrlX[ ,c("x1", "x2")]) %*% coef(mod.updated))
  y1hat.updated <- plogis(qlogis(y0hat.updated) + deltahatlp)
  cal.mod.update <- cal(y0hat.updated,y1hat.updated,y,ind.B)

  ext1.total.cal$oe <- oe(y0hat.updated,y1hat.updated,y,ind.A,ind.B)
  ext1.total.cal$empirical <- coef(cal.mod.update)
  ext1.total.cal$true.oe <- oe(y0hat.updated, y1hat.updated, p01, ind.A, ind.B)
  y1hat.updated <- plogis(qlogis(p0) + deltahatlp)
  ext1.total.cal$true.cal <- coef(cal(p0,y1hat.updated,p1,ind.B))

  # population (total.update)
  y0hat.ref <- plogis(cbind(1, val1.ctrlX[ ,c("x1", "x2")]) %*% coef(mod.updated))
  y1hat.ref <- plogis(qlogis(y0hat.ref) + deltahatlp.ref)
  deltahat.ref <- y1hat.ref-y0hat.ref
  deltahatlp.ref <- qlogis(y1hat.ref) - qlogis(y0hat.ref)
  ext1.total.cal$pop.oe <- oe(y0hat.ref, y1hat.ref, val1.out$y, val1.ind.A, val1.ind.B)
  ext1.total.cal$pop.empirical <- coef(cal(y0hat.ref,y1hat.ref,val1.out$y,val1.ind.B))
  y1.pop.updated <- plogis(qlogis(val1.out$p0) + deltahatlp.ref)
  ext1.total.cal$pop.cal <- coef(cal(val1.out$p0,y1.pop.updated,val1.out$p1))

  # risk.cal val DGM 1
  risk.cal$ctrl.val1 <- val.prob(y0hat[ind.A], y[ind.A], pl=FALSE)[val.prob.cols]
  risk.cal$trt.val1 <- val.prob(y1hat[ind.B], y[ind.B], pl=FALSE)[val.prob.cols]
  risk.cal$overall.val1 <- val.prob(yhat, y, pl=FALSE)[val.prob.cols]

  # rmspe
  rmspe$rmspe.p01.val1   <- sqrt(mean((p01 - yhat)^2))
  rmspe$rmspe.delta.val1 <- sqrt(mean((delta - deltahat)^2))

  ############################
  ## new data different DGM ##
  ############################
  Sample <- c(sample(which(val2.out$trt==0), n.val/2, replace=FALSE),
              sample(which(val2.out$trt==1), n.val/2, replace=FALSE))

  ctrlX <- val2.ctrlX[Sample, ]
  trtX <- val2.trtX[Sample, ]
  X <- val2.X[Sample, ]
  trt <- val2.out$trt[Sample]
  ind.A <- which(trt == 0)
  ind.B <- which(trt == 1)

  p0 <- plogis(ctrlX %*% beta.alt)
  p1 <- plogis(trtX %*% beta.alt)
  p01 <- ifelse(trt==0, p0, p1)
  delta <- p1 - p0

  y0 <- rbinom(length(p0),1,p0)
  y1 <- rbinom(length(p1),1,p1)
  ben <- y1-y0 < 0

  y <- ifelse(trt==0, y0, y1)

  y0hat <- plogis(ctrlX %*% coef(mod))
  y1hat <- plogis(trtX %*% coef(mod))
  deltahat <- y1hat - y0hat
  deltahatlp <- qlogis(y1hat) - qlogis(y0hat)
  yhat <- ifelse(trt==0, y0hat, y1hat)

  # apparent.discr
  ext2.app.discr$cforbenefit <- cbendelta(deltahat, y, ind.A, ind.B, nresample=1000)
  ext2.app.discr$cforbenefit.y0hat <- cbeny0(y0hat, y1hat, y, ind.A, ind.B, match="y0hat", nresample=1000)
  ext2.app.discr$mbcforbenefit <- mbcb(y0hat, y1hat)
  ext2.app.discr$true.diff.y0y1.mbcben <- mbcb(y0hat, y1hat, p0, p1)

  # population
  y0hat.ref <- plogis(val2.ctrlX %*% coef(mod))
  y1hat.ref <- plogis(val2.trtX %*% coef(mod))
  deltahat.ref <- y1hat.ref-y0hat.ref
  ext2.app.discr$true.diff.y0y1.pop <- rcorr.cens(deltahat.ref, I(val2.out$y1-val2.out$y0))["C Index"]

  # total.update
  mod.updated <- glm(formula = y ~ -1 + X, family = "binomial")
  y0hat.updated <- plogis(ctrlX %*% coef(mod.updated))
  y1hat.updated <- plogis(trtX %*% coef(mod.updated))

  ext2.total.discr$cforbenefit <- NA
  ext2.total.discr$cforbenefit.y0hat <- cbeny0(y0hat.updated, I(y0hat.updated + deltahat), y, ind.A, ind.B, match="y0hat", nresample=1000)
  ext2.total.discr$mbcforbenefit <- mbcb(y0hat, y1hat, y0hat.updated, y1hat.updated)
  ext2.total.discr$true.diff.y0y1.mbcben <- ext2.app.discr$true.diff.y0y1.mbcben
  ext2.total.discr$true.diff.y0y1.pop <- ext2.app.discr$true.diff.y0y1.pop

  # apparent calibration
  ext2.app.cal$oe <- oe(y0hat,y1hat,y,ind.A,ind.B)
  ext2.app.cal$empirical <- coef(cal(y0hat,y1hat,y,ind.B))
  ext2.app.cal$true.oe <- oe(y0hat, y1hat, p01, ind.A, ind.B)
  y1hat.updated <- plogis(qlogis(p0) + deltahatlp)
  ext2.app.cal$true.cal <- coef(cal(p0,y1hat.updated,p1,ind.B))

  # population
  y0hat.ref <- plogis(val2.ctrlX %*% coef(mod))
  y1hat.ref <- plogis(val2.trtX %*% coef(mod))
  deltahat.ref <- y1hat.ref-y0hat.ref
  deltahatlp.ref <- (val2.trtX - val2.ctrlX) %*% coef(mod)

  ext2.app.cal$pop.oe <- oe(y0hat.ref, y1hat.ref, val2.out$y, val2.ind.A, val2.ind.B)
  ext2.app.cal$pop.empirical <- coef(cal(y0hat.ref,y1hat.ref,val2.out$y,val2.ind.B))
  y1.pop.updated <- plogis(qlogis(val2.out$p0) + deltahatlp.ref)
  ext2.app.cal$pop.cal  <- coef(cal(val2.out$p0,y1.pop.updated,val2.out$p1))

  # total.update calibration
  mod.updated <- glm(formula = y[ind.A] ~ 1 + X[ind.A, c("x1", "x2") ], family = "binomial")
  y0hat.updated <- plogis(cbind(1, ctrlX[ ,c("x1", "x2")]) %*% coef(mod.updated))
  y1hat.updated <- plogis(qlogis(y0hat.updated) + deltahatlp)
  cal.mod.update <- cal(y0hat.updated,y1hat.updated,y,ind.B)

  ext2.total.cal$oe <- oe(y0hat.updated,y1hat.updated,y,ind.A,ind.B)
  ext2.total.cal$empirical <- coef(cal.mod.update)
  ext2.total.cal$true.oe <- oe(y0hat.updated, y1hat.updated, p01, ind.A, ind.B)
  y1hat.updated <- plogis(qlogis(p0) + deltahatlp)
  ext2.total.cal$true.cal <- coef(cal(p0,y1hat.updated,p1,ind.B))

  # population (total.update)
  y0hat.ref <- plogis(cbind(1, val2.ctrlX[ ,c("x1", "x2")]) %*% coef(mod.updated))
  y1hat.ref <- plogis(qlogis(y0hat.ref) + deltahatlp.ref)
  deltahat.ref <- y1hat.ref-y0hat.ref
  deltahatlp.ref <- qlogis(y1hat.ref) - qlogis(y0hat.ref)
  ext2.total.cal$pop.oe <- oe(y0hat.ref, y1hat.ref, val2.out$y, val2.ind.A, val2.ind.B)
  ext2.total.cal$pop.empirical <- coef(cal(y0hat.ref,y1hat.ref,val2.out$y,val2.ind.B))
  y1.pop.updated <- plogis(qlogis(val2.out$p0) + deltahatlp.ref)
  ext2.total.cal$pop.cal <- coef(cal(val2.out$p0,y1.pop.updated,val2.out$p1))

  # risk.cal val DGM 2
  risk.cal$ctrl.val2 <- val.prob(y0hat[ind.A], y[ind.A], pl=FALSE)[val.prob.cols]
  risk.cal$trt.val2 <- val.prob(y1hat[ind.B], y[ind.B], pl=FALSE)[val.prob.cols]
  risk.cal$overall.val2 <- val.prob(yhat, y, pl=FALSE)[val.prob.cols]

  # rmspe
  rmspe$rmspe.p01.val2   <- sqrt(mean((p01 - yhat)^2))
  rmspe$rmspe.delta.val2 <- sqrt(mean((delta - deltahat)^2))

  results[[g]] <- list(apparent.discr, ext1.app.discr, ext1.total.discr, ext2.app.discr, ext2.total.discr,
                       boot0.632.discr, boot.opt.discr,
                       apparent.cal, ext1.app.cal, ext1.total.cal, ext2.app.cal, ext2.total.cal,
                      boot0.632.cal, boot.opt.cal,
                      risk.cal, rmspe,
                      "time.elapsed.sec"=(proc.time() - ptm)["elapsed"])
  names(results[[g]]) <- c("apparent.discr", "ext1.app.discr", "ext1.total.discr", "ext2.app.discr", "ext2.total.discr",
                           "boot0.632.discr", "boot.opt.discr",
                           "apparent.cal", "ext1.app.cal", "ext1.total.cal", "ext2.app.cal", "ext2.total.cal",
                           "boot0.632.cal", "boot.opt.cal",
                           "risk.cal", "rmspe",
                           "time.elapsed.sec")
}

names(results)[1:3] <- paste("n", sample.sizes, sep="")

# # Save the results for each iteration (iter) with suffix ".RData"
# save(results, file=paste0("ite", iter, suffix))

