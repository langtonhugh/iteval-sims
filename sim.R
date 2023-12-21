library(doParallel)
library(pbapply) # adds progress bar
numCores <- parallel::detectCores()
print(numCores)
cl <- makeCluster(numCores-1)

info <- sessionInfo()

nsim=500

todo <- c()
for(ss in 1:nsim){
  if(!file.exists(file=paste0("results/", "ite", ss, ".Rdata"))){
    todo <- c(todo,ss)
  }
}

ss=1
pbres <- pblapply(todo, function(ss){

  set.seed(ss)

  library(MatchIt)
  library(Hmisc)
  library(rms)
  library(iteval)

  source("support_functions.R")

  load("population.RData")

  results <- list()

  app.model.coefs <- list()

  list.of.ite.discr.stats <- c("cforbenefit", "cforbenefit.y0hat",
                               "mbcforbenefit",
                               "cforbenefit.new",
                               "true.diff.y0y1.mbcb")

  list.of.ite.cal.stats <- c("oe", "empirical.trt", "empirical.both",
                             "true.oe", "true.trt", "true.both")

  list.of.risk.cal.stats <- c("ctrl.app", "trt.app", "overall.app",
                              "ctrl.val1", "trt.val1", "overall.val1",
                              "ctrl.val2", "trt.val2", "overall.val2")

  list.of.true.rmspe.stats <- c("rmspe.p01.app", "rmspe.p01.val1", "rmspe.p01.val2",
                                "rmspe.delta.app", "rmspe.delta.val1", "rmspe.delta.val2")

  sample.sizes <- c(500, 750, 1000)

  g=1
  for(g in 1:length(sample.sizes)){

    ptm <- proc.time()

    apparent.discr <- vector("list", length = length(list.of.ite.discr.stats))
    names(apparent.discr) <- list.of.ite.discr.stats
    ext1.app.discr <- ext1.total.discr <- ext2.app.discr <- ext2.total.discr <-
      boot0.632.discr <- boot.opt.discr <- boot.oos.90ci.discr <- boot.oos.95ci.discr <-
      boot0.632.discr.90ci <- boot0.632.discr.95ci <- apparent.discr

    apparent.cal <- vector("list", length = length(list.of.ite.cal.stats))
    names(apparent.cal) <- list.of.ite.cal.stats
    ext1.app.cal <- ext1.total.cal <- ext2.app.cal <- ext2.total.cal <-
      boot0.632.cal <- boot.opt.cal <- boot.oos.90ci.cal <- boot.oos.95ci.cal <-
      boot0.632.cal.90ci <- boot0.632.cal.95ci <- apparent.cal

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
    trt <- X[ ,"trt"]
    ind.A <- which(trt == 0)
    ind.B <- which(trt == 1)

    p0 <- val1.out$p0[Sample] # == plogis(ctrlX %*% beta)
    p1 <- val1.out$p1[Sample] # == plogis(trtX %*% beta)
    p01 <- ifelse(trt==0, p0, p1)
    delta <- p1 - p0

    y0 <- rbinom(length(p0),1,p0) # take a new sample
    y1 <- rbinom(length(p1),1,p1) # take a new sample
    ben <- y1-y0 < 0

    # y <- val1.out$y[Sample]
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
    apparent.discr$cforbenefit.new <- cben.new(form=trt ~ x1 + x2, X=X, y=y, y0hat=y0hat, y1hat=y1hat) # Maas et al. 2023; van Klaveren et al. 2023
    apparent.discr$true.diff.y0y1.mbcb <- mbcb(y0hat, y1hat, p0, p1) # estimand!

    ## Calibration ##
    apparent.cal$oe <- oe(y0hat, y1hat, y, ind.A, ind.B)
    apparent.cal$empirical.trt <- coef(cal(y0hat,y1hat,y,trt=trt,type="treated"))
    apparent.cal$empirical.both <- coef(cal(y0hat,y1hat,y,trt=trt,type="both"))
    apparent.cal$true.oe <- oe(y0hat, y1hat, p01, ind.A, ind.B) # so against true probabilities
    apparent.cal$true.trt <- coef(cal(y0hat,y1hat,y,trt=trt,p0=p0,p1=p1,type="treated"))
    apparent.cal$true.both <- coef(cal(y0hat,y1hat,y,trt=trt,p0=p0,p1=p1,type="both"))

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
    cben.boot <- cbeny0.boot <- ben.boot <- cben.new.boot <- list()
    # for opt
    cben.boot.opt.discr <- cbeny0.boot.opt.discr <- ben.boot.opt.discr <- cben.new.boot.opt.discr <- list()
    cben.orig.opt.discr <- cbeny0.orig.opt.discr <- ben.orig.opt.discr <- cben.new.orig.opt.discr <- list()

    ## cal
    # for 0.632
    oe.boot.0.632 <- rep(NA, B)
    empirical.trt.0.632 <- empirical.both.0.632 <- matrix(NA, B, 2)
    # for opt
    oe.boot <- oe.boot.orig <- rep(NA, B)
    empirical.trt.boot <- empirical.both.boot <- empirical.trt.orig <- empirical.both.orig <- matrix(NA, B, 2)

    b=1
    for(b in 1:B){
      # sample
      ind <- c(sample(ind.A, n.dev/2, replace=TRUE), sample(ind.B, n.dev/2, replace=TRUE))

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
      deltahatlp.boot.oos <- (trtX.out - ctrlX.out) %*% coef(mod.boot) # new (to check)

      cben.boot[[b]] <- cbendelta(deltahat.boot.oos, y.out,
                             ind.A.out, ind.B.out, nresample=1000)
      cbeny0.boot[[b]] <- cbeny0(y0hat.boot.oos, y1hat.boot.oos, y.out,
                                     ind.A.out, ind.B.out, match="y0hat", nresample=1000)
      cben.new.boot[[b]] <- cben.new(form=trt ~ x1 + x2, X=X.out, y=y.out, y0hat=y0hat.boot.oos, y1hat=y1hat.boot.oos)
      mod.oos <- glm(y.out ~ -1 + X.out, family = "binomial")
      y0hat.true.oos <- plogis(ctrlX.out %*% coef(mod.oos))
      y1hat.true.oos <- plogis(trtX.out %*% coef(mod.oos))
      ben.boot[[b]] <-  mbcb(y0hat.boot.oos, y1hat.boot.oos, y0hat.true.oos, y1hat.true.oos)

      # opt corrected discr
      y0hat.boot <- plogis(ctrlX.in %*% coef(mod.boot))
      y1hat.boot <- plogis(trtX.in %*% coef(mod.boot))
      deltahat.boot <- y1hat.boot - y0hat.boot

      y0hat.boot.orig <- plogis(ctrlX %*% coef(mod.boot))
      y1hat.boot.orig <- plogis(trtX %*% coef(mod.boot))
      deltahat.boot.orig <- y1hat.boot.orig - y0hat.boot.orig

      cben.boot.opt.discr[[b]] <- cbendelta(deltahat.boot, y.in,
                                       ind.A.boot, ind.B.boot, nresample=1000)
      cben.orig.opt.discr[[b]] <- cbendelta(deltahat.boot.orig, y,
                                       ind.A, ind.B, nresample=1000)
      cbeny0.boot.opt.discr[[b]] <- cbeny0(y0hat.boot, y1hat.boot, y.in,
                                               ind.A.boot, ind.B.boot, match="y0hat", nresample=1000)
      cbeny0.orig.opt.discr[[b]] <- cbeny0(y0hat.boot.orig, y1hat.boot.orig, y,
                                               ind.A, ind.B, match="y0hat", nresample=1000)
      ben.boot.opt.discr[[b]] <- mbcb(y0hat.boot, y1hat.boot)
      ben.orig.opt.discr[[b]] <- mbcb(y0hat.boot.orig, y1hat.boot.orig, y0hat, y1hat) # original full sample model as reference
      cben.new.boot.opt.discr[[b]] <- cben.new(form=trt ~ x1 + x2, X=X.in, y=y.in, y0hat=y0hat.boot, y1hat=y1hat.boot)
      cben.new.orig.opt.discr[[b]] <- cben.new(form=trt ~ x1 + x2, X=X, y=y, y0hat=y0hat.boot.orig, y1hat=y1hat.boot.orig)


      ## 0.632 calibration
      oe.boot.0.632[b] <- oe(y0hat.boot.oos, y1hat.boot.oos, y.out, ind.A.out, ind.B.out)
      empirical.trt.0.632[b,] <- coef(cal(y0hat.boot.oos,y1hat.boot.oos,y.out,trt=X.out[ ,"trt"],
                                           p0=y0hat.true.oos, p1=y1hat.true.oos, type="treated"))
      empirical.both.0.632[b,] <- coef(cal(y0hat.boot.oos,y1hat.boot.oos,y.out,trt=X.out[ ,"trt"],
                                            p0=y0hat.true.oos, p1=y1hat.true.oos, type="both"))

      # optimism corrected cal
      oe.boot[b] <- oe(y0hat.boot,y1hat.boot, y.in, ind.A.boot,ind.B.boot) # always 0 for maxlik
      empirical.trt.boot[b,] <- coef(cal(y0hat.boot,y1hat.boot,y.in,trt=X.in[ ,"trt"], type="treated")) # always 0,1 for maxlik
      empirical.both.boot[b,] <- coef(cal(y0hat.boot,y1hat.boot,y.in,trt=X.in[ ,"trt"], type="both")) # always 0,1 for maxlik

      oe.boot.orig[b] <- oe(y0hat.boot.orig,y1hat.boot.orig,y,ind.A,ind.B)
      empirical.trt.orig[b,] <- coef(cal(y0hat.boot.orig,y1hat.boot.orig,y,trt=trt, type="treated"))
      empirical.both.orig[b,] <- coef(cal(y0hat.boot.orig,y1hat.boot.orig,y,trt=trt, type="both"))

      print(b)
    }


    ## Discrimination ####

    # CIs based on oos bootstrap evaluations
    boot.oos.90ci.discr$true.diff.y0y1.mbcb <- NULL
    boot.oos.90ci.discr$cforbenefit <- quantile(unlist(cben.boot), probs=c(.05,.5,.95))
    boot.oos.90ci.discr$cforbenefit.y0hat <- quantile(unlist(cbeny0.boot), probs=c(.05,.5,.95))
    boot.oos.90ci.discr$mbcforbenefit <- quantile(unlist(ben.boot), probs=c(.05,.5,.95))
    boot.oos.90ci.discr$cforbenefit.new <- quantile(unlist(cben.new.boot), probs=c(.05,.5,.95))

    boot.oos.95ci.discr$true.diff.y0y1.mbcb <- NULL
    boot.oos.95ci.discr$cforbenefit <- quantile(unlist(cben.boot), probs=c(.025,.5,.975))
    boot.oos.95ci.discr$cforbenefit.y0hat <- quantile(unlist(cbeny0.boot), probs=c(.025,.5,.975))
    boot.oos.95ci.discr$mbcforbenefit <- quantile(unlist(ben.boot), probs=c(.025,.5,.975))
    boot.oos.95ci.discr$cforbenefit.new <- quantile(unlist(cben.new.boot), probs=c(.025,.5,.975))

    # 0.632+ discr
    app <- apparent.discr$cforbenefit
    boot.oos <- mean(unlist(cben.boot))
    boot0.632.discr$cforbenefit <- boot0.632plus_cstat(app, boot.oos)
    # naive CI (with just 1 oos evaluation in the inner loop, corresponding to boot0.632plus_cstat with just 1 random element of cben.boot)
    boot0.632.discr.90ci$cforbenefit <- quantile(sapply(unlist(cben.boot), function(x) boot0.632plus_cstat(app, x)), probs=c(.05,.5,.95))
    boot0.632.discr.95ci$cforbenefit <- quantile(sapply(unlist(cben.boot), function(x) boot0.632plus_cstat(app, x)), probs=c(.025,.5,.975))

    app <- apparent.discr$cforbenefit.y0hat
    boot.oos <- mean(unlist(cbeny0.boot))
    boot0.632.discr$cforbenefit.y0hat <- boot0.632plus_cstat(app, boot.oos)
    # naive CI (with just 1 oos evaluation in the inner loop)
    boot0.632.discr.90ci$cforbenefit.y0hat <- quantile(sapply(unlist(cbeny0.boot), function(x) boot0.632plus_cstat(app, x)), probs=c(.05,.5,.95))
    boot0.632.discr.95ci$cforbenefit.y0hat <- quantile(sapply(unlist(cbeny0.boot), function(x) boot0.632plus_cstat(app, x)), probs=c(.025,.5,.975))

    app <- apparent.discr$mbcforbenefit
    boot.oos <- mean(unlist(ben.boot))
    boot0.632.discr$mbcforbenefit <- boot0.632plus_cstat(app, boot.oos)
    # naive CI (with just 1 oos evaluation in the inner loop)
    boot0.632.discr.90ci$mbcforbenefit <- quantile(sapply(unlist(ben.boot), function(x) boot0.632plus_cstat(app, x)), probs=c(.05,.5,.95))
    boot0.632.discr.95ci$mbcforbenefit <- quantile(sapply(unlist(ben.boot), function(x) boot0.632plus_cstat(app, x)), probs=c(.025,.5,.975))

    app <- apparent.discr$cforbenefit.new
    boot.oos <- mean(unlist(cben.new.boot))
    boot0.632.discr$cforbenefit.new <- boot0.632plus_cstat(app, boot.oos)
    # naive CI (with just 1 oos evaluation in the inner loop)
    boot0.632.discr.90ci$cforbenefit.new <- quantile(sapply(unlist(cben.new.boot), function(x) boot0.632plus_cstat(app, x)), probs=c(.05,.5,.95))
    boot0.632.discr.95ci$cforbenefit.new <- quantile(sapply(unlist(cben.new.boot), function(x) boot0.632plus_cstat(app, x)), probs=c(.025,.5,.975))

    boot0.632.discr$true.diff.y0y1.mbcb <- apparent.discr$true.diff.y0y1.mbcb
    boot0.632.discr.90ci$true.diff.y0y1.mbcb <- NULL
    boot0.632.discr.95ci$true.diff.y0y1.mbcb <- NULL

    # optimism corrected discr
    boot.opt.discr$cforbenefit <- apparent.discr$cforbenefit - mean(unlist(cben.boot.opt.discr) - unlist(cben.orig.opt.discr))
    boot.opt.discr$cforbenefit.y0hat <- apparent.discr$cforbenefit.y0hat - mean(unlist(cbeny0.boot.opt.discr) - unlist(cbeny0.orig.opt.discr))
    boot.opt.discr$mbcforbenefit <- apparent.discr$mbcforbenefit - mean(unlist(ben.boot.opt.discr) - unlist(ben.orig.opt.discr))
    boot.opt.discr$cforbenefit.new <- apparent.discr$cforbenefit.new - mean(unlist(cben.new.boot.opt.discr) - unlist(cben.new.orig.opt.discr))
    boot.opt.discr$true.diff.y0y1.mbcb <- apparent.discr$true.diff.y0y1.mbcb


    ## Calibration ####

    # CIs based on oos bootstrap evaluations
    boot.oos.90ci.cal$true.oe <- NULL
    boot.oos.90ci.cal$true.trt <- NULL
    boot.oos.90ci.cal$true.both <- NULL
    boot.oos.90ci.cal$oe <- quantile(unlist(oe.boot.0.632), probs=c(.05,.5,.95))
    boot.oos.90ci.cal$empirical.trt <- apply(empirical.trt.0.632, 2, function(x) quantile(x, probs=c(.05,.5,.95)))
    boot.oos.90ci.cal$empirical.both <- apply(empirical.both.0.632, 2, function(x) quantile(x, probs=c(.05,.5,.95)))

    boot.oos.95ci.cal$true.oe <- NULL
    boot.oos.95ci.cal$true.trt <- NULL
    boot.oos.95ci.cal$true.both <- NULL
    boot.oos.95ci.cal$oe <- quantile(unlist(oe.boot.0.632), probs=c(.025,.5,.975))
    boot.oos.95ci.cal$empirical.trt <- apply(empirical.trt.0.632, 2, function(x) quantile(x, probs=c(.025,.5,.975)))
    boot.oos.95ci.cal$empirical.both <- apply(empirical.both.0.632, 2, function(x) quantile(x, probs=c(.025,.5,.975)))

    ## 0.632+ calibration
    app <- apparent.cal$empirical.trt[2]
    boot.oos <- mean(empirical.trt.0.632[ ,2])
    boot0.632.cal$empirical.trt <- c(NA, boot0.632plus_cal_slope(app, boot.oos))
    # naive CIs based on just B=1 inner loop (so not representative for the mean 0.632)
    boot0.632.cal.90ci$empirical.trt <- quantile(sapply(unlist(empirical.trt.0.632[ ,2]), function(x) boot0.632plus_cal_slope(app, x)), probs=c(.05,.5,.95))
    boot0.632.cal.95ci$empirical.trt <- quantile(sapply(unlist(empirical.trt.0.632[ ,2]), function(x) boot0.632plus_cal_slope(app, x)), probs=c(.025,.5,.975))

    app <- apparent.cal$empirical.both[2]
    boot.oos <- mean(empirical.both.0.632[ ,2])
    boot0.632.cal$empirical.both <- c(NA, boot0.632plus_cal_slope(app, boot.oos))
    # naive CIs based on just B=1 inner loop (so not representative for the mean 0.632)
    boot0.632.cal.90ci$empirical.both <- quantile(sapply(unlist(empirical.both.0.632[ ,2]), function(x) boot0.632plus_cal_slope(app, x)), probs=c(.05,.5,.95))
    boot0.632.cal.95ci$empirical.both <- quantile(sapply(unlist(empirical.both.0.632[ ,2]), function(x) boot0.632plus_cal_slope(app, x)), probs=c(.025,.5,.975))

    boot0.632.cal.90ci$oe <- boot0.632.cal.95ci$oe <- NULL
    boot0.632.cal.90ci$true.oe <- boot0.632.cal.95ci$true.oe <- NULL
    boot0.632.cal.90ci$true.trt <- boot0.632.cal.95ci$true.trt <- NULL
    boot0.632.cal.90ci$true.both <- boot0.632.cal.95ci$true.both  <- NULL

    # optimism corrected
    boot.opt.cal$oe <- apparent.cal$oe + mean(oe.boot.orig - oe.boot)
    boot.opt.cal$empirical.trt <- apparent.cal$empirical.trt + apply(empirical.trt.orig - empirical.trt.boot, 2, function(x) mean(x))
    boot.opt.cal$empirical.both <- apparent.cal$empirical.both + apply(empirical.both.orig - empirical.both.boot, 2, function(x) mean(x))

    # boot cal references
    boot0.632.cal$true.oe <- boot.opt.cal$true.oe <- apparent.cal$true.oe
    boot0.632.cal$true.trt <- boot.opt.cal$true.trt  <- apparent.cal$true.trt
    boot0.632.cal$true.both <- boot.opt.cal$true.both  <- apparent.cal$true.both


    #######################
    ## new data same DGM ##
    #######################
    Sample <- c(sample(which(val1.out$trt==0), n.val/2, replace=FALSE),
                sample(which(val1.out$trt==1), n.val/2, replace=FALSE))

    ctrlX <- val1.ctrlX[Sample, ]
    trtX <- val1.trtX[Sample, ]
    X <- val1.X[Sample, ]
    trt <- X[ ,"trt"]
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
    ext1.app.discr$cforbenefit.new <- cben.new(form=trt ~ x1 + x2, X=X, y=y, y0hat=y0hat, y1hat=y1hat)
    ext1.app.discr$true.diff.y0y1.mbcb <- mbcb(y0hat, y1hat, p0, p1)

    # total.update discr
    mod.updated <- glm(formula = y ~ -1 + X, family = "binomial")
    y0hat.updated <- plogis(ctrlX %*% coef(mod.updated))
    y1hat.updated <- plogis(trtX %*% coef(mod.updated))

    ext1.total.discr$cforbenefit <- NA
    ext1.total.discr$cforbenefit.y0hat <- cbeny0(y0hat.updated, I(y0hat.updated + deltahat), y, ind.A, ind.B, match="y0hat", nresample=1000)
    ext1.total.discr$mbcforbenefit <- mbcb(y0hat, y1hat, y0hat.updated, y1hat.updated)
    ext1.total.discr$cforbenefit.new <- NA
    ext1.total.discr$true.diff.y0y1.mbcb <- ext1.app.discr$true.diff.y0y1.mbcb

    # apparent cal
    ext1.app.cal$oe <- oe(y0hat,y1hat,y,ind.A,ind.B)
    ext1.app.cal$empirical.trt <- coef(cal(y0hat,y1hat,y,trt=trt,type="treated"))
    ext1.app.cal$empirical.both <- coef(cal(y0hat,y1hat,y,trt=trt,type="both"))
    ext1.app.cal$true.oe <- oe(y0hat,y1hat,p01,ind.A,ind.B)
    ext1.app.cal$true.trt <- coef(cal(y0hat,y1hat,y,trt=trt,p0=p0,p1=p1,type="treated"))
    ext1.app.cal$true.both <- coef(cal(y0hat,y1hat,y,trt=trt,p0=p0,p1=p1,type="both"))

    # total.update
    ext1.total.cal$oe <- oe(y0hat.updated,y1hat.updated,y,ind.A,ind.B)
    ext1.total.cal$empirical.trt <- coef(cal(y0hat,y1hat,y,trt=trt,p0=y0hat.updated,p1=y1hat.updated,type="treated"))
    ext1.total.cal$empirical.both <- coef(cal(y0hat,y1hat,y,trt=trt,p0=y0hat.updated,p1=y1hat.updated,type="both"))
    ext1.total.cal$true.oe <- oe(y0hat.updated, y1hat.updated, p01, ind.A, ind.B)
    ext1.total.cal$true.trt <- coef(cal(y0hat,y1hat,y,trt=trt,p0=p0,p1=p1,type="treated"))
    ext1.total.cal$true.both <- coef(cal(y0hat,y1hat,y,trt=trt,p0=p0,p1=p1,type="both"))

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
    trt <- X[ ,"trt"]
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
    ext2.app.discr$cforbenefit.new <- cben.new(form=trt ~ x1 + x2, X=X, y=y, y0hat=y0hat, y1hat=y1hat)
    ext2.app.discr$true.diff.y0y1.mbcb <- mbcb(y0hat, y1hat, p0, p1)

    # total.update
    mod.updated <- glm(formula = y ~ -1 + X, family = "binomial")
    y0hat.updated <- plogis(ctrlX %*% coef(mod.updated))
    y1hat.updated <- plogis(trtX %*% coef(mod.updated))

    ext2.total.discr$cforbenefit <- NA
    ext2.total.discr$cforbenefit.y0hat <- cbeny0(y0hat.updated, I(y0hat.updated + deltahat), y, ind.A, ind.B, match="y0hat", nresample=1000)
    ext2.total.discr$mbcforbenefit <- mbcb(y0hat, y1hat, y0hat.updated, y1hat.updated)
    ext2.total.discr$cforbenefit.new <- NA
    ext2.total.discr$true.diff.y0y1.mbcb <- ext2.app.discr$true.diff.y0y1.mbcb

    # apparent calibration
    ext2.app.cal$oe <- oe(y0hat,y1hat,y,ind.A,ind.B)
    ext2.app.cal$empirical.trt <- coef(cal(y0hat,y1hat,y,trt=trt,type="treated"))
    ext2.app.cal$empirical.both <- coef(cal(y0hat,y1hat,y,trt=trt,type="both"))
    ext2.app.cal$true.oe <- oe(y0hat,y1hat,p01,ind.A,ind.B)
    ext2.app.cal$true.trt <- coef(cal(y0hat,y1hat,y,trt=trt,p0=p0,p1=p1,type="treated"))
    ext2.app.cal$true.both <- coef(cal(y0hat,y1hat,y,trt=trt,p0=p0,p1=p1,type="both"))

    ext2.total.cal$oe <- oe(y0hat.updated,y1hat.updated,y,ind.A,ind.B)
    ext2.total.cal$empirical.trt <- coef(cal(y0hat,y1hat,y,trt=trt,p0=y0hat.updated,p1=y1hat.updated,type="treated"))
    ext2.total.cal$empirical.both <- coef(cal(y0hat,y1hat,y,trt=trt,p0=y0hat.updated,p1=y1hat.updated,type="both"))
    ext2.total.cal$true.oe <- oe(y0hat.updated, y1hat.updated, p01, ind.A, ind.B)
    ext2.total.cal$true.trt <- coef(cal(y0hat,y1hat,y,trt=trt,p0=p0,p1=p1,type="treated"))
    ext2.total.cal$true.both <- coef(cal(y0hat,y1hat,y,trt=trt,p0=p0,p1=p1,type="both"))

    # risk.cal val DGM 2
    risk.cal$ctrl.val2 <- val.prob(y0hat[ind.A], y[ind.A], pl=FALSE)[val.prob.cols]
    risk.cal$trt.val2 <- val.prob(y1hat[ind.B], y[ind.B], pl=FALSE)[val.prob.cols]
    risk.cal$overall.val2 <- val.prob(yhat, y, pl=FALSE)[val.prob.cols]

    # rmspe
    rmspe$rmspe.p01.val2   <- sqrt(mean((p01 - yhat)^2))
    rmspe$rmspe.delta.val2 <- sqrt(mean((delta - deltahat)^2))

    results[[g]] <- list(apparent.discr, ext1.app.discr, ext1.total.discr, ext2.app.discr, ext2.total.discr,
                         boot0.632.discr, boot.opt.discr, boot0.632.discr.90ci, boot0.632.discr.95ci,
                         boot.oos.90ci.discr, boot.oos.95ci.discr,
                         apparent.cal, ext1.app.cal, ext1.total.cal, ext2.app.cal, ext2.total.cal,
                         boot0.632.cal, boot.opt.cal, boot0.632.cal.90ci, boot0.632.cal.95ci,
                         boot.oos.90ci.cal, boot.oos.90ci.cal,
                         risk.cal, rmspe,
                         "time.elapsed.sec"=(proc.time() - ptm)["elapsed"])
    names(results[[g]]) <- c("apparent.discr", "ext1.app.discr", "ext1.total.discr", "ext2.app.discr", "ext2.total.discr",
                             "boot0.632.discr", "boot.opt.discr", "boot0.632.discr.90ci", "boot0.632.discr.95ci",
                             "boot.oos.90ci.discr", "boot.oos.95ci.discr",
                             "apparent.cal", "ext1.app.cal", "ext1.total.cal", "ext2.app.cal", "ext2.total.cal",
                             "boot0.632.cal", "boot.opt.cal", "boot0.632.cal.90ci", "boot0.632.cal.95ci",
                             "boot.oos.90ci.cal", "boot.oos.90ci.cal",
                             "risk.cal", "rmspe",
                             "time.elapsed.sec")
  }

  names(results)[1:length(sample.sizes)] <- paste("n", sample.sizes, sep="")

  save(results, file=paste0("results/", "ite", ss, ".Rdata"))

  results

}, cl=cl)

stopCluster(cl)

save(info, file="results/info.RData")
