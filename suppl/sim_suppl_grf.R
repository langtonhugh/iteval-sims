## supplementary sims on misspecification
library(doParallel)
library(pbapply) # adds progress bar
numCores <- parallel::detectCores()
print(numCores)
cl <- makeCluster(8)
registerDoParallel(cl)

info <- sessionInfo()

nsim <- 500

todo <- c()
for(ss in 1:nsim){
  if(!file.exists(file=paste0("suppl/suppl_sim_results_grf/", "ite_supp_grf", ss, ".Rdata"))){ # supp sim (new folder)
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
    library(grf)

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
      ext1.app.discr <- ext1.total.discr <- ext2.app.discr <-
        ext2.total.discr <- apparent.discr

      apparent.cal <- vector("list", length = length(list.of.ite.cal.stats))
      names(apparent.cal) <- list.of.ite.cal.stats
      ext1.app.cal <- ext1.total.cal <- ext2.app.cal <-
        ext2.total.cal <- apparent.cal

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

      y <- ifelse(trt==0, y0, y1)

      # robust model y0
      Xrf <- X[ ,c("x1", "x2")]
      rfy0 <- grf::regression_forest(X = Xrf[trt==0, ], Y=y[trt==0])
      # robust model y1
      rfy1 <- grf::regression_forest(X = Xrf[trt==1, ], Y=y[trt==1])

      y0hat.rob <- ifelse(trt==0,
                          predict(rfy0)$predictions,
                          predict(rfy0, newdata = ctrlX[ ,c("x1", "x2")])$predictions)
      y1hat.rob <- ifelse(trt==1,
                          predict(rfy1)$predictions,
                          predict(rfy1, newdata = trtX[ ,c("x1", "x2")])$predictions)
      deltahat.rob <- y1hat.rob - y0hat.rob
      yhat.rob <- ifelse(trt==0, y0hat.rob, y1hat.rob)

      # misspecified model
      Xmiss <- X[ ,-which(colnames(X) == "trt:x2")]
      mod.miss <- glm(y ~ -1 + Xmiss, family = "binomial")  # supp sim omit x1:trt interaction
      mod.miss.coefs <- as.numeric(c(coef(mod.miss)[1:5], 0)) # supp sim (insert the required 0)

      # supp sim (new code block for predictions according to the misspecified model)
      y0hat.miss <- plogis(ctrlX %*% mod.miss.coefs)
      y1hat.miss <- plogis(trtX %*% mod.miss.coefs)
      deltahat.miss <- y1hat.miss - y0hat.miss
      yhat.miss <- ifelse(trt==0, y0hat.miss, y1hat.miss)

      #######################
      ## Apparent          ##
      #######################

      # supp sim: .miss predictions introduced throughout

      ## Discrimination ##
      apparent.discr$cforbenefit <- cbendelta(deltahat.miss, y, ind.A, ind.B, nresample=1000)
      apparent.discr$cforbenefit.y0hat <- cbeny0(y0hat.miss, y1hat.miss, y, ind.A, ind.B, match="y0hat", nresample=1000)
      apparent.discr$mbcforbenefit <- mbcb(y0hat.miss, y1hat.miss, y0hat.rob, y1hat.rob)
      apparent.discr$cforbenefit.new <- cben.new(form=trt ~ x1 + x2, X=X, y=y, y0hat=y0hat.miss, y1hat=y1hat.miss) # Maas et al. 2023; van Klaveren et al. 2023
      apparent.discr$true.diff.y0y1.mbcb <- mbcb(y0hat.miss, y1hat.miss, p0, p1) # estimand!

      ## Calibration ##
      apparent.cal$oe <- oe(y0hat.miss, y1hat.miss, y, ind.A, ind.B)
      apparent.cal$empirical.trt <- coef(cal(y0hat.miss,y1hat.miss,y,trt=trt,p0=y0hat.rob,p1=y1hat.rob,type="treated"))
      apparent.cal$empirical.both <- coef(cal(y0hat.miss,y1hat.miss,y,trt=trt,p0=y0hat.rob,p1=y1hat.rob,type="both"))
      apparent.cal$true.oe <- oe(y0hat.miss, y1hat.miss, p01, ind.A, ind.B) # so against true probabilities
      apparent.cal$true.trt <- coef(cal(y0hat.miss,y1hat.miss,y,trt=trt,p0=p0,p1=p1,type="treated"))
      apparent.cal$true.both <- coef(cal(y0hat.miss,y1hat.miss,y,trt=trt,p0=p0,p1=p1,type="both"))

      # risk.cal app
      val.prob.cols <- c("C (ROC)", "R2", "Brier", "Intercept", "Slope", "Emax", "E90", "Eavg")
      risk.cal$ctrl.app <- val.prob(y0hat.miss[ind.A], y[ind.A], pl=FALSE)[val.prob.cols]
      risk.cal$trt.app <- val.prob(y1hat.miss[ind.B], y[ind.B], pl=FALSE)[val.prob.cols]
      risk.cal$overall.app <- val.prob(yhat.miss, y, pl=FALSE)[val.prob.cols]

      # rmspe app
      rmspe$rmspe.p01.app <- sqrt(mean((p01 - yhat.miss)^2))
      rmspe$rmspe.delta.app <- sqrt(mean((delta - deltahat.miss)^2))

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

      # apply development rf models
      Xrf <- X[ ,c("x1", "x2")]
      y0hat.rob <- ifelse(trt==0,
                          predict(rfy0)$predictions,
                          predict(rfy0, newdata = ctrlX[ ,c("x1", "x2")])$predictions)
      y1hat.rob <- ifelse(trt==1,
                          predict(rfy1)$predictions,
                          predict(rfy1, newdata = trtX[ ,c("x1", "x2")])$predictions)
      deltahat.rob <- y1hat.rob - y0hat.rob

      # supp sim (new code block for predictions according to the misspecified model)
      y0hat.miss <- plogis(ctrlX %*% mod.miss.coefs)
      y1hat.miss <- plogis(trtX %*% mod.miss.coefs)
      deltahat.miss <- y1hat.miss - y0hat.miss
      deltahatlp.miss <- qlogis(y1hat.miss) - qlogis(y0hat.miss)
      yhat.miss <- ifelse(trt==0, y0hat.miss, y1hat.miss)

      # apparent discr
      ext1.app.discr$cforbenefit <- cbendelta(deltahat.miss, y, ind.A, ind.B, nresample=1000)
      ext1.app.discr$cforbenefit.y0hat <- cbeny0(y0hat.miss, y1hat.miss, y, ind.A, ind.B, match="y0hat", nresample=1000)
      ext1.app.discr$mbcforbenefit <- mbcb(y0hat.miss, y1hat.miss, y0hat.rob, y1hat.rob)
      ext1.app.discr$cforbenefit.new <- cben.new(form=trt ~ x1 + x2, X=X, y=y, y0hat=y0hat.miss, y1hat=y1hat.miss)
      ext1.app.discr$true.diff.y0y1.mbcb <- mbcb(y0hat.miss, y1hat.miss, p0, p1)

      # total.update discr
      # robust models y0,y1
      rfy0.updated <- grf::regression_forest(X = Xrf[trt==0, ], Y=y[trt==0])
      rfy1.updated <- grf::regression_forest(X = Xrf[trt==1, ], Y=y[trt==1])
      y0hat.updated <- ifelse(trt==0,
                          predict(rfy0.updated)$predictions,
                          predict(rfy0.updated, newdata = ctrlX[ ,c("x1", "x2")])$predictions)
      y1hat.updated <- ifelse(trt==1,
                              predict(rfy1.updated)$predictions,
                              predict(rfy1.updated, newdata = trtX[ ,c("x1", "x2")])$predictions)
      deltahat.updated <- y1hat.updated - y0hat.updated

      ext1.total.discr$cforbenefit <- ext1.app.discr$cforbenefit
      ext1.total.discr$cforbenefit.y0hat <- cbeny0(y0hat.updated, I(y0hat.updated + deltahat.miss), y, ind.A, ind.B, match="y0hat", nresample=1000)
      ext1.total.discr$mbcforbenefit <- mbcb(y0hat.miss, y1hat.miss, y0hat.updated, y1hat.updated)
      ext1.total.discr$cforbenefit.new <- ext1.app.discr$cforbenefit.new
      ext1.total.discr$true.diff.y0y1.mbcb <- ext1.app.discr$true.diff.y0y1.mbcb

      # apparent cal
      ext1.app.cal$oe <- oe(y0hat.miss,y1hat.miss,y,ind.A,ind.B)
      ext1.app.cal$empirical.trt <- coef(cal(y0hat.miss,y1hat.miss,y,trt=trt,p0=y0hat.rob,p1=y1hat.rob,type="treated"))
      ext1.app.cal$empirical.both <- coef(cal(y0hat.miss,y1hat.miss,y,trt=trt,p0=y0hat.rob,p1=y1hat.rob,type="both"))
      ext1.app.cal$true.oe <- oe(y0hat.miss,y1hat.miss,p01,ind.A,ind.B)
      ext1.app.cal$true.trt <- coef(cal(y0hat.miss,y1hat.miss,y,trt=trt,p0=p0,p1=p1,type="treated"))
      ext1.app.cal$true.both <- coef(cal(y0hat.miss,y1hat.miss,y,trt=trt,p0=p0,p1=p1,type="both"))

      # total.update
      ext1.total.cal$oe <- NA
      ext1.total.cal$empirical.trt <- coef(cal(y0hat.miss,y1hat.miss,y,trt=trt,p0=y0hat.updated,p1=y1hat.updated,type="treated"))
      ext1.total.cal$empirical.both <- coef(cal(y0hat.miss,y1hat.miss,y,trt=trt,p0=y0hat.updated,p1=y1hat.updated,type="both"))
      ext1.total.cal$true.oe <- NA
      ext1.total.cal$true.trt <- coef(cal(y0hat.miss,y1hat.miss,y,trt=trt,p0=p0,p1=p1,type="treated"))
      ext1.total.cal$true.both <- coef(cal(y0hat.miss,y1hat.miss,y,trt=trt,p0=p0,p1=p1,type="both"))

      # risk.cal val DGM 1
      risk.cal$ctrl.val1 <- val.prob(y0hat.miss[ind.A], y[ind.A], pl=FALSE)[val.prob.cols]
      risk.cal$trt.val1 <- val.prob(y1hat.miss[ind.B], y[ind.B], pl=FALSE)[val.prob.cols]
      risk.cal$overall.val1 <- val.prob(yhat.miss, y, pl=FALSE)[val.prob.cols]

      # rmspe
      rmspe$rmspe.p01.val1   <- sqrt(mean((p01 - yhat.miss)^2))
      rmspe$rmspe.delta.val1 <- sqrt(mean((delta - deltahat.miss)^2))

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

      # apply development rf models
      Xrf <- X[ ,c("x1", "x2")]
      y0hat.rob <- ifelse(trt==0,
                          predict(rfy0)$predictions,
                          predict(rfy0, newdata = ctrlX[ ,c("x1", "x2")])$predictions)
      y1hat.rob <- ifelse(trt==1,
                          predict(rfy1)$predictions,
                          predict(rfy1, newdata = trtX[ ,c("x1", "x2")])$predictions)
      deltahat.rob <- y1hat.rob - y0hat.rob

      # supp sim (new code block for predictions according to the misspecified model)
      y0hat.miss <- plogis(ctrlX %*% mod.miss.coefs)
      y1hat.miss <- plogis(trtX %*% mod.miss.coefs)
      deltahat.miss <- y1hat.miss - y0hat.miss
      deltahatlp.miss <- qlogis(y1hat.miss) - qlogis(y0hat.miss)
      yhat.miss <- ifelse(trt==0, y0hat.miss, y1hat.miss)

      # apparent.discr
      ext2.app.discr$cforbenefit <- cbendelta(deltahat.miss, y, ind.A, ind.B, nresample=1000)
      ext2.app.discr$cforbenefit.y0hat <- cbeny0(y0hat.miss, y1hat.miss, y, ind.A, ind.B, match="y0hat", nresample=1000)
      ext2.app.discr$mbcforbenefit <- mbcb(y0hat.miss, y1hat.miss, y0hat.rob, y1hat.rob)
      ext2.app.discr$cforbenefit.new <- cben.new(form=trt ~ x1 + x2, X=X, y=y, y0hat=y0hat.miss, y1hat=y1hat.miss)
      ext2.app.discr$true.diff.y0y1.mbcb <- mbcb(y0hat.miss, y1hat.miss, p0, p1)

      # total.update discr
      # robust models y0,y1
      rfy0.updated <- grf::regression_forest(X = Xrf[trt==0, ], Y=y[trt==0])
      rfy1.updated <- grf::regression_forest(X = Xrf[trt==1, ], Y=y[trt==1])
      y0hat.updated <- ifelse(trt==0,
                              predict(rfy0.updated)$predictions,
                              predict(rfy0.updated, newdata = ctrlX[ ,c("x1", "x2")])$predictions)
      y1hat.updated <- ifelse(trt==1,
                              predict(rfy1.updated)$predictions,
                              predict(rfy1.updated, newdata = trtX[ ,c("x1", "x2")])$predictions)
      deltahat.updated <- y1hat.updated - y0hat.updated

      ext2.total.discr$cforbenefit <- NA
      ext2.total.discr$cforbenefit.y0hat <- cbeny0(y0hat.updated, I(y0hat.updated + deltahat.miss), y, ind.A, ind.B, match="y0hat", nresample=1000)
      ext2.total.discr$mbcforbenefit <- mbcb(y0hat.miss, y1hat.miss, y0hat.updated, y1hat.updated)
      ext2.total.discr$cforbenefit.new <- NA
      ext2.total.discr$true.diff.y0y1.mbcb <- ext2.app.discr$true.diff.y0y1.mbcb

      # apparent calibration
      ext2.app.cal$oe <- oe(y0hat.miss,y1hat.miss,y,ind.A,ind.B)
      ext2.app.cal$empirical.trt <- coef(cal(y0hat.miss,y1hat.miss,y,trt=trt,p0=y0hat.rob,p1=y1hat.rob,type="treated"))
      ext2.app.cal$empirical.both <- coef(cal(y0hat.miss,y1hat.miss,y,trt=trt,p0=y0hat.rob,p1=y1hat.rob,type="both"))
      ext2.app.cal$true.oe <- oe(y0hat.miss,y1hat.miss,p01,ind.A,ind.B)
      ext2.app.cal$true.trt <- coef(cal(y0hat.miss,y1hat.miss,y,trt=trt,p0=p0,p1=p1,type="treated"))
      ext2.app.cal$true.both <- coef(cal(y0hat.miss,y1hat.miss,y,trt=trt,p0=p0,p1=p1,type="both"))

      ext2.total.cal$oe <- NA
      ext2.total.cal$empirical.trt <- coef(cal(y0hat.miss,y1hat.miss,y,trt=trt,p0=y0hat.updated,p1=y1hat.updated,type="treated"))
      ext2.total.cal$empirical.both <- coef(cal(y0hat.miss,y1hat.miss,y,trt=trt,p0=y0hat.updated,p1=y1hat.updated,type="both"))
      ext2.total.cal$true.oe <- NA
      ext2.total.cal$true.trt <- coef(cal(y0hat.miss,y1hat.miss,y,trt=trt,p0=p0,p1=p1,type="treated"))
      ext2.total.cal$true.both <- coef(cal(y0hat.miss,y1hat.miss,y,trt=trt,p0=p0,p1=p1,type="both"))

      # risk.cal val DGM 2
      risk.cal$ctrl.val2 <- val.prob(y0hat.miss[ind.A], y[ind.A], pl=FALSE)[val.prob.cols]
      risk.cal$trt.val2 <- val.prob(y1hat.miss[ind.B], y[ind.B], pl=FALSE)[val.prob.cols]
      risk.cal$overall.val2 <- val.prob(yhat.miss, y, pl=FALSE)[val.prob.cols]

      # rmspe
      rmspe$rmspe.p01.val2   <- sqrt(mean((p01 - yhat.miss)^2))
      rmspe$rmspe.delta.val2 <- sqrt(mean((delta - deltahat.miss)^2))

      results[[g]] <- list(apparent.discr, ext1.app.discr, ext1.total.discr, ext2.app.discr, ext2.total.discr,
                           apparent.cal, ext1.app.cal, ext1.total.cal, ext2.app.cal, ext2.total.cal,
                           risk.cal, rmspe,
                           "time.elapsed.sec"=(proc.time() - ptm)["elapsed"])
      names(results[[g]]) <- c("apparent.discr", "ext1.app.discr", "ext1.total.discr", "ext2.app.discr", "ext2.total.discr",
                               "apparent.cal", "ext1.app.cal", "ext1.total.cal", "ext2.app.cal", "ext2.total.cal",
                               "risk.cal", "rmspe", "time.elapsed.sec")
    }

    names(results)[1:length(sample.sizes)] <- paste("n", sample.sizes, sep="")

    save(results, file=paste0("suppl/suppl_sim_results_grf/", "ite_supp_grf", ss, ".Rdata"))

    results

}, cl=cl)

stopCluster(cl)

save(info, file="suppl/suppl_sim_results_grf/info.RData")
