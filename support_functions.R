## functions
cben.new <- function(form, X, y, y0hat, y1hat){ # according to Maas et al. 2023
  form <- as.formula(form)
  X <- as.data.frame(X)
  rownames(X) <- NULL

  tab <- table(X[ ,"trt"])
  m <- matchit(form, data=X, distance="mahalanobis",
               method = "nearest", estimand = ifelse(tab[1] < tab[2], "ATC", "ATT"))
  M <- get_matches(m, data=X)
  ind.A <- as.numeric(M[which(M$trt %in% 0), "id"])
  ind.B <- as.numeric(M[which(M$trt %in% 1), "id"])
  deltahat_ij <- y1hat[ind.B] - y0hat[ind.A]
  O_ij <- y[ind.B] - y[ind.A]
  return(rcorr.cens(deltahat_ij, O_ij)["C Index"])
}

boot0.632plus_cstat <- function(app, boot){
  boot.prime <- ifelse(app > 0.5, pmax(boot, 0.5), pmin(boot, 0.5))
  R <- ifelse(abs(boot.prime - 0.5) < abs(app - 0.5),
              abs(app - boot.prime) / abs(app - 0.5),
              0)
  weight <- (1-exp(-1)) / (1 - exp(-1)*R)
  app*(1-weight) + weight*boot.prime
}

boot0.632plus_cal_slope <- function(app, boot){
  boot.prime <- ifelse(app > 0, pmax(boot, 0), pmin(boot, 0))
  R <- ifelse(abs(boot.prime - 0) < abs(app - 0),
              abs(app - boot.prime) / abs(app - 0),
              0)
  weight <- (1-exp(-1)) / (1 - exp(-1)*R)
  # app*(1-weight) + weight*boot.prime
  app*(1-weight) + weight*boot.prime
}

oe <- function(y0hat,y1hat,y,ind.A,ind.B){
  delta.obs <- mean(y[ind.B]) - mean(y[ind.A])
  delta.pred <- mean(y1hat[ind.B]) - mean(y0hat[ind.A])
  delta.obs - delta.pred # O-E
}
