## functions
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
