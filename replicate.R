library(ggplot2)
library(patchwork)
library(stringr)
library(xtable)

# folder where the simulation results are stored
dir <- "replicate/"
files <- list.files(dir)
files <- grep("ite*.*RData", files, value = TRUE)
nsim <- length(files)
nss <- 3

numextract <- function(string){
  as.numeric(str_extract(string, "\\-*\\d+"))
}

files <- files[order(numextract(files))]

sr <- plyr::llply(paste0(dir, files), function(x){
  load(x)
  return(list=(r=results))}, .progress = "text")

names(sr) <- paste0("sim", numextract(files))

sample.sizes <- c(500, 750, 1000)

# List structure: sr, nsim, sample.size, measure

names(sr$sim1$n500)

## Discrimination ####

# Monte Carlo error for the individual measures
discr.list <- c("apparent.discr", "ext1.app.discr", "ext2.app.discr",
                "boot0.632.discr", "boot.opt.discr", "ext1.total.discr",
                "ext2.total.discr")

plot.discr <- function(sr, m, columns, print.msd=TRUE, ylim=NULL, ...){

  ellipsis <- list(...)

  plot.data <- array(as.numeric(sapply(sr, function(x) sapply(x, function(xx) unlist(xx[m])))),
        dim=c(5, nss, nsim))

  msd <- cbind(
    apply(cbind(t(plot.data[ ,1, ]), t(plot.data[ ,2, ]), t(plot.data[ ,3, ]))[,columns],
          2, mean),
    apply(cbind(t(plot.data[ ,1, ]), t(plot.data[ ,2, ]), t(plot.data[ ,3, ]))[,columns],
          2, sd))

  msd <- data.frame("cstat"=msd[,1], "se"=msd[,2],
                    upp=msd[,1]+1*msd[,2], low=msd[,1]-1*msd[,2],
                    "statistic"=factor(rep(c("cforbendelta", "cforbeny0", "mbcb", "sample_ref", "population_ref"), 3),
                                       levels = c("cforbendelta", "cforbeny0", "mbcb", "sample_ref", "population_ref")),
                    "sample.size"=factor(rep(sample.sizes, each=5)))

  if(print.msd){
    print(msd)
  }

  if(!is.null(ellipsis$ylab)) ylab <- ellipsis$ylab else ylab <- ""
  if(!is.null(ellipsis$title)) title <- ellipsis$title else title <- ""

  # https://astrostatistics.psu.edu/su07/R/html/grDevices/html/plotmath.html
  ggplot(msd, aes(x=statistic, y=cstat, group=sample.size, colour=sample.size)) +
    geom_point(position=position_dodge(width=0.7)) +
    geom_errorbar(aes(x=statistic, ymax=upp, ymin=low), na.rm=TRUE, position=position_dodge(width=0.7)) +
    theme_bw() + ylim(ylim) +
    guides(x =  guide_axis(angle = 45)) +
    labs(x = "",y=ylab, title=title, col="Sample size") +
    scale_x_discrete(labels = c(expression("cben-"*hat(delta)),
                                expression("cben-"*hat(y)[0]),
                                "mbcm",
                                "sample ref.",
                                "population ref.")) +
    theme(text = element_text(size=14),
          axis.text.x=element_text(size=rel(1.3)),
          axis.text.y=element_text(size=rel(1.3)))
}

columns <- 1:15
ylim <- c(0.49, 0.65)
layout <- "
#AAABBBCCC###
DDDEEEFFFGGGG
"
p1 <- plot.discr(sr, "apparent.discr", columns, ylim=ylim, title="Apparent (DGM1)", ylab="C-statistic")
p2 <- plot.discr(sr, "ext1.app.discr", columns, ylim=ylim, title="External data DGM1")
p3 <- plot.discr(sr, "ext2.app.discr", columns, ylim=ylim, title="External data DGM2")
p4 <- plot.discr(sr, "boot0.632.discr", columns, ylim=ylim, title="Bootstrap\n0.632+", ylab="C-statistic")
p5 <- plot.discr(sr, "boot.opt.discr", columns, ylim=ylim, title="Bootstrap\nOptimism corrected")
p6 <- plot.discr(sr, "ext1.total.discr", columns, ylim=ylim, title="External data DGM1*")
p7 <- plot.discr(sr, "ext2.total.discr", columns, ylim=ylim, title="External data DGM2*")
pl <- p1 + p2 + p3 + p4 + p5 + p6 + p7 +
  plot_layout(guides = "collect", design=layout) &
  theme(legend.position = "bottom")
pl

# Mean table
means <- matrix(NA,  15, length(discr.list))
i <- 1
for(m in discr.list[c(1,4,5,2,3,6,7)]){
  plot.data <- array(as.numeric(sapply(sr, function(x) sapply(x, function(xx) unlist(xx[m])))),
                     dim=c(5, nss, nsim))
  means[ ,i] <- apply(cbind(t(plot.data[ ,1, ]), t(plot.data[ ,2, ]), t(plot.data[ ,3, ]))[,columns],
                 2, mean)
  i <- i+1
}
means <- rbind(t(means[1:5, ]), t(means[1:5+5, ]), t(means[1:5+10, ]))
xtable(means, digits=3)


# rmse
rmse.discr <- function(sr, m, ref=c("sample", "population"), pl=FALSE, digits=NULL){
  rmse.data <- array(as.numeric(sapply(sr, function(x) sapply(x, function(xx) unlist(xx[m])))),
                     dim=c(5, nss, nsim))
  if(ref=="sample"){
    main <- c("cbendelta", "cbeny0hat", "mbcb")
    err <- lapply(1:3, function(x) rmse.data[x,,] - rmse.data[4,,])
    out <- t(sapply(err, function(x) apply(x, 1, function(xx) sqrt(mean(xx^2)))))
    if(pl){
      sapply(1:3, function(p)
        if(sum(is.na(err[[p]])) == 0){
          matplot(apply(err[[p]], 1, function(x) sqrt(cumsum(x^2) / (1:length(x)))), type="l", main=paste(m, main[p]))
        } else {
          plot.new()
        })
    }
  }
  if(ref=="population"){
    main <- c("cbendelta", "cbeny0hat", "mbcb")
    err <- lapply(1:3, function(x) rmse.data[x,,] - rmse.data[5,,])
    out <- t(sapply(err, function(x) apply(x, 1, function(xx) sqrt(mean(xx^2)))))
    if(pl){
      sapply(1:3, function(p)
        if(sum(is.na(err[[p]])) == 0){
          matplot(apply(err[[p]], 1, function(x) sqrt(cumsum(x^2) / (1:length(x)))), type="l", main=paste(m, main[p]))
        } else {
          plot.new()
        })
    }
  }
  colnames(out) <- paste0("n", sample.sizes)
  rownames(out) <- c("cbendelta", "cbeny0hat", "mbcb")
  if(!is.null(digits)) out <- round(out, digits)
  return(out)
}
xtable(t(sapply(discr.list, function(m){
  xx <- t(rmse.discr(sr, m, ref="population", pl=FALSE, digits=3))
  c(paste(xx[1], xx[2], xx[3], sep=","),
    paste(xx[4], xx[5], xx[6], sep=","),
    paste(xx[7], xx[8], xx[9], sep=","))
})))


## Calibration ####
cal.list <- c("apparent.cal", "ext1.app.cal", "ext2.app.cal",
              "boot0.632.cal", "boot.opt.cal", "ext1.total.cal",
              "ext2.total.cal")

plot.cal <- function(sr, m, columns, param, print.msd=TRUE, ylim=NULL, ...){

  ellipsis <- list(...)

  empirical <- ifelse(!is.null(sr$sim1$n500[[m]]$empirical), TRUE, FALSE)
  true.cal <- ifelse(!is.null(sr$sim1$n500[[m]]$true.cal), TRUE, FALSE)
  pop.empirical <- ifelse(!is.null(sr$sim1$n500[[m]]$pop.empirical), TRUE, FALSE)
  pop.cal <- ifelse(!is.null(sr$sim1$n500[[m]]$pop.cal), TRUE, FALSE)

  plot.data <- cbind(
    if(empirical) t(array(as.numeric(sapply(sr, function(x) sapply(x, function(xx) unlist(xx[[m]]$empirical[param])))),
                   dim=c(nss, nsim))) else matrix(NA, nsim, nss),
    if(true.cal) t(array(as.numeric(sapply(sr, function(x) sapply(x, function(xx) unlist(xx[[m]]$true.cal[param])))),
                        dim=c(nss, nsim))) else matrix(NA, nsim, nss),
    if(pop.empirical) t(array(as.numeric(sapply(sr, function(x) sapply(x, function(xx) unlist(xx[[m]]$pop.empirical[param])))),
                       dim=c(nss, nsim))) else matrix(NA, nsim, nss),
    if(pop.cal) t(array(as.numeric(sapply(sr, function(x) sapply(x, function(xx) unlist(xx[[m]]$pop.cal[param])))),
                       dim=c(nss, nsim))) else matrix(NA, nsim, nss))

  plot.data <- plot.data[ ,as.numeric(matrix(1:12, 4, 3, byrow = TRUE))]

  msd <- cbind(
    apply(plot.data[,columns], 2, mean),
    apply(plot.data[,columns], 2, sd))

  msd <- data.frame("mean"=msd[,1], "se"=msd[,2],
                    upp=msd[,1]+1*msd[,2], low=msd[,1]-1*msd[,2],
                    "statistic"=factor(rep(c("beta0hat", "sample_ref", "population_emp", "population_ref"), 3),
                                       levels = c("beta0hat", "population_emp", "sample_ref", "population_ref")),
                    "sample.size"=factor(rep(sample.sizes, each=4)))

  if(print.msd){
    print(msd)
  }

  if(!is.null(ellipsis$ylab)) ylab <- ellipsis$ylab else ylab <- ""
  if(!is.null(ellipsis$title)) title <- ellipsis$title else title <- ""

  # https://astrostatistics.psu.edu/su07/R/html/grDevices/html/plotmath.html
  ggplot(msd, aes(x=statistic, y=mean, group=sample.size, colour=sample.size)) +
    geom_point(position=position_dodge(width=0.7)) +
    geom_errorbar(aes(x=statistic, ymax=upp, ymin=low), na.rm=TRUE, position=position_dodge(width=0.7)) +
    theme_bw() + ylim(ylim) +
    guides(x =  guide_axis(angle = 45)) +
    labs(x = "",y=ylab, title=title, col="Sample size") +
    scale_x_discrete(labels = c(#expression(hat(beta)[param]),
                                bquote(hat(beta)[.(param-1)]),
                                "population naive",
                                "sample ref.",
                                "population ref.")) +
    theme(text = element_text(size=14),
          axis.text.x=element_text(size=rel(1.3)),
          axis.text.y=element_text(size=rel(1.3)))
}
# empirical: coef(cal(y0hat,y1hat,y,ind.B))
# true.cal: coef(cal(p0,y1hat.updated,p1,ind.B))
# pop.naive: coef(cal(y0hat.ref,y1hat.ref,val1.out$y,val1.ind.B))
# pop.cal: coef(cal(val1.out$p0,y1.pop.updated,val1.out$p1))

# Intercepts
columns <- 1:12
ylim <- c(-1,1.5)
p1 <- plot.cal(sr, "apparent.cal", columns, param=1, ylim=ylim, title="Apparent (DGM1)", ylab="Intercept")
p2 <- plot.cal(sr, "ext1.app.cal", columns, param=1, ylim=ylim, title="External data DGM1")
p3 <- plot.cal(sr, "ext2.app.cal", columns, param=1, ylim=ylim, title="External data DGM2")
# p4 <- plot.cal(sr, "boot0.632.cal", columns, param=1, ylim=ylim, title="Bootstrap\n0.632+", ylab="Intercept")
p5 <- plot.cal(sr, "boot.opt.cal", columns, param=1, ylim=ylim, title="Bootstrap\nOptimism corrected", ylab="Intercept")
p6 <- plot.cal(sr, "ext1.total.cal", columns, param=1, ylim=ylim, title="External data DGM1*")
p7 <- plot.cal(sr, "ext2.total.cal", columns, param=1, ylim=ylim, title="External data DGM2*")
layout <- "
AAABBBCCC
DDDEEEFFF
"
pl <- p1 + p2 + p3 + p5 + p6 + p7 +
  plot_layout(guides = "collect", design=layout) &
  theme(legend.position = "bottom")
pl

# Slopes
columns <- 1:12
ylim <- c(-0.5,2)
p1 <- plot.cal(sr, "apparent.cal", columns, param=2, ylim=ylim, title="Apparent (DGM1)", ylab="Slope")
p2 <- plot.cal(sr, "ext1.app.cal", columns, param=2, ylim=ylim, title="External data DGM1")
p3 <- plot.cal(sr, "ext2.app.cal", columns, param=2, ylim=ylim, title="External data DGM2")
p4 <- plot.cal(sr, "boot0.632.cal", columns, param=2, ylim=ylim, title="Bootstrap\n0.632+", ylab="Slope")
p5 <- plot.cal(sr, "boot.opt.cal", columns, param=2, ylim=ylim, title="Bootstrap\nOptimism corrected")
p6 <- plot.cal(sr, "ext1.total.cal", columns, param=2, ylim=ylim, title="External data DGM1*")
p7 <- plot.cal(sr, "ext2.total.cal", columns, param=2, ylim=ylim, title="External data DGM2*")
layout <- "
#AAABBBCCC###
DDDEEEFFFGGGG
"
pl <- p1 + p2 + p3 + p4 + p5 + p6 + p7 +
  plot_layout(guides = "collect", design=layout) &
  theme(legend.position = "bottom")
pl

# means tables
param=1 # param=1 for intercepts and param=2 for slopes
means <- matrix(NA, 12, length(cal.list))
i=1
for(m in cal.list[c(1,4,5,2,3,6,7)]){
  plot.data <- cbind(
    t(array(as.numeric(sapply(sr, function(x) sapply(x, function(xx) unlist(xx[[m]]$empirical[param])))),
            dim=c(nss, nsim))),
    t(array(as.numeric(sapply(sr, function(x) sapply(x, function(xx) unlist(xx[[m]]$pop.empirical[param])))),
            dim=c(nss, nsim))),
    t(array(as.numeric(sapply(sr, function(x) sapply(x, function(xx) unlist(xx[[m]]$true.cal[param])))),
            dim=c(nss, nsim))),
    t(array(as.numeric(sapply(sr, function(x) sapply(x, function(xx) unlist(xx[[m]]$pop.cal[param])))),
            dim=c(nss, nsim))))
  plot.data <- plot.data[ ,as.numeric(matrix(1:12, 4, 3, byrow = TRUE))]
  means[,i] <- apply(plot.data[,columns], 2, mean)
  i <- i+1
}
means <- rbind(t(means[1:4, ]), t(means[1:4+4, ]), t(means[1:4+8, ]))
xtable(means, digits=3)

# rmse for cal function
rmse.cal <- function(sr, m, param, ref=c("sample.p0p1", "population.empirical", "population.p0p1"), pl=FALSE, digits=NULL){
  empirical <- ifelse(!is.null(sr$sim1$n500[[m]]$empirical), TRUE, FALSE)
  true.cal <- ifelse(!is.null(sr$sim1$n500[[m]]$true.cal), TRUE, FALSE)
  pop.empirical <- ifelse(!is.null(sr$sim1$n500[[m]]$pop.empirical), TRUE, FALSE)
  pop.cal <- ifelse(!is.null(sr$sim1$n500[[m]]$pop.cal), TRUE, FALSE)

  rmse.data <- cbind(
    if(empirical) t(array(as.numeric(sapply(sr, function(x) sapply(x, function(xx) unlist(xx[[m]]$empirical[param])))),
                          dim=c(nss, nsim))) else matrix(NA, nsim, nss),
    if(true.cal) t(array(as.numeric(sapply(sr, function(x) sapply(x, function(xx) unlist(xx[[m]]$true.cal[param])))),
                         dim=c(nss, nsim))) else matrix(NA, nsim, nss),
    if(pop.empirical) t(array(as.numeric(sapply(sr, function(x) sapply(x, function(xx) unlist(xx[[m]]$pop.empirical[param])))),
                              dim=c(nss, nsim))) else matrix(NA, nsim, nss),
    if(pop.cal) t(array(as.numeric(sapply(sr, function(x) sapply(x, function(xx) unlist(xx[[m]]$pop.cal[param])))),
                        dim=c(nss, nsim))) else matrix(NA, nsim, nss))

  if(!(empirical & true.cal & pop.empirical & pop.cal)){
    warning("not all measures available; manually check columns")
  }

  rmse.data <- rmse.data[ ,as.numeric(matrix(1:12, 4, 3, byrow = TRUE))]

  colnames(rmse.data) <- paste0(rep(c("beta0hat", "sample_ref", "population_emp", "population_ref"), 3),
                             paste0("_n", rep(sample.sizes, each=4)))

  if(ref=="sample.p0p1"){
    err <- rmse.data[,c(1,5,9)] - rmse.data[,c(1,5,9) + 1]
    if(pl){
      if(sum(is.na(err)) == 0){
        matplot(apply(err, 2, function(x) sqrt(cumsum(x^2) / (1:length(x)))), type="l", main=paste(m, param))
      } else {
        plot.new()
      }
    }
  }
  if(ref=="population.empirical"){
    err <- rmse.data[,c(1,5,9)] - rmse.data[,c(1,5,9) + 2]
    if(pl){
      if(sum(is.na(err)) == 0){
        matplot(apply(err, 2, function(x) sqrt(cumsum(x^2) / (1:length(x)))), type="l", main=paste(m, param))
      } else {
        plot.new()
      }
    }
  }
  if(ref=="population.p0p1"){
    err <- rmse.data[,c(1,5,9)] - rmse.data[,c(1,5,9) + 3]
    if(pl){
      if(sum(is.na(err)) == 0){
        matplot(apply(err, 2, function(x) sqrt(cumsum(x^2) / (1:length(x)))), type="l", main=paste(m, param))
      } else {
        plot.new()
      }
    }
  }
  out <- apply(err, 2, function(x) sqrt(mean(x^2)))
  names(out) <- paste0("n", sample.sizes)
  return(out)
}

# rmse cal estimates
xtable(t(sapply(cal.list, function(m){
  xx <- t(cbind(rmse.cal(sr, m, param=1, ref="population.p0p1", pl=F, digits=3),
                rmse.cal(sr, m, param=2, ref="population.p0p1", pl=F, digits=3)))
  xx <- round(xx, 3)
  c(paste(xx[1,1], xx[1,2], xx[1,3], sep=","),
    paste(xx[2,1], xx[2,2], xx[2,3], sep=","))
})))

