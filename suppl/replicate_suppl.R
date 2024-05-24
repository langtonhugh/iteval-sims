library(ggplot2)
library(patchwork)
library(stringr)

# grf == TRUE for grf results and grf==FALSE for result estimating g0 and g1 according to model (19).
grf <- TRUE

if(grf){
  dir <- "suppl/suppl_sim_results_grf/"
} else {
  dir <- "suppl_sim_results"
}

files <- list.files(dir)
files <- grep("ite*", files, value = TRUE)
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

# Structure: sr, nsim, sample.size, measure


## Discrimination ####

discr.list <- c("apparent.discr", "ext1.app.discr", "ext2.app.discr",
                "ext1.total.discr", "ext2.total.discr")

## absolute values version
plot.discr <- function(sr, m, columns, print.msd=TRUE, ylim=NULL, ...){

  ellipsis <- list(...)

  plot.data <- array(as.numeric(sapply(sr, function(x) sapply(x, function(xx) unlist(xx[m])))),
                     dim=c(5, nss, nsim))

  msd <- cbind(
    apply(cbind(t(plot.data[ ,1, ]), t(plot.data[ ,2, ]), t(plot.data[ ,3, ]))[,columns],
          2, mean),
    apply(cbind(t(plot.data[ ,1, ]), t(plot.data[ ,2, ]), t(plot.data[ ,3, ]))[,columns],
          2, sd))

  msd <- data.frame("cstat"=msd[,1], "sd"=msd[,2],
                    upp=msd[,1]+1*msd[,2], low=msd[,1]-1*msd[,2],
                    "statistic"=factor(rep(c("cforbendelta", "cforbeny0", "mbcb", "cforbennew", "theta_d"), 3),
                                       levels = c("cforbendelta", "cforbeny0", "mbcb", "cforbennew", "theta_d")),
                    "sample.size"=factor(rep(sample.sizes, each=5)))

  # order as in manuscript
  msd$statistic <- factor(msd$statistic, levels=levels(msd$statistic)[c(1,4,2,3,5)])

  if(print.msd){
    print(msd)
  }

  if(!is.null(ellipsis$ylab)) ylab <- ellipsis$ylab else ylab <- ""
  if(!is.null(ellipsis$title)) title <- ellipsis$title else title <- ""

  ggplot(msd, aes(x=statistic, y=cstat, group=sample.size, colour=sample.size)) +
    geom_point(position=position_dodge(width=0.7)) +
    geom_errorbar(aes(x=statistic, ymax=upp, ymin=low), na.rm=TRUE, position=position_dodge(width=0.7)) +
    theme_bw() + ylim(ylim) +
    guides(x =  guide_axis(angle = 45)) +
    labs(x = "",y=ylab, title=title, col="Sample size") +
    scale_x_discrete(labels = c(expression("cben-"*hat(delta)),
                                expression("cben"[ppte]),
                                expression("cben-"*hat(y)[0]),
                                "mbcm",
                                bquote(theta["d"]))) +

    theme(text = element_text(size=14),
          axis.text.x=element_text(size=rel(1.3)),
          axis.text.y=element_text(size=rel(1.3)))
}

# make sure the ext1.app.discr results for cforbendelta and cforbennew are available in ext1.total.discr
for(i in 1:nsim){
  for(j in 1:length(sample.sizes)){
    sr[[i]][[j]]$ext1.total.discr$cforbenefit <- sr[[i]][[j]]$ext1.app.discr$cforbenefit
    sr[[i]][[j]]$ext1.total.discr$cforbenefit.new <- sr[[i]][[j]]$ext1.app.discr$cforbenefit.new

    sr[[i]][[j]]$ext2.total.discr$cforbenefit <- sr[[i]][[j]]$ext2.app.discr$cforbenefit
    sr[[i]][[j]]$ext2.total.discr$cforbenefit.new <- sr[[i]][[j]]$ext2.app.discr$cforbenefit.new
  }
}
columns <- 1:15
ylim <- c(0.44, 0.62)
layout <- "
AAABBBCCC
"
p1 <- plot.discr(sr, "apparent.discr", columns, ylim=ylim, title="Apparent (DGM1)", ylab="C-statistic")
p6 <- plot.discr(sr, "ext1.total.discr", columns, ylim=ylim, title="External data DGM1*", ylab="C-statistic")
p7 <- plot.discr(sr, "ext2.total.discr", columns, ylim=ylim, title="External data DGM2*")
pl <- p1 + p6 + p7 +
  plot_layout(guides = "collect", design=layout) &
  theme(legend.position = "bottom")
pl

## plot deviation from estimand
columns <- 1:12 # 4 methods * 3 sample sizes
plot.discr.estimand <- function(sr, m, columns, print.msd=TRUE, ylim=NULL, ...){

  ellipsis <- list(...)

  plot.data <- array(as.numeric(sapply(sr, function(x) sapply(x, function(xx){
    zz <- unlist(xx[m])
    zz[1:4] - zz[5]
  }))), dim=c(4, nss, nsim))

  msd <- cbind(
    apply(cbind(t(plot.data[ ,1, ]), t(plot.data[ ,2, ]), t(plot.data[ ,3, ]))[,columns],
          2, mean),
    apply(cbind(t(plot.data[ ,1, ]), t(plot.data[ ,2, ]), t(plot.data[ ,3, ]))[,columns],
          2, sd))

  msd <- data.frame("cstat"=msd[,1], "sd"=msd[,2],
                    upp=msd[,1]+1*msd[,2], low=msd[,1]-1*msd[,2],
                    "statistic"=factor(rep(c("cforbendelta", "cforbeny0", "mbcb", "cforbennew"), 3),
                                       levels = c("cforbendelta", "cforbeny0", "mbcb", "cforbennew")),
                    "sample.size"=factor(rep(sample.sizes, each=4)))

  # order as in manuscript
  msd$statistic <- factor(msd$statistic, levels=levels(msd$statistic)[c(1,4,2,3)])

  if(print.msd){
    print(msd)
  }

  if(!is.null(ellipsis$ylab)) ylab <- ellipsis$ylab else ylab <- ""
  if(!is.null(ellipsis$title)) title <- ellipsis$title else title <- ""

  ggplot(msd, aes(x=statistic, y=cstat, group=sample.size, colour=sample.size)) +
    geom_hline(yintercept=0, linetype="dashed",
               color = "grey", size=1) +
    geom_point(position=position_dodge(width=0.7)) +
    geom_errorbar(aes(x=statistic, ymax=upp, ymin=low), na.rm=TRUE, position=position_dodge(width=0.7)) +
    theme_bw() + ylim(ylim) +
    guides(x =  guide_axis(angle = 45)) +
    labs(x = "",y=ylab, title=title, col="Sample size") +
    scale_x_discrete(labels = c(expression("cben-"*hat(delta)),
                                expression("cben"[ppte]),
                                expression("cben-"*hat(y)[0]),
                                "mbcm")) +
    theme(text = element_text(size=14),
          axis.text.x=element_text(size=rel(1.3)),
          axis.text.y=element_text(size=rel(1.3)))
}
columns <- 1:12
ylim <- c(-0.05, 0.085)
layout <- "
AAABBBCCC
"
p1 <- plot.discr.estimand(sr, "apparent.discr", columns, ylim=ylim, title="Apparent (DGM1)", ylab="Deviation from estimand")
p6 <- plot.discr.estimand(sr, "ext1.total.discr", columns, ylim=ylim, title="External data DGM1",ylab="Deviation from estimand")
p7 <- plot.discr.estimand(sr, "ext2.total.discr", columns, ylim=ylim, title="External data DGM2")
pl <- p1 + p6 + p7 +
  plot_layout(guides = "collect", design=layout) &
  theme(legend.position = "bottom")
pl


## Calibration ####

library(DescTools)
plot.cal <- function(sr, m, columns, print.msd=TRUE, ylim=NULL, ...){

  ellipsis <- list(...)

  plot.data <- cbind(
    t(array(as.numeric(sapply(sr, function(x) sapply(x, function(xx) unlist(xx[[m]]$empirical.both[1])))),
                               dim=c(nss, nsim))),
    t(array(as.numeric(sapply(sr, function(x) sapply(x, function(xx) unlist(xx[[m]]$true.both[1])))),
                                     dim=c(nss, nsim))),
    t(array(as.numeric(sapply(sr, function(x) sapply(x, function(xx) unlist(xx[[m]]$empirical.both[2])))),
                               dim=c(nss, nsim))),
    t(array(as.numeric(sapply(sr, function(x) sapply(x, function(xx) unlist(xx[[m]]$true.both[2])))),
                          dim=c(nss, nsim))))

  plot.data <- plot.data[ ,as.numeric(matrix(1:12, 4, 3, byrow = TRUE))]

  msd <- cbind(
    apply(plot.data[,columns], 2, function(x) mean(x, trim=.1)),
    apply(plot.data[,columns], 2, function(x) sd(Trim(x, trim=.1))))

  msd <- data.frame("mean"=msd[,1], "sd"=msd[,2],
                    upp=msd[,1]+1*msd[,2], low=msd[,1]-1*msd[,2],
                    "statistic"=factor(rep(c("emp.int", "true.int", "emp.slope", "true.slope"), 3),
                                       levels=c("emp.int", "true.int", "emp.slope", "true.slope")),
                    "sample.size"=factor(rep(sample.sizes, each=4)))

  if(print.msd){
    print(msd)
  }

  if(!is.null(ellipsis$ylab)) ylab <- ellipsis$ylab else ylab <- ""
  if(!is.null(ellipsis$title)) title <- ellipsis$title else title <- ""

  ggplot(msd, aes(x=statistic, y=mean, group=sample.size, colour=sample.size)) +
    geom_point(position=position_dodge(width=0.7)) +
    geom_errorbar(aes(x=statistic, ymax=upp, ymin=low), na.rm=TRUE, position=position_dodge(width=0.7)) +
    theme_bw() + ylim(ylim) +
    guides(x =  guide_axis(angle = 45)) +
    labs(x = "",y=ylab, title=title, col="Sample size") +
    scale_x_discrete(labels = c(
      bquote(hat(beta)[.(0)]),
      bquote(beta[.(0)]),
      bquote(hat(beta)[.(1)]),
      bquote(beta[.(1)]))) +
    theme(text = element_text(size=14),
          axis.text.x=element_text(size=rel(1.3)),
          axis.text.y=element_text(size=rel(1.3)))
}
columns <- 1:12
layout <- "
AAABBBCCC
"
ylim <- c(-2.5, 2.5)
p1 <- plot.cal(sr, "apparent.cal", columns, ylim=ylim, title="Apparent (DGM1)", ylab="Coefficient")
p6 <- plot.cal(sr, "ext1.total.cal", columns, ylim=ylim, title="External data DGM1", ylab="Coefficient")
p7 <- plot.cal(sr, "ext2.total.cal", columns, ylim=ylim, title="External data DGM2")
pl <- p1 + p6 + p7 +
  plot_layout(guides = "collect", design=layout) &
  theme(legend.position = "bottom")
pl

## Deviation from estimands
## Plot int and slope in same plot
plot.cal <- function(sr, m, columns, print.msd=TRUE, ylim=NULL, ...){

  ellipsis <- list(...)

  plot.data <- cbind(
    t(array(as.numeric(sapply(sr, function(x) sapply(x, function(xx) unlist(xx[[m]]$empirical.both[1])))),
                               dim=c(nss, nsim))),
    t(array(as.numeric(sapply(sr, function(x) sapply(x, function(xx) unlist(xx[[m]]$true.both[1])))),
                                     dim=c(nss, nsim))),
    t(array(as.numeric(sapply(sr, function(x) sapply(x, function(xx) unlist(xx[[m]]$empirical.both[2])))),
                               dim=c(nss, nsim))),
    t(array(as.numeric(sapply(sr, function(x) sapply(x, function(xx) unlist(xx[[m]]$true.both[2])))),
                          dim=c(nss, nsim))))

  plot.data <- plot.data[ ,as.numeric(matrix(1:12, 4, 3, byrow = TRUE))]

  plot.data <- plot.data[ ,c(1,3,5,7,9,11)] - plot.data[ ,c(2,4,6,8,10,12)]

  msd <- cbind(
    apply(plot.data, 2, function(x) mean(x, trim=.1)),
    apply(plot.data, 2, function(x) sd(Trim(x, trim=.1))))

  msd <- data.frame("mean"=msd[,1], "sd"=msd[,2],
                    upp=msd[,1]+1*msd[,2], low=msd[,1]-1*msd[,2],
                    "statistic"=factor(rep(c("int", "slope"), 3),
                                       levels=c("int", "slope")),
                    "sample.size"=factor(rep(sample.sizes, each=2)))

  if(print.msd){
    print(msd)
  }

  if(!is.null(ellipsis$ylab)) ylab <- ellipsis$ylab else ylab <- ""
  if(!is.null(ellipsis$title)) title <- ellipsis$title else title <- ""

  ggplot(msd, aes(x=statistic, y=mean, group=sample.size, colour=sample.size)) +
    geom_hline(yintercept=0, linetype="dashed",
               color = "grey", size=1) +
    geom_point(position=position_dodge(width=0.7)) +
    geom_errorbar(aes(x=statistic, ymax=upp, ymin=low), na.rm=TRUE, position=position_dodge(width=0.7)) +
    theme_bw() + ylim(ylim) +
    guides(x =  guide_axis(angle = 45)) +
    labs(x = "",y=ylab, title=title, col="Sample size") +
    scale_x_discrete(labels = c(
      bquote(hat(beta)[.(0)]),
      bquote(hat(beta)[.(1)]))) +
    theme(text = element_text(size=14),
          axis.text.x=element_text(size=rel(1.3)),
          axis.text.y=element_text(size=rel(1.3)))
}
layout <- "
AAABBBCCC
"
ylim <- c(-2.5, 2.5)
p1 <- plot.cal(sr, "apparent.cal", columns, ylim=ylim, title="Apparent (DGM1)", ylab="Deviation from estimand")
p6 <- plot.cal(sr, "ext1.total.cal", columns, ylim=ylim, title="External data DGM1", ylab="Deviation from estimand")
p7 <- plot.cal(sr, "ext2.total.cal", columns, ylim=ylim, title="External data DGM2")
pl <- p1 + p6 + p7 +
  plot_layout(guides = "collect", design=layout) &
  theme(legend.position = "bottom")
pl
