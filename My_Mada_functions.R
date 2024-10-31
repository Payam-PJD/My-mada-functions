library(scales)
library(mada)
library(meta)
library(dmetatools)
library(MASS)
library(doParallel)
library(foreach)
library(SciViews)
library(MVPBT)
library(metafor)
library(dplyr)
library(yarrr)
library(ggplot2)
library(patchwork)
library(grid)
library(gridExtra)
library(gtable)
library(gridGraphics)
library(ggplotify) 



rename <- function(dat, option){
  
  for (metric in c('TP', 'TN', 'FP', 'FN')){
    oldvar <- paste(metric, option)
    dat[[metric]] <- round(dat[[oldvar]])
  }
  return(dat)
}

rename.first <- function(dat, option){
  
  for (metric in c('TP', 'TN', 'FP', 'FN')){
    oldvar <- paste(option, metric)
    dat[[metric]] <- round(dat[[oldvar]])
  }
  return(dat)
}


# het.string <- function(reitsma.fit){
#   summary.reitsma <- summary(reitsma.fit)
#   Zhou <- percent(summary.reitsma$i2$Zhou, accuracy = 0.1)
#   Holling.ua <- paste("[",
#                       percent(summary.reitsma$i2$HollingUnadjusted3, accuracy = 0.1),
#                       "-", 
#                       percent(summary.reitsma$i2$HollingUnadjusted2, accuracy = 0.1),
#                       "]",
#                       sep = "")
#   Holling.a <- paste("[", 
#                      percent(summary.reitsma$i2$HollingAdjusted3, accuracy = 0.1),
#                      "-", 
#                      percent(summary.reitsma$i2$HollingAdjusted2, accuracy = 0.1),
#                      "]",
#                      sep = "")
#   combined.string <- paste(
#     "I^2: ",
#     "ZD: ",
#     Zhou,
#     ", Holling(ua): ",
#     Holling.ua,
#     ", Holling(a): ",
#     Holling.a,
#     sep = ""
#   )
#   return(combined.string)
# }

het.string <- function(reitsma.fit){
  summary.reitsma <- summary(reitsma.fit)
  Holling.ua <- paste("[",
                      percent(summary.reitsma$i2$HollingUnadjusted3, accuracy = 0.1),
                      "-", 
                      percent(summary.reitsma$i2$HollingUnadjusted2, accuracy = 0.1),
                      "]",
                      sep = "")
  combined.string <- paste(
    "I^2: ",
    Holling.ua,
    sep = ""
  )
  return(combined.string)
}








forest.diag <- function(dat,
                        lcols1 =NULL,
                        llab1 = NULL,
                        lcols2 =NULL,
                        llab2 = NULL,
                        object.return = F,
                        sens.forest = T,
                        spec.forest = F,
                        plot.het = T,
                        calcwidth.addline.opt = T,
                        xlim = c(50,100))
  {
  dat[["lcols1"]] <- lcols1
  dat[["lcols2"]] <- lcols2
  dat <- dat[which(!is.na(dat[["TP"]]) & !is.na(dat[["TN"]]) & !is.na(dat[["FP"]]) & !is.na(dat[["FN"]])), ]
  reitsma.fit <- reitsma(dat, method ='ml')
  summary.reitsma <- summary(reitsma.fit)
  metaprop.sens <- metaprop(data = dat,
                            studlab = dat[["names"]],
                            event = dat[["TP"]],
                            n = dat[["TP"]] + dat[["FN"]], method.tau =  "ml",
                            sm = "PRAW",
                            outclab = "sensitivity",
                            method.ci = "WS",
                            common = FALSE)
  metaprop.spec <- metaprop(data = dat,
                            studlab = dat[["names"]],
                            event = dat[["TN"]],
                            n = dat[["TN"]] + dat[["FP"]], method.tau =  "ml",
                            sm = "PRAW",
                            outclab = "specificity",
                            method.ci = "WS",
                            common = FALSE)
  
  
  madad.dat <- madad(dat)
  metaprop.sens$upper <- madad.dat$sens$sens.ci[,2]
  metaprop.sens$lower <- madad.dat$sens$sens.ci[,1]
  metaprop.sens$TE.random <- summary.reitsma$coefficients["sensitivity", "Estimate"]
  metaprop.sens$lower.random <- summary.reitsma$coefficients["sensitivity", "95%ci.lb"]
  metaprop.sens$upper.random <- summary.reitsma$coefficients["sensitivity", "95%ci.ub"]
  # metaprop.sens$lcols <- dat[["lcols"]]
  metaprop.spec$upper <- madad.dat$spec$spec.ci[,2]
  metaprop.spec$lower <- madad.dat$spec$spec.ci[,1]
  metaprop.spec$TE.random <- 1- summary.reitsma$coefficients[4,1]
  metaprop.spec$upper.random <- 1- summary.reitsma$coefficients[4,5]
  metaprop.spec$lower.random <- 1- summary.reitsma$coefficients[4,6]
  # metaprop.spec$lcols <- dat[["lcols"]]
  combined.set <- list()
  combined.set$reitsma <- reitsma.fit
  combined.set$summary.reitsma <- summary.reitsma
  combined.set$metaprop.sens <- metaprop.sens
  combined.set$metaprop.spec <- metaprop.spec
  if (plot.het){
    plot.string.het.overall <- het.string(reitsma.fit)
  }else{
    plot.string.het.overall <- " "
  }
  
  
  if (is.null(lcols1)){
    lcols1.string <- NULL
  }else{
    lcols1.string <- "lcols1"
  }
  if (is.null(lcols2)){
    lcols2.string <- NULL
  }else{
    lcols2.string <- "lcols2"
  }
  if (sens.forest){
    meta::forest(metaprop.sens,
           xlim = xlim,
           pscale = 100,
           just.addcols.right = "center",
           rightcols = c("effect", "ci"),
           rightlabs = c("Sensitivity %", "95% C.I.     "),
           leftcols = c("studlab", lcols1.string, lcols2.string),
           leftlabs = c("Study", llab1, llab2),
           xlab = "Sensitivity", smlab = "",
           weight.study = "random", squaresize = 0.7, col.square = "navy",
           col.square.lines = "navy",
           col.diamond = "maroon", col.diamond.lines = "maroon",
           pooled.totals = FALSE,
           comb.fixed = FALSE,
           fs.hetstat = 10,
           print.tau2 = FALSE,
           print.Q = FALSE,
           print.pval.Q = FALSE,
           print.I2 = FALSE,
           digits = 1,
           hetstat = FALSE,
           text.addline1 = plot.string.het.overall,
           ref = 100 * combined.set$metaprop.sens$TE.random,
           calcwidth.addline = calcwidth.addline.opt,
           just = "center"
    )
  }
  if (spec.forest){
    meta::forest(metaprop.spec,
           xlim = xlim,
           pscale = 100,
           just.addcols.right = "center",
           rightcols = c("effect", "ci"),
           rightlabs = c("Specificty %", "95% C.I.     "),
           leftcols = c("studlab", lcols1.string, lcols2.string),
           leftlabs = c("Study", llab1, llab2),
           xlab = "Specifcity", smlab = "",
           weight.study = "random", squaresize = 0.7, col.square = "navy",
           col.square.lines = "navy",
           col.diamond = "maroon", col.diamond.lines = "maroon",
           pooled.totals = FALSE,
           comb.fixed = FALSE,
           fs.hetstat = 10,
           print.tau2 = FALSE,
           print.Q = FALSE,
           print.pval.Q = FALSE,
           print.I2 = FALSE,
           digits = 1,
           hetstat = FALSE,
           text.addline1 = plot.string.het.overall,
           ref = 100 * combined.set$metaprop.spec$TE.random,
           calcwidth.addline = calcwidth.addline.opt,
           just = "center"
    )
  }
    # Calculate sensitivity estimates and confidence intervals
  sens_estimate <- 100 * combined.set$metaprop.sens$TE.random
  sens_lower <- 100 * combined.set$metaprop.sens$lower.random
  sens_upper <- 100 * combined.set$metaprop.sens$upper.random
  
  # Format to 2 decimal places
  sens_estimate_formatted <- sprintf("%.2f", sens_estimate)
  sens_lower_formatted <- sprintf("%.2f", sens_lower)
  sens_upper_formatted <- sprintf("%.2f", sens_upper)
  
  # Print sensitivity
  cat("Sensitivity:", sens_estimate_formatted, "% (95% CI:", sens_lower_formatted, "-", sens_upper_formatted, "%)\n")
  
  # Calculate specificity estimates and confidence intervals
  spec_estimate <- 100 * combined.set$metaprop.spec$TE.random
  spec_lower <- 100 * combined.set$metaprop.spec$lower.random
  spec_upper <- 100 * combined.set$metaprop.spec$upper.random
  
  # Format to 2 decimal places
  spec_estimate_formatted <- sprintf("%.2f", spec_estimate)
  spec_lower_formatted <- sprintf("%.2f", spec_lower)
  spec_upper_formatted <- sprintf("%.2f", spec_upper)
  
  # Print specificity
  cat("Specificity:", spec_estimate_formatted, "% (95% CI:", spec_lower_formatted, "-", spec_upper_formatted, "%)\n")
  if (object.return){
    return(combined.set)
  }
  }


forest.diag.combined <- function(dat,
                                 lcols1 = NULL,
                                 llab1 = NULL,
                                 lcols2 = NULL,
                                 llab2 = NULL,
                                 object.return = FALSE,
                                 plot.het = TRUE,
                                 calcwidth.addline.opt = TRUE,
                                 xlim = c(50, 100),
                                 leftspace = "    ",
                                 rightspace = "        ",
                                 ratio = c(.7,.3),
                                 ...
                                 ) {
  # Prepare data
  dat[["lcols1"]] <- lcols1
  dat[["lcols2"]] <- lcols2
  dat <- dat[complete.cases(dat[, c("TP", "TN", "FP", "FN")]), ]
  
  # Fit Reitsma model
  reitsma.fit <- reitsma(dat, method = 'ml')
  summary.reitsma <- summary(reitsma.fit)
  
  # Meta-analysis for sensitivity
  metaprop.sens <- metaprop(
    data = dat,
    studlab = dat[["names"]],
    event = dat[["TP"]],
    n = dat[["TP"]] + dat[["FN"]],
    method.tau = "ML",
    sm = "PRAW",
    outclab = "Sensitivity",
    method.ci = "WS",
    common = FALSE
  )
  
  # Meta-analysis for specificity
  metaprop.spec <- metaprop(
    data = dat,
    studlab = dat[["names"]],
    event = dat[["TN"]],
    n = dat[["TN"]] + dat[["FP"]],
    method.tau = "ML",
    sm = "PRAW",
    outclab = "Specificity",
    method.ci = "WS",
    common = FALSE
  )
  
  # Update metaprop objects with confidence intervals from mada
  madad.dat <- madad(dat)
  
  # Sensitivity
  metaprop.sens$upper <- madad.dat$sens$sens.ci[, 2]
  metaprop.sens$lower <- madad.dat$sens$sens.ci[, 1]
  metaprop.sens$TE.random <- summary.reitsma$coefficients["sensitivity", "Estimate"]
  metaprop.sens$lower.random <- summary.reitsma$coefficients["sensitivity", "95%ci.lb"]
  metaprop.sens$upper.random <- summary.reitsma$coefficients["sensitivity", "95%ci.ub"]
  
  # Specificity (using numeric indices as per your instructions)
  metaprop.spec$upper <- madad.dat$spec$spec.ci[, 2]
  metaprop.spec$lower <- madad.dat$spec$spec.ci[, 1]
  metaprop.spec$TE.random <- 1 - summary.reitsma$coefficients[4, 1]
  metaprop.spec$lower.random <- 1 - summary.reitsma$coefficients[4, 5]
  metaprop.spec$upper.random <- 1 - summary.reitsma$coefficients[4, 6]
  
  # Prepare left column labels
  if (is.null(lcols1)) {
    lcols1.string <- NULL
  } else {
    lcols1.string <- "lcols1"
  }
  if (is.null(lcols2)) {
    lcols2.string <- NULL
  } else {
    lcols2.string <- "lcols2"
  }
  
  leftcols_sens <- c("studlab", lcols1.string, lcols2.string)
  leftlabs_sens <- c("Study", llab1, llab2)
  
  # Remove NULL values from leftcols and leftlabs
  leftcols_sens <- leftcols_sens[!sapply(leftcols_sens, is.null)]
  leftlabs_sens <- leftlabs_sens[!sapply(leftlabs_sens, is.null)]
  
  # Generate the sensitivity forest plot and capture it as a grob
  metaprop.sens$text.random <- "Bivariate model"
  plot_sens <- function() {
    meta::forest(
      metaprop.sens,
      xlim = xlim,
      pscale = 100,
      just.addcols.right = "center",
      rightcols = c("effect", "ci"),
      rightlabs = c("Sensitivity %", "95% C.I."),
      leftcols = leftcols_sens,
      leftlabs = leftlabs_sens,
      xlab = "Sensitivity",
      smlab = "",
      weight.study = "random",
      squaresize = 0.7,
      col.square = "navy",
      col.square.lines = "navy",
      col.diamond = "maroon",
      col.diamond.lines = "maroon",
      pooled.totals = FALSE,
      comb.fixed = FALSE,
      fs.hetstat = 10,
      print.tau2 = FALSE,
      print.Q = FALSE,
      print.pval.Q = FALSE,
      print.I2 = FALSE,
      digits = 1,
      hetstat = FALSE,
      text.addline1 = if (plot.het) het.string(reitsma.fit) else " ",
      ref = 100 * metaprop.sens$TE.random,
      calcwidth.addline = calcwidth.addline.opt,
      just = "center"
    )
  }
  grob.sens <- ggplotify::as.grob(plot_sens)
  
  # Generate the specificity forest plot and capture it as a grob
  
  metaprop.spec$leftspace <- rep(leftspace, length(metaprop.spec$studlab))
  metaprop.spec$rightspace <- rep(rightspace, length(metaprop.spec$studlab))
  metaprop.spec$text.random <- "        "
  plot_spec <- function() {
    meta::forest(
      metaprop.spec,
      xlim = xlim,
      pscale = 100,
      just.addcols.right = "center",
      rightcols = c("effect", "ci", "rightspace"),
      rightlabs = c("Specificity %", "95% C.I.", "  "),
      leftcols = "leftspace", # Exclude left columns
      leftlabs = NULL,
      studlab = F,
      xlab = "Specificity",
      smlab = "",
      weight.study = "random",
      squaresize = 0.7,
      col.square = "navy",
      col.square.lines = "navy",
      col.diamond = "maroon",
      col.diamond.lines = "maroon",
      pooled.totals = FALSE,
      comb.fixed = FALSE,
      fs.hetstat = 10,
      print.tau2 = FALSE,
      print.Q = FALSE,
      print.pval.Q = FALSE,
      print.I2 = FALSE,
      digits = 1,
      hetstat = FALSE,
      text.addline1 =  " ",
      ref = 100 * metaprop.spec$TE.random,
      calcwidth.addline = calcwidth.addline.opt,
      just = "center"
    )
  }
  grob.spec <- ggplotify::as.grob(plot_spec)
  
  # Adjust grob.spec to remove left padding
  grob.spec$widths[1] <- unit(0, "cm")
  # Combine the two grobs
  combined <- grid.arrange(grob.sens, grob.spec, ncol = 2, widths = ratio)
  
  # Print overall estimates
  sens_estimate <- 100 * metaprop.sens$TE.random
  sens_lower <- 100 * metaprop.sens$lower.random
  sens_upper <- 100 * metaprop.sens$upper.random
  
  spec_estimate <- 100 * metaprop.spec$TE.random
  spec_lower <- 100 * metaprop.spec$lower.random
  spec_upper <- 100 * metaprop.spec$upper.random
  
  cat("Sensitivity:", sprintf("%.2f", sens_estimate), "% (95% CI:", sprintf("%.2f", sens_lower), "-", sprintf("%.2f", sens_upper), "%)\n")
  cat("Specificity:", sprintf("%.2f", spec_estimate), "% (95% CI:", sprintf("%.2f", spec_lower), "-", sprintf("%.2f", spec_upper), "%)\n")
  
  # Return objects if requested
  if (object.return) {
    combined.set <- list(
      reitsma = reitsma.fit,
      summary.reitsma = summary.reitsma,
      metaprop.sens = metaprop.sens,
      metaprop.spec = metaprop.spec
    )
    return(combined.set)
  }
}



forest.diag.subgroup <- function (dat, #dataframe with TP, TN, FP, and FN and a variable named "names"
                                  subgrouping.variable, #dat$subgrouping.variable
                                  sortvar = NULL ,
                                  sglabel = "subgroup" , #string for subgroup labe;
                                  lcols1 = NULL, #Variable to be shown in left side dat$lcols
                                  llab1 = NULL, #a string for the label of lcols
                                  lcols2 = NULL, #Variable to be shown in left side dat$lcols
                                  llab2 = NULL, #a string for the label of lcols
                                  object.return = T, # IF True it returns an object any object used in the function
                                  sens.forest = T, #If true, it draws  sensitivity forest plot
                                  spec.forest = F,# If true it draws specificity forrest plot
                                  only.subgroups.bigger.than.3 = T, # If true it excludes subgroups smaller than 3 subjets.
                                  plot.het.overall =T,
                                  plot.het.subgroup =T,
                                  plot.overall = T,
                                  calcwidth.shet.opt =F,
                                  forest.xlim = c(50,100)
                                  )
{
  dat[["subgrouping.variable"]] <- subgrouping.variable
  dat[["lcols1"]] <- lcols1
  dat[["lcols2"]] <- lcols2
  if (is.null(sortvar)) {
    dat <- dat[order(dat[["subgrouping.variable"]]), ]
  } else {
    dat <- dat[order(dat[["subgrouping.variable"]], sortvar), ]
  }
  dat <- dat[which(!is.na(dat[["TP"]]) & !is.na(dat[["TN"]]) & !is.na(dat[["FP"]]) & !is.na(dat[["FN"]])), ]
  subgroup.list <- unique(dat[["subgrouping.variable"]])
  counts <- table(dat[["subgrouping.variable"]])
  if (only.subgroups.bigger.than.3){
    subgroup.list <- names(counts[counts >= 3])
  }
  subgroup.list <- sort(subgroup.list)
  valid.subgroup.list <- make.names(subgroup.list, unique = T)
  dat <- subset(dat, dat[["subgrouping.variable"]] %in% subgroup.list)
  reitsmas <- list()
  reitsmas$reitsma.overall <- reitsma(data = dat, method = "ml")
  reitsmas$reitsma.reg.fit <- reitsma(data = dat, method = "ml", formula = cbind(tsens , tfpr) ~ dat[["subgrouping.variable"]])
  reitsmas$reitsma.intercept <- reitsma(data = dat, method = "ml", formula = cbind(tsens , tfpr) ~ 1)
  reitsmas$anova.reitsma <- anova(reitsmas$reitsma.reg.fit, reitsmas$reitsma.intercept)
  reitsmas$subgroups <- list()
  summaries <- list()
  summaries$reitsma.overalll <- summary(reitsmas$reitsma.overall)
  summaries$reitsma.reg.fit <- summary(reitsmas$reitsma.reg.fit)
  summaries$subgroups <- list() 
  madads <- list()
  madads$overall <- madad(dat)
  madads$subgroups <- list()
  props <- list()
  props$sens <- list()
  props$spec <- list()
  props$sens$overall <-  metaprop(data = dat,
                                  studlab = dat[["names"]],
                                  event = dat[["TP"]],
                                  n = dat[["TP"]] + dat[["FN"]], method.tau =  "ml",
                                  sm = "PRAW",
                                  outclab = "sensitivity",
                                  method.ci = "WS",
                                  common = FALSE,
                                  subgroup = dat[["subgrouping.variable"]],
                                  subgroup.name = sglabel)
  props$sens$overall$upper <- madads$overall$sens$sens.ci[,2]
  props$sens$overall$lower <- madads$overall$sens$sens.ci[,1]
  props$sens$overall$TE.random <- summaries$reitsma.overall$coefficients[3,1]
  props$sens$overall$lower.random <- summaries$reitsma.overall$coefficients[3,5]
  props$sens$overall$upper.random <- summaries$reitsma.overall$coefficients[3,6]
  props$spec$overall <- metaprop(data = dat,
                                 studlab = dat[["names"]],
                                 event = dat[["TN"]],
                                 n = dat[["FP"]] + dat[["TN"]], method.tau =  "ml",
                                 sm = "PRAW",
                                 outclab = "specificity",
                                 method.ci = "WS",
                                 common = FALSE,
                                 subgroup = dat[["subgrouping.variable"]],
                                 subgroup.name = sglabel,
  )
  props$spec$overall$upper <- madads$overall$spec$spec.ci[,2]
  props$spec$overall$lower <- madads$overall$spec$spec.ci[,1]
  props$spec$overall$TE.random <- 1 - summaries$reitsma.overall$coefficients[4,1]
  props$spec$overall$upper.random <- 1 - summaries$reitsma.overall$coefficients[4,5]
  props$spec$overall$lower.random <- 1 - summaries$reitsma.overall$coefficients[4,6]
  hets.list <- list()
  for (sg in 1:length(subgroup.list)){
    reitsma.sg <- reitsma(data = dat[which(dat[["subgrouping.variable"]] == subgroup.list[sg]), ], method = 'ml')
    madad.sg <- madad(dat[which(dat[["subgrouping.variable"]] == subgroup.list[sg]), ])
    name <- valid.subgroup.list[sg]
    reitsmas$subgroups[[name]] <- reitsma.sg 
    madads$subgroups[[name]] <- madad.sg
    summary.sg <- summary(reitsma.sg)
    summaries$subgroups[[name]] <- summary.sg
    props$sens$overall$TE.random.w[subgroup.list[sg]] <- summary.sg$coefficients[3,1]
    props$sens$overall$upper.random.w[subgroup.list[sg]] <- summary.sg$coefficients[3,6]
    props$sens$overall$lower.random.w[subgroup.list[sg]] <- summary.sg$coefficients[3,5]
    props$spec$overall$TE.random.w[subgroup.list[sg]] <- 1 - summary.sg$coefficients[4,1]
    props$spec$overall$upper.random.w[subgroup.list[sg]] <- 1 - summary.sg$coefficients[4,5]
    props$spec$overall$lower.random.w[subgroup.list[sg]] <- 1 - summary.sg$coefficients[4,6]
    if (plot.het.subgroup) {
      hets.list[sg] <- het.string(reitsma.sg)
    }
  }
  if (plot.het.overall){
    overall.hetstring <- het.string(reitsmas$reitsma.overall)
  } else {
    overall.hetstring <- " "
  }
  if (!plot.het.subgroup){
    hets.list <- "Subgroup pooled effect"
  }
  if (is.null(lcols1)){
    lcols1.string <- NULL
  }else{
    lcols1.string <- "lcols1"
  }
  if (is.null(lcols2)){
    lcols2.string <- NULL
  }else{
    lcols2.string <- "lcols2"
  }
  if (sens.forest){
    meta::forest(props$sens$overall,
           xlim = forest.xlim,
           pscale = 100,
           just.addcols.right = "center",
           rightcols = c("effect", "ci"),
           rightlabs = c("Sensitivity %", "95% C.I."),
           leftcols = c("studlab", lcols1.string, lcols2.string),
           leftlabs = c("Study                                                   ", llab1 , llab2),
           xlab = "Sensitivity", smlab = "",
           weight.study = "random", squaresize = 0.7, col.square = "navy",
           col.square.lines = "navy",
           col.diamond = "maroon", col.diamond.lines = "maroon",
           pooled.totals = FALSE,
           comb.fixed = FALSE,
           fs.hetstat = 10,
           print.tau2 = FALSE,
           print.Q = FALSE,
           print.pval.Q = FALSE,
           print.I2 = FALSE,
           # overall = T,
           lty.random = 0,
           digits = 1,
           subgroup.hetstat = F,
           hetstat = FALSE,
           test.subgroup = FALSE,
           text.addline2 =  if (plot.overall) het.string(reitsmas$reitsma.overall) else " ",
           text.addline1 = paste("Between-group difference (p): ", round(reitsmas$anova.reitsma$statistic[3], digits = 3)),
           text.random.w = hets.list,
           text.random = "Random effects bivariate model",
           ref = if (plot.overall) 100 * props$sens$overall$TE.random else NA,
           calcwidth.random = calcwidth.shet.opt,
           just = "center",
           overall=plot.overall
    
    )
  }
  
  if (spec.forest){
    meta::forest(props$spec$overall,
           xlim = forest.xlim,
           pscale = 100,
           just.addcols.right = "center",
           rightcols = c("effect", "ci"),
           rightlabs = c("Specificity %", "95% C.I."),
           leftcols = c("studlab", lcols1.string, lcols2.string),
           leftlabs = c("Study                                                   ", llab1, llab2),
           xlab = "Specifcity", smlab = "",
           weight.study = "random", squaresize = 0.7, col.square = "navy",
           col.square.lines = "navy",
           col.diamond = "maroon", col.diamond.lines = "maroon",
           pooled.totals = FALSE,
           comb.fixed = FALSE,
           fs.hetstat = 10,
           print.tau2 = FALSE,
           print.Q = FALSE,
           print.pval.Q = FALSE,
           print.I2 = FALSE,
           # overall = T,
           lty.random = 0,
           digits = 1,
           subgroup.hetstat = F,
           hetstat = FALSE,
           test.subgroup = FALSE,
           text.addline2 =  if (plot.overall) het.string(reitsmas$reitsma.overall) else " ",
           text.addline1 = paste("Between-group difference: p = ", round(reitsmas$anova.reitsma$statistic[3], digits = 3)),
           text.random.w = hets.list,
           text.random = "Random effects bivariate model",
           ref = if (plot.overall) 100 * props$spec$overall$TE.random else NA,
           calcwidth.random = calcwidth.shet.opt,
           just = "center",
           overall = plot.overall
    )
  }



  # Loop over subgroups
for (sg in 1:length(subgroup.list)) {
    subgroup_name <- subgroup.list[sg]
    valid_name <- valid.subgroup.list[sg]
    
    # Extract sensitivity estimates and confidence intervals
    sens_estimate <- 100 * props$sens$overall$TE.random.w[subgroup_name]
    sens_lower <- 100 * props$sens$overall$lower.random.w[subgroup_name]
    sens_upper <- 100 * props$sens$overall$upper.random.w[subgroup_name]
    
    # Format to 2 decimal places
    sens_estimate_formatted <- sprintf("%.2f", sens_estimate)
    sens_lower_formatted <- sprintf("%.2f", sens_lower)
    sens_upper_formatted <- sprintf("%.2f", sens_upper)
    
    # Extract specificity estimates and confidence intervals
    spec_estimate <- 100 * props$spec$overall$TE.random.w[subgroup_name]
    spec_lower <- 100 * props$spec$overall$lower.random.w[subgroup_name]
    spec_upper <- 100 * props$spec$overall$upper.random.w[subgroup_name]
    
    # Format to 2 decimal places
    spec_estimate_formatted <- sprintf("%.2f", spec_estimate)
    spec_lower_formatted <- sprintf("%.2f", spec_lower)
    spec_upper_formatted <- sprintf("%.2f", spec_upper)
    
    # Print subgroup name
    cat("\nSubgroup:", subgroup_name, "\n")
    # Print sensitivity
    cat("  Sensitivity:", sens_estimate_formatted, "% (95% CI:", sens_lower_formatted, "-", sens_upper_formatted, "%)\n")
    # Print specificity
    cat("  Specificity:", spec_estimate_formatted, "% (95% CI:", spec_lower_formatted, "-", spec_upper_formatted, "%)\n")
}
  
  print(summary(reitsmas$reitsma.reg.fit))
          
  if (object.return){
    returned.object <- list()
    returned.object$reitsmas <- reitsmas
    returned.object$madads <- madads
    returned.object$metaprops <- props
    returned.object$summaries <- summaries
    returned.object$valid.subgroup.names <- valid.subgroup.list
    returned.object$subgroup.names <- subgroup.list
    return(returned.object)
  }
}




forest.diag.subgroup.combined <- function(dat,                # Dataframe with TP, TN, FP, FN, and a variable named "names"
                                          subgrouping.variable, # dat$subgrouping.variable
                                          sortvar = NULL,
                                          sglabel = "Subgroup", # String for subgroup label
                                          lcols1 = NULL,        # Variable to be shown in left side dat$lcols1
                                          llab1 = NULL,         # Label for lcols1
                                          lcols2 = NULL,        # Variable to be shown in left side dat$lcols2
                                          llab2 = NULL,         # Label for lcols2
                                          object.return = TRUE, # If TRUE, returns an object used in the function
                                          plot.het.overall = TRUE,
                                          plot.het.subgroup = TRUE,
                                          plot.overall = TRUE,
                                          calcwidth.shet.opt = FALSE,
                                          forest.xlim = c(50, 100),
                                          leftspace = "    ",
                                          rightspace = "        ",
                                          ratio = c(0.7, 0.3),
                                          ...) {
# Prepare data
  dat[["subgrouping.variable"]] <- subgrouping.variable
  dat[["lcols1"]] <- lcols1
  dat[["lcols2"]] <- lcols2
  if (is.null(sortvar)) {
    dat <- dat[order(dat[["subgrouping.variable"]]), ]
  } else {
    dat <- dat[order(dat[["subgrouping.variable"]], sortvar), ]
  }
  dat <- dat[complete.cases(dat[, c("TP", "TN", "FP", "FN")]), ]
  
  # Get list of subgroups
  subgroup.list <- unique(dat[["subgrouping.variable"]])
  counts <- table(dat[["subgrouping.variable"]])
  # Exclude subgroups with less than 3 studies
  subgroup.list <- names(counts[counts >= 3])
  subgroup.list <- sort(subgroup.list)
  valid.subgroup.list <- make.names(subgroup.list, unique = TRUE)
  dat <- subset(dat, dat[["subgrouping.variable"]] %in% subgroup.list)
  
  # Initialize lists to store results
  reitsmas <- list()
  reitsmas$subgroups <- list()
  summaries <- list()
  summaries$subgroups <- list()
  madads <- list()
  madads$subgroups <- list()
  props <- list()
  props$sens <- list()
  props$spec <- list()
  
  # Overall Reitsma model
  reitsmas$reitsma.overall <- reitsma(data = dat, method = "ml")
  summaries$reitsma.overall <- summary(reitsmas$reitsma.overall)
  madads$overall <- madad(dat)
  
  # Meta-analysis for sensitivity and specificity (overall)
  props$sens$overall <- metaprop(
    data = dat,
    studlab = dat[["names"]],
    event = dat[["TP"]],
    n = dat[["TP"]] + dat[["FN"]],
    method.tau = "ML",
    sm = "PRAW",
    outclab = "Sensitivity",
    method.ci = "WS",
    common = FALSE,
    subgroup = dat[["subgrouping.variable"]],
    subgroup.name = sglabel
  )

  props$spec$overall <- metaprop(
    data = dat,
    studlab = dat[["names"]],
    event = dat[["TN"]],
    n = dat[["TN"]] + dat[["FP"]],
    method.tau = "ML",
    sm = "PRAW",
    outclab = "Specificity",
    method.ci = "WS",
    common = FALSE,
    subgroup = dat[["subgrouping.variable"]],
    subgroup.name = sglabel
  )
  
  # Update confidence intervals from mada
  props$sens$overall$upper <- madads$overall$sens$sens.ci[, 2]
  props$sens$overall$lower <- madads$overall$sens$sens.ci[, 1]
  props$sens$overall$TE.random <- summaries$reitsma.overall$coefficients[3, 1]
  props$sens$overall$lower.random <- summaries$reitsma.overall$coefficients[3, 5]
  props$sens$overall$upper.random <- summaries$reitsma.overall$coefficients[3, 6]

  props$spec$overall$upper <- madads$overall$spec$spec.ci[, 2]
  props$spec$overall$lower <- madads$overall$spec$spec.ci[, 1]
  props$spec$overall$TE.random <- 1 - summaries$reitsma.overall$coefficients[4, 1]
  props$spec$overall$lower.random <- 1 - summaries$reitsma.overall$coefficients[4, 6]
  props$spec$overall$upper.random <- 1 - summaries$reitsma.overall$coefficients[4, 5]

  
  hets.list <- list()
  
  # Loop over subgroups
  for (sg in seq_along(subgroup.list)) {
    subgroup_name <- subgroup.list[sg]
    valid_name <- valid.subgroup.list[sg]
    subgroup_data <- dat[dat[["subgrouping.variable"]] == subgroup_name, ]
    
    # Reitsma model for subgroup
    reitsma.sg <- reitsma(data = subgroup_data, method = 'ml')
    reitsmas$subgroups[[valid_name]] <- reitsma.sg
    summaries$subgroups[[valid_name]] <- summary(reitsma.sg)
    madads$subgroups[[valid_name]] <- madad(subgroup_data)
    
    # Sensitivity estimates for subgroup
    props$sens$overall$TE.random.w[[subgroup_name]] <- summaries$subgroups[[valid_name]]$coefficients[3, 1]
    props$sens$overall$lower.random.w[[subgroup_name]] <- summaries$subgroups[[valid_name]]$coefficients[3, 5]
    props$sens$overall$upper.random.w[[subgroup_name]] <- summaries$subgroups[[valid_name]]$coefficients[3, 6]
    
    # Specificity estimates for subgroup
    props$spec$overall$TE.random.w[[subgroup_name]] <- 1 - summaries$subgroups[[valid_name]]$coefficients[4, 1]
    props$spec$overall$lower.random.w[[subgroup_name]] <- 1 - summaries$subgroups[[valid_name]]$coefficients[4, 6]
    props$spec$overall$upper.random.w[[subgroup_name]] <- 1 - summaries$subgroups[[valid_name]]$coefficients[4, 5]
    
    if (plot.het.subgroup) {
      hets.list[[subgroup_name]] <- het.string(reitsma.sg)
    }
  }

  # Prepare left column labels
  if (is.null(lcols1)) {
    lcols1.string <- NULL
  } else {
    lcols1.string <- "lcols1"
  }
  if (is.null(lcols2)) {
    lcols2.string <- NULL
  } else {
    lcols2.string <- "lcols2"
  }
  
  leftcols <- c("studlab", lcols1.string, lcols2.string)
  leftlabs <- c("Study", llab1, llab2)
  
  # Remove NULL values from leftcols and leftlabs
  leftcols <- leftcols[!sapply(leftcols, is.null)]
  leftlabs <- leftlabs[!sapply(leftlabs, is.null)]
  
  # Create dummy variables for spacing
  props$spec$overall$leftspace <- rep(leftspace, length(props$spec$overall$studlab))
  props$spec$overall$rightspace <- rep(rightspace, length(props$spec$overall$studlab))
  
  # Generate the sensitivity forest plot and capture it as a grob
  plot_sens <- function() {
    meta::forest(
      props$sens$overall,
      xlim = forest.xlim,
      pscale = 100,
      just.addcols.right = "center",
      rightcols = c("effect", "ci"),
      rightlabs = c("Sensitivity %", "95% C.I."),
      leftcols = leftcols,
      leftlabs = leftlabs,
      xlab = "Sensitivity",
      smlab = "",
      weight.study = "random",
      squaresize = 0.7,
      col.square = "navy",
      col.square.lines = "navy",
      col.diamond = "maroon",
      col.diamond.lines = "maroon",
      pooled.totals = FALSE,
      comb.fixed = FALSE,
      fs.hetstat = 10,
      print.tau2 = FALSE,
      print.Q = FALSE,
      print.pval.Q = FALSE,
      print.I2 = FALSE,
      digits = 1,
      subgroup.hetstat = FALSE,
      hetstat = FALSE,
      test.subgroup = FALSE,
      text.random.w = if (plot.het.subgroup) hets.list else "Subgroup pooled effect",
      text.random = if (plot.overall) "Random effects bivariate model" else "",
      ref = if (plot.overall) 100 * props$sens$overall$TE.random else NA,
      calcwidth.random = calcwidth.shet.opt,
      just = "center",
      overall = plot.overall
    )
  }
  grob.sens <- ggplotify::as.grob(plot_sens)
  
  # Generate the specificity forest plot and capture it as a grob
  props$spec$overall$subgroup.name <- ""
  props$spec$overall$subgroup.levels <- strrep(" ", seq_along(props$spec$overall$subgroup.levels))
  plot_spec <- function() {
    meta::forest(
      props$spec$overall,
      xlim = forest.xlim,
      pscale = 100,
      just.addcols.right = "center",
      rightcols = c("effect", "ci", "rightspace"),
      rightlabs = c("Specificity %", "95% C.I.", " "),
      leftcols = c("leftspace"),  # Use dummy variable for spacing
      leftlabs = c(" "),
      studlab = FALSE,
      xlab = "Specificity",
      smlab = "",
      weight.study = "random",
      squaresize = 0.7,
      col.square = "navy",
      col.square.lines = "navy",
      col.diamond = "maroon",
      col.diamond.lines = "maroon",
      pooled.totals = FALSE,
      comb.fixed = FALSE,
      fs.hetstat = 10,
      print.tau2 = FALSE,
      print.Q = FALSE,
      print.pval.Q = FALSE,
      print.I2 = FALSE,
      digits = 1,
      subgroup.hetstat = FALSE,
      hetstat = FALSE,
      test.subgroup = FALSE,
      text.random.w = "",
      text.random = "",
      ref = if (plot.overall) 100 * props$spec$overall$TE.random else NA,
      calcwidth.random = calcwidth.shet.opt,
      just = "center",
      overall = plot.overall
    )
  }
  grob.spec <- ggplotify::as.grob(plot_spec)
  
  # Adjust grob.spec to remove left padding
  grob.spec$widths[1] <- unit(0, "cm")
  
  # Combine the two grobs
  combined <- grid.arrange(grob.sens, grob.spec, ncol = 2, widths = ratio)
  
  # Print overall estimates
  sens_estimate <- 100 * props$sens$overall$TE.random
  sens_lower <- 100 * props$sens$overall$lower.random
  sens_upper <- 100 * props$sens$overall$upper.random
  
  spec_estimate <- 100 * props$spec$overall$TE.random
  spec_lower <- 100 * props$spec$overall$lower.random
  spec_upper <- 100 * props$spec$overall$upper.random
  
  cat("Overall Sensitivity:", sprintf("%.2f", sens_estimate), "% (95% CI:", sprintf("%.2f", sens_lower), "-", sprintf("%.2f", sens_upper), "%)\n")
  cat("Overall Specificity:", sprintf("%.2f", spec_estimate), "% (95% CI:", sprintf("%.2f", spec_lower), "-", sprintf("%.2f", spec_upper), "%)\n")
  
  # Loop over subgroups to print estimates
  for (sg in seq_along(subgroup.list)) {
    subgroup_name <- subgroup.list[sg]
    valid_name <- valid.subgroup.list[sg]
    
    # Extract sensitivity estimates and confidence intervals
    sens_estimate <- 100 * props$sens$overall$TE.random.w[[subgroup_name]]
    sens_lower <- 100 * props$sens$overall$lower.random.w[[subgroup_name]]
    sens_upper <- 100 * props$sens$overall$upper.random.w[[subgroup_name]]
    
    # Format to 2 decimal places
    sens_estimate_formatted <- sprintf("%.2f", sens_estimate)
    sens_lower_formatted <- sprintf("%.2f", sens_lower)
    sens_upper_formatted <- sprintf("%.2f", sens_upper)
    
    # Extract specificity estimates and confidence intervals
    spec_estimate <- 100 * props$spec$overall$TE.random.w[[subgroup_name]]
    spec_lower <- 100 * props$spec$overall$lower.random.w[[subgroup_name]]
    spec_upper <- 100 * props$spec$overall$upper.random.w[[subgroup_name]]
    
    # Format to 2 decimal places
    spec_estimate_formatted <- sprintf("%.2f", spec_estimate)
    spec_lower_formatted <- sprintf("%.2f", spec_lower)
    spec_upper_formatted <- sprintf("%.2f", spec_upper)
    
    # Print subgroup name
    cat("\nSubgroup:", subgroup_name, "\n")
    # Print sensitivity
    cat("  Sensitivity:", sens_estimate_formatted, "% (95% CI:", sens_lower_formatted, "-", sens_upper_formatted, "%)\n")
    # Print specificity
    cat("  Specificity:", spec_estimate_formatted, "% (95% CI:", spec_lower_formatted, "-", spec_upper_formatted, "%)\n")
  }
  
  # Print summary of the Reitsma regression
  print(summary(reitsmas$reitsma.overall))
  
  # Return objects if requested
  if (object.return) {
    returned.object <- list()
    returned.object$reitsmas <- reitsmas
    returned.object$madads <- madads
    returned.object$metaprops <- props
    returned.object$summaries <- summaries
    returned.object$valid.subgroup.names <- valid.subgroup.list
    returned.object$subgroup.names <- subgroup.list
    return(returned.object)
  }
}






AUC_boot_paralell <- function(TP, FP, FN, TN, B=2000, alpha=0.95)
  {
  N <- length(TP)
  p <- 2
  n1 <- TP + FN
  n2 <- TN + FP
  
  expit <- function(x) exp(x) / (1 + exp(x))
  
  dt1 <- data.frame(TP, FP, FN, TN)
  fit0 <- mada::reitsma(dt1)
  auc <- summary(fit0)$AUC$AUC
  
  mu1 <- as.numeric(fit0$coefficients)
  G1 <- fit0$Psi
  
  auc.pb <- foreach(b = 1:B, .combine = c, .packages = c("MASS", "mada")) %dopar% {
    
    t.pb <- expit(MASS::mvrnorm(N, mu1, G1))
    TPb <- rbinom(N, prob = t.pb[, 1], size = n1)
    FPb <- rbinom(N, prob = t.pb[, 2], size = n2)
    boot.data <- list()
    boot.data[["TP"]] <- TPb
    boot.data[["FP"]] <- FPb
    boot.data[["FN"]] <- n1 - TPb
    boot.data[["TN"]] <- n2 - FPb
    # Explicitly naming columns
    fit.pb <- mada::reitsma(boot.data)
    return(summary(fit.pb)$AUC$AUC)
  }
  
  
  Q1 <- quantile(auc.pb, c(.5 * (1 - alpha), 1 - .5 * (1 - alpha)))
  
  return(list(AUC = auc, CI = Q1))
} #this is a version of dmeta tools AUC_CI with parallel bootstrapping


multiple.srocs <- function(dat, # a dataset with TP, TN, FP, FN
                           subgrouping.variable = NULL, #dat$subgrouping variable (max is 6 subgroups but you can change it with making longer lists for the next 6 arguments)
                           sroc.colors = c("blue", "maroon", "black", "skyblue", "#20cb20", "red"), #colors for the sroc and summary estimates and ellipse 
                           points.colors = c("#0000FF20", "#A52A2A20","#1d1c1c20" ,"#87ceeb30", "#20cb2020", "#ff000350"),#colors for the point estimates
                           pch.list = c(16, 15, 17, 18, 15, 13), #shape for point estimate default is : c(1, 0, 2, 5, 7, 13)
                           summary.pch.list = c(16, 15, 17, 18, 15, 13), #shape for summary estimate
                           plot.ellipse = T,
                           plot.points = T,
                           plot.legend =T,
                           main.title = "SROC curves for subgroups",
                           object.return = T, #if true, it returns an object with all used objects
                           legend.AUC = T,
                           AUC.CI = F, # If true it also returns confidence interval for AUC. Be cautious since it takes ~ 2 minutes for each subgroup
                           n.boots = 500, # number of bootstraps for AUC CI
                           AUC.CI.object = NULL, # if you have already done auc ci calculation and have the results  as an object, place it here for faster implementation
                           magnify = 2 # argument to magnify weight size of point estimates in the sroc
){
  on.exit(eval(quote(closeAllConnections()), envir = .GlobalEnv))
  n.cores <- detectCores() - 1
  registerDoParallel(cores = n.cores)
  ellipse.lty <- 0
  if (plot.ellipse){
    ellipse.lty <- 2
  }
  if (is.null(subgrouping.variable)){
    dat[["subgrouping.variable"]] <- "Pooled"
  }else{
    dat[["subgrouping.variable"]] <- subgrouping.variable  
  }
  dat <- dat[which(!is.na(dat[["TP"]]) & !is.na(dat[["TN"]]) & !is.na(dat[["FP"]]) & !is.na(dat[["FN"]])), ]
  subgroup.list <- unique(dat[["subgrouping.variable"]])
  counts <- table(dat[["subgrouping.variable"]])
  subgroup.list <- names(counts[counts >= 3])
  subgroup.list <- sort(subgroup.list)
  valid.subgroup.list <- make.names(subgroup.list, unique = T)
  dat <- subset(dat, dat[["subgrouping.variable"]] %in% subgroup.list)
  reitsmas <- list()
  reitsmas$reitsma.overall <- reitsma(data = dat, method = "ml")
  reitsmas$subgroups <- list()
  summaries <- list()
  if (length(subgroup.list)>1){
    summaries$reitsma.overalll <- summary(reitsmas$reitssma.overall)
    summaries$reitsma.reg.fit <- summary (reitsmas$reitsma.reg.fit)
    reitsmas$reitsma.reg.fit <- reitsma(data = dat, method = "ml", formula = cbind(tsens , tfpr) ~ dat[["subgrouping.variable"]])
    reitsmas$reitsma.intercept <- reitsma(data = dat, method = "ml", formula = cbind(tsens , tfpr) ~ 1)
    reitsmas$anova.reitsma <- anova(reitsmas$reitsma.reg.fit, reitsmas$reitsma.intercept)
  }
  summaries$subgroups <- list()
  madauni.main <- madauni(dat, type = "DOR", method = "DSL")
  dat[["weights"]] <- magnify * madauni.main$weights
  datasets <- list()
  datasets$main <- dat
  datasets$subgroups <- list()
  for (sg in 1:length(subgroup.list)){
    reitsma.sg <- reitsma(data = dat[which(dat[["subgrouping.variable"]] == subgroup.list[sg]), ], method = 'ml')
    name <- valid.subgroup.list[sg]
    reitsmas$subgroups[[name]] <- reitsma.sg 
    summaries$subgroups[[name]] <- summary(reitsma.sg)
    datasets$subgroups[[name]] <- dat[which(dat[["subgrouping.variable"]] == subgroup.list[sg]), ]
  }
  if (length(subgroup.list)>1) {
    plot(reitsmas$subgroups[[valid.subgroup.list[1]]],
         extrapolate = TRUE,
         sroclwd = 2,
         predict = FALSE,
         pch = summary.pch.list[1],
         plotsumm = FALSE,
         main = main.title,
         sub = if (reitsmas$anova.reitsma$statistic[3]<0.01) "Difference between subgroups (bivariate model): p < 0.01" else paste("Difference between subgroups (bivariate model): p =" , round(reitsmas$anova.reitsma$statistic[3], digits = 3))
    )
  }else{
    plot(reitsmas$subgroups[[valid.subgroup.list[1]]],
         extrapolate = TRUE,
         sroclwd = 2,
         predict = FALSE,
         pch = summary.pch.list[1],
         plotsumm = FALSE,
         main = main.title,
    )
  }
  
  
  ROCellipse(reitsmas$subgroups[[valid.subgroup.list[1]]],
             lty = ellipse.lty,
             pch = summary.pch.list[1],
             add = TRUE,
             col = sroc.colors[1]
  )
  lines(sroc(reitsmas$subgroups[[valid.subgroup.list[1]]],
             extrapolate = TRUE
  ),
  lty = 1,
  col = sroc.colors[1],
  lwd = 2
  )
  if (plot.points){
    points(fpr(dat[which(dat[["subgrouping.variable"]] == subgroup.list[1]), ]),
           sens(dat[which(dat[["subgrouping.variable"]] == subgroup.list[1]), ]),
           pch = pch.list[1],
           cex = dat[which(dat[["subgrouping.variable"]] == subgroup.list[1]), ][["weights"]],
           col = points.colors[1]
    )
  }
  
  if (length(subgroup.list)>1){
    for (sg in 2:length(subgroup.list)){
      lines(sroc(reitsmas$subgroups[[valid.subgroup.list[sg]]],
                 extrapolate = TRUE
      ),
      lty = 1,
      col = sroc.colors[sg],
      lwd = 2
      )
      ROCellipse(reitsmas$subgroups[[valid.subgroup.list[sg]]],
                 lty = ellipse.lty,
                 pch = summary.pch.list[sg],
                 add = TRUE,
                 col = sroc.colors[sg],
                 cex = 1.5
      )
      if(plot.points){
        points(fpr(dat[which(dat[["subgrouping.variable"]] == subgroup.list[sg]), ]),
               sens(dat[which(dat[["subgrouping.variable"]] == subgroup.list[sg]), ]),
               pch = pch.list[sg],
               cex = dat[which(dat[["subgrouping.variable"]] == subgroup.list[sg]), ][["weights"]],
               col = points.colors[sg]
        )
      }
      
    }
  }
  
  
  
  
  
  
  no.subgroups <- length(subgroup.list)
  legend.pch <- summary.pch.list[1:no.subgroups]
  legend.col <- sroc.colors[1:no.subgroups]
  subgroups.auc.list <- list()
  AUC_CIs <- list()
  if (AUC.CI & is.null(AUC.CI.object)){
    registerDoParallel(cores = detectCores())
  }
  for (sg in 1:length(subgroup.list)){
    summary.sg <- summaries$subgroups[[valid.subgroup.list[sg]]]
    if(!AUC.CI){
      subgroups.auc.list[sg] <-  if (legend.AUC) paste(subgroup.list[sg],
                                      " (AUC: ",
                                      round(summary.sg$AUC$AUC, digits = 2),
                                      ")",
                                      sep = ""
      ) else subgroup.list[sg]
      
    } else {
      if (!is.null(AUC.CI.object)){
        AUC_CIs <- AUC.CI.object
      } else {
        sg.dat <- dat[which(dat[["subgrouping.variable"]] == subgroup.list[sg]), ]
        AUC_CIs[[valid.subgroup.list[sg]]] <- AUC_boot_paralell(TP = sg.dat[["TP"]],
                                                                FP = sg.dat[["FP"]],
                                                                FN = sg.dat[["FN"]],
                                                                TN = sg.dat[["TN"]],
                                                                B = n.boots)
      }
      AUC_CIs_sg <- AUC_CIs[[valid.subgroup.list[sg]]]
      subgroups.auc.list[sg] <- paste(subgroup.list[sg],
                                      " - AUC: ",
                                      round(AUC_CIs_sg$AUC, digits = 2),
                                      " [",
                                      round(AUC_CIs_sg$CI[1], digits = 2),
                                      " - ",
                                      round(AUC_CIs_sg$CI[2], digits = 2),
                                      "]",
                                      sep = ""
      )
    }
    
  }
  if (plot.legend){
    legend("bottomright",
           legend = subgroups.auc.list,
           pch = legend.pch,
           lty = 1,
           col = legend.col
           , lwd =1.5
    )
  }

  # Check if there are multiple subgroups
if (length(subgroup.list) > 1) {
    # Loop over each subgroup
    for (sg in 1:length(subgroup.list)) {
        subgroup_name <- subgroup.list[sg]
        valid_name <- valid.subgroup.list[sg]
        summary.sg <- summaries$subgroups[[valid_name]]
        
        # Extract AUC estimate
        auc_estimate <- summary.sg$AUC$AUC
        
        if (AUC.CI) {
            # Confidence intervals are stored in AUC_CIs
            auc_ci <- AUC_CIs[[valid_name]]$CI
            auc_lower <- auc_ci[1]
            auc_upper <- auc_ci[2]
            
            # Format to 2 decimal places
            auc_estimate_formatted <- sprintf("%.2f", auc_estimate)
            auc_lower_formatted <- sprintf("%.2f", auc_lower)
            auc_upper_formatted <- sprintf("%.2f", auc_upper)
            
            # Print AUC with confidence intervals
            cat("\nSubgroup:", subgroup_name, "\n")
            cat("  AUC:", auc_estimate_formatted, "(95% CI:", auc_lower_formatted, "-", auc_upper_formatted, ")\n")
        } else {
            # Format AUC estimate to 2 decimal places
            auc_estimate_formatted <- sprintf("%.2f", auc_estimate)
            
            # Print AUC without confidence intervals
            cat("\nSubgroup:", subgroup_name, "\n")
            cat("  AUC:", auc_estimate_formatted, "\n")
        }
    }
} else {
    # Only one subgroup or overall data
    summary.overall <- summaries$subgroups[[valid.subgroup.list[1]]]
    auc_estimate <- summary.overall$AUC$AUC
    
    if (AUC.CI) {
        # Confidence intervals are stored in AUC_CIs
        auc_ci <- AUC_CIs[[valid.subgroup.list[1]]]$CI
        auc_lower <- auc_ci[1]
        auc_upper <- auc_ci[2]
        
        # Format to 2 decimal places
        auc_estimate_formatted <- sprintf("%.2f", auc_estimate)
        auc_lower_formatted <- sprintf("%.2f", auc_lower)
        auc_upper_formatted <- sprintf("%.2f", auc_upper)
        
        # Print AUC with confidence intervals
        cat("AUC:", auc_estimate_formatted, "(95% CI:", auc_lower_formatted, "-", auc_upper_formatted, ")\n")
    } else {
        # Format AUC estimate to 2 decimal places
        auc_estimate_formatted <- sprintf("%.2f", auc_estimate)
        
        # Print AUC without confidence intervals
        cat("AUC:", auc_estimate_formatted, "\n")
    }
}       
  if (object.return){
    returned.object <- list()
    returned.object$reitsmas <- reitsmas
    returned.object$summaries <- summaries
    returned.object$madauni <- madauni.main
    if(length(subgroup.list)>1){
      returned.object$anova <- reitsmas$anova.reitsma
    }
    returned.object$datasets <- datasets
    returned.object$AUC_CIs <-AUC_CIs
    return(returned.object)
  }
  for (sg  in 1:length(subgroup.list)){
    print(paste(subgroup.list[sg], het.string(reitsmas$subgroups[[valid.subgroup.list[sg]]]), sep = " : "))
    
    if(plot.points){
      points(fpr(dat[which(dat[["subgrouping.variable"]] == subgroup.list[sg]), ]),
             sens(dat[which(dat[["subgrouping.variable"]] == subgroup.list[sg]), ]),
             pch = pch.list[sg],
             cex = dat[which(dat[["subgrouping.variable"]] == subgroup.list[sg]), ][["weights"]],
             col = points.colors[sg]
      )
    }
  }
}

dta.outliers <- function(dat, object.return = FALSE) {
  dat <- dat[which(!is.na(dat[["TP"]]) & !is.na(dat[["TN"]]) & 
                   !is.na(dat[["FP"]]) & !is.na(dat[["FN"]])), ]
  madauni.dat <- madauni(x = dat, type = "DOR", method = "DSL")
  metafor.object <- rma(yi = ln(madauni.dat$descr$DOR$DOR),
                        sei = madauni.dat$descr$DOR$se.lnDOR,
                        weights = madauni.dat$weights,
                        method = "DL")
  inf.object <- influence(metafor.object)
  plot(inf.object)
  
  outlier_indices <- which(abs(inf.object$inf$rstudent) > 2)
  
  for (std in outlier_indices) {
    print(paste("Model No. ", 
                dat[["Model No."]][std],
                '.   ',
                dat[["names"]][std],
                " is outlier",
                sep = ""))
  }
  
  if (object.return) {
    returned.object <- list()
    returned.object$madauni <- madauni.dat
    returned.object$inf <- inf.object
    returned.object$outlier_indices <- outlier_indices
    return(returned.object)
  }
}

forest.diag.no <- function(dat, combined = T, ...) {
  # Capture additional arguments
  args_list <- list(...)
  
  # First, run dta.outliers and capture the result
  outliers_result <- dta.outliers(dat, object.return = TRUE)
  
  # The dta.outliers function already prints outlier information
  # Now, extract the indices of outlier studies
  outlier_indices <- outliers_result$outlier_indices
  
  # Remove outlier studies from the dataframe
  if (length(outlier_indices) > 0) {
    dat_no_outliers <- dat[-outlier_indices, ]
  } else {
    dat_no_outliers <- dat
  }
  
  # Adjust any data-dependent arguments to match the modified dataframe
  adjust_args <- function(arg_value) {
    if (is.vector(arg_value) && length(arg_value) == nrow(dat)) {
      return(arg_value[-outlier_indices])
    } else {
      return(arg_value)
    }
  }
  
  # Apply the adjustment to all arguments
  args_list <- lapply(args_list, adjust_args)
  
  # Now, run forest.diag with the modified dataframe and adjusted arguments
  if (combined){
do.call(forest.diag.combined, c(list(dat = dat_no_outliers), args_list))
    } else {
  do.call(forest.diag, c(list(dat = dat_no_outliers), args_list))
    }
}

         

dta.outliers.multi <- function(dat,
                               subgrouping.variable,
                               object.return = FALSE) {
  dat[["subgrouping.variable"]] <- subgrouping.variable
  dat <- dat[which(!is.na(dat[["TP"]]) & !is.na(dat[["TN"]]) & 
                   !is.na(dat[["FP"]]) & !is.na(dat[["FN"]])), ]
  subgroup.list <- unique(dat[["subgrouping.variable"]])
  counts <- table(dat[["subgrouping.variable"]])
  subgroup.list <- names(counts[counts >= 3])
  subgroup.list <- sort(subgroup.list)
  valid.subgroup.list <- make.names(subgroup.list, unique = TRUE)
  dat <- subset(dat, dat[["subgrouping.variable"]] %in% subgroup.list)
  
  metafors <- list()
  madaunis <- list()
  infs <- list()
  outlier_indices <- list()  # To store outlier indices
  
  for (sg in 1:length(subgroup.list)) {
    datsg <- dat[which(dat[["subgrouping.variable"]] == subgroup.list[sg]), ]
    madauni.sg <- madauni(datsg, type = "DOR", method = "DSL")
    metafor.object.sg <- rma(yi = ln(madauni.sg$descr$DOR$DOR),
                             sei = madauni.sg$descr$DOR$se.lnDOR,
                             weights = madauni.sg$weights,
                             method = "DL")
    inf.object.sg <- influence(metafor.object.sg)
    name.sg <- valid.subgroup.list[sg]
    metafors[[name.sg]] <- metafor.object.sg
    madaunis[[name.sg]] <- madauni.sg
    infs[[name.sg]] <- inf.object.sg
    cat(paste0("*** In subgroup [", subgroup.list[sg], "] :\n"))
    
    # Identify outliers within the subgroup
    sg_outlier_indices <- which(abs(inf.object.sg$inf$rstudent) > 2)
    number_of_outliers <- length(sg_outlier_indices)
    global_outlier_indices <- which(dat[["subgrouping.variable"]] == subgroup.list[sg])[sg_outlier_indices]
    
    if (number_of_outliers > 0) {
      for (std in sg_outlier_indices) {
        cat(paste("Model No. ", 
                  datsg[["Model No."]][std],
                  '.   ',
                  datsg[["names"]][std],
                  " is outlier",
                  sep = ""), "\n")
      }
    }
    cat(paste("Number of Outliers within subgroup:", number_of_outliers, "\n"))
    
    # Store the global indices of outliers
    outlier_indices[[name.sg]] <- global_outlier_indices
  }
  
  if (object.return) {
    returned.object <- list()
    returned.object$metafors <- metafors
    returned.object$madaunis <- madaunis
    returned.object$infs <- infs
    returned.object$outlier_indices <- outlier_indices
    return(returned.object)
  }
}

forest.diag.subgroup.no <- function(dat, 
                                    subgrouping.variable,
                                    combined = T,
                                    ..., 
                                    only.subgroups.bigger.than.3 = TRUE) {
  # Capture additional arguments
  args_list <- list(...)
  
  # Run dta.outliers.multi and capture the result
  outliers_result <- dta.outliers.multi(dat, 
                                        subgrouping.variable = subgrouping.variable, 
                                        object.return = TRUE)
  
  # The dta.outliers.multi function already prints outlier information
  # Extract all outlier indices across subgroups
  outlier_indices <- unlist(outliers_result$outlier_indices)
  
  # Remove outlier studies from the dataframe
  if (length(outlier_indices) > 0) {
    dat_no_outliers <- dat[-outlier_indices, ]
    subgrouping.variable_no_outliers <- subgrouping.variable[-outlier_indices]
  } else {
    dat_no_outliers <- dat
    subgrouping.variable_no_outliers <- subgrouping.variable
  }
  
  # Adjust any data-dependent arguments to match the modified dataframe
  adjust_args <- function(arg_value) {
    if (is.vector(arg_value) && length(arg_value) == nrow(dat)) {
      return(arg_value[-outlier_indices])
    } else if (is.list(arg_value) && length(arg_value) == nrow(dat)) {
      return(arg_value[-outlier_indices])
    } else {
      return(arg_value)
    }
  }
  
  # Apply the adjustment to all arguments except subgrouping.variable
  args_list <- lapply(args_list, adjust_args)
  
  # Now, run forest.diag.subgroup with the modified dataframe and adjusted arguments
  if (combined) {
  do.call(forest.diag.subgroup.combined, c(list(dat = dat_no_outliers, 
                                       subgrouping.variable = subgrouping.variable_no_outliers, 
                                       only.subgroups.bigger.than.3 = only.subgroups.bigger.than.3), 
                                  args_list))
    } else {
  do.call(forest.diag.subgroup, c(list(dat = dat_no_outliers, 
                                       subgrouping.variable = subgrouping.variable_no_outliers, 
                                       only.subgroups.bigger.than.3 = only.subgroups.bigger.than.3), 
                                  args_list))
    }
}


multiple.srocs.no <- function(dat, 
                              subgrouping.variable = NULL, 
                              ..., 
                              object.return = TRUE,
                              AUC.CI.object = NULL) {
  # Capture additional arguments
  args_list <- list(...)
  
  # Add default values for data-dependent arguments if not provided
  if (is.null(args_list$sroc.colors)) {
    args_list$sroc.colors <- c("blue", "maroon", "black", "skyblue", "#20cb20", "red")
  }
  if (is.null(args_list$points.colors)) {
    args_list$points.colors <- c("#0000FF20", "#A52A2A20","#1d1c1c20" ,"#87ceeb30", "#20cb2020", "#ff000350")
  }
  if (is.null(args_list$pch.list)) {
    args_list$pch.list <- c(16, 15, 17, 18, 15, 13)
  }
  if (is.null(args_list$summary.pch.list)) {
    args_list$summary.pch.list <- c(16, 15, 17, 18, 15, 13)
  }
  
  # Determine the number of subgroups
  if (is.null(subgrouping.variable)) {
    subgrouping.variable <- rep("Pooled", nrow(dat))
  }
  dat[["subgrouping.variable"]] <- subgrouping.variable
  subgroup.list <- unique(subgrouping.variable)
  
  # Run appropriate outlier detection function
  if (length(subgroup.list) <= 1) {
    # Only one subgroup, use dta.outliers
    outliers_result <- dta.outliers(dat, object.return = TRUE)
    outlier_indices <- outliers_result$outlier_indices
  } else {
    # Multiple subgroups, use dta.outliers.multi
    outliers_result <- dta.outliers.multi(dat, 
                                          subgrouping.variable = subgrouping.variable, 
                                          object.return = TRUE)
    # Collect all outlier indices across subgroups
    outlier_indices <- unlist(outliers_result$outlier_indices)
  }
  
  # Remove outlier studies from the dataframe
  if (length(outlier_indices) > 0) {
    dat_no_outliers <- dat[-outlier_indices, ]
    subgrouping.variable_no_outliers <- subgrouping.variable[-outlier_indices]
  } else {
    dat_no_outliers <- dat
    subgrouping.variable_no_outliers <- subgrouping.variable
  }
  
  # Adjust any data-dependent arguments to match the modified dataframe
  adjust_args <- function(arg_value) {
    if (is.vector(arg_value) && length(arg_value) == nrow(dat)) {
      return(arg_value[-outlier_indices])
    } else if (is.list(arg_value) && length(arg_value) == nrow(dat)) {
      return(arg_value[-outlier_indices])
    } else {
      return(arg_value)
    }
  }
  
  # Apply the adjustment to all arguments except subgrouping.variable
  args_list <- lapply(args_list, adjust_args)
  
  # Now, run multiple.srocs with the modified dataframe and adjusted arguments
  do.call(multiple.srocs, c(list(dat = dat_no_outliers, 
                                 subgrouping.variable = subgrouping.variable_no_outliers, 
                                 object.return = object.return, 
                                 AUC.CI.object = AUC.CI.object), 
                            args_list))
}


         
PBS3 <- function(y,S,b0,V0){
  
  N <- dim(y)[1]
  p <- dim(y)[2]
  
  y.pb <- matrix(numeric(N*p),N)
  
  for(i in 1:N){
    
    yi <- y[i,]
    Si <- matrix(c(S[i,1],S[i,2],S[i,2],S[i,3]),p)
    
    Vi <- Si + V0
    
    Xi <- diag( sqrt(diag(Si) + diag(V0))^-1 )
    Psii <- Xi %*% Vi %*% t(Xi)
    
    mui <- Xi %*% b0
    
    y.pb[i,] <- MASS::ginv(Xi) %*% MASS::mvrnorm(1, mui, Psii)
    
  }
  
  return(y.pb)
  
}


MVPBT_boot <- function(y, S, B = 2000) {

  V0 <- mvmeta::mvmeta(y, S)$Psi
  Q0 <- MVPBT::MVPBT2(y, S)
  
  # Parallelize the loop
  results <- foreach(b=1:B, .packages=c("mvmeta", "MVPBT"), .export="PBS3", .combine=rbind) %dopar% {
    y.pb <- PBS3(y, S, Q0$b0, V0)
    mm.pb <- mvmeta::mvmeta(y.pb, S)
    Q.b <- MVPBT::MVPBT2(y.pb, S)
    T.b <- Q.b$T
    return(T.b)
  }
  
  T.b <- unlist(results)
  
  QT <- function(x, x0) {
    x1 <- sort(c(x, x0))
    w1 <- which(x1 == as.numeric(x0))
    qt <- 1 - w1/(length(x) + 1)
    return(qt)
  }
  
  P <- QT(T.b, Q0$T)
  R1 <- list(T.b = T.b, T = Q0$T, P = P)
  return(R1)
}


pubbias.diag <- function(dat, n.boots = 1000){
  on.exit(eval(quote(closeAllConnections()), envir = .GlobalEnv))
  n.cores <- detectCores() - 1
  registerDoParallel(cores = n.cores)
  dat <- dat[which(!is.na(dat[["TP"]]) & !is.na(dat[["TN"]]) & !is.na(dat[["FP"]]) & !is.na(dat[["FN"]])), ]
  dta.dat <- edta(TP = dat[["TP"]], FN = dat[["FN"]], TN = dat[["TN"]], FP = dat[["FP"]])
  oldpar <- par(mfrow=c(1,1))
  par(mfrow=c(1,2))
  attach(dta.dat)
  res1 <- rma(y[,1], S[,1])
  funnel(res1,main="(b) Funnel plot for logit(Se)")
  res2 <- rma(y[,2], S[,3])
  
  
  funnel(res2,main="(c) Funnel plot for logit(FPR)")
  registerDoParallel(cores = detectCores())
  MVPBT3.dat <- MVPBT_boot(y,S,B=n.boots)
  print (MVPBT3.dat[["P"]])
  detach(dta.dat)
  par(oldpar)
  return(MVPBT3.dat)
}



find_repeated_studies <- function(dat, unique_var = "Model No.", repeated_var = "Study No.") {
  
  # Identify Study Nos that have been repeated
  repeated_studies <- dat %>%
    group_by(!!sym(repeated_var)) %>%
    tally() %>%
    filter(n > 1) %>%
    pull(!!sym(repeated_var))
  
  # For each repeated Study No., get all associated Model Nos.
  if (length(repeated_studies) > 0) {
    for(study in repeated_studies) {
      associated_models <- dat %>%
        filter(!!sym(repeated_var) == study) %>%
        pull(!!sym(unique_var))
      
      cat("Study number", study, "has the following model numbers:\n")
      print(unique(associated_models))
      cat("\n")
    }
  } else {
    print("No repeated studies found.")
  }
  
  # Return the repeated study numbers
  return(repeated_studies)
}

find_repeated_studies_by_subgroup <- function(dat, unique_var = "Model No.", repeated_var = "Study No.", subgroup_var) {
  
  
  
  # Unique subgroups
  subgroups <- unique(dat[[subgroup_var]])
  
  for (subgroup in subgroups) {
    cat("For subgroup:", subgroup, "\n")
    cat("-------------------------------\n")
    
    subgroup_data <- dat %>% filter(!!sym(subgroup_var) == subgroup)
    
    # Identify Study Nos that have been repeated within the subgroup
    repeated_studies <- subgroup_data %>%
      group_by(!!sym(repeated_var)) %>%
      tally() %>%
      filter(n > 1) %>%
      pull(!!sym(repeated_var))
    
    # For each repeated Study No., get all associated Model Nos.
    if (length(repeated_studies) > 0) {
      for(study in repeated_studies) {
        associated_models <- subgroup_data %>%
          filter(!!sym(repeated_var) == study) %>%
          pull(!!sym(unique_var))
        
        cat("Study number", study, "has the following model numbers:\n")
        print(unique(associated_models))
        cat("\n")
      }
    } else {
      print(paste("No repeated studies found in subgroup", subgroup))
      cat("\n")
    }
    cat("-------------------------------\n\n")
  }
}

multiple.LRmats <- function(dat,
                            subgrouping.variable = NULL,
                            sum.colors = c("blue", "maroon", "black", "#20cb20", "skyblue", "red"),
                            points.colors = c("#0000FF20", "#A52A2A20","#1d1c1c20" ,"#20cb2020" ,"#87ceeb30", "#ff000350"),
                            inset_var = -.4
){
  if (is.null(subgrouping.variable)){
    dat[["subgrouping.variable"]] <- "Pooled"
  }else{
    dat[["subgrouping.variable"]] <- subgrouping.variable  
  }
  dat <- dat[which(!is.na(dat[["TP"]]) & !is.na(dat[["TN"]]) & !is.na(dat[["FP"]]) & !is.na(dat[["FN"]])), ]
  subgroup.list <- unique(dat[["subgrouping.variable"]])
  counts <- table(dat[["subgrouping.variable"]])
  subgroup.list <- names(counts[counts >= 3])
  subgroup.list <- sort(subgroup.list)
  valid.subgroup.list <- make.names(subgroup.list, unique = T)
  dat <- subset(dat, dat[["subgrouping.variable"]] %in% subgroup.list)
  reitsmas <- list()
  reitsmas$subgroups <- list()
  pLRs <- list()
  pLRs$subgroups <- list()
  nLRs <- list()
  nLRs$subgroups <- list()
  summaries <- list()
  summaries$subgroups <- list()
  datasets <- list()
  datasets$subgroups <- list()
  for (sg in 1:length(subgroup.list)){
    reitsma.sg <- reitsma(data = dat[which(dat[["subgrouping.variable"]] == subgroup.list[sg]), ], method = 'ml')
    name <- valid.subgroup.list[sg]
    reitsmas$subgroups[[name]] <- reitsma.sg 
    summaries$subgroups[[name]] <- summary(SummaryPts(reitsma.sg))
    pLRs$subgroups[[name]] <- madad(dat[which(dat[["subgrouping.variable"]] == subgroup.list[sg]), ])$posLR$posLR
    nLRs$subgroups[[name]] <- madad(dat[which(dat[["subgrouping.variable"]] == subgroup.list[sg]), ])$negLR$negLR
    datasets$subgroups[[name]] <- dat[which(dat[["subgrouping.variable"]] == subgroup.list[sg]), ]
  }
  
  par(mar = c(5, 4, 4, 12) + 0.1, xpd = TRUE)
  plot(nLRs$subgroups[[valid.subgroup.list[1]]], pLRs$subgroups[[valid.subgroup.list[1]]], 
       log = "xy", 
       type = "n",  # 'n' for no plotting of points or lines
       xlim = c(0.01, 1), 
       ylim = c(1,100),
         # if (max(madad(dat)$posLR$posLR)<=100) c(1, 100) else if (max(madad(dat)$posLR$posLR)<=200) c(1, 200) else max(madad(dat)$posLR$posLR),
       xlab = "Negative Likelihood Ratio", ylab = "Positive Likelihood Ratio",
       cex.lab = 1.2, cex.axis = 0.7, las = 1
       
  )
  
  usr <- par("usr")
  
  # Add horizontal pink dashed line at y=10 within plot region
  segments(x0 = 10^usr[1], x1 = 10^usr[2], y0 = 10, y1 = 10, col = "darkgray", lty = 2, lwd = 4)
  
  
  # Add vertical pink dashed line at x=0.1 within plot region
  segments(x0 = 0.1, x1 = 0.1, y0 = 10^usr[3], y1 = 10^usr[4], col = "darkgray", lty = 2, lwd = 4)
  
  legend_text <- c("LUQ: Exclusion & Confirmation\nLRP > 10, LRN < 0.1",
                   "RUQ: Confirmation Only\nLRP > 10, LRN > 0.1",
                   "LLQ: Exclusion Only\nLRP < 10, LRN < 0.1",
                   "RLQ: No Exclusion or Confirmation\nLRP < 10, LRN > 0.1",
                   "\n")
  # ,
  # "Summary LRP & LRN with 95% CI", 
  # "Study-specific LRP & LRN")
  var_pch = c(NA, NA, NA, NA, NA)
  var_col = c("black", "black", "black", "black", "black")
  # sum.colors[1], points.colors[1]
  var_pt.cex = c(1,1,1,1,1)
  
  for (sg in 1:length(subgroup.list)){
    
    # Plot the data points for subgroups on top of the gray background
    points(nLRs$subgroups[[valid.subgroup.list[sg]]], pLRs$subgroups[[valid.subgroup.list[sg]]], pch = 19, col = points.colors[sg], cex = 1.5)
    points(summaries$subgroups[[valid.subgroup.list[sg]]]["negLR", 1], summaries$subgroups[[valid.subgroup.list[sg]]]["posLR", 1], pch = 18, col = sum.colors[sg], cex = 3)
    
    arrows(summaries$subgroups[[valid.subgroup.list[sg]]]["negLR", 1], summaries$subgroups[[valid.subgroup.list[sg]]]["posLR", 3], summaries$subgroups[[valid.subgroup.list[sg]]]["negLR", 1], summaries$subgroups[[valid.subgroup.list[sg]]]["posLR", 4], length = 0.05, col = yarrr::transparent(sum.colors[sg], 0.6), angle = 90, code = 3, lwd = 2)
    
    
    arrows(summaries$subgroups[[valid.subgroup.list[sg]]]["negLR", 3], summaries$subgroups[[valid.subgroup.list[sg]]]["posLR", 1], summaries$subgroups[[valid.subgroup.list[sg]]]["negLR", 4], summaries$subgroups[[valid.subgroup.list[sg]]]["posLR", 1], length = 0.05, col = yarrr::transparent(sum.colors[sg], 0.6), angle = 90, code = 3, lwd = 2)
    var_pch <- c(var_pch, c(NA, 18,16))
    var_pt.cex <- c(var_pt.cex, c(1,2,1))
    
    if (length(subgroup.list) ==1){
      var_pch <- c(NA,NA,NA,NA,NA,18,16)
      var_pt.cex <- c(1,1,1,1,1,2,1)
      var_col <- c(var_col, c(sum.colors[sg], points.colors[sg]))
      legend_text <- c(legend_text, c("Summary LRP & LRN with 95% CI" , "Study-specific LRP & LRN"))
    }else{
      var_pch <- c(var_pch, c(NA, 18,16))
      var_pt.cex <- c(var_pt.cex, c(1, 2,1))
      var_col <- c(var_col, c( "black", sum.colors[sg], points.colors[sg]))
      sg_string <- paste("Subgroup:", subgroup.list[sg])
      legend_text <- c(legend_text, c(sg_string, "Summary LRP & LRN with 95% CI" , "Study-specific LRP & LRN"))
    }
    
  }
  legend("topright", inset = c(inset_var, 0), legend = legend_text, 
         pch = var_pch, col = var_col, 
         pt.cex = var_pt.cex, cex = 0.7, bty = "o", box.lwd = 1, box.col = "black", y.intersp = 1.2)
  par(mar = c(5, 4, 4, 2) + 0.1, xpd = FALSE)
  returned.obj <- list()
  returned.obj$reitsmas <- reitsmas
  returned.obj$pLRs <- pLRs
  returned.obj$nLRs <- nLRs
  returned.obj$summaries <- summaries
  returned.obj$datasets <- datasets
  return(returned.obj)
  # return(list(reitsmas, pLRs, nLRs, summaries, datasets))
}













nomogrammer <- function(Prevalence,
                        Sens = NULL,
                        Spec = NULL,
                        Plr = NULL,
                        Nlr = NULL,
                        Detail = T,
                        NullLine = T,
                        LabelSize = (12/5),
                        Verbose = FALSE,
                        x_var = .75,
                        y_var = 2){
  
  ## Function inputs:
  # Prevalence (prior probability) as a number between 0 and 1
  # Either
  # Sens & Spec
  # model sensitivity and specificity as a number between 0 and 1
  # Or
  # Likelihood ratios
  # Positive and Negative LRs (numeric)
  
  ## Function options:
  # Detail: If true, will overlay key statistics onto the plot
  # NullLine: If true, will add a line from prior prob through LR = 1
  # LabelSize: Tweak this number to change the label sizes
  # Verbose: Print out relevant metrics in the conso
  
  
  
  
  
  ## Helper functions
  ##   (defined inside nomogrammer, so remain local only & wont clutter user env)
  odds         <- function(p){
    # Function converts probability into odds
    o <- p/(1-p)
    return(o)
  }
  
  logodds      <- function(p){
    # Function returns logodds for a probability
    lo <- log10(p/(1-p))
    return(lo)
  }
  
  logodds_to_p <- function(lo){
    # Function goes from logodds back to a probability
    o <- 10^lo
    p <- o/(1+o)
    return(p)
  }
  
  p2percent <- function(p){
    # Function turns numeric probability into string percentage
    # e.g. 0.6346111 -> 63.5% 
    # scales::percent(signif(p, digits = 4))
    round(p*100, digits = 2)
    }
  
  
  ######################################
  ########## Calculations     ##########
  ######################################
  
  ## Checking inputs
  
  ## Prevalence
  # needs to exist
  if(missing(Prevalence)){
    stop("Prevalence is missing")
  }
  # needs to be numeric
  if(!is.numeric(Prevalence)){stop("Prevalence should be numeric")}
  # needs to be a prob not a percent
  if((Prevalence > 1) | (Prevalence <= 0)){stop("Prevalence should be a probability (did you give a %?)")}
  
  # Did user give sens & spec?
  if(missing(Sens) | missing(Spec)){
    sensspec <- FALSE
  } else{ sensspec <- TRUE}
  # if yes, make sure they are numbers
  if(sensspec == TRUE){
    if(!is.numeric(Sens)){stop("Sensitivity should be numeric")}
    if(!is.numeric(Spec)){stop("Specificity should be numeric")}
    # numbers that are probabilities not percentages
    if((Sens > 1) | (Sens <= 0)){stop("Sensitivity should be a probability (did you give a %?)")}
    if((Spec > 1) | (Spec <= 0)){stop("Specificity should be a probability (did you give a %?)")}
  }
  
  
  # Did user give PLR & NLR?
  if(missing(Plr) | missing(Nlr)){
    plrnlr <- FALSE
  } else{plrnlr <- TRUE}
  # if yes, make sure they are numbers
  if(plrnlr == TRUE){
    if(!is.numeric(Plr)){stop("PLR should be numeric")}
    if(!is.numeric(Nlr)){stop("NLR should be numeric")}
    # numbers that vaguely make sense
    if(Plr < 1){stop("PLR shouldn't be less than 1")}
    if(Nlr < 0){stop("NLR shouldn't be below zero")}
    if(Nlr > 1){stop("NLR shouldn't be more than 1")}
  }
  
  # Did they give a valid sensspec and plrnlr? If yes, ignore the LRs and tell them
  if((sensspec == TRUE) && (plrnlr == TRUE) ){
    warning("You provided sens/spec as well as likelihood ratios-- I ignored the LRs!")
  }
  
  
  ## If sens/spec provided, we calculate posterior probabilities & odds using sens & spec
  ##  otherwise, if plr and nlr provided, we calculate posteriors using them
  ##  if neither exist, then return an error
  if(sensspec == TRUE){
    prior_prob  <- Prevalence
    prior_odds  <- odds(prior_prob)
    sensitivity <- Sens
    specificity <- Spec
    PLR <- sensitivity/(1-specificity)
    NLR <- (1-sensitivity)/specificity
    post_odds_pos  <- prior_odds * PLR
    post_odds_neg  <- prior_odds * NLR
    post_prob_pos  <- post_odds_pos/(1+post_odds_pos)
    post_prob_neg  <- post_odds_neg/(1+post_odds_neg)
  } else if(plrnlr == TRUE){
    prior_prob  <- Prevalence
    prior_odds  <- odds(prior_prob)
    PLR <- Plr
    NLR <- Nlr
    sensitivity <- (PLR*(1-NLR))/(PLR-NLR)    ## TODO: check Adam's math! 
    specificity <- (1-PLR)/(NLR-PLR)          ## TODO: check Adam's math! 
    post_odds_pos  <- prior_odds * PLR
    post_odds_neg  <- prior_odds * NLR
    post_prob_pos  <- post_odds_pos/(1+post_odds_pos)
    post_prob_neg  <- post_odds_neg/(1+post_odds_neg)
  } else{
    stop("Couldn't find sens & spec, or positive & negative likelihood ratios")
  }
  
  
  
  ######################################
  ########## Plotting (prep)  ##########
  ######################################
  
  
  ## Set common theme preferences up front
  theme_set(theme_bw() +
              theme(axis.text.x = element_blank(),
                    axis.ticks.x = element_blank(),
                    axis.title.x = element_blank(),
                    axis.title.y = element_text(angle = 90, face = "bold"),
                    axis.title.y.right = element_text(angle = 90, face = "bold"),
                    axis.line = element_blank(),
                    panel.grid = element_blank(),
                    legend.position = "none",
                    panel.background = element_rect(fill = "gray90"),
                    plot.title = element_text(hjust = 0.5)
              )
  )
  
  ## Setting up the points of interest along the y-axes
  
  # Select probabilities of interest (nb as percentages)
  ticks_prob    <- c(0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 30,
                     40, 50, 60, 70, 80, 90, 95, 99, 99.5, 99.9)
  # Convert % to odds
  ticks_odds    <- odds(ticks_prob/100)
  # Convert % to logodds 
  ticks_logodds <- logodds(ticks_prob/100)
  
  # Select the likelihood ratios of interest (for the middle y-axis)
  ticks_lrs     <- sort(c(10^(-3:3), 2*(10^(-3:2)), 5*(10^(-3:2))))
  # Log10 them since plot is in logodds space
  ticks_log_lrs <- log10(ticks_lrs)
  
  
  
  
  ## Fixing particular x-coordinates
  left     <- 0
  right    <- 1
  middle   <- 0.5
  midright <- 0.75
  
  ## Lay out the four key plot points
  ##  (the start and finish of the positive and negative lines)
  
  # Initially these are expressed as probabilities
  df <- data.frame(x=c(left, right, left, right), 
                   y=c(prior_prob, post_prob_pos, prior_prob, post_prob_neg), 
                   line = c("pos", "pos", "neg", "neg"))
  
  adj_min      <- range(ticks_logodds)[1]
  adj_max      <- range(ticks_logodds)[2]
  adj_diff     <- adj_max - adj_min
  scale_factor <- abs(adj_min) - adj_diff/2
  #df$lo_y <- ifelse(df$x==left,(10/adj_diff)*logodds(1-df$y)-1,logodds(df$y))
  
  # Convert probabilities to logodds for plotting
  df$lo_y  <- ifelse(df$x==left,logodds(1-df$y)-scale_factor,logodds(df$y))
  # zero         <- data.frame(x = c(left,right),
  #                            y = c(0,0),
  #                            line = c('pos','pos'),
  #                            lo_y = c(-scale_factor,0))
  
  df$lo_y <- ifelse(is.na(df$lo_y)|is.infinite(df$lo_y)|is.null(df$lo_y), 0,df$lo_y)
  
  
  
  rescale   <- range(ticks_logodds) + abs(adj_min) - adj_diff/2
  rescale_x_breaks  <- ticks_logodds + abs(adj_min) - adj_diff/2  
  
  
  
  ######################################
  ########## Plot             ##########
  ######################################
  
  
  p <- ggplot(df) +
    geom_line(aes(x = x, y = lo_y, color = line), size = 1) +
    geom_vline(xintercept = middle) +
    annotate(geom = "text",
             x = rep(middle+.075, length(ticks_log_lrs)),
             y = (ticks_log_lrs-scale_factor)/2,
             label = ticks_lrs,
             size = rel(LabelSize)) +
    annotate(geom="point",
             x = rep(middle, length(ticks_log_lrs)),
             y = (ticks_log_lrs-scale_factor)/2,
             size = 1) +
    scale_x_continuous(expand = c(0,0)) + 
    scale_y_continuous(expand = c(0,0),
                       limits = rescale,
                       breaks = -rescale_x_breaks,
                       labels = ticks_prob,
                       name = "Pre-test probablity (%)",
                       sec.axis = sec_axis(trans = ~.,
                                           name = "Post-test probablity (%)",
                                           labels = ticks_prob,
                                           breaks = ticks_logodds))
  
  ## Optional overlay text: prevalence, PLR/NLR, and posterior probabilities
  detailedAnnotation <- paste(
    paste("prevalence = ", p2percent(prior_prob), "%", sep= ""),
    paste("PLR =", signif(PLR, 3),", NLR =", signif(NLR, 3)),
    paste("post. pos = ", p2percent(post_prob_pos),
          "% , neg = ", p2percent(post_prob_neg), "%", sep = ""),
    sep = "\n")
  
  
  ## Optional amendments to the plot
  
  ## Do we add the null line i.e. LR = 1, illustrating an uninformative model
  if(NullLine == TRUE){
    ## If yes, first calculate the start and end points
    uninformative <- data.frame(
      x = c(left,right),
      lo_y = c( (logodds(1-prior_prob) - scale_factor) , logodds(prior_prob))
    ) 
    
    p <- p + geom_line(aes(x = x, y = lo_y), data = uninformative,
                       color = "gray20", 
                       lty = 2,
                       inherit.aes = FALSE)
  }
  
  
  ## Do we add the detailed stats to the top right?
  if(Detail == TRUE){
    p <- p + annotate(geom = "text",
                      x = x_var,
                      y = y_var,
                      label = detailedAnnotation,
                      size = rel(LabelSize))
  }
  
  if(Verbose == TRUE){
    writeLines(
      text = c(
        paste0("prevalence = ", p2percent(prior_prob)),
        paste("PLR =", signif(PLR, 3)),
        paste("NLR =", signif(NLR, 3)),
        paste("posterior probability (positive) =", p2percent(post_prob_pos)),
        paste("posterior probability (negative) =", p2percent(post_prob_neg)),
        paste("sensitivity =", p2percent(sensitivity)),
        paste("specificity =", p2percent(specificity))
        # sep = "\n"
      )
    )
  }
  
  
  return(p)
  
}



nomogrammer_plus <- function(dat, prevalence, x_var = .75, y_var = 2, return.list = F, alphabet = T) {
  posLR <- summary(SummaryPts(reitsma(dat, method = "ml")))["posLR", 1]
  negLR <- summary(SummaryPts(reitsma(dat, method = "ml")))["negLR", 1]
  if (alphabet){
      alphabet <- c("A", "B", "C", "D", "E", "F")
    } else {
      alphabet <- c("", "", "", "", "", "")
    }
  
  plots <- list() # List to store ggplot objects
  for(PR in seq_along(prevalence)) {
    plot <- nomogrammer(Prevalence = prevalence[PR], 
                        Plr = posLR,
                        Nlr = negLR,
                        x_var = if (length(x_var) == 1) x_var else x_var[PR] ,
                        y_var = if(length(y_var) == 1) y_var else y_var[PR]) + 
      labs(title = paste(alphabet[PR],". Pre-test probablity", ": ", prevalence[PR]*100, "%", sep = ""))
    plots[[PR]] <- plot
   }
  if (return.list){
    return(plots)
   } else {
  # Combine plots using patchwork
  combined_plot <- wrap_plots(plots, ncol = length(prevalence), )
  return(combined_plot)
  print(combined_plot) # Print the combined plot
   }
  }
  


nomogrammer_cutoffs <- function(DM_object, cutoffs = DM_object$optcut, x_var = .75, y_var = 2, prevalence = 0.1, return.list = F, alphabet = T){
  LRs <- list()
  LRs$pos <- list()
  LRs$neg <- list()
  diagstat_object <- diagstats(DM_object, cutoff = cutoffs)
  if (alphabet){
      alphabet <- c("A", "B", "C", "D", "E", "F")
    } else {
      alphabet <- c("", "", "", "", "", "")
    }
  
  plots <- list()
  for (co in 1:length(cutoffs)){
    LRs$pos[co] <- diagstat_object[co, 2]/(1 - diagstat_object[co, 6])
    LRs$neg[co] <- (1 - diagstat_object[co, 2])/ diagstat_object[co, 6]
    plots[[co]] <- nomogrammer(Prevalence = prevalence,
                             Sens = diagstat_object[co, 2],
                             Spec = diagstat_object[co, 6],
                             x_var = if (length(x_var) == 1) x_var else x_var[co] ,
                             y_var = if(length(y_var) == 1) y_var else y_var[co]) + 
      labs(title = paste(alphabet[co],"Cutoff", ": ", cutoffs[co], sep = ""))
  }

  if (length(cutoffs) == 1){
    return(plots[[1]])
    print(plots[[1]])
  } else {
    if(return.list){
      return(plots)
    }
    else{
      combined_plot <- wrap_plots(plots, ncol = length(cutoffs), )
      return(combined_plot)
      print(combined_plot) # Print the combined plot
    }
  }
}


nomogrammer_cutoffs_DSO <- function(diagstat_object, x_var = .75, y_var = 2, prevalence = 0.1, return.list = F){
  LRs <- list()
  LRs$pos <- list()
  LRs$neg <- list()
  alphabet <- c("", "", "", "", "", "")
  plots <- list()
  for (co in 1:length(diagstat_object$cutoff)){
    LRs$pos[co] <- diagstat_object[co, 2]/(1 - diagstat_object[co, 6])
    LRs$neg[co] <- (1 - diagstat_object[co, 2])/ diagstat_object[co, 6]
    plots[[co]] <- nomogrammer(Prevalence = prevalence,
                               Sens = diagstat_object[co, 2],
                               Spec = diagstat_object[co, 6],
                               x_var = if (length(x_var) == 1) x_var else x_var[co] ,
                               y_var = if (length(y_var) == 1) y_var else y_var[co]) + 
      labs(title = paste(alphabet[co],"Cutoff", ": ", diagstat_object$cutoff[co], sep = ""))
  }
  
  if (length(diagstat_object$cutoff) == 1){
    return(plots[[1]])
    print(plots[[1]])
  } else {

    if(return.list){
      return(plots)
    }
    else{
      combined_plot <- wrap_plots(plots, ncol = length(diagstat_object$cutoff), )
      return(combined_plot)
      print(combined_plot) # Print the combined plot
    }
  }
}




nomogrammer_subgroups <- function(dat,
                                  subgrouping.variable,
                                  prevalence,
                                  x_var = .75,
                                  y_var = 2,
                                  return.list = F
                                  ){
dat[["subgrouping.variable"]] <- subgrouping.variable  
dat <- dat[which(!is.na(dat[["TP"]]) & !is.na(dat[["TN"]]) & !is.na(dat[["FP"]]) & !is.na(dat[["FN"]])), ]
subgroup.list <- unique(dat[["subgrouping.variable"]])
counts <- table(dat[["subgrouping.variable"]])
subgroup.list <- names(counts[counts >= 3])
subgroup.list <- sort(subgroup.list)
valid.subgroup.list <- make.names(subgroup.list, unique = T)
dat <- subset(dat, dat[["subgrouping.variable"]] %in% subgroup.list)
plots <- list()
for (sg in 1:length(subgroup.list)){
  sgdat <- dat[which(dat$subgrouping.variable == subgroup.list[sg]), ] 
  plr <- summary(SummaryPts(reitsma(sgdat, method = "ml")))["posLR", 1]
  nlr <- summary(SummaryPts(reitsma(sgdat, method = "ml")))["negLR", 1]
  plot <- nomogrammer(Prevalence = prevalence, 
                      Plr = plr,
                      Nlr = nlr,
                      x_var = if (length(x_var) == 1) x_var else x_var[sg] ,
                      y_var = if (length(y_var) == 1) y_var else y_var[sg]) + 
    labs(title = paste("Test: ",subgroup.list[sg], sep = ""))
  plots[[sg]] <- plot
 }
if(return.list){
 return(plots)
 } else {
  combined_plot <- wrap_plots(plots, ncol = length(subgroup.list), )
  return(combined_plot)
  print(combined_plot)
  }
}










forest.cutoffs <- function(dat,
                           cutoff,
                           lcols1 = cutoff,
                           lcols2 = NULL,
                           llab1 = "Cutoff",
                           llab2 = NULL,
                           sens.forest = T,
                           spec.forest = F
    ){
  dat[["cutoff"]] <- cutoff
  # dat <- dat[order(dat$cutoff), ]
  dat[["lcols1"]] <- lcols1
  dat[["lcols2"]] <- lcols2
  dat <- dat[which(!is.na(dat[["TP"]]) & !is.na(dat[["TN"]]) & !is.na(dat[["FP"]]) & !is.na(dat[["FN"]])), ]
  metaprop.sens <- metaprop(data = dat,
                            studlab = dat[["names"]],
                            event = dat[["TP"]],
                            n = dat[["TP"]] + dat[["FN"]], method.tau =  "ml",
                            sm = "PRAW",
                            outclab = "sensitivity",
                            method.ci = "WS",
                            common = FALSE,
                            )
  metaprop.spec <- metaprop(data = dat,
                            studlab = dat[["names"]],
                            event = dat[["TN"]],
                            n = dat[["TN"]] + dat[["FP"]], method.tau =  "ml",
                            sm = "PRAW",
                            outclab = "specificity",
                            method.ci = "WS",
                            common = FALSE)
  madad.dat <- madad(dat)
metaprop.sens$upper <- madad.dat$sens$sens.ci[,2]
metaprop.sens$lower <- madad.dat$sens$sens.ci[,1]
# metaprop.sens$lcols <- dat[["lcols"]]
metaprop.spec$upper <- madad.dat$spec$spec.ci[,2]
metaprop.spec$lower <- madad.dat$spec$spec.ci[,1]
# metaprop.spec$lcols <- dat[["lcols"]]



if (is.null(lcols1)){
  lcols1.string <- NULL
}else{
  lcols1.string <- "lcols1"
}
if (is.null(lcols2)){
  lcols2.string <- NULL
}else{
  lcols2.string <- "lcols2"
}
if (sens.forest){
  meta::forest(metaprop.sens,
               xlim = c(0,100),
               pscale = 100,
               just.addcols.right = "center",
               rightcols = c("effect", "ci"),
               rightlabs = c("Sensitivity %", "95% C.I.     "),
               leftcols = c("studlab", lcols1.string, lcols2.string),
               leftlabs = c("Study", llab1, llab2),
               xlab = "Sensitivity", smlab = "",
               weight.study = "random", squaresize = 0.7, col.square = "navy",
               col.square.lines = "navy",
               col.diamond = "maroon",
               col.diamond.lines = "maroon",
               pooled.totals = FALSE,
               comb.fixed = FALSE,
               fs.hetstat = 10,
               print.tau2 = FALSE,
               print.Q = FALSE,
               print.pval.Q = FALSE,
               print.I2 = FALSE,
               digits = 1,
               hetstat = FALSE,
               just = "center",
               random = F,
               sortvar = dat$cutoff
  )
}
if (spec.forest){
  meta::forest(metaprop.spec,
               xlim = c(70,100),
               pscale = 100,
               just.addcols.right = "center",
               rightcols = c("effect", "ci"),
               rightlabs = c("Specificty %", "95% C.I.     "),
               leftcols = c("studlab", lcols1.string, lcols2.string),
               leftlabs = c("Study", llab1, llab2),
               xlab = "Specifcity", smlab = "",
               weight.study = "random", squaresize = 0.7, col.square = "navy",
               col.square.lines = "navy",
               col.diamond = "maroon", col.diamond.lines = "maroon",
               pooled.totals = FALSE,
               comb.fixed = FALSE,
               fs.hetstat = 10,
               print.tau2 = FALSE,
               print.Q = FALSE,
               print.pval.Q = FALSE,
               print.I2 = FALSE,
               digits = 1,
               hetstat = FALSE,
               just = "center",
               random=FALSE,
               sortvar = dat$cutoff
  )
}
}



forest.diag.wrapper <- function(dat,
                                lcols1 =NULL,
                                llab1 = NULL,
                                lcols2 =NULL,
                                llab2 = NULL,
                                object.return = F,
                                sens.forest = T,
                                spec.forest = F,
                                plot.het = F,
                                calcwidth.addline.opt = T,
                                dso)  
{
  dat[["lcols1"]] <- lcols1
  dat[["lcols2"]] <- lcols2
  dat <- dat[which(!is.na(dat[["TP"]]) & !is.na(dat[["TN"]]) & !is.na(dat[["FP"]]) & !is.na(dat[["FN"]])), ]
  reitsma.fit <- reitsma(dat, method ='ml')
  summary.reitsma <- summary(reitsma.fit)
  metaprop.sens <- metaprop(data = dat,
                            studlab = dat[["names"]],
                            event = dat[["TP"]],
                            n = dat[["TP"]] + dat[["FN"]], method.tau =  "ml",
                            sm = "PRAW",
                            outclab = "sensitivity",
                            method.ci = "WS",
                            common = FALSE)
  metaprop.spec <- metaprop(data = dat,
                            studlab = dat[["names"]],
                            event = dat[["TN"]],
                            n = dat[["TN"]] + dat[["FP"]], method.tau =  "ml",
                            sm = "PRAW",
                            outclab = "specificity",
                            method.ci = "WS",
                            common = FALSE)
  
  
  madad.dat <- madad(dat)
  metaprop.sens$upper <- madad.dat$sens$sens.ci[,2]
  metaprop.sens$lower <- madad.dat$sens$sens.ci[,1]
  metaprop.sens$TE.random <- dso$Sens
  metaprop.sens$lower.random <- dso$lower.Sens
  metaprop.sens$upper.random <- dso$upper.Sens
  # metaprop.sens$lcols <- dat[["lcols"]]
  metaprop.spec$upper <- madad.dat$spec$spec.ci[,2]
  metaprop.spec$lower <- madad.dat$spec$spec.ci[,1]
  metaprop.spec$TE.random <- dso$Spec
  metaprop.spec$upper.random <- dso$upper.Spec
  metaprop.spec$lower.random <- dso$lower.Spec
  # metaprop.spec$lcols <- dat[["lcols"]]
  combined.set <- list()
  combined.set$reitsma <- reitsma.fit
  combined.set$summary.reitsma <- summary.reitsma
  combined.set$metaprop.sens <- metaprop.sens
  combined.set$metaprop.spec <- metaprop.spec
  if (plot.het){
    plot.string.het.overall <- het.string(reitsma.fit)
  }else{
    plot.string.het.overall <- " "
  }
  
  
  if (is.null(lcols1)){
    lcols1.string <- NULL
  }else{
    lcols1.string <- "lcols1"
  }
  if (is.null(lcols2)){
    lcols2.string <- NULL
  }else{
    lcols2.string <- "lcols2"
  }
  if (sens.forest){
    meta::forest(metaprop.sens,
                 xlim = c(0,100),
                 pscale = 100,
                 just.addcols.right = "center",
                 rightcols = c("effect", "ci"),
                 rightlabs = c("Sensitivity %", "95% C.I.     "),
                 leftcols = c("studlab", lcols1.string, lcols2.string),
                 leftlabs = c("Study", llab1, llab2),
                 xlab = "Sensitivity", smlab = "",
                 weight.study = "random", squaresize = 0.7, col.square = "navy",
                 col.square.lines = "navy",
                 col.diamond = "maroon", col.diamond.lines = "maroon",
                 pooled.totals = FALSE,
                 comb.fixed = FALSE,
                 fs.hetstat = 10,
                 print.tau2 = FALSE,
                 print.Q = FALSE,
                 print.pval.Q = FALSE,
                 print.I2 = FALSE,
                 digits = 1,
                 hetstat = FALSE,
                 text.addline1 = plot.string.het.overall,
                 ref = 100 * combined.set$metaprop.sens$TE.random,
                 calcwidth.addline = calcwidth.addline.opt,
                 just = "center"
    )
  }
  if (spec.forest){
    meta::forest(metaprop.spec,
                 xlim = c(50,100),
                 pscale = 100,
                 just.addcols.right = "center",
                 rightcols = c("effect", "ci"),
                 rightlabs = c("Specificty %", "95% C.I.     "),
                 leftcols = c("studlab", lcols1.string, lcols2.string),
                 leftlabs = c("Study", llab1, llab2),
                 xlab = "Specifcity", smlab = "",
                 weight.study = "random", squaresize = 0.7, col.square = "navy",
                 col.square.lines = "navy",
                 col.diamond = "maroon", col.diamond.lines = "maroon",
                 pooled.totals = FALSE,
                 comb.fixed = FALSE,
                 fs.hetstat = 10,
                 print.tau2 = FALSE,
                 print.Q = FALSE,
                 print.pval.Q = FALSE,
                 print.I2 = FALSE,
                 digits = 1,
                 hetstat = FALSE,
                 text.addline1 = plot.string.het.overall,
                 ref = 100 * combined.set$metaprop.spec$TE.random,
                 calcwidth.addline = calcwidth.addline.opt,
                 just = "center"
    )
  }
  if (object.return){
    return(combined.set)
  }
}










df.mixer <- function (training, intval, extval){
  training[["Cohort"]] <- "Training"
  intval[["Cohort"]] <- "Internal Validation"
  extval[["Cohort"]] <- "External Validation"
  overall <- rbind(training,intval,extval)
  return(overall)
}



