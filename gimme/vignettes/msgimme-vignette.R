## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ------------------------------------------------------------------------
library(gimme)
data(ms.fit, package="gimme")

## ---- eval=FALSE---------------------------------------------------------
#  ms.fit <- gimme(data        = ms.sim,
#                  out         = '~/outputfolder/',
#                  ar          = FALSE,
#                  subgroup    = FALSE,
#                  standardize = TRUE,
#                  ms_allow    = TRUE,
#                  ms_tol      = 1e-5)

## ------------------------------------------------------------------------
solution.tree(ms.fit, level = "group", cols="stage")

## ------------------------------------------------------------------------
solution.tree(ms.fit, level = "group", plot.tree = TRUE)

## ------------------------------------------------------------------------
ms.fit$tables$summaryFit[1:10,]

