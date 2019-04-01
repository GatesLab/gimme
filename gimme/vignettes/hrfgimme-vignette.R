## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ------------------------------------------------------------------------
library(gimme)
data(HRFsim, package="gimme")

## ---- results='hide',message=FALSE---------------------------------------
?gimme::gimmeSEM 

## ----results='hide',message=FALSE----------------------------------------
HRF.fit <- gimme(data = data,
                 out = '~/outputfolder/',
                 ar = TRUE,
                 standardize = TRUE,
                 exogenous      = "V5",
                 conv_vars = "V5", 
                 conv_length = 16,
                 conv_interval = 1,
                 mult_vars      = "V3*V5",
                 mean_center_mult = TRUE, 
                 groupcutoff = .75)

## ---- fig.show='hold', results='hide', message=FALSE---------------------
plot(HRF.fit)

