# MIIVefa

A quick tutotial for using the MIIVefa package in R.

MIIVefa uses Model Implied Instrumental Variables (MIIVs) to perform Exploratory Factor Analysis (EFA). 

# The Basics

MIIVefa is data-driven algorithm for Exploratory Factor Analysis (EFA) that uses Model Implied Instrumental Variables (MIIVs). The method starts with a one factor model and arrives at a suggested model with enhanced interpretability that allows cross-loadings and correlated errors.

# Running MIIVefa

1, Prepare your data.

- The input data frame should be in a wide format: columns being different observations and rows being the specific data entries.

- Column names should be clearly labeled.

2, Installing MIIVefa.

- In the R console, enter and execute 'install.packages("MIIVefa")' or 'devtools::install_github("https://github.com/lluo0/MIIVefa")' after installing the "devtools" package.

- Load the MIIVefa by executing 'library(MIIVefa)' after installing. 

3, Running miivefa.

- The only necessarily required input is the raw data matrix.

- All 4 arguments are shown below.

- 'sigLevel' is the significance level with a default of 0.05. 'scalingCrit' is the specified criterion for selecting the scaling indicator whenever a new latent factor is created and the default is 'sargan+factorloading_R2.' And 'CorrelatedErrors' is a vector containing correlated error relations between observed variables with a default of NULL.


                    EFAmiiv <- function(data,
 
                    sigLevel = .05,
                    
                    scalingCrit = "sargan+factorloading_R2",
                    
                    correlatedErrors = NULL)
                    
# Output of MIIVefa

- The output of a miivefa object contains 2 parts:

- 1, a suggested model, of which the syntax is identical to a 'lavaan' model. Accessible via output$model.

- 2, a miivsem model fit of the suggested model. The suggested model is run and evaluated using 'MIIvsem' and all miivsem attributes can be accessed. Accessible via output$fit.
# Examples of MIIVefa

Please refer to the package vignette.
  

  
