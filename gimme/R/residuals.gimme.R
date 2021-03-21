#' @name residuals.gimme
#' @title GIMME Residuals.
#' @description This function calculates the unstandardized and standardized
#'   residuals of a fitted gimme model.
#' @usage residuals.gimme(x)
#' @param x A fitted gimme object.
#' @return List of two lists of data frames. \describe{ \item{residuals}{List of
#'   the unstandardized residuals per subject.}
#'   \item{standardized.residuals}{List of the standardized residuals per
#'   subject.}}
#' @author Sebastian Castro-Alvarez
#' @examples
#'  \dontrun{
#' paths <- 'V2 ~ V1
#'           V3 ~ V4lag'
#'
#' fit <- gimmeSEM(data     = simData,
#'                 out      = "C:/simData_out",
#'                 subgroup = TRUE,
#'                 paths    = paths)
#'
#' residuals <- residuals.gimme(fit)
#'  }
#' @keywords internal

residuals.gimme <- function(x) {
  # Verify that x is a 'gimmep' object. If not, abort.
  if (!inherits(x, "gimmep")) {
    stop("Object ", x, " is not of class 'gimmep'. Aborted.")
  }
  
  # Extract number of subjects and variables from the 'gimmep' object.
  # number of subjects
  nsub <- length(x$data)
  # number of lag variables
  nlag <- x$n_lagged
  # number of endogenous variables
  nend <- x$n_endog
  
  # Compute predicted values.
  pred <- predict.gimme(x)
  
  # Compute unstandardized and standardized residuals.
  res       <- list()
  stand.res <- list() 
  
  for (i in 1:nsub) {
    observed   <- x$data[[i]][, (nlag + 1):(nlag + nend)]
    prediction <- pred[[i]]
    res[[i]]   <- observed - prediction
    
    stand.res[[i]] <- as.data.frame(apply(res[[i]], 2, scale))
  }
  rm(i, observed, prediction, pred)
  
  names(res) <- names(x$data)
  names(stand.res) <- names(x$data)
  
  # Return residuals in a list.
  out <- list(residuals = res, standardized.residuals = stand.res)
  return(out)
}
