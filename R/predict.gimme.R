#' @name predict.gimme
#' @title GIMME Predicted Values.
#' @description This function calculates the predicted values of a fitted gimme
#'   model.
#' @usage predict.gimme(x)
#' @param x A fitted gimme object.
#' @return List of data frames. Each data frame contains the predicted values of
#'   a subject in the data.
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
#' predictions <- predict.gimme(fit)
#'  }
#' @keywords internal

predict.gimme <- function(x) {
  # Verify that x is a 'gimmep' object. If not, abort.
  if (!inherits(x, "gimmep")) {
    stop("Object ", x, " is not of class 'gimmep'. Aborted.")
  }
  
  # Extract number of subjects and variables from the 'gimmep' object.
  # number of subjects
  nsub <- length(x$data)
  # number of endogenous variables
  nend <- x$n_endog
  
  # Compute predicted values.
  predictions <- list()
  
  for (i in 1:nsub) {
    dataX <- as.matrix(x$data[[i]])
    coeff <- t(x$path_est_mats[[i]])[, 1:nend]
    predictions[[i]] <- as.data.frame(dataX %*% coeff)
  }
  rm(dataX, coeff, i)
  
  names(predictions) <- names(x$data)
  
  # Return predicted values in a list.
  return(predictions)
}
