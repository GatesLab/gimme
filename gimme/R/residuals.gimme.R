#' @name residuals.gimme
#' @title GIMME Residuals.
#' @description This function calculates the unstandardized and standardized
#'   residuals of a fitted gimme model.
#' @usage residuals.gimme(x, lag)
#' @param x A fitted gimme object.
#' @param lag The number of lags tested in the Box-Pierce and Ljung-Box tests of the residuals.
#' If user does not specify a value, default is the smaller of 10 or the length of the time series divided by 5.
#' @return List of four lists of data frames: \describe{ \item{residuals}{List of
#'   the unstandardized residuals per subject.}
#'   \item{standardized.residuals}{List of the standardized residuals per
#'   subject.}
#'   \item{Box.Pierce.test}{List of the results of the Box-Pierce test for each subject's residuals.}
#'   \item{Ljung.Box.test}{List of the results of the Ljung-Box test for each subject's residuals.}}
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
#' residuals <- residuals.gimme(fit, lag = 5)
#'  }
#' @keywords internal

residuals.gimme <- function(x, lag = NULL) {
  # Verify that x is a 'gimmep' object. If not, abort.
  if (!inherits(x, "gimmep")) {
    stop("Object ", x, " is not of class 'gimmep'. Please input output object from running 'gimme'.")
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
  
  # Create function for running Box-Pierce and Ljung-Box tests
  residual.tests <- function(x, lag = NULL, test.name, AR){
    # If user does not specify a lag, do the smaller of 10 or the length of the time series divided by 5
    if(is.null(lag)){
      lag = min(10, nrow(x))
    }
    
    # Adjust df for test if AR has been specified 
    if(!AR){
      fitdf = 0
    }else{
      fitdf = 1
    }
    
    # Check if lag is df is greater than df (i.e. that df is not 0 or negative)
    if(lag <= fitdf){
      stop('Degrees of freedom for ', test.name, ' is 0 or negative. Specify a longer lag.')
    }
    
    # Set up output structure
    col.names <- c('X-squared','df','p-value')
    out <- matrix(NA, nrow = ncol(x), ncol = length(col.names))
    colnames(out) <- col.names; row.names(out) <- names(x) 
    
    # Run test
    for(i in 1:ncol(x)){
      test <- stats::Box.test(x[,i], lag = lag, type = test.name, fitdf = fitdf)
      out[i,'X-squared'] <- test$statistic
      out[i, 'df'] <- test$parameter
      out[i, 'p-value'] <- test$p.value
    }
    
    return(out)
  }
  
  # Specify whether AR is true for test function
  AR <- NULL
  if(nlag == 0){AR = FALSE}else{AR = TRUE}
  
  # Run Box-Pierce Test
  Box.Pierce.test <- lapply(res, function(x) residual.tests(x, lag = lag, test.name = "Box-Pierce", AR = AR))

  # Run Ljung-Box Test
  Ljung.Box.test <- lapply(res, function(x) residual.tests(x, lag = lag, test.name = "Ljung-Box", AR = AR))
  
  # Return residuals in a list along with test results
  out <- list(residuals = res, standardized.residuals = stand.res,
              Box.Pierce.test = Box.Pierce.test, Ljung.Box.test = Ljung.Box.test)
  return(out)
}
