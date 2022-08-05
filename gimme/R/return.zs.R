#' Returns z values from lavaan fit object.
#' @param fit An object from lavaan.
#' @param hybrid logical indicating if hybrid gimme rules apply 
#' @return If successful, returns z values for an individual. If unsuccessful, 
#' returns NA.
#' @keywords internal 
return.zs <- function(fit, hybrid){
  
  op  = NULL # appease CRAN check
  
  error   <- any(grepl("error", class(fit)))
  
  if (!error) {
    converge <- lavInspect(fit, "converged")
    zero_se  <- sum(lavInspect(fit, "se")$beta, na.rm = TRUE) == 0
  } else {
    converge <- FALSE
    zero_se <- TRUE
  }
  
  if (!error & !zero_se & converge & !hybrid){
    zs <- tryCatch(subset(standardizedSolution(fit), 
                          op == "~" | op == "~~"),
                   error = function(e) e)
    error <- any(grepl("error", class(zs)))
    if (error) zs <- NA 
  } else if (!error & !zero_se & converge & hybrid){
    zs <- tryCatch(subset(standardizedSolution(fit), 
                          op == "~" | op == "~~"),
                   error = function(e) e)
    error <- any(grepl("error", class(zs)))
    if (error) zs <- NA else {
      zs <- NA} 
    } else {
    zs <- NA
  }
  return(zs)
}
