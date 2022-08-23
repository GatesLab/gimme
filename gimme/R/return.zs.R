#' Returns z values from lavaan fit object.
#' @param fit An object from lavaan.
#' @param elig_paths eligable paths at this stage. For subgrouping, group and fixed paths. 
#' For pruning, only group paths. 
#' @return If successful, returns z values for an individual. If unsuccessful, 
#' returns NA.
#' @keywords internal 
return.zs <- function(fit, elig_paths){
  
  op  = NULL # appease CRAN check
  
  error   <- any(grepl("error", class(fit)))
  
  if (!error) {
    converge <- lavInspect(fit, "converged")
    zero_se  <- sum(lavInspect(fit, "se")$beta, na.rm = TRUE) == 0
  } else {
    converge <- FALSE
    zero_se <- TRUE
  }
  
if (!error & !zero_se & converge){

    zs0 <- tryCatch(standardizedSolution(fit),
                    error = function(e) e)
    zs0_idx <- paste0(zs0$lhs,zs0$op,zs0$rhs)
    zs <- zs0[zs0_idx %in% elig_paths,]
} else {
  zs <- NA
}

  return(zs)
}
