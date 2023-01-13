#' Returns MIs from lavaan fit object.
#' @param fit An object from lavaan.
#' @return If successful, returns MIs for an individual. If unsuccessful, 
#' returns NA.
#' @keywords internal 
return.mis <- function(fit, elig_paths){
  zero_se  <- FALSE
  no_paths <- FALSE
  error    <- inherits(fit, "try-error")
  if (!error){
    no_paths <- sum(lavInspect(fit, "free")$beta, na.rm = TRUE) == 0 
  }
  if (!error & !no_paths){
    zero_se <- sum(lavInspect(fit, "se")$beta, na.rm = TRUE) == 0
  } 
  
  if (!error & !zero_se){
    #commented out by lan 4.11.2019
    #mis   <- tryCatch(modindices(fit, op = "~", 
    #                             standardized = FALSE,
    #                             sort. = FALSE), 
    #                  error = function(e) e)
    # 
    lanMod <- function(fit, elig_paths){
      mis0 <- modindices(fit, standardized = FALSE, sort. = FALSE)
      mis0_idx <- paste0(mis0$lhs,mis0$op,mis0$rhs)
      mis <- mis0[mis0_idx %in% elig_paths,]
      return(mis)
    }
    
    mis   <- try(lanMod(fit, elig_paths))

    error <- inherits(mis, "try-error")
    if (error) mis <- NA 
  } else {
    mis <- NA
  }
  return(mis)
}
