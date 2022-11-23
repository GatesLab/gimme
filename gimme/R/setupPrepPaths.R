#' Allows user to open and close certain paths.
#' @param paths \code{lavaan}-style syntax containing paths with which
#' to begin model estimation (optional). That is, Y~X indicates that Y
#' is regressed on X, or X predicts Y. If no header is used,
#' then variables should be referred to with V followed (with no separation)
#' by the column number. If a header is used, variables should be referred to using 
#' variable names. To reference lag variables, "lag" should be added to the 
#' end of the variable name with no separation. Defaults to NULL.
#' @keywords internal
setupPrepPaths  <- function(paths, varLabels, ctrlOpts){
  
  # Satisfy CRAN checks
  varnames = NULL
  lvarnames = NULL
  
  table   <- lavaan::lavParTable(paths)
  
  # include paths in the syntax which are specified free by the user
  # allows the possibility for the user to fix certain paths to zero (or now specific values, done below)
  tableFree    <- table[table$op == "~" & table$free != 0, ]
  
  dvsFree      <- tableFree$lhs
  ivsFree      <- tableFree$rhs
  
  # check if any exogenous variables have been incorrectly specified
  # for free paths
  if(!is.null(varLabels$uexo)){
    for (exog in varLabels$uexo){
      if (exog %in% dvsFree){
        stop(paste0('gimme ERROR: an exogenous variable was treated as endogenous in 
                    specified paths.  Please remove variable from exogenous list or 
                    correct path specification'))
      }
    }
  }
  
  # check if any outcome variables have been incorrectly specified
  # for free paths
  if(!is.null(varLabels$outc)){
    for (out in varLabels$outc){
      if (out %in% ivsFree){
        stop(paste0('gimme ERROR: an outcome variable was included as a predictor in 
                    user-specified paths.  Please remove variable from outcome list or 
                    correct path specification'))
      }
    }
  }
  
  
  if (nrow(tableFree) != 0){
    
    vsFree <- paste0(dvsFree, "~", ivsFree)
    
  } else {
    
    vsFree <- NULL
    
  }
  
  # 7.16.22 kad: ALSO include paths in the syntax which are specified to a certain non-zero value by the user
  # Paths fixed to zero are separated below, as they will be removed, whereas these should be included in base syntax
  tableSpecific    <- table[table$op == "~" & table$free == 0 & table$ustart != 0, ]
  
  dvsSpecific      <- tableSpecific$lhs
  ivsSpecific      <- tableSpecific$rhs
  ustartSpecific   <- tableSpecific$ustart
  
  # check if any exogenous variables have been incorrectly specified
  # for specified value paths
  if(!is.null(varLabels$uexo)){
    for (exog in varLabels$uexo){
      if (exog %in% dvsSpecific){
        stop(paste0('gimme ERROR: an exogenous variable was treated as endogenous in 
                    specified paths.  Please remove variable from exogenous list or 
                    correct path specification'))
      }
    }
  }
  
  # check if any outcome variables have been incorrectly specified
  # for specified value paths
  if(!is.null(varLabels$outc)){
    for (out in varLabels$outc){
      if (out %in% ivsSpecific){
        stop(paste0('gimme ERROR: an outcome variable was treated as predictor in 
                    specified paths.  Please remove variable from outcome list or 
                    correct path specification'))
      }
    }
  }
  
  if (nrow(tableSpecific) != 0){
    
    vsSpecific <- paste0(dvsSpecific, "~", ustartSpecific, "*", ivsSpecific)
    
  } else {
    
    vsSpecific <- NULL
    
  }
  
  # 7.16.22 kad: will also need to remove these paths from candidate paths
  if (nrow(tableSpecific) != 0){
    
    vsSpecificRemove <- paste0(dvsSpecific, "~", ivsSpecific)
    
  } else {
    
    vsSpecificRemove <- NULL
    
  }
  
  # 7.16.22 kad: table up the paths which are fixed to ZERO by the user (and thus will be removed)
  tableFixed   <- table[table$op == "~" & table$free == 0 & table$ustart == 0,]
  
  if (nrow(tableFixed) > 0){
    
    dvsFixed     <- tableFixed$lhs
    
    # check if any exogenous variables have been incorrectly specified
    # for fixed to zero paths
    if(!is.null(varLabels$uexo)){
      for (exog in varLabels$uexo){
        if (exog %in% dvsFixed){
          stop(paste0('gimme ERROR: an exogenous variable was treated as endogenous in 
                      specified paths.  Please remove variable from exogenous list or 
                      correct path specification'))
        }
      }
    }
    
    # check if any outcome variables have been incorrectly specified
    # for fixed paths
    if(!is.null(varLabels$outc)){
      for (out in varLabels$outc){
        if (out %in% ivsFixed){
          stop(paste0('gimme ERROR: an outcome variable was included as a predictor in 
                    user-specified paths.  Please remove variable from outcome list or 
                    correct path specification'))
        }
      }
    }
    
    
    ivsFixed     <- recode.vars(tableFixed$rhs, varnames, lvarnames)
    vsFixed      <- paste0(dvsFixed, "~", ivsFixed)
    
  } else {
        
    vsFixed <- NULL
      
  }
  
  list = list(
    paths  = c(vsFree,vsSpecific), # 7.16.22 kad: now include free and specific (non-zero) value paths in base syntax
    remove = c(vsFixed,vsSpecificRemove), # 7.16.22 kad: also remove the specific value paths from candidate paths
    zero.paths = vsFixed # 8.13.22 kad: also return the paths set to 0 by user, so this info can be returned in output
  )
  
  return(list)
}


