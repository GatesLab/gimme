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

  table   <- lavaan::lavParTable(paths)
  
  # include paths in the syntax which are specified free by the user
  # allows the possibility for the user to fix certain paths to zero (or now specific values, done below)
  
  # check if any exogenous variables have been incorrectly specified
  # for free paths
  if(!is.null(varLabels$uexo)){
    for (exog in varLabels$uexo){
      if (exog %in% table[which(table$op == "~"),]$lhs){
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
      if (out %in% table[which(table$op == "~"),]$rhs){
        stop(paste0('gimme ERROR: an outcome variable was included as a predictor in 
                    user-specified paths.  Please remove variable from outcome list or 
                    correct path specification'))
      }
    }
  }
  
  
  if (any(is.na(table$ustart))){
    
    free <- which(is.na(table$ustart))
                  
    vsFree <- paste0(table[free,]$lhs, 
                     table[free,]$op, 
                     table[free,]$rhs)
    
  } else {
    
    vsFree <- NULL
    
  }
  
  # 7.16.22 kad: ALSO include paths in the syntax which are specified to a certain non-zero value by the user
  # Paths fixed to zero are separated below, as they will be removed, whereas these should be included in base syntax
  # 8.30.23 kmg: updated to include ~~ paths
  
  if (nrow(table[table$free == 0 & table$ustart != 0, ])>0) {
    
    specific <- which(table$ustart != 0)
    
    vsSpecific <- paste0(table[specific,]$lhs, 
                         table[specific,]$op, 
                         table[specific,]$ustart, 
                         "*", 
                        table[specific,]$rhs)
    

  
  # 7.16.22 kad: will also need to remove these paths from candidate paths
    
    vsSpecificRemove <- paste0(table[specific,]$lhs, table[specific,]$op, table[specific,]$rhs)
     
    if (any(table[specific,]$op == '~~')) {
      specificBi <- which(table[specific,]$op == '~~')
      vsSpecificRemoveBi <- paste0(table[specificBi,]$rhs, table[specificBi,]$op, table[specificBi,]$lhs)
      vsSpecificRemove <- c(vsSpecificRemove, vsSpecificRemoveBi)
    }
      
    
    
  } else {
    
    vsSpecific <- NULL
 
    vsSpecificRemove <- NULL
    
  }
  
  # 7.16.22 kad: table up the paths which are fixed to ZERO by the user (and thus will be removed)

  if (any(table$user == 1 & table$ustart == 0 & !is.na(table$ustart))){
    
    fixed        <- which(table$user == 1 & table$ustart == 0 & !is.na(table$ustart))
    
    vsFixed      <- paste0(table[fixed,]$lhs, table[fixed,]$op, table[fixed,]$rhs)
    
    if (any(table[fixed,]$op == '~~')) {
      fixedBi <- which(table[specific,]$op == '~~')
      vsFixedBi <- paste0(table[fixedBi,]$rhs, table[fixedBi,]$op, table[fixedBi,]$lhs)
      vsSpecificRemove <- c(vsSpecificRemove, vsFixedBi)
    }
    
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


