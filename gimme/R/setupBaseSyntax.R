#' Set up base syntax file.
#' @keywords internal
setupBaseSyntax  <- function(paths, varLabels, ctrlOpts){
  
    #-------------------------------------------------------------#
    # NULL MODEL CONTAINS
    #-------------------------------------------------------------#
      # Cov & Var among variables than cannot be predicted.
      # Var of variables than can be predicted.
      # Means of exogenous variables
      # Intercepts of endogenous variables
      # Nonsense paths (fixed to zero)
      # Any fixed paths
    #-------------------------------------------------------------#
  
    # Var among variables than can be predicted.
    var.endo <- paste0(varLabels$endo, "~~", varLabels$endo)
  
    # Intercepts of endogenous variables
    int.endo  <- paste0(varLabels$endo, "~1")
  
    # Cov & Var among variables than cannot be predicted.
    cov.exog <- outer(varLabels$exog, varLabels$exog, function(x, y) paste0(x, "~~", y))
    cov.exog <- cov.exog[lower.tri(cov.exog, diag = TRUE)]
    
    # Means of exogenous variables
    mean.exog <- paste0(varLabels$exog, "~1")
    
    # Nonsense paths (fixed to zero)
    nons.reg <- outer(varLabels$exog, varLabels$endo, function(x, y) paste0(x, "~0*", y))
    nons.reg <- c(nons.reg[lower.tri(nons.reg, diag = FALSE)],nons.reg[upper.tri(nons.reg, diag = FALSE)])
    
    # Nonsense paths (not fixed to zero)
    nons.paths <- outer(varLabels$exog, varLabels$endo, function(x, y) paste0(x, "~", y))
    nons.paths <- c(nons.paths[lower.tri(nons.paths, diag = FALSE)],nons.paths[upper.tri(nons.paths, diag = FALSE)])
    
    # Any fixed paths
    fixed.paths <- paths
    
    if(ctrlOpts$ar) {
      ar.paths    <- paste0(
        setdiff(varLabels$orig, varLabels$uexo), 
        "~", paste0(setdiff(varLabels$orig, varLabels$uexo),"lag")
      )
      fixed.paths <- c(fixed.paths, ar.paths)
    }
    
    # All Possible Paths
    all.poss <- outer(varLabels$endo, c(varLabels$endo, varLabels$exog), function(x, y) paste0(x, "~", y))
    all.poss <- c(all.poss[lower.tri(all.poss, diag = FALSE)], all.poss[upper.tri(all.poss, diag = FALSE)])
    
    
    base.syntax <- c(
      var.endo,
      int.endo,
      cov.exog,
      mean.exog,
      nons.reg,
      fixed.paths
    )
    
    return(list(
      syntax = base.syntax,
      fixed.paths = fixed.paths,
      candidate.paths = setdiff(all.poss, fixed.paths),
      nonsense.paths =  nons.paths 
    ))
  
   
}


