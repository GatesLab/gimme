#' Identifies highest MI from list of MIs.
#' @param mi_list A list of MIs across individuals
#' @param indices A list of fit indices. Only relevant at the individual-level.
#' @param elig_paths A character vector containing eligible paths that
#' gimme is allowed to add to a model (e.g., no nonsense paths).
#' @param prop_cutoff The proportion of individuals for whom a path must be
#' significant in order for it to be added to the models. NULL if used 
#' at the individual-level.
#' @param n_subj The number of subjects in a given stage of the search. If
#' in the group stage, n_subj equals the number of subjects. If in the subgroup
#' stage, n_subj equals the number of individuals in a given subgroup. At the 
#' individual stage, n_subj = 1.
#' @param chisq_cutoff Cutoff used in order for MI to be considered significant.
#' Value varies depending on stage of search (e.g., group, subgroup, 
#' individual).
#' @return Returns name of parameter associated with highest MI. If no MI meets 
#' the criteria, returns NA.
#' @keywords internal 
highest.mi <- function(mi_list, 
                       indices,
                       elig_paths, 
                       prop_cutoff, 
                       n_subj,
                       chisq_cutoff,
                       allow.mult,
                       ms_tol){
  

  
  mi  = NULL # appease CRAN check
  sig = NULL # appease CRAN check
  param  = NULL # appease CRAN check
  
  mi_list       <- mi_list[!is.na(mi_list)]    
  n_converge    <- length(mi_list)
  mi_list       <- do.call("rbind", mi_list)
  
  mi_list$param <- paste0(mi_list$lhs, mi_list$op, mi_list$rhs)
  
  mi_list <- subset(mi_list, param %in% elig_paths, 
                    select = c("param", "mi", "epc"))
  
  mi_list$sig <- ifelse(mi_list$mi >= chisq_cutoff, 1, 0)
  
  mi_list <- transform(mi_list, 
                       sum = ave(mi, param, FUN = sum),
                       count = ave(sig, param, FUN = sum),
                       mean = ave(mi, param, FUN = mean))
  
  mi_list   <- subset(mi_list, !duplicated(param))
  mi_list   <- mi_list[order(-mi_list$count, -mi_list$sum), ]
  
  
  #------------------------------------------------------#
  # Group search ongoing...
  #------------------------------------------------------#
  if (!is.null(prop_cutoff)){
    
    if(allow.mult){
      
      # if there are good solutions
      if (mi_list$count[1] > (prop_cutoff*n_converge)){
        
        # we need to look at the means rather than the sum
        mi_list_ms <- mi_list[order(-mi_list$count, -mi_list$mean), ]
        red_mi     <- mi_list_ms[mi_list_ms$mean >= (mi_list_ms$mean[1]-ms_tol) & mi_list_ms$count == mi_list$count[1], , drop = FALSE]
        add_param  <- unlist(lapply(seq_along(1:nrow(red_mi)), function(i){
          red_mi$param[i]
        }))
        
        
      } else {
        
        add_param <- NA
        
      }
       
  
    } else {
      
      add_param <- ifelse(
        mi_list$count[1] > (prop_cutoff*n_converge),
        mi_list$param[1], 
        NA
      )
      
      if (n_converge <= (n_subj/2)) { 
        
        add_param <- NA
        
      }
      
    }

    
  #------------------------------------------------------#
  # Individual search ongoing...
  #------------------------------------------------------#
  } else {
    
    
    if(allow.mult){
      
      # we need to look at the means rather than the sum
      mi_list_ms <- mi_list[order(-mi_list$count, -mi_list$mean), ]
      red_mi     <- mi_list_ms[mi_list_ms$mean >= (mi_list_ms$mean[1]-ms_tol), , drop = FALSE]
      add_param  <- unlist(lapply(seq_along(1:nrow(red_mi)), function(i){
        red_mi$param[i]
      }))
        
    } else {
        
        add_param <- mi_list_ms$param[1L]
        
    }
       

    if (count.excellent(indices) >= 2) {
      
      add_param <- NA
      
    }
    
  }
  
  return(list(add_param=add_param,mi_list=mi_list))
}