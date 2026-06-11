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
#' @param stop_crit Stopping criterion for the individual-level search.
#' "standard" stops when fit is adequate or no significant paths remain.
#' "model fit" continues adding the highest-MI path until fit is adequate,
#' even if the path is not significant. "significance" continues adding
#' significant paths even after fit is adequate.
#' @param n_excellent Number of fit indices needed to surpass their cutoffs for an
#' individual model to be considered excellent. Default is 2. Max is 4.
#' @param alpha Alpha level for significance testing. Default is .05.
#' @param correction Multiple comparison correction method. "Bonferroni" (default)
#' divides alpha by the number of candidate paths. "fdr" applies the
#' Benjamini-Hochberg false discovery rate procedure.
#' @inheritParams count.excellent
#' @return Returns name of parameter associated with highest MI. If no MI meets
#' the criteria, returns NA.
#' @keywords internal
highest.mi <- function(mi_list,
                       indices,
                       elig_paths,
                       prop_cutoff,
                       n_subj,
                       chisq_cutoff,
                       stop_crit = "standard",
                       allow.mult,
                       ms_tol,
                       hybrid,
                       dir_prop_cutoff,
                       rmsea_cutoff = .05,
                       srmr_cutoff = .05,
                       nnfi_cutoff = .95,
                       cfi_cutoff = .95,
                       n_excellent = 2,
                       alpha = .05,
                       correction = "Bonferroni"){
  
  mi  = NULL # appease CRAN check
  sig = NULL # appease CRAN check
  param  = NULL # appease CRAN check
  augmented_prohibited = NULL # appease CRAN check
  goodfit = NULL
  
  mi_list_na       <- mi_list[!is.na(mi_list)]    # retain list form for later use
  mi_list_na <- lapply(mi_list_na, function(x) {
    if (!"augmented_prohibited" %in% names(x)) {
      x$augmented_prohibited <- FALSE
    }
    x
  })
  n_converge    <- length(mi_list_na)
  mi_list       <- do.call("rbind", mi_list_na)
  
  mi_list$param <- paste0(mi_list$lhs, mi_list$op, mi_list$rhs)
  
  mi_list <- subset(mi_list, param %in% elig_paths, 
                    select = c("param", "mi", "epc", "augmented_prohibited"))

  # Augmented-variable constraints set prohibited MIs to zero. Exclude these
  # rows so an individual "model fit" search cannot select a forbidden path.
  mi_list <- subset(mi_list, !augmented_prohibited)
  if (nrow(mi_list) == 0) {
    if (n_subj == 1 && !is.null(indices)) {
      goodfit <- count.excellent(indices,
                                 rmsea_cutoff,
                                 srmr_cutoff,
                                 nnfi_cutoff,
                                 cfi_cutoff) >= n_excellent
    }
    return(list(add_param = NA, mi_list = mi_list, goodfit = goodfit))
  }
  
  # Assess per-subject significance of each MI.
  # For FDR, apply Benjamini-Hochberg across all individual MI p-values so that
  # count = ave(sig, param, FUN = sum) below correctly reflects FDR-adjusted
  # significance counts per path.
  if (identical(correction, "fdr") && nrow(mi_list) > 0) {
    p_vals    <- 1 - stats::pchisq(mi_list$mi, df = 1)
    m         <- length(p_vals)
    sorted_p  <- sort(p_vals)
    bh_thresh <- alpha * seq_len(m) / m
    k_star    <- max(c(0L, which(sorted_p <= bh_thresh)))
    bh_cutoff <- if (k_star > 0L) sorted_p[k_star] else -Inf
    mi_list$sig <- as.integer(p_vals <= bh_cutoff)
  } else {
    mi_list$sig <- ifelse(mi_list$mi >= chisq_cutoff, 1, 0)
  }

  mi_list <- transform(mi_list,
                       sum = ave(mi, param, FUN = sum),
                       count = ave(sig, param, FUN = sum),
                       mean = ave(mi, param, FUN = mean))

  mi_list   <- subset(mi_list, !duplicated(param))

  mi_list   <- mi_list[order(-mi_list$count, -mi_list$sum), ]

  # we need to look at the means rather than the sum
  mi_list_ms <- mi_list[order(-mi_list$count, -mi_list$mean), ]
  
  #------------------------------------------------------#
  # Group search ongoing...
  #------------------------------------------------------#
  if (!is.null(prop_cutoff)){
    if(allow.mult){
      # if there are good solutions
      if (mi_list$count[1] > (prop_cutoff*n_converge)){
        
        red_mi     <- mi_list_ms[mi_list_ms$mean >= (mi_list_ms$mean[1]-ms_tol) & mi_list_ms$count == mi_list$count[1], , drop = FALSE]
        add_param  <- unlist(lapply(seq_along(1:nrow(red_mi)), function(i){
          red_mi$param[i]
        }))
        
        
      } else {
        
        add_param <- NA
        
      }
      
    } else {
      whichone <- 1
      go <- 0
      while(go == 0) { # for directionality search implemented after convo w Peter, Sy-Miin, & Jonathon
        
        mi_list_ms <- mi_list[order(-mi_list$count, -mi_list$mean), ]
        red_mi     <- mi_list_ms[whichone,]
        
        # if the opposite direction is significant for > proportion, go indiv by indiv and ensure its preferred in one direction > groupcutoff in one direction
        red_lhs <- strsplit(mi_list_ms$param[whichone], "~")[[1]][1]
        red_rhs <- strsplit(mi_list_ms$param[whichone], "~")[[1]][2]
        opposite <- paste0(red_rhs, "~", red_lhs)
        count_red <- 0
        count_opp <- 0 
        if (length(which(mi_list_ms$param == opposite))>0 && dir_prop_cutoff>0 && hybrid==FALSE) {
          if (!grepl("lag", mi_list_ms$param[whichone]) && mi_list[which(mi_list_ms$param == opposite), 6] >= (prop_cutoff*n_converge))
          {
            for (p in 1:length(mi_list_na)){
              grab_red_mi  <- mi_list_na[[p]][which(mi_list_na[[p]]$lhs == red_lhs),]
              red_mi_value <- grab_red_mi[which(grab_red_mi$rhs == red_rhs),]$mi
              grab_opp_mi  <- mi_list_na[[p]][which(mi_list_na[[p]]$lhs == red_rhs),]
              opp_mi_value <- grab_opp_mi[which(grab_opp_mi$rhs == red_lhs),]$mi
              ifelse(red_mi_value>opp_mi_value, (count_red = count_red+1), (count_opp = count_opp+1))
            }
          }
          if (((count_opp-count_red)/length(mi_list_na)) > dir_prop_cutoff){ 
            add_param<- opposite
            go <- 1} else if (((count_red-count_opp)/length(mi_list_na)) >= dir_prop_cutoff){
              add_param<- red_mi$param
              go <- 1} else{
                go <- 0 # directionality could not be determined for the group level, try next MI
                whichone <- whichone + 1}
        } else {
          go <- 1
          add_param <- ifelse(
            mi_list$count[whichone] > (prop_cutoff*n_converge),
            mi_list$param[whichone], 
            NA)
        }
        
      }
    }
    
    if (n_converge <= (n_subj/2)) {
      add_param <- NA
    }
    #------------------------------------------------------#
    # Individual search ongoing...
    #------------------------------------------------------#
  } else {
    if (!stop_crit %in% c("standard", "model fit", "significance")) {
      stop("gimme ERROR: stop_crit must be one of 'standard', 'model fit', or 'significance'.")
    }
    
    if (n_subj == 1) {
      goodfit <- count.excellent(indices,
                                 rmsea_cutoff,
                                 srmr_cutoff,
                                 nnfi_cutoff,
                                 cfi_cutoff) >= n_excellent
    } else {
      goodfit <- NULL
    }
    
    
    if(allow.mult){
      
      # we need to look at the means rather than the sum
      mi_list_ms <- mi_list[order(-mi_list$count, -mi_list$mean), ]
      red_mi     <- mi_list_ms[mi_list_ms$mean >= (mi_list_ms$mean[1]-ms_tol), , drop = FALSE]
      add_param  <- unlist(lapply(seq_along(1:nrow(red_mi)), function(i){
        red_mi$param[i]
      }))
      
    } else {
      if (identical(stop_crit, "model fit") && (is.null(goodfit) || !goodfit)) {
        add_param <- mi_list_ms$param[1L]
      } else {
        add_param <- ifelse(mi_list_ms$sig[1L] == 1, mi_list_ms$param[1L], NA)
      }
      
    }
    
    
    if (n_subj == 1 && identical(stop_crit, "standard") && goodfit) {
      add_param <- NA
    }
    
    if (n_subj == 1 && identical(stop_crit, "model fit") && goodfit) {
      add_param <- NA
    }
    
  }
  
  return(list(add_param=add_param,mi_list=mi_list, goodfit = goodfit ))
}
