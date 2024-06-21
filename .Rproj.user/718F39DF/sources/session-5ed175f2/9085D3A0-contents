########NOTES: ############
# test eigen() of each matrix separately 
# keep testing while this is not stable and >0 ind paths; remove remove etc. until it works
# Start from scratch and create outline 

#' Grabs final coefficients for each individual.
#' @param dat A list containing information created in setup().
#' @param grp A list containing group-level information. NULL in aggSEM and
#' indSEM.
#' @param ind A list containing individual- and (potentially) subgroup-level
#' information.
#' @param k The counter indicating the individual.
#' @return Individual-level information on fit, coefficients, and plots.
#' @keywords internal
get.params <- function(dat, grp, ind, k, ms.print = TRUE){
  
  op  = NULL # appease CRAN check
  ind_plot = NA
  ind_plot_psi = NA
  
  ###### FUNCTION FOR TESTING FOR STABILITY ######
  testWeights <- function(fit, dat){
    ind_betas <- round(lavInspect(fit, "std")$beta, digits = 4)
    #added to ensure correct ordering in matrices
    ind_betas <- ind_betas[dat$varLabels$endo,]
    ind_betas <- ind_betas[,dat$varLabels$coln]
    test      <- any(Re(eigen(ind_betas[,1:dat$n_endog])$values)>=1) | 
      any(Re(eigen(ind_betas[,(dat$n_endog+1):(dat$n_endog*2)])$values)>=1)
    return(test)
  }
  
  ###### RUN INITIAL FIT ######
  if (!dat$agg){
    fit <- fit.model(syntax    = c(dat$syntax, 
                                   grp$group_paths, 
                                   ind$sub_paths[[k]], 
                                   ind$ind_paths[[k]]), 
                     data_file = dat$ts_list[[k]])
  } else {
    data_file <- do.call("rbind", dat$ts_list)
    fit        <- fit.model(syntax    = c(dat$syntax, 
                                          grp$group_paths,
                                          ind$sub_paths[[k]], 
                                          ind$ind_paths[[k]]), 
                            data_file = data_file)
  }
  
  error   <- inherits(fit, "try-error")

    if (!error) {
    converge <- lavInspect(fit, "converged")
    zero_se  <- sum(lavInspect(fit, "se")$beta, na.rm = TRUE) == 0
  } else {
    converge <- FALSE
    zero_se  <- TRUE
  }
  
  ###### IF NON CONVERGE, ROLL BACK ######
  
  if (converge & !zero_se & !testWeights(fit, dat)){
    status   <- "converged normally" } else {
  
  # if no convergence or unstable, roll back one path at individual level, try again 
  if ((!converge | zero_se | testWeights(fit, dat)) & (length(ind$ind_paths[[k]])!= 0)){
    keepgoing = TRUE } else if ((!converge | zero_se | testWeights(fit, dat)) & (length(ind$ind_paths[[k]])== 0)) {
      if (converge & (!testWeights(fit, dat))) # if stable, testWeights(fit) = FALSE
        status <- "last known convergence"
      if (testWeights(fit, dat))
        status <- "unstable solution"
      if (!converge)
        status <- "nonconvergence"
      keepgoing = FALSE
    }
  
    while(keepgoing){
      ind$ind_paths[[k]] <- ind$ind_paths[[k]][-length(ind$ind_paths[[k]])]
      if (!dat$agg){
        fit <- fit.model(syntax    = c(dat$syntax, 
                                       grp$group_paths, 
                                       ind$sub_paths[[k]], 
                                       ind$ind_paths[[k]]), 
                         data_file = dat$ts_list[[k]])
      } else {
        data_file <- do.call("rbind", dat$ts_list)
        fit        <- fit.model(syntax    = c(dat$syntax, 
                                              grp$group_paths,
                                              ind$sub_paths[[k]], 
                                              ind$ind_paths[[k]]), 
                                data_file = data_file)
      }
      error   <- inherits(fit, "try-error")
      if (!error) {
        converge <- lavInspect(fit, "converged")
        zero_se  <- sum(lavInspect(fit, "se")$beta, na.rm = TRUE) == 0
      } else {
        converge <- FALSE
        zero_se  <- TRUE
      }
      if ((!converge | zero_se | testWeights(fit, dat)) & (length(ind$ind_paths[[k]])!= 0)){
        keepgoing = TRUE } else if (length(ind$ind_paths[[k]])== 0) {
          if (converge & (!testWeights(fit, dat))) # if stable, testWeights(fit) = FALSE
            status <- "last known convergence"
          if (testWeights(fit, dat))
            status <- "unstable solution"
          keepgoing = FALSE
        }
    }
}
  
  ###### COMPILE RESULTS IF CONVERGED ######
  
  if (converge & !zero_se){ #& (ind$n_ind_paths[k] >0) ){
    ind_fit    <- fitMeasures(fit, c("chisq", "df", "npar", "pvalue", "rmsea", 
                                     "srmr", "nnfi", "cfi", "bic", "aic", "logl"))
    ind_fit    <- round(ind_fit, digits = 4)
    ind_fit[2] <- round(ind_fit[2], digits = 0)
    
    r2         <- inspect(fit, "rsquare")
    r2         <- r2[dat$varLabels$endo]
    
    ind_fit    <- c(ind_fit, round(r2, digits = 4))
    
    ind_vcov_full <- lavInspect(fit, "vcov.std.all")
    keep          <- rownames(ind_vcov_full) %in% dat$candidate_paths
    ind_vcov      <- ind_vcov_full[keep, keep]
    
    ind_coefs_unst0 <- parameterEstimates(fit)
    ind_coefs_unst_idx <- paste0(ind_coefs_unst0$lhs,ind_coefs_unst0$op,ind_coefs_unst0$rhs)
    ind_coefs_unst <- ind_coefs_unst0[ind_coefs_unst0$op == "~" |
                                        ind_coefs_unst_idx %in% c(dat$candidate_paths, dat$candidate_corr),]
    
    ind_coefs0 <- standardizedSolution(fit)
    ind_coefs_idx <- paste0(ind_coefs0$lhs,ind_coefs0$op,ind_coefs0$rhs)
    ind_coefs <- ind_coefs0[ind_coefs0$op == "~" |
                              ind_coefs_idx %in% c(dat$candidate_paths, dat$candidate_corr),]
    
    ind_coefs <- cbind(ind_coefs[,1:3], ind_coefs_unst$est, ind_coefs[,4:9])
    colnames(ind_coefs) <- c("lhs", "op", "rhs", "est", "est.std", "se", "z", "pvalue", "ci.lower", "ci.upper")
    
    #ind_coefs <- subset(standardizedSolution(fit), op == "~")
    
    ind_betas <- round(lavInspect(fit, "std")$beta, digits = 4)
    ind_ses   <- round(lavInspect(fit, "se")$beta, digits = 4)
    
    #added to ensure correct ordering in matrices
    ind_betas <- ind_betas[dat$varLabels$endo,]
    ind_betas <- ind_betas[,dat$varLabels$coln]
    
    ind_ses <- ind_ses[dat$varLabels$endo,]
    ind_ses <- ind_ses[,dat$varLabels$coln]
    
    # zf added 2019-01-23
    ind_psi <- round(lavInspect(fit, "std")$psi, digits = 4)
    ind_psi_unstd <- round(lavInspect(fit, "estimates")$psi, digits = 4)
    
    ind_psi <- ind_psi[dat$varLabels$endo,]
    ind_psi <- ind_psi[,dat$varLabels$coln]
    ind_psi_unstd <- ind_psi_unstd[dat$varLabels$endo,]
    ind_psi_unstd <- ind_psi_unstd[,dat$varLabels$coln]
    #rownames(ind_betas) <- rownames(ind_ses) <- dat$varnames[(dat$n_lagged+1):(dat$n_vars_total)]
    #colnames(ind_betas) <- colnames(ind_ses) <- dat$varnames
    #   } # stl comment out 11.20.17 
    
    if(ms.print){
      if (dat$agg & !is.null(dat$out)){
        
        write.csv(ind_betas, file.path(dat$out, "allBetas.csv"), 
                  row.names = TRUE)
        
        # write.csv(ind_vcov_full, file.path(dat$out, "allvcov.csv"), 
        #           row.names = TRUE)
        
        write.csv(ind_ses, file.path(dat$out, "allStdErrors.csv"), 
                  row.names = TRUE)
        
        # zf added 2019-01-23
        write.csv(ind_psi, file.path(dat$out, "allPsi.csv"),row.names = TRUE)
        write.csv(ind_psi_unstd, file.path(dat$out, "allPsiUnstd.csv"),row.names = TRUE)
        
      } else if (!dat$agg & !is.null(dat$out)) { # & ind$n_ind_paths[k]>0)
        write.csv(ind_betas, file.path(dat$ind_dir, 
                                       paste0(dat$file_order[k,2], 
                                              "BetasStd.csv")), row.names = TRUE)
        
        # write.csv(ind_vcov_full, file.path(dat$ind_dir, 
        #                                paste0(dat$file_order[k,2], 
        #                                       "vcov.csv")), row.names = TRUE)
        # zf added 2019-01-23
        write.csv(ind_psi, file.path(dat$ind_dir, 
                                     paste0(dat$file_order[k,2], 
                                            "Psi.csv")), row.names = TRUE)
        write.csv(ind_psi_unstd, file.path(dat$ind_dir, 
                                           paste0(dat$file_order[k,2], 
                                                  "PsiUnstd.csv")), row.names = TRUE)
        write.csv(ind_ses, file.path(dat$ind_dir,
                                     paste0(dat$file_order[k,2], 
                                            "StdErrors.csv")), row.names = TRUE)
      }
    }
    
    if (dat$plot){
      ind_betas_t <- t(ind_betas)
      lagged      <- ind_betas_t[1:dat$n_lagged, ]
      contemp     <- ind_betas_t[(dat$n_lagged+1):(dat$n_vars_total), ]
      plot_vals   <- rbind(w2e(lagged), w2e(contemp))
      is_lagged   <- c(rep(TRUE, sum(lagged != 0)), 
                       rep(FALSE, sum(contemp != 0)))
      
      plot_file   <- ifelse(dat$agg, 
                            file.path(dat$out, "summaryPathsPlot.pdf"),
                            file.path(dat$ind_dir, 
                                      paste0(dat$file_order[k,2], "Plot.pdf")))
      
      ind_plot <- try(qgraph(plot_vals,
                             layout       = "circle",
                             lty          = ifelse(is_lagged, 2, 1),
                             edge.labels  = FALSE,
                             curve        = FALSE,
                             parallelEdge = TRUE,
                             fade         = FALSE,
                             posCol       = "red",
                             negCol       = "blue",
                             labels       = 
                               dat$varnames[(dat$n_lagged+1):(dat$n_vars_total)],
                             label.cex    = 2,
                             DoNotPlot    = TRUE))
      
      if (!is.null(dat$out) & !inherits(ind_plot, "try-error") & ms.print){
        pdf(plot_file)
        plot(ind_plot)
        dev.off()
      }
      if(dat$hybrid){
        covpsi     <- ind_psi[, (dat$n_lagged+1):(dat$n_vars_total)]
        covpsi[lower.tri(covpsi)] <- 0 # so we don't get duplicates
        diag(covpsi)    <- 0
        plot_vals_psi   <- w2e(covpsi)
        
        plot_file_psi   <- ifelse(dat$agg, 
                                  file.path(dat$out, "summaryCovPlot.pdf"),
                                  file.path(dat$ind_dir, 
                                            paste0(dat$file_order[k,2], "PlotCov.pdf")))
        
        ind_plot_psi <- try(qgraph(plot_vals_psi,
                                   layout       = "circle",
                                   lty          = 1,
                                   edge.labels  = FALSE,
                                   curve        = FALSE,
                                   parallelEdge = TRUE,
                                   fade         = FALSE,
                                   posCol       = "red",
                                   negCol       = "blue",
                                   arrows       = FALSE,
                                   labels       = 
                                     dat$varnames[(dat$n_lagged+1):(dat$n_vars_total)],
                                   label.cex    = 2,
                                   DoNotPlot    = TRUE))
        
        if (!is.null(dat$out) & !inherits(ind_plot_psi, "try-error") & ms.print){
          pdf(plot_file_psi)
          plot(ind_plot_psi)
          dev.off()
        }
        
      } else {
        ind_plot_psi <- NA
      }
      
    }
  } 
  
  ###### COMPILE RESULTS IF NOT CONVERGED ######
  
  if (!converge | zero_se){
    if (!converge) status <- "nonconvergence"
    if (sum(lavInspect(fit, "se")$beta, na.rm = TRUE) == 0) status <- "computationally singular"
    ind_fit   <- rep(NA, 11)
    ind_coefs <- matrix(NA, nrow = 1, ncol = 10)
    colnames(ind_coefs) <- c("lhs", "op", "rhs", "est", "est.std", "se", "z", "pvalue", "ci.lower", "ci.upper")
    ind_betas <- NA
    ind_vcov  <- NA
    ind_plot  <- NA
    ind_plot_psi <- NA
    ind_psi   <- NA
    ind_psi_unstd   <- NA
    ind_vcov_full <- NA
  }
  
  res <- list("status"    = status, 
              "ind_fit"   = ind_fit, 
              "ind_coefs" = ind_coefs, 
              "ind_betas" = ind_betas, 
              "ind_psi"   = ind_psi, 
              "ind_psi_unstd"  = ind_psi_unstd, 
              "ind_vcov"       = ind_vcov,
              "ind_vcov_full"  = ind_vcov_full,
              "ind_plot"       = ind_plot,
              "ind_plot_cov"   = ind_plot_psi,
              "ind_syntax" = c(dat$syntax, grp$group_paths,ind$sub_paths[[k]], ind$ind_paths[[k]])
  )
  return(res)
}
