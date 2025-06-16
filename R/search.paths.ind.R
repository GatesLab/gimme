#' Searches for paths. Ties together highest.mi and return.mis functions.
#' @param base_syntax A character vector containing syntax that never changes.
#' @param fixed_syntax A character vector containing syntax that does not change
#' in a given stage of searching.
#' @param data_list A list of datasets to be used in a given stage of the 
#' search. Varies based on group, subgroup, or individual-level stage.
#' @param n_subj The number of subjects in a given stage of the search. If
#' in the group stage, n_subj equals the number of subjects. If in the subgroup
#' stage, n_subj equals the number of individuals in a given subgroup. At the 
#' individual stage, n_subj = 1.
#' @param prop_cutoff The proportion of individuals for whom a path must be
#' nonsignificant in order for it to be dropped from the models. NULL if used 
#' at the individual-level.
#' @param elig_paths A character vector containing eligible paths that
#' gimme is allowed to add to the model at a given stage.
#' @param chisq_cutoff Cutoff used in order for MI to be considered significant.
#' Value varies depending on stage of search (e.g., group, subgroup, 
#' individual).
#' @param subgroup_stage Logical. Only present in order to instruct gimme
#' what message to print to console using writeLines.
#' @return Returns updated values of n_paths and add_syntax.
#' @keywords internal 
search.paths.ind <- function(data_list,
                            base_syntax, 
                            fixed_syntax,
                            elig_paths, 
                            prop_cutoff, 
                            n_subj, 
                            chisq_cutoff,
                            subgroup_stage,
                            hybrid,
                            dir_prop_cutoff,
                            ind_z_cutoff){
  
  obj <- replicate(1, 
                   list(
                     add_syntax     = character(),
                     n_paths        = 0, 
                     final.sol      = FALSE
                   ), simplify = FALSE
  )
  
  cnt <- 1
  history <- list()
  
  #-----------------------------------------------#
  # Determine the stage                           #
  #-----------------------------------------------#
    
    stage <- "individual search"
 
  #-----------------------------------------------#

        search    <- TRUE
        
        while(search){ # begin search
          
          mi_list <- list() 
          
          indices <- NULL

              
              # individual level search
              fit <- fit.model(
                syntax = c(base_syntax, fixed_syntax, obj[[1]]$add_syntax),
                data_file = data_list
              )
              
              #------------------------------------------------------#
              # Check to see if model converged.
              #------------------------------------------------------#
              
              if (!inherits(fit, "try-error")){
                # stl 2018/08/16 separated convergence check from error check
                # can't inspect convergence of an error object
                if (lavaan::lavInspect(fit, "converged") & !any(is.na(lavInspect(fit, what = "list")$se))){ 
                  indices    <- fitMeasures(fit, c("chisq", "df", "pvalue", "rmsea", 
                                                   "srmr", "nnfi", "cfi"))
                  converge <- lavInspect(fit, "converged")
                  zero_se  <- sum(lavInspect(fit, "se")$beta, na.rm = TRUE) == 0
                } else {
                  indices <- NULL
                  converge <-  lavInspect(fit, "converged")
                  zero_se  <- sum(lavInspect(fit, "se")$beta, na.rm = TRUE) == 0
                }
              } else {
                indices <- NULL
                converge <- FALSE
                zero_se  <- TRUE
              }
              mi_list[[1]] <- return.mis(fit, elig_paths)
              
            
    
          
          
          #------------------------------------------------------#
          # Add the parameters with the largest MI
          #------------------------------------------------------#
          if (!all(is.na(unlist(mi_list)))){
            if (ms_allow)
              add_p     <- highest.mi(mi_list      = mi_list,
                                      indices      = indices,
                                      elig_paths   = elig_paths,
                                      prop_cutoff  = prop_cutoff, 
                                      n_subj       = n_subj,
                                      chisq_cutoff = chisq_cutoff,
                                      allow.mult   = TRUE,
                                      ms_tol       = ms_tol,
                                      hybrid       = hybrid,
                                      dir_prop_cutoff = dir_prop_cutoff)
            if (!ms_allow)
              add_p     <- highest.mi(mi_list      = mi_list,
                                      indices      = indices,
                                      elig_paths   = elig_paths,
                                      prop_cutoff  = prop_cutoff, 
                                      n_subj       = n_subj,
                                      chisq_cutoff = chisq_cutoff,
                                      allow.mult   = FALSE,
                                      ms_tol       = ms_tol,
                                      hybrid       = hybrid, 
                                      dir_prop_cutoff = dir_prop_cutoff)
            
            add_param <- add_p$add_param
            mi_info   <- add_p$mi_list
            
          } else {
            
            add_param            <- NA
            mi_info              <- NA
            
          }
          #------------------------------------------------------#
          
          
          #------------------------------------------------------#
          # If there are no paths to add.
          #------------------------------------------------------#
          if(all(is.na(add_param))){
            
            search               <- FALSE
            obj[[1]]$final.sol   <- TRUE
            
            res      <- list()
            res[[1]] <- list(
              add_syntax     = obj[[1]]$add_syntax,
              n_paths        = obj[[1]]$n_paths,
              final.sol      = obj[[1]]$final.sol
            )
            
            #------------------------------------------------------#
            # If there is only 1 path to add, 
            #  still searching...
            #------------------------------------------------------#
          } else {
            
            
            obj[[1]]$n_paths     <- obj[[1]]$n_paths + 1
            obj[[1]]$add_syntax  <- append(obj[[1]]$add_syntax, add_param[1])
            
            #------------------------------------------------------#
            #  If there is > 1 path to add, stop searching.
            #------------------------------------------------------#  
          } 
          
        } # end search
      
        
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
        
        #### start get.params ####
        # prune is integrated here #

          ###### IF NON CONVERGE, ROLL BACK ######
          
          if (converge & !zero_se & !testWeights(fit, dat)){
            status1   <- "converged normally"
            keepgoing = FALSE} else {
              
              # if no convergence or unstable, roll back one path at individual level, try again 
              if ((!converge | zero_se | testWeights(fit, dat)) & (length(obj[[1]]$add_syntax)!= 0)){
                keepgoing = TRUE } else if ((!converge | zero_se | testWeights(fit, dat)) & (length(obj[[1]]$add_syntax)== 0)) {
                  if (converge & (!testWeights(fit, dat))) # if stable, testWeights(fit) = FALSE
                    status1 <- "last known convergence"
                  if (testWeights(fit, dat))
                    status1 <- "unstable solution"
                  if (!converge)
                    status1 <- "nonconvergence"
                  keepgoing = FALSE
                }
              
              while(keepgoing){
                obj[[1]]$add_syntax <- obj[[1]]$add_syntax[-length(obj[[1]]$add_syntax)]
                fit <- fit.model(syntax    = c(base_syntax, fixed_syntax, obj[[1]]$add_syntax), 
                                 data_file = dat$ts_list[[k]])
                error   <- inherits(fit, "try-error")
                if (!error) {
                  converge <- lavInspect(fit, "converged")
                  zero_se  <- sum(lavInspect(fit, "se")$beta, na.rm = TRUE) == 0
                } else {
                  converge <- FALSE
                  zero_se  <- TRUE
                }
                if ((!converge | zero_se | testWeights(fit, dat)) & (length(obj[[1]]$add_syntax)!= 0)){
                  keepgoing = TRUE } else if ((!converge | zero_se | testWeights(fit, dat)) & (length(obj[[1]]$add_syntax)== 0)) {
                    if (converge & (!testWeights(fit, dat))) # if stable, testWeights(fit) = FALSE
                      status1 <- "last known convergence"
                    if (testWeights(fit, dat))
                      status1 <- "unstable solution"
                    if (!converge)
                      status1 <- "nonconvergence"
                    keepgoing = FALSE
                  } else {
                    if (converge & (!testWeights(fit, dat))) # if stable, testWeights(fit) = FALSE
                      status1 <- "last known convergence"
 
                    keepgoing = FALSE
                    
                  }
              }
              res[[1]] <- list(
                add_syntax     = obj[[1]]$add_syntax,
                n_paths        = length(obj[[1]]$add_syntax),
                final.sol      = obj[[1]]$final.sol
              )
            } # end non converge roll back
          
          
          ##                                ##
          ########### PRUNE ##################
          ##                                ##
          
          prune <- TRUE
      
          
          if(converge){
            z_list <- list()
            while(prune){ 
              if (lavInspect(fit, "converged")){
                z_list <- list()
              z_list[[1]] <- return.zs(fit, elig_paths)
              if(!all(is.na(z_list))){
                drop_param <- lowest.z(z_list,
                                       elig_paths  =   obj[[1]]$add_syntax,
                                       prop_cutoff = prop_cutoff, 
                                       n_subj      = n_subj, 
                                       test_cutoff = ind_z_cutoff)
              } else {
                drop_param <- NA
              }
              
              if (!is.na(drop_param)){
                obj[[1]]$add_syntax <- obj[[1]]$add_syntax[!obj[[1]]$add_syntax %in% drop_param] 
                
                fit <- fit.model(
                  syntax = c(base_syntax, fixed_syntax, obj[[1]]$add_syntax),
                  data_file = data_list
                )
                
              } else {
                
                prune <- FALSE
                
              }
              } else {
                status1 <- "nonconvergence"
                prune <- FALSE
              }
              if (testWeights(fit, dat))
                status1 <- "unstable solution"
            }
          }
          
          ###### COMPILE RESULTS IF CONVERGED ######
        op  = NULL # appease CRAN check
        ind_plot = NA
        ind_plot_psi = NA
        
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
              
              if (!is.null(dat$out) & !inherits(ind_plot, "try-error")){
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
                
                if (!is.null(dat$out) & !inherits(ind_plot_psi, "try-error")){
                  pdf(plot_file_psi)
                  plot(ind_plot_psi)
                  dev.off()
                }
                
              } else {
                ind_plot_psi <- NA
              }
              
            }
          } # end compile results if converged
          
          ###### COMPILE RESULTS IF NOT CONVERGED ######
          
          if (!converge | zero_se){
            if (!converge) status1 <- "nonconvergence"
            #if (sum(lavInspect(fit, "se")$beta, na.rm = TRUE) == 0) status <- "computationally singular"
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
          
          ######## WRAP UP ###########
          
          syntax<- c(base_syntax, fixed_syntax, obj[[1]]$add_syntax)
          name <- names(dat$ts_list)[k]
       
          new.obj <- list(status1, 
                          ind_fit, 
                          ind_coefs,
                          ind_betas,
                          ind_vcov,
                          ind_plot,
                          ind_plot_psi, 
                          ind_psi, 
                          ind_psi_unstd, 
                          ind_vcov_full,
                          syntax, 
                          name)
          names(new.obj) <- c("status", "ind_fit", "ind_coefs", "ind_betas", "ind_vcov", "ind_plot", 
                              "ind_plot_psi", "ind_psi", "ind_psi_unstd", "ind_vcov_full", "syntax")
         ## end get.params if n_sub == 1
  
  return(new.obj)
  
}
