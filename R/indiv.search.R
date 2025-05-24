#' Individual-level search. Used in gimmeSEM, aggSEM, indSEM.
#' @param dat A list containing information created in setup().
#' @param grp A list containing group-level information. NULL in aggSEM and
#' indSEM.
#' @param ind A list containing individual- and (potentially) subgroup-level
#' information.
#' @param ind_cutoff Chi square cutoff, .05 level adjusted for multiple tests.
#' @param ind_z_cutoff  Z score cutoff, .05 level adjusted for multiple tests. 
#' @return Lists associated with coefficients, fit indices, etc.
#' @keywords internal 
indiv.search <- function(dat, grp, ind, ind_cutoff = NULL, ind_z_cutoff = 1.96){
  
  if (!dat$agg){
    ind$ind_paths   <-  vector("list", dat$n_subj)
    ind$n_ind_paths <- 0
  } else {
    ind <- NULL
  }
  
  if(!dat$hybrid){
    elig_paths   = dat$candidate_paths
  } else{
    elig_paths   = c(dat$candidate_paths, dat$candidate_corr)
  }
  
  status   <- list()
  fits     <- list()
  coefs    <- list()
  betas    <- list()
  vcov     <- list()
  vcovfull <- list()
  plots    <- list()
  syntax   <- list()
  psi      <- list()
  psiunstd <- list()
  plots_cov <- list()
  
  n_ind    <- ifelse(dat$agg, 1, dat$n_subj) 
  name     <- matrix(, nrow = n_ind, ncol = 1)
  
  
  if (dat$agg){
    data_all <- do.call(rbind, dat$ts_list)
    colnames(data_all) <- dat$varnames
  }
  
  for (k in 1:n_ind){
    
    if (dat$agg){
      data_list <- data_all
    } else {
      data_list <- dat$ts_list[[k]]
    }
    
    if (!dat$agg){
      writeLines(paste0("individual-level search, subject ", k, " (", names(dat$ts_list)[k],")"))
    }
    
    k_ind <- which(ind$names == names(dat$ts_list)[k])
    
    ind_spec <- search.paths(base_syntax  = dat$syntax, 
                             fixed_syntax = c(grp$group_paths, 
                                              ind$sub_paths[[k_ind]]),
                             add_syntax   = character(),
                             n_paths      = 0,
                             data_list    = data_list,
                             elig_paths   = elig_paths,
                             prop_cutoff  = NULL,
                             n_subj       = 1,
                             chisq_cutoff = ind_cutoff
            )
    
    temp_ind_spec <- ind_spec
    
    ind_spec <- prune.paths(base_syntax  = dat$syntax,
                            fixed_syntax = c(grp$group_paths, 
                                             ind$sub_paths[[k_ind]]),
                            add_syntax   = ind_spec[[1]][[1]]$add_syntax,
                            data_list    = data_list,
                            n_paths      = ind_spec[[1]][[1]]$n_paths,
                            n_subj       = 1,
                            prop_cutoff  = NULL,
                            elig_paths   = ind_spec[[1]][[1]]$add_syntax,
                            test_cutoff = ind_z_cutoff)
    
    if (!identical(temp_ind_spec[[1]][[1]]$add_syntax, ind_spec$add_syntax)){
      ind_spec <- search.paths(base_syntax  = dat$syntax, 
                               fixed_syntax = c(grp$group_paths, 
                                                ind$sub_paths[[k_ind]]),
                               add_syntax   = ind_spec$add_syntax,
                               n_paths      = ind_spec$n_paths,
                               data_list    = data_list,
                               elig_paths   = elig_paths,
                               prop_cutoff  = NULL,
                               n_subj       = 1,
                               chisq_cutoff = 0)
      ind$ind_paths[[k_ind]] <- ind_spec[[1]][[1]]$add_syntax
      ind$n_ind_paths[k_ind] <- ind_spec[[1]][[1]]$n_paths
      
    } else {
      
      ind$ind_paths[[k_ind]] <- ind_spec$add_syntax
      ind$n_ind_paths[k_ind] <- ind_spec$n_paths
      
    }
    
    
    s10 <- get.params(dat, 
                      grp, 
                      ind,
                      k)
    
    status[[k]] <- s10$status
    fits[[k]]   <- s10$ind_fit
    coefs[[k]]  <- s10$ind_coefs
    betas[[k]]  <- s10$ind_betas
    vcov[[k]]   <- s10$ind_vcov
    vcovfull[[k]]   <- s10$ind_vcov_full
    plots[[k]]  <- s10$ind_plot
    plots_cov[[k]] <- s10$ind_plot_cov
    syntax[[k]] <- c(dat$syntax,  grp$group_paths, ind$sub_paths[[k_ind]], ind$ind_paths[[k_ind]])
    psi[[k]]      <- s10$ind_psi
    psiunstd[[k]] <- s10$ind_psi_unstd
    name[k] <- names(dat$ts_list)[k]
  }
  
  if (dat$agg){
    names(status) <- names(fits) <- names(coefs) <- 
      names(betas) <- names(vcov) <- names(plots) <- names(psi) <- names(psiunstd) <-"all"
    # } else if (ind$n_ind_paths[k] > 0 & !dat$agg){
    #   names(status) <- names(fits) <- names(coefs) <- 
    #     names(betas) <- names(vcov) <- names(plots) <- names(dat$ts_list)
  } else {
    names(status) <- names(fits) <- names(coefs) <- 
      names(betas) <- names(vcov) <- names(vcovfull) <- names(plots) <- 
      names(psi) <- names(psiunstd) <- names(dat$ts_list)
  }
  
  res <- list("status" = status,
              "fits"   = fits,
              "coefs"  = coefs,
              "betas"  = betas,
              "psi"     = psi,
              "psiunstd" = psiunstd,
              "vcov"   = vcov,
              "vcovfull"   = vcovfull,
              "plots"  = plots,
              "plots_cov" = plots_cov,
              "ind"    = ind,
              "syntax" = syntax)
  return(res)
}