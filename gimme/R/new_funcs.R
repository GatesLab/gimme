#' Recode variable names.
#' @param data The vector of variable names to be recoded
#' @param oldvalue A vector containing the latent variable names used internally.
#' @param newvalue A vector containing the observed variable names, either
#' provided by the user (as a header) or provided by R (e.g., V1, V2).
#' @return Recoded vector of variable names.
#' @keywords internal 
recode.vars <- function(data,
                        oldvalue,
                        newvalue){
  if (is.factor(data))     data     <- as.character(data)
  if (is.factor(oldvalue)) oldvalue <- as.character(oldvalue)
  if (is.factor(newvalue)) newvalue <- as.character(newvalue)
  newvec <- data
  for (i in unique(oldvalue)) newvec[data == i] <- newvalue[oldvalue == i]
  newvec
}

#' Create edge list from weight matrix.
#' @param x The coefficient matrix from an individual
#' @return A list of all non-zero edges to feed to qgraph
#' @keywords internal
w2e <- function(x) cbind(which(x != 0, arr.ind = TRUE), x[x != 0])


#' Attempt to fit lavaan model.
#' @param syntax A character vector containing syntax.
#' @param data_file A data frame containing individual data set.
#' @return If successful, returns fitted lavaan object. If not successful,
#' catches and returns error message.
#' @keywords internal 
fit.model <- function (syntax,
                       data_file) {
  
  fit <- tryCatch(lavaan::lavaan(syntax,
                                 data            = data_file,
                                 model.type      = "sem",
                                 missing         = "fiml",
                                 estimator       = "ml",
                                 int.ov.free     = FALSE,
                                 int.lv.free     = TRUE,
                                 auto.fix.first  = TRUE,
                                 auto.var        = TRUE,
                                 auto.cov.lv.x   = TRUE,
                                 auto.th         = TRUE,
                                 auto.delta      = TRUE,
                                 auto.cov.y      = FALSE,
                                 auto.fix.single = TRUE,
                                 warn            = FALSE),
                  error = function(e) e)
  return(fit)
}

#' Counts number of excellent fit indices
#' @param indices A vector of fit indices from lavaan.
#' @return The number of fit indices that are excellent.
#' @keywords internal 
count.excellent <- function(indices){
  rmseaE    <- ifelse(indices[4] < .05, 1, 0)
  srmrE     <- ifelse(indices[5] < .05, 1, 0)
  nnfiE     <- ifelse(indices[6] > .95, 1, 0)
  cfiE      <- ifelse(indices[7] > .95, 1, 0)
  excellent <- sum(rmseaE, srmrE, nnfiE, cfiE, na.rm = TRUE)
  return(excellent)
}

#' Provides unique combinations of two vectors.
#' @param x A character vector containing variable names.
#' @param y A character vector containing variable names.
#' @param incl.eq Logical. TRUE means that combinations are kept where
#' a variable appears twice.
#' @return The unique combinations of the variable names. Used in syntax
#' creation.
#' @keywords internal 
expand.grid.unique <- function(x, y, incl.eq = TRUE){
  g <- function(i){
    z <- setdiff(y, x[seq_len(i - incl.eq)])
    if(length(z)) cbind(x[i], z, deparse.level = 0)
  }
  do.call(rbind, lapply(seq_along(x), g))
}


#' Returns MIs from lavaan fit object.
#' @param fit An object from lavaan.
#' @return If successful, returns MIs for an individual. If unsuccessful, 
#' returns NA.
#' @keywords internal 
return.mis <- function(fit){
  zero_se  <- FALSE
  no_paths <- FALSE
  error    <- any(grepl("error", class(fit)))
  if (!error){
    no_paths <- sum(lavInspect(fit, "free")$beta, na.rm = TRUE) == 0 
  }
  if (!error & !no_paths){
    zero_se <- sum(lavInspect(fit, "se")$beta, na.rm = TRUE) == 0
  } 
  
  if (!error & !zero_se){
    mis   <- tryCatch(modindices(fit, op = "~", 
                                 standardized = FALSE,
                                 sort. = FALSE), 
                      error = function(e) e)
    error <- any(grepl("error", class(mis)))
    if (error) mis <- NA 
  } else {
    mis <- NA
  }
  return(mis)
}


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
                       chisq_cutoff){
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
                       count = ave(sig, param, FUN = sum))
  
  mi_list   <- subset(mi_list, !duplicated(param))
  mi_list   <- mi_list[order(-mi_list$count, -mi_list$sum), ]
  
  if (!is.null(prop_cutoff)){
    add_param <- ifelse(mi_list$count[1] > (prop_cutoff*n_converge),
                        mi_list$param[1], NA)
    if (n_converge <= (n_subj/2)) add_param <- NA
  } else {
    add_param <- mi_list$param[1L]
    if (count.excellent(indices) >= 2) add_param <- NA
  }
  
  return(add_param)
}

#' Returns z values from lavaan fit object.
#' @param fit An object from lavaan.
#' @return If successful, returns z values for an individual. If unsuccessful, 
#' returns NA.
#' @keywords internal 
return.zs <- function(fit){
  
  op  = NULL # appease CRAN check
  
  converge <- lavInspect(fit, "converged")
  error   <- any(grepl("error", class(fit)))
  zero_se <- sum(lavInspect(fit, "se")$beta, na.rm = TRUE) == 0
  if (!error & !zero_se & converge){
    zs <- tryCatch(subset(standardizedSolution(fit), 
                          op == "~"),
                   error = function(e) e)
    error <- any(grepl("error", class(zs)))
    if (error) zs <- NA 
  } else {
    zs <- NA
  }
  return(zs)
}


#' Identifies lowest z value from list of z values.
#' @param z_list A list of z values across individuals.
#' @param elig_paths A character vector containing eligible paths that
#' gimme is allowed to drop from the model at a given stage.
#' @param prop_cutoff The proportion of individuals for whom a path must be
#' nonsignificant in order for it to be dropped from the models. NULL if used 
#' at the individual-level.
#' @param n_subj The number of subjects in a given stage of the search. If
#' in the group stage, n_subj equals the number of subjects. If in the subgroup
#' stage, n_subj equals the number of individuals in a given subgroup. At the 
#' individual stage, n_subj = 1.
#' @return Returns name of parameter associated with lowest z. If no z meets 
#' the criteria, returns NA.
#' @keywords internal 
lowest.z <- function(z_list, 
                     elig_paths, 
                     prop_cutoff,
                     n_subj){
  
  param  = NULL # appease CRAN check
  z      = NULL # appease CRAN check
  sig    = NULL # appease CRAN check
  
  z_list     <- z_list[!is.na(z_list)]
  n_converge <- length(z_list)
  z_list     <- do.call("rbind", z_list)
  z_list$param <- paste0(z_list$lhs, z_list$op, z_list$rhs)
  z_list       <- subset(z_list, param %in% elig_paths,
                         select = c("param", "z"))
  z_list$sig   <- ifelse(abs(z_list$z) > abs(qnorm(.025/n_subj)), 1, 0)
  z_list <- transform(z_list,
                      sum   = ave(abs(z), param, FUN = sum),
                      count = ave(sig, param, FUN = sum))
  z_list <- subset(z_list, !duplicated(param))
  z_list <- z_list[order(z_list$count, z_list$sum), ]
  if (!is.null(prop_cutoff)){
    drop_param <- ifelse(z_list$count[1L] <= (prop_cutoff*n_converge),
                         z_list$param[1L], NA)
    if (n_converge <= (n_subj/2)) drop_param <- NA
  } else {
    drop_param <- ifelse(z_list$sig[1L] == 0, z_list$param[1L], NA)
  }
  return(drop_param)
}

#' Prunes paths. Ties together lowest.z and return.zs functions. 
#' @param base_syntax A character vector containing syntax that never changes.
#' @param fixed_syntax A character vector containing syntax that does not change
#' in a given stage of pruning.
#' @param add_syntax A character vector containing the syntax that is allowed
#' to change in a given stage of pruning.
#' @param data_list A list of datasets to be used in a given stage of the 
#' search. Varies based on group, subgroup, or individual-level stage.
#' @param n_paths The number of paths that are eligible for pruning. Equal
#' to the number of paths in add_syntax.
#' @param n_subj The number of subjects in a given stage of the search. If
#' in the group stage, n_subj equals the number of subjects. If in the subgroup
#' stage, n_subj equals the number of individuals in a given subgroup. At the 
#' individual stage, n_subj = 1.
#' @param prop_cutoff The proportion of individuals for whom a path must be
#' nonsignificant in order for it to be dropped from the models. NULL if used 
#' at the individual-level.
#' @param elig_paths A character vector containing eligible paths that
#' gimme is allowed to drop from the model at a given stage.
#' @param subgroup_stage Logical. Only present in order to instruct gimme
#' what message to print to console using writeLines.
#' @return Returns updated values of n_paths and add_syntax.
#' @keywords internal 
prune.paths <- function(base_syntax, 
                        fixed_syntax,
                        add_syntax,
                        data_list, 
                        n_paths, 
                        n_subj, 
                        prop_cutoff, 
                        elig_paths, 
                        subgroup_stage){
  prune <- TRUE
  while(prune){      
    z_list <- list()
    for (k in 1:n_subj){
      if (!is.null(prop_cutoff)){
        if (is.list(fixed_syntax)){
          fit         <- fit.model(syntax    = c(base_syntax, 
                                                 fixed_syntax[[k]], 
                                                 add_syntax),
                                   data_file = data_list[[k]])
        } else {
          fit         <- fit.model(syntax    = c(base_syntax, 
                                                 fixed_syntax, 
                                                 add_syntax),
                                   data_file = data_list[[k]]) 
        } 
        if (subgroup_stage){
          writeLines(paste("subgroup-level pruning, subject", k))
        } else {
          writeLines(paste("group-level pruning, subject", k))
        }
      } else{
        fit         <- fit.model(syntax    = c(base_syntax, 
                                               fixed_syntax, 
                                               add_syntax),
                                 data_file = data_list) 
      }
      z_list[[k]] <- return.zs(fit)
    }
    
    if(!all(is.na(z_list))){
      drop_param <- lowest.z(z_list,
                             elig_paths  = elig_paths,
                             prop_cutoff = prop_cutoff, 
                             n_subj      = n_subj)
    } else {
      drop_param <- NA
    }
    if (!is.na(drop_param)){
      n_paths <- n_paths - 1
      add_syntax <- add_syntax[!add_syntax %in% drop_param]
    } else prune <- FALSE
  }
  res <- list(n_paths, add_syntax)
  return(res)
}

#' Searches for paths. Ties together highest.mi and return.mis functions.
#' @param base_syntax A character vector containing syntax that never changes.
#' @param fixed_syntax A character vector containing syntax that does not change
#' in a given stage of searching.
#' @param add_syntax A character vector containing the syntax that is allowed
#' to change in a given stage of searching.
#' @param data_list A list of datasets to be used in a given stage of the 
#' search. Varies based on group, subgroup, or individual-level stage.
#' @param n_paths The number of paths present in a given stage of searching.
#' Equal to the number of paths in add_syntax.
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
search.paths <- function(base_syntax, 
                         fixed_syntax,
                         add_syntax,
                         n_paths, 
                         data_list, 
                         elig_paths, 
                         prop_cutoff, 
                         n_subj, 
                         chisq_cutoff,
                         subgroup_stage){
  search <- TRUE
  while(search){
    mi_list <- list() 
    indices <- NULL
    for (k in 1:n_subj){
      if (!is.null(prop_cutoff)){
        fit        <- fit.model(syntax    = c(base_syntax, 
                                              fixed_syntax, 
                                              add_syntax),
                                data_file = data_list[[k]])
        if(!subgroup_stage){
          writeLines(paste("group-level search, subject", k))
        } else {
          writeLines(paste("subgroup-level search, subject", k))
        }
      } else 
        {
        fit        <- fit.model(syntax    = c(base_syntax, 
                                              fixed_syntax, 
                                              add_syntax),
                                data_file = data_list)
        if (!"error" %in% class(fit) & lavInspect(fit, "converged")){
          indices    <- fitMeasures(fit, c("chisq", "df", "pvalue", "rmsea", 
                                           "srmr", "nnfi", "cfi"))
        } else indices <- NULL
      }
      mi_list[[k]] <- return.mis(fit)
    }
    
    if (!all(is.na(mi_list))){
      add_param <- highest.mi(mi_list      = mi_list,
                              indices      = indices,
                              elig_paths   = elig_paths,
                              prop_cutoff  = prop_cutoff, 
                              n_subj       = n_subj,
                              chisq_cutoff = chisq_cutoff)
    } else {
      add_param <- NA
    }
    if (!is.na(add_param)){
      n_paths     <- n_paths + 1
      add_syntax  <- append(add_syntax, add_param)
    } else search <- FALSE
  }  
  res <- list(n_paths, add_syntax)
  return(res)
}


#' Determines subgroups.
#' @param base_syntax A character vector containing syntax that never changes.
#' @param data_list A list of all datasets.
#' @param n_subj The number of subjects in the sample.
#' @param file_order A data frame containing the order of the files and 
#' the names of the files. Used to merge in subgroup assignment and preserve
#' order. 
#' @param chisq_cutoff Cutoff used in order for MI to be considered significant.
#' @param elig_paths A character vector containing eligible paths that
#' gimme is allowed to add to the model. Ensures only EPCs from allowable paths
#' are considered in the creation of the similarity matrix.
#' @param confirm_subgroup A dataframe with the first column a string vector of data file names
#' without extensions and the second vector a integer vector of subgroup labels.
#' @return Returns sub object containing similarity matrix, the number of
#' subgroups, the modularity associated with the subgroup memberships, 
#' and a data frame containing the file names and subgroup memberships.
#' @keywords internal 
determine.subgroups <- function(data_list,
                                base_syntax, 
                                n_subj,
                                chisq_cutoff,
                                file_order,
                                elig_paths,
                                confirm_subgroup, 
                                out_path = NULL,
                                sub_feature = sub_feature){
  
  sub_membership  = NULL # appease CRAN check
  
  sub     <- list()
  z_list  <- list()
  mi_list <- list()
  
  for (k in 1:n_subj){
    writeLines(paste("subgroup search, subject", k))
    fit          <- fit.model(syntax    = base_syntax,
                              data_file = data_list[[k]])
    z_list[[k]]  <- return.zs(fit)
    mi_list[[k]] <- return.mis(fit)
  }
  
  names(z_list) <- names(mi_list) <- names(data_list)
  
  # drop individuals who did not converge
  drop    <- unique(c(which(is.na(z_list)), which(is.na(mi_list))))
  
  if (length(drop) != 0){
    mi_list <- mi_list[-drop]
    z_list  <- z_list[-drop]
  }
  
  mi_list_temp <- lapply(mi_list, 
                         function(x){x$param <- paste0(x$lhs, x$op, x$rhs)
                         x$sig   <- ifelse(x$mi > chisq_cutoff, 1, 0)
                         return(x)})
  
  mi_list <- lapply(mi_list_temp, 
                    function(x){subset(x, x$param %in% elig_paths)})
  

  z_list <- lapply(z_list, 
                   function(x){x$sig <- ifelse(x$p < .05/n_subj, 1, 0)
                   return(x)})
  
  #subgroup based only on contemporaneous paths kmg
  if(sub_feature == "contemp"){
    mi_list <- lapply(mi_list, 
                      function(x){x[-grep('lag',mi_list[[1]]$rhs),]})
    z_list <- lapply(z_list, 
                     function(x){x[-grep('lag',z_list[[1]]$rhs),]})
  }
  
  #subgroup based only on lagged paths kmg
  if(sub_feature == "lagged"){
    mi_list <- lapply(mi_list, 
                      function(x){x[grep('lag',mi_list[[1]]$rhs),]})
    z_list <- lapply(z_list, 
                     function(x){x[grep('lag',z_list[[1]]$rhs),]})
  }
  
  # remove lines that have "NA" from z_list (occurs for exog and ar=FALSE)
  # commented out by stl april 2018 - likely to have unintended consequences
  # because it will cause off-by-one errors in the creation of the similarity matrix
  # these NA issues should instead by captured in the na.rm = T arguments
  # added to the similarity matrix. if not, we can revisit
  # z_list <- lapply(z_list, na.exclude) 
  
  sim_mi <- matrix(0, ncol = length(mi_list), nrow = length(mi_list))
  sim_z  <- sim_mi
  
  ## march 2018 stl - na.rm = TRUE added in creation of similarity matrix. 
  ## This was added to address cases where the standard errors for certain paths,
  ## particularly bidirectional paths, were NA. This caused NAs throughout
  ## the adjacency matrix. We should consider whether this is the permanent 
  ## solution we want, as it means that individuals contribute different numbers
  ## of paths to the adjacency matrix (i.e., those individuals with paths
  ## that have NA standard errors contribute fewer paths to the matrix)
  
  for (i in 1:length(mi_list)){
    for (j in 1:length(mi_list)){
      sim_mi[i,j] <- sum(mi_list[[i]]$sig == 1 & mi_list[[j]]$sig == 1 & 
                           sign(mi_list[[i]]$epc) == sign(mi_list[[j]]$epc), na.rm = TRUE)
      sim_z[i,j]  <- sum(z_list[[i]]$sig == 1 & z_list[[j]]$sig == 1 &
                           sign(z_list[[i]]$z) == sign(z_list[[j]]$z), na.rm = TRUE)
    }
  }
  
  sim           <- sim_mi + sim_z
  sim           <- sim - min(sim, na.rm = TRUE)
  diag(sim)     <- 0
  colnames(sim) <- rownames(sim) <- names(mi_list)
  if(is.null(confirm_subgroup)){
    res        <- walktrap.community(graph.adjacency(sim, mode = "undirected"), 
                                     steps = 4)
    sub_mem    <- data.frame(names      = names(membership(res)), 
                             sub_membership = as.numeric(membership(res)))
    sub$sim         <- sim
    sub$n_subgroups <- length(unique(na.omit(sub_mem$sub_membership))) 
    sub$modularity  <- modularity(res)
    sub$sub_mem     <- merge(file_order, sub_mem, by = "names", all.x = TRUE)
  } else {
    sub_mem         <- confirm_subgroup
    names(sub_mem)  <- c("names", "sub_membership")
    sub$sim         <- sim
    sub$n_subgroups <- length(unique(na.omit(sub_mem$sub_membership))) 
    sub$sub_mem     <- merge(file_order, sub_mem, by = "names", all.x = TRUE)
    sub$modularity  <- modularity(graph.adjacency(sim, mode = "undirected"), (sub$sub_mem)$sub_membership)
    
    
  }
  return(sub)
}

#' Individual-level search. Used in gimmeSEM, aggSEM, indSEM.
#' @param dat A list containing information created in setup().
#' @param grp A list containing group-level information. NULL in aggSEM and
#' indSEM.
#' @param ind A list containing individual- and (potentially) subgroup-level
#' information.
#' @return Lists associated with coefficients, fit indices, etc.
#' @keywords internal 
indiv.search <- function(dat, grp, ind){
  
  if (!dat$agg){
    ind$ind_paths   <-  vector("list", dat$n_subj)
    ind$n_ind_paths <- 0
  } else {
    ind <- NULL
  }
  
  status  <- list()
  fits    <- list()
  coefs   <- list()
  betas   <- list()
  vcov    <- list()
  plots   <- list()
  
  n_ind    <- ifelse(dat$agg, 1, dat$n_subj) 
  
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
      writeLines(paste("individual-level search, subject", k))
    }
    
    ind_spec <- search.paths(base_syntax  = dat$syntax, 
                             fixed_syntax = c(grp$group_paths, 
                                              ind$sub_paths[[k]]),
                             add_syntax   = character(),
                             n_paths      = 0,
                             data_list    = data_list,
                             elig_paths   = dat$candidate_paths,
                             prop_cutoff  = NULL,
                             n_subj       = 1,
                             chisq_cutoff = qchisq(.99, 1))
    
    temp_ind_spec <- ind_spec
    
    ind_spec <- prune.paths(base_syntax  = dat$syntax,
                            fixed_syntax = c(grp$group_paths, 
                                             ind$sub_paths[[k]]),
                            add_syntax   = ind_spec[[2]],
                            data_list    = data_list,
                            n_paths      = ind_spec[[1]],
                            n_subj       = 1,
                            prop_cutoff  = NULL,
                            elig_paths   = ind_spec[[2]])
    
    if (!identical(temp_ind_spec, ind_spec)){
      ind_spec <- search.paths(base_syntax  = dat$syntax, 
                               fixed_syntax = c(grp$group_paths, 
                                                ind$sub_paths[[k]]),
                               add_syntax   = ind_spec[[2]],
                               n_paths      = ind_spec[[1]],
                               data_list    = data_list,
                               elig_paths   = dat$candidate_paths,
                               prop_cutoff  = NULL,
                               n_subj       = 1,
                               chisq_cutoff = 0)
    }
    
    ind$ind_paths[[k]] <- ind_spec[[2]]
    ind$n_ind_paths[k] <- ind_spec[[1]]
    
    s10 <- get.params(dat, 
                      grp, 
                      ind,
                      k)
    
    status[[k]] <- s10$status
    fits[[k]]   <- s10$ind_fit
    coefs[[k]]  <- s10$ind_coefs
    betas[[k]]  <- s10$ind_betas
    vcov[[k]]   <- s10$ind_vcov
    plots[[k]]  <- s10$ind_plot
  }
  
  if (dat$agg){
    names(status) <- names(fits) <- names(coefs) <- 
      names(betas) <- names(vcov) <- names(plots) <- "all"
  # } else if (ind$n_ind_paths[k] > 0 & !dat$agg){
  #   names(status) <- names(fits) <- names(coefs) <- 
  #     names(betas) <- names(vcov) <- names(plots) <- names(dat$ts_list)
  } else {
    names(status) <- names(fits) <- names(coefs) <- 
      names(betas) <- names(vcov) <- names(plots) <- names(dat$ts_list)
  }
  
  res <- list("status" = status,
              "fits"   = fits,
              "coefs"  = coefs,
              "betas"  = betas,
              "vcov"   = vcov,
              "plots"  = plots,
              "ind"    = ind)
  return(res)
}


#' Grabs final coefficients for each individual.
#' @param dat A list containing information created in setup().
#' @param grp A list containing group-level information. NULL in aggSEM and
#' indSEM.
#' @param ind A list containing individual- and (potentially) subgroup-level
#' information.
#' @param k The counter indicating the individual.
#' @return Individual-level information on fit, coefficients, and plots.
#' @keywords internal
get.params <- function(dat, grp, ind, k){
  
  op  = NULL # appease CRAN check
  ind_plot = NA
  
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
  converge <- lavInspect(fit, "converged")
  # if (ind$n_ind_paths[k] > 0){ commented out on 11.20.17 by stl
  # potentially insert some other check for an empty model
  zero_se  <- sum(lavInspect(fit, "se")$beta, na.rm = TRUE) == 0
  # } else{ zero_se <- FALSE} commented out on 11.20.17 by stl
  
  # if no convergence, roll back one path at individual level, try again 
  if (!converge | zero_se){
    status <- "nonconvergence"
    if (length(ind$ind_paths[[k]]!= 0)){
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
    }
    converge <- lavInspect(fit, "converged")
    ind_coefs <- subset(standardizedSolution(fit), op == "~") # if betas = 0, no SEs
    if (length(ind_coefs[,1]) > 0){
      zero_se  <- sum(lavInspect(fit, "se")$beta, na.rm = TRUE) == 0}
    else
    {zero_se <- FALSE}
    if (converge){
      status <- "last known convergence"
    }
  }
  
  if (converge & !zero_se){#& (ind$n_ind_paths[k] >0) ){
    status   <- "converged normally"
    
    ind_fit    <- fitMeasures(fit, c("chisq", "df", "npar", "pvalue", "rmsea", 
                                     "srmr", "nnfi", "cfi", "bic", "aic", "logl"))
    ind_fit    <- round(ind_fit, digits = 4)
    ind_fit[2] <- round(ind_fit[2], digits = 0)
    
    ind_vcov  <- lavInspect(fit, "vcov.std.all")
    keep      <- rownames(ind_vcov) %in% dat$candidate_paths
    ind_vcov  <- ind_vcov[keep, keep]
    
    ind_coefs <- subset(standardizedSolution(fit), op == "~")
    
   # if (length(ind_coefs[,1]) > 0){ # stl comment out 11.20.17
      ind_betas <- round(lavInspect(fit, "std")$beta, digits = 4)
      ind_ses   <- round(lavInspect(fit, "se")$beta, digits = 4)
      
      ind_betas <- ind_betas[(dat$n_lagged+1):(dat$n_lagged + dat$n_lagged), ]
      ind_ses   <- ind_ses[(dat$n_lagged+1):(dat$n_lagged + dat$n_lagged), ]
      
      rownames(ind_betas) <- rownames(ind_ses) <- dat$varnames[(dat$n_lagged+1):(dat$n_lagged + dat$n_lagged)]
      colnames(ind_betas) <- colnames(ind_ses) <- dat$varnames
 #   } # stl comment out 11.20.17 
    
    if (dat$agg & !is.null(dat$out)){
      write.csv(ind_betas, file.path(dat$out, "allBetas.csv"), 
                row.names = TRUE)
      write.csv(ind_ses, file.path(dat$out, "allStdErrors.csv"), 
                row.names = TRUE)
    } else if (!dat$agg & !is.null(dat$out)) { # & ind$n_ind_paths[k]>0)
      write.csv(ind_betas, file.path(dat$ind_dir, 
                                     paste0(dat$file_order[k,2], 
                                            "Betas.csv")), row.names = TRUE)
      write.csv(ind_ses, file.path(dat$ind_dir,
                                   paste0(dat$file_order[k,2], 
                                          "StdErrors.csv")), row.names = TRUE)
    }
    
    ind_plot  <- NA
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
      
      ind_plot <- tryCatch(qgraph(plot_vals,
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
                                  DoNotPlot    = TRUE), 
                           error = function(e) e)
      
      if (!is.null(dat$out) & !"error" %in% class(ind_plot)){
        pdf(plot_file)
        plot(ind_plot)
        dev.off()
      }
    }
  } 
  
  # commented out on 11.20.17 by stl 
  # if (ind$n_ind_paths[k] ==0 & converge) {
  #   status     <- "no paths added"
  #   ind_fit    <- fitMeasures(fit, c("chisq", "df", "npar", "pvalue", "rmsea", 
  #                                    "srmr", "nnfi", "cfi", "bic", "aic", "logl"))
  #   ind_fit    <- round(ind_fit, digits = 4)
  #   ind_fit[2] <- round(ind_fit[2], digits = 0)
  #   
  #   ind_vcov  <- lavInspect(fit, "vcov.std.all")
  #   keep      <- rownames(ind_vcov) %in% dat$candidate_paths
  #   ind_vcov  <- ind_vcov[keep, keep]
  #   
  #   ind_betas <- NULL
  #   ind_coefs <- subset(standardizedSolution(fit), op == "~")
  #   
  # } 
  
  if (!converge | zero_se){
    if (!converge) status <- "nonconvergence"
    if (zero_se)   status <- "computationally singular"
    ind_fit   <- rep(NA, 8)
    ind_coefs <- matrix(NA, nrow = 1, ncol = 7)
    colnames(ind_coefs) <- c("lhs", "op", "rhs", "est.std", "se", "z", "pvalue")
    ind_betas <- NA
    ind_vcov  <- NA
    ind_plot  <- NA
  }
  
  if (!dat$plot)
    ind_plot  <- NA
  
  res <- list("status"    = status, 
              "ind_fit"   = ind_fit, 
              "ind_coefs" = ind_coefs, 
              "ind_betas" = ind_betas, 
              "ind_vcov"  = ind_vcov,
              "ind_plot"  = ind_plot)
  return(res)
}


#' Wrapup, create output files.
#' @param dat A list containing information created in setup().
#' @param grp A list containing group-level information. NULL in aggSEM and
#' indSEM.
#' @param ind A list containing individual- and (potentially) subgroup-level
#' information.
#' @param sub A list containing subgroup information.
#' @param sub_spec A list containing information specific to each subgroup.
#' @param store A list containing output from indiv.search().
#' @return Aggregated information, such as counts, levels, and plots.
#' @keywords internal
final.org <- function(dat, grp, ind, sub, sub_spec, store){
  
  sub_coefs  <- list()
  sub_summ   <- list()
  sub_plots  <- list()
  sub_counts <- list()
  
  param  = NULL # appease CRAN check
  est.std = NULL # appease CRAN check
  count  = NULL # appease CRAN check
  
  if (!dat$agg){
    
    coefs       <- do.call("rbind", store$coefs)
    
    if(length(coefs[,1])>0){
      coefs$id    <- rep(names(store$coefs), sapply(store$coefs, nrow))
      coefs$param <- paste0(coefs$lhs, coefs$op, coefs$rhs)
      
      coefs$level[coefs$param %in% c(grp$group_paths, dat$syntax)] <- "group"
      coefs$level[coefs$param %in% unique(unlist(ind$ind_paths))]  <- "ind"
      coefs$color[coefs$level == "group"] <- "black"
      coefs$color[coefs$level == "ind"]   <- "gray50"}
    
    indiv_paths <- NULL
    samp_plot <- NULL
    sample_counts <- NULL
    # if (length(coefs[,1])>0){ # commented out stl 11.20.17
    if (dat$subgroup) {
      if (sub$n_subgroups != dat$n_subj){ # ensure everyone isn't in their own subgroup
        
        sub_paths_count <- table(unlist(
          lapply(sub_spec, FUN = function(x) c(x$sub_paths))))
        sub_to_group    <- names(
          sub_paths_count[sub_paths_count == sub$n_subgroups])
        
        for (s in 1:sub$n_subgroups){
          sub_s_mat_counts <- matrix(0, nrow = (dat$n_vars_total), 
                                     ncol = (dat$n_vars_total))
          sub_s_mat_means  <- sub_s_mat_counts
          sub_s_mat_colors <- matrix(NA, nrow = (dat$n_vars_total), 
                                     ncol = (dat$n_vars_total))
          
          sub_s_coefs <- coefs[coefs$id %in% sub_spec[[s]]$sub_s_subjids, ]
          sub_s_coefs$level[sub_s_coefs$param %in% sub_spec[[s]]$sub_paths] <- "sub"
          sub_s_coefs$level[sub_s_coefs$param %in% sub_to_group] <- "group"
          sub_s_coefs$level[sub_s_coefs$param %in% unique(
            unlist(ind[ind$sub_membership == s, ]$ind_paths))] <- "ind"
          sub_s_coefs$color[sub_s_coefs$level == "group"] <- "black"
          sub_s_coefs$color[sub_s_coefs$level == "sub"]   <- "green3"
          sub_s_coefs$color[sub_s_coefs$level == "ind"]   <- "gray50"
          
          ## march 2018 stl - fix to remove error caused where lhs and rhs 
          ## values are NA. there's no deeper trouble here - it was just due to an 
          ## rbind where individuals with no paths (e.g., entirely NA) were included
          ## in the full rbind, which led to variable names of "NA" 
          sub_s_coefs <- sub_s_coefs[!is.na(sub_s_coefs$lhs), ]
          sub_s_coefs <- sub_s_coefs[!is.na(sub_s_coefs$rhs), ]
          
          sub_s_summ <- transform(sub_s_coefs, 
                                  count = as.numeric(
                                    ave(param, param, FUN = length)),
                                  mean  = ave(est.std, param, FUN = sum)/sub_spec[[s]]$n_sub_subj)
          sub_s_summ <- subset(sub_s_summ, !duplicated(param))
          sub_s_summ$row <- match(sub_s_summ$lhs, dat$lvarnames)
          sub_s_summ$col <- match(sub_s_summ$rhs, dat$lvarnames)
          sub_s_summ$mem <- s
          
          sub_s_mat_counts[cbind(sub_s_summ$row, sub_s_summ$col)] <- 
            as.numeric(as.character(sub_s_summ$count))
          sub_s_mat_counts <- sub_s_mat_counts[(dat$n_lagged+1):(dat$n_vars_total), ]
          colnames(sub_s_mat_counts) <- dat$varnames
          rownames(sub_s_mat_counts) <- dat$varnames[(dat$n_lagged+1):(dat$n_vars_total)]
          
          sub_s_mat_means[cbind(sub_s_summ$row, sub_s_summ$col)]  <- sub_s_summ$mean
          sub_s_mat_colors[cbind(sub_s_summ$row, sub_s_summ$col)] <- sub_s_summ$color
          sub_s_mat_colors <- sub_s_mat_colors[(dat$n_lagged+1):(dat$n_vars_total), ]
          
          if (dat$plot & sub_spec[[s]]$n_sub_subj != 1){ #plot subgroup plot if >1 nodes in subgroup
            
            sub_s_counts <- t(sub_s_mat_counts/sub_spec[[s]]$n_sub_subj)
            lagged     <- sub_s_counts[1:(dat$n_lagged), ]
            contemp    <- sub_s_counts[(dat$n_lagged+1):(dat$n_vars_total), ]
            plot_vals  <- rbind(w2e(lagged), w2e(contemp))
            is_lagged  <- c(rep(TRUE, sum(lagged != 0)), rep(FALSE, sum(contemp != 0)))
            
            sub_colors <- t(sub_s_mat_colors)
            colors     <- c(sub_colors[1:(dat$n_lagged), ],
                            sub_colors[(dat$n_lagged+1):(dat$n_vars_total), ])
            colors     <- colors[!is.na(colors)]
            
            sub_plot <- tryCatch(qgraph(plot_vals,
                                        layout       = "circle",
                                        lty          = ifelse(is_lagged, 2, 1),
                                        edge.labels  = FALSE,
                                        edge.color   = colors,
                                        parallelEdge = TRUE,
                                        fade         = FALSE,
                                        labels       = 
                                          dat$varnames[(dat$n_lagged+1):(dat$n_vars_total)],
                                        label.cex    = 2,
                                        DoNotPlot    = TRUE), 
                                 error = function(e) e)
            
            if (!is.null(dat$out)){
              pdf(file.path(dat$subgroup_dir, 
                            paste0("subgroup", s, "Plot.pdf")))
              plot(sub_plot)
              dev.off()
            }
            
          } else {
            sub_plot         <- NULL
            sub_s_mat_counts <- NULL
          }
          
          if (sub_spec[[s]]$n_sub_subj != 1 & !is.null(dat$out)){
            write.csv(sub_s_mat_counts, 
                      file = file.path(dat$subgroup_dir, 
                                       paste0("subgroup", s, 
                                              "PathCountsMatrix.csv")), 
                      row.names = TRUE)
          }
          sub_coefs[[s]] <- sub_s_coefs
          sub_summ[[s]]  <- sub_s_summ
          sub_plots[[s]] <- sub_plots
          sub_counts[[s]] <- sub_s_mat_counts
        }
        
        summ <- do.call("rbind", sub_summ)
        coefs <- do.call("rbind", sub_coefs)
        
      } else {
        sub_coefs <- NULL
        sub_plots <- NULL
        sub_paths <- NULL
        summ <- transform(coefs, count = as.numeric(
          ave(param, param, FUN = length)))
        summ <- subset(summ, !duplicated(param)) 
      }
    }
    else {
      sub_coefs <- NULL
      sub_plots <- NULL
      sub_paths <- NULL
      summ <- transform(coefs, count = as.numeric(
        ave(param, param, FUN = length)))
      summ <- subset(summ, !duplicated(param)) 
    }
    
    # combining and creating wide summaryPathCounts -------------------------- #
    summ$label <- ifelse(summ$level == "sub", 
                         paste0("subgroup", summ$mem),
                         summ$level)
    a <- aggregate(count ~ lhs + rhs + label, data = summ, sum)
    
    a <- a[order(-a$count, a$label),]
    a <- reshape(a, timevar = "label", idvar = c("lhs", "rhs"), 
                 direction = "wide")
    a[is.na(a)] <- 0
    a$lhs <- recode.vars(a$lhs, dat$lvarnames, dat$varnames)
    a$rhs <- recode.vars(a$rhs, dat$lvarnames, dat$varnames)
    colnames(a)[1:2] <- c("dv", "iv")
    
    if (!is.null(dat$out)){
      write.csv(a, file.path(dat$out, "summaryPathCounts.csv"), 
                row.names = FALSE)
    }
    
    # end creating wide summaryPathCounts ------------------------------------ #
    
    b <- aggregate(count ~ lhs + rhs + color + label + param, data = summ, sum)
    b <- transform(b, xcount = ave(count, param, FUN = sum))
    # sorting by count and then dropping duplicated parameters
    # ensures that subgroup paths will be displayed as green
    # and individual paths that appear in another subgroup
    # will not cause subgroup paths to all display as individual
    b <- b[order(-b$count), ]
    c <- b[!duplicated(b$param), c("lhs", "rhs", "color", "xcount")] 
    
    c$row <- match(c$lhs, dat$lvarnames) - dat$n_lagged
    c$col <- match(c$rhs, dat$lvarnames)
    
    sample_counts <- matrix(0, ncol = (dat$n_vars_total), nrow = dat$n_lagged)
    sample_counts[cbind(c$row, c$col)] <- c$xcount
    colnames(sample_counts) <- dat$varnames
    rownames(sample_counts) <- dat$varnames[(dat$n_lagged+1):(dat$n_lagged+dat$n_lagged)]
    
    if (dat$plot){
      
      sample_colors <- matrix(NA, ncol = (dat$n_vars_total), nrow = dat$n_lagged)
      sample_colors[cbind(c$row, c$col)] <- c$color
      
      sample_paths  <- t(sample_counts)/dat$n_subj
      
      lagged     <- sample_paths[1:(dat$n_lagged), ]
      contemp    <- sample_paths[(dat$n_lagged+1):(dat$n_vars_total), ]
      plot_vals  <- rbind(w2e(lagged), w2e(contemp))
      is_lagged  <- c(rep(TRUE, sum(lagged != 0)),
                      rep(FALSE, sum(contemp != 0)))
      
      samp_colors <- t(sample_colors)
      colors      <- c(samp_colors[1:(dat$n_lagged), ],
                       samp_colors[(dat$n_lagged+1):(dat$n_vars_total), ])
      colors      <- colors[!is.na(colors)]
      
      samp_plot <- tryCatch(qgraph(plot_vals,
                                   layout       = "circle",
                                   lty          = ifelse(is_lagged, 2, 1),
                                   edge.labels  = FALSE,
                                   edge.color   = colors,
                                   parallelEdge = TRUE,
                                   fade         = FALSE,
                                   labels       = 
                                     dat$varnames[(dat$n_lagged+1):(dat$n_vars_total)],
                                   label.cex    = 2,
                                   DoNotPlot    = TRUE), 
                            error = function(e) e)
      
      if (!is.null(dat$out)){
        pdf(file.path(dat$out, "summaryPathsPlot.pdf"))
        plot(samp_plot)
        dev.off()
      }
      
    } else samp_plot <- NULL
    
    indiv_paths     <- coefs[, c("id", "lhs", "rhs", "est.std", 
                                 "se", "z", "pvalue", "level")]
    indiv_paths$lhs <- recode.vars(indiv_paths$lhs, dat$lvarnames, dat$varnames)
    indiv_paths$rhs <- recode.vars(indiv_paths$rhs, dat$lvarnames, dat$varnames)
    indiv_paths     <- indiv_paths[order(indiv_paths$id, indiv_paths$level), ]
    colnames(indiv_paths) <- c("file", "dv", "iv", "beta", "se", 
                               "z", "pval", "level")
    # } # end "if no coefficients" commented out stl 11.20.17
    # combine fit information for summaryFit.csv
    
    fits        <- as.data.frame(do.call(rbind, store$fits))
    fits$file   <- rownames(fits)
    fits$status <- do.call(rbind, store$status)
    fits        <- fits[ ,c(12, 1:11, 13)]
    
    if (dat$subgroup){
      fits <- merge(fits, sub$sub_mem[ ,c(1,3)], by.x = "file", by.y = "names")  
      fits$modularity <- c(round(sub$modularity, digits = 4), 
                           rep("", (nrow(fits) - 1)))
      indiv_paths <- merge(indiv_paths, sub$sub_mem[ ,c(1,3)], 
                           by.x = "file", by.y = "names")
    }
    
    if (!is.null(dat$out)){ #& length(coefs[,1]) > 0){ # commented out stl 11.20.17
      write.csv(indiv_paths, file.path(dat$out, "indivPathEstimates.csv"),
                row.names = FALSE)
      write.csv(sample_counts, file.path(dat$out,
                                         "summaryPathCountsMatrix.csv"),
                row.names = FALSE)
      write.csv(fits, file.path(dat$out, "summaryFit.csv"), row.names = FALSE)
      write.csv(sub$sim, file.path(dat$out, "similarityMatrix.csv"), row.names = FALSE)
    }
    
  } else {
    indiv_paths <- store$coefs[[1L]]
    indiv_paths$file <- "all"
    indiv_paths$lhs  <- recode.vars(indiv_paths$lhs, dat$lvarnames, dat$varnames)
    indiv_paths$rhs  <- recode.vars(indiv_paths$rhs, dat$lvarnames, dat$varnames)
    indiv_paths      <- indiv_paths[order(indiv_paths$file), ]
    indiv_paths      <- indiv_paths[ ,c("file", "lhs", "rhs", "est.std", 
                                        "se", "z", "pvalue")]
    colnames(indiv_paths) <- c("file", "dv", "iv", "beta", "se", "z", "pval")
    
    fits          <- store$fits[[1L]]
    file          <- c("all")
    names(file)   <- "file"
    status        <- store$status[[1L]]
    names(status) <- "status"
    fits          <- c(file, fits, status)
    fits <- t(fits)
    
    if (!is.null(dat$out)){
      write.csv(indiv_paths, file.path(dat$out, "allPathEstimates.csv"), 
                row.names = FALSE)
      write.csv(fits, file.path(dat$out, "summaryFit.csv"), row.names = FALSE)
    }
    
    sample_counts <- NULL
    samp_plot     <- NULL
  }
  
  res <- list(fit           = fits,
              param_est     = indiv_paths,
              samp_plot     = samp_plot,
              sub_plots     = sub_plots,
              sample_counts = sample_counts,
              sub_counts    = sub_counts)  
  return(res)
  
}

