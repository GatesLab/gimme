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
                                sub_feature){
  #######################
  # base_syntax  = c(dat$syntax, grp[[i]]$group_paths)
  # data_list    = dat$ts_list
  # n_subj       = dat$n_subj
  # chisq_cutoff = dat$chisq_cutoff_mi_epc
  # file_order   = dat$file_order
  # elig_paths   = dat$candidate_paths
  # confirm_subgroup = confirm_subgroup
  # out_path     = dat$out
  # sub_feature  = sub_feature
  #######################
  
  sub_membership  = NULL # appease CRAN check
  
  sub     <- list()
  z_list  <- list()
  mi_list <- list()
  converge <- matrix(,n_subj,1)

  for (k in 1:n_subj){
    writeLines(paste0("subgroup search, subject ", k, " (",names(data_list)[k],")"))
    fit          <- fit.model(syntax    = base_syntax,
                              data_file = data_list[[k]])
    z_list[[k]]  <- return.zs(fit)
    mi_list[[k]] <- return.mis(fit)
    converge[k]  <- lavInspect(fit, "converged")
  }
  
  names(z_list) <- names(mi_list) <- names(data_list)
  
  # drop individuals who did not converge
  # kmg: this is an imperfect approach because the z_list is NA if no paths 
  # were added. Caused errors, commenting out, added "converge" to return.zs output. 
  # drop    <- unique(c(which(is.na(z_list)), which(is.na(mi_list))))
   drop <- which(converge==FALSE)
   
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
  
  # if no group-level paths added, don't consider z_list
  if (length(which(is.na(z_list)))==0)
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
      if (length(which(is.na(z_list)))==0)
      sim_z[i,j]  <- sum(z_list[[i]]$sig == 1 & z_list[[j]]$sig == 1 &
                           sign(z_list[[i]]$z) == sign(z_list[[j]]$z), na.rm = TRUE)
    }
  }
  
  sim           <- sim_mi + sim_z
  sim           <- sim - min(sim, na.rm = TRUE)
  diag(sim)     <- 0
  colnames(sim) <- rownames(sim) <- names(mi_list)
  if(is.null(confirm_subgroup)){
    res        <- igraph::cluster_walktrap(graph.adjacency(sim, mode = "undirected"), 
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