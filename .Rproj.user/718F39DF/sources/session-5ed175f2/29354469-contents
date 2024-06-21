#' Create summary matrix of path counts and subgroup plots
#' @param dat A list containing information created in setup().
#' @param grp A list containing group-level information. NULL in aggSEM and
#' indSEM.
#' @param sub A list containing subgroup information.
#' @param sub_spec A list containing information specific to each subgroup.
#' @param store A list containing output from indiv.search().
#' @return Aggregated information, such as counts, levels, and plots.
#' @keywords internal
#' 
summaryPathsCounts <- function(dat, grp, store, sub, sub_spec){
 
  param  = NULL # appease CRAN check
  est.std = NULL # appease CRAN check
  sub_plots  <- list()
  sub_plots_cov <- list()
  sub_coefs  <- list()
  sub_counts <- list()
  sub_counts_cov <- list()
  
  coefs       <- do.call("rbind", store$coefs)
  sub_summ   <- list()
  
  if(length(coefs[,1])>0){
    coefs$id    <- rep(names(store$coefs), sapply(store$coefs, nrow))
    coefs$param <- paste0(coefs$lhs, coefs$op, coefs$rhs)
    coefs <- coefs[!coefs$param %in% dat$nonsense_paths,] # Removes non-sense paths that occur when ar = FALSE or mult_vars is not null from output 
    
    ### kad 7.26.22: make sure paths set to a specific value (e.g. V2 ~ 0.5*V1) are included in group output
    ### with "level" specified as "group"
    # Set default to null
    specificValuePaths <- NULL
    # Check if any paths set to a specific value exist in fixed_paths [note paths set to 0 are already removed]
    if(any(grepl("\\*",dat$fixed_paths))){
      # Separate at both "~" and "*", then paste together var names, skipping the specific multiplier value
      specificPaths <- dat$fixed_paths[grep("\\*",dat$fixed_paths)]
      specificPathsSplit <- strsplit(specificPaths,"\\*|~")
      specificPathsList <- lapply(specificPathsSplit, function(x) paste0(x[1],"~",x[3]))
      # Return paths so they can be included in group-level "coefs$level" list
      specificValuePaths <- unlist(specificPathsList)
    }
    
    coefs$level[coefs$param %in% c(grp$group_paths, dat$syntax, specificValuePaths)] <- "group" # kad 7.26.22 added specificValuePaths created above
    coefs$level[coefs$param %in% unique(unlist(store$ind$ind_paths))]  <- "ind"
    coefs$color[coefs$level == "group"] <- "black"
      coefs$color[coefs$level == "ind"]   <- "gray50"
  }
  
  indiv_paths <- NULL
  samp_plot <- NULL
  samp_plot_cov <- NULL
  sample_counts <- NULL
  sample_counts_corr <- NULL
  # if (length(coefs[,1])>0){ # commented out stl 11.20.17
  if (dat$subgroup) {
    if (sub$n_subgroups != dat$n_subj){ # ensure everyone isn't in their own subgroup
      
      sub_paths_count <- table(unlist(
        lapply(sub_spec, FUN = function(x) c(x$sub_paths))))
      # if path exists in all subgroups, push it up to group level
      sub_to_group    <- names(
        sub_paths_count[sub_paths_count == sub$n_subgroups])
      
      
      for (s in 1:sub$n_subgroups){
        sub_s_mat_counts <- matrix(0, nrow = (dat$n_vars_total), 
                                   ncol = (dat$n_vars_total))
        sub_s_mat_counts_cov <- sub_s_mat_counts
        sub_s_mat_means  <- sub_s_mat_counts
        sub_s_mat_means_cov  <- sub_s_mat_counts
        sub_s_mat_colors <- matrix(NA, nrow = (dat$n_vars_total), 
                                   ncol = (dat$n_vars_total))
        sub_s_mat_colors_cov <- matrix(NA, nrow = (dat$n_vars_total), 
                                       ncol = (dat$n_vars_total))
        
        sub_s_coefs <- coefs[coefs$id %in% sub_spec[[s]]$sub_s_subjids, ]
        sub_s_coefs$level[sub_s_coefs$param %in% sub_spec[[s]]$sub_paths] <- "sub"
        sub_s_coefs$level[sub_s_coefs$param %in% sub_to_group] <- "group"
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
            
            regressions <- 
              sub_s_summ[which(sub_s_summ$op == "~"),]
            sub_s_mat_counts[cbind(regressions$row, regressions$col)] <- 
              as.numeric(as.character(regressions$count))
            sub_s_mat_counts <- sub_s_mat_counts[(dat$n_lagged+1):(dat$n_vars_total), ]
            colnames(sub_s_mat_counts) <- dat$varnames
            rownames(sub_s_mat_counts) <- dat$varnames[(dat$n_lagged+1):(dat$n_vars_total)]
            
            sub_s_mat_means[cbind(regressions$row, regressions$col)]  <- regressions$mean
            sub_s_mat_colors[cbind(regressions$row, regressions$col)] <- regressions$color
            sub_s_mat_colors <- sub_s_mat_colors[(dat$n_lagged+1):(dat$n_vars_total), ]
            
            cov <- 
              sub_s_summ[(which(sub_s_summ$op == "~~")),]
            sub_s_mat_counts_cov[cbind(cov$row, cov$col)] <- 
              as.numeric(as.character(cov$count))
            sub_s_mat_counts_cov <- sub_s_mat_counts_cov[(dat$n_lagged+1):(dat$n_vars_total), ]
            colnames(sub_s_mat_counts_cov) <- dat$varnames
            rownames(sub_s_mat_counts_cov) <- dat$varnames[(dat$n_lagged+1):(dat$n_vars_total)]
            
            sub_s_mat_means_cov[cbind(cov$row, cov$col)]  <- cov$mean
            sub_s_mat_colors_cov[cbind(cov$row, cov$col)] <- cov$color
            sub_s_mat_colors_cov <- sub_s_mat_colors_cov[(dat$n_lagged+1):(dat$n_vars_total), ]
            
            if (dat$plot & sub_spec[[s]]$n_sub_subj != 1){ #plot subgroup plot if >1 nodes in subgroup
              
              sub_s_counts <- t(sub_s_mat_counts/sub_spec[[s]]$n_sub_subj)
              sub_s_counts_cov <- t(sub_s_mat_counts_cov/sub_spec[[s]]$n_sub_subj)
              contemp_cov    <- sub_s_counts_cov[(dat$n_lagged+1):(dat$n_vars_total), ]
              lagged     <- sub_s_counts[1:(dat$n_lagged), ]
              
              contemp    <- sub_s_counts[(dat$n_lagged+1):(dat$n_vars_total), ]
              plot_vals  <- rbind(w2e(lagged), w2e(contemp))
              is_lagged  <- c(rep(TRUE, sum(lagged != 0)), rep(FALSE, sum(contemp != 0)))
              
              sub_colors <- t(sub_s_mat_colors)
              colors     <- c(sub_colors[1:(dat$n_lagged), ],
                              sub_colors[(dat$n_lagged+1):(dat$n_vars_total), ])
              colors     <- colors[!is.na(colors)]
              
              sub_plot <- try(qgraph(plot_vals,
                                     layout       = "circle",
                                     lty          = ifelse(is_lagged, 2, 1),
                                     edge.labels  = FALSE,
                                     edge.color   = colors,
                                     parallelEdge = TRUE,
                                     fade         = FALSE,
                                     labels       = 
                                       dat$varnames[(dat$n_lagged+1):(dat$n_vars_total)],
                                     label.cex    = 2,
                                     DoNotPlot    = TRUE))
              
              sub_plot$graphAttributes$Edges$width <- (plot_vals[,3])*7.137138 
              
              if (sum(contemp_cov)>0){
                
                plot_vals_cov  <- w2e(contemp_cov)
                sub_colors_cov <- t(sub_s_mat_colors_cov)
                #commented out by lan 2.10.2020
                #colors     <- c(sub_colors_cov[1:(dat$n_lagged), ],
                #sub_colors_cov[(dat$n_lagged+1):(dat$n_vars_total), ])
                colors    <- sub_colors_cov[(dat$n_lagged+1):(dat$n_vars_total), ]
                colors     <- colors[!is.na(colors)]
                sub_plot_cov <- try(qgraph(plot_vals_cov,
                                           layout       = "circle",
                                           edge.labels  = FALSE,
                                           edge.color   = colors,
                                           parallelEdge = TRUE,
                                           fade         = FALSE,
                                           arrows       = FALSE,
                                           labels       = 
                                             dat$varnames[(dat$n_lagged+1):(dat$n_vars_total)],
                                           label.cex    = 2,
                                           DoNotPlot    = TRUE))
                
                sub_plot_cov$graphAttributes$Edges$width <- (plot_vals_cov[,3])*7.137138 
              }
              if (!is.null(dat$out) & !inherits(sub_plot, "try-error")){
                pdf(file.path(dat$subgroup_dir, 
                              paste0("subgroup", s, "Plot.pdf")))
                plot(sub_plot)
                dev.off()
                if(sum(sub_s_counts_cov)>0){
                  pdf(file.path(dat$subgroup_dir, 
                                paste0("subgroup", s, "Plot_cov.pdf")))
                  plot(sub_plot_cov)
                  dev.off()}
              }
              
              sub_plots[[s]] <- sub_plot
              if (sum(contemp_cov)>0) sub_plots_cov[[s]] <- sub_plot_cov
              
            } else {
              sub_plot         <- NULL
              sub_s_mat_counts <- NULL
            }
            
            ### Write subgroup path counts matrix to output
            if (sub_spec[[s]]$n_sub_subj != 1 & !is.null(dat$out)){
              write.csv(sub_s_mat_counts, 
                        file = file.path(dat$subgroup_dir, 
                                         paste0("subgroup", s, 
                                                "PathCountsMatrix.csv")), 
                        row.names = TRUE)
              
              ### If hybrid=TRUE or VAR=TRUE, also output covariance counts matrix
              if(dat$hybrid|dat$VAR){
                write.csv(sub_s_mat_counts_cov, 
                          file = file.path(dat$subgroup_dir, 
                                           paste0("subgroup", s, 
                                                  "CovCountsMatrix.csv")), 
                          row.names = TRUE)
              }
            }
            #if (dat$plot & sub_spec[[s]]$n_sub_subj != 1){ ##add by lan 021220: store the sub_plot & sub_plot_cov to the plots when n>1
              sub_summ[[s]]  <- sub_s_summ
              sub_counts[[s]] <- sub_s_mat_counts
              sub_counts_cov[[s]] <- sub_s_mat_counts_cov
             
            #} removed by Katie - prevented singletons from being included in path counts 01302023
            sub_coefs[[s]] <- sub_s_coefs      
      }
      
      summ <- do.call("rbind", sub_summ)
      coefs <- do.call("rbind", sub_coefs)
      
    } else {
      sub_coefs <- NULL
      sub_plots <- NULL
      sub_paths <- NULL
      sub_counts_cov <- NULL
      summ <- transform(coefs, count = as.numeric(
        ave(param, param, FUN = length)))
      summ <- subset(summ, !duplicated(param)) 
    }
    
  } else {
    sub_coefs <- NULL
    sub_plots <- NULL
    sub_paths <- NULL
    sub_plots_cov <- NULL
    sub_counts_cov <- NULL
    summ <- transform(coefs, count = as.numeric(
      ave(param, param, FUN = length)))
    summ <- subset(summ, !duplicated(param)) 
  }
  
  # combining and creating wide summaryPathCounts -------------------------- #
  summ$label <- ifelse(summ$level == "sub", 
                       paste0("subgroup", summ$mem),
                       summ$level)
  a <- aggregate(count ~ lhs + op + rhs + label, data = summ, sum)
  
  a <- a[order(-a$count, a$label),]
  a <- reshape(a, timevar = "label", idvar = c("lhs", "op", "rhs"), 
               direction = "wide")
  a[is.na(a)] <- 0
  a$lhs <- recode.vars(a$lhs, dat$lvarnames, dat$varnames)
  a$rhs <- recode.vars(a$rhs, dat$lvarnames, dat$varnames)
  
  res <- list(a = a, 
              summ = summ,
              coefs = coefs,
              sub_plots = sub_plots,
              sub_counts = sub_counts, 
              sub_plots_cov = sub_plots_cov, 
              sub_counts_cov = sub_counts_cov)
return(res)
}
