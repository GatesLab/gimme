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
final.org <- function(dat, grp, sub, sub_spec, diagnos=FALSE, store){
  
  ind = store$ind
  
  param  = NULL # appease CRAN check
  est.std = NULL # appease CRAN check
  count  = NULL # appease CRAN check
  
  if (!dat$agg){
    
    summarize <- summaryPathsCounts(dat, grp, store, sub, sub_spec)
  
    ### If path now exists for >= groupcutoff, rerun individual search with it estimated for all
    if(any(summarize$a$count.ind/dat$n_subj >= dat$groupcutoff)){
      loc <- which(summarize$a$count.ind/dat$n_subj >= dat$groupcutoff)
      grp$group_paths <-c(grp$group_paths, paste0(summarize$a$lhs[loc],summarize$a$op[loc], summarize$a$rhs[loc]))
      if(dat$subgroup){
        store <- indiv.search(dat, grp, ind[[1]])
      } else {
        store <- indiv.search(dat, grp, ind[1])
      }
      
      summarize <- summaryPathsCounts(dat, grp, store, sub, sub_spec)
    }
    
    if (!is.null(dat$out)){
      write.csv(summarize$a, file.path(dat$out, "summaryPathCounts.csv"), 
                row.names = FALSE)
    }
    
    # end creating wide summaryPathCounts ------------------------------------ #
    
    b <- aggregate(count ~ lhs + op + rhs + color + label + param, data = summarize$summ, sum)
    b <- transform(b, xcount = ave(count, param, FUN = sum))
    # sorting by count and then dropping duplicated parameters
    # ensures that subgroup paths will be displayed as green
    # and individual paths that appear in another subgroup
    # will not cause subgroup paths to all display as individual
    # CA 10.5.18 created variable to order by label.  Some individual paths were 
    # being selected over subgroup paths in the duplicated function.
    
    b$labelnum[b$label=='group'] <- 1
    b$labelnum[b$label=='ind'] <- 3
    b$labelnum[is.na(b$labelnum)] <-2
    
    b <- b[order(b$labelnum), ]
    d <- b[!duplicated(b$param), c("lhs", "op", "rhs", "color", "xcount")] 
    
    c_direct <- d[which(d$op == "~"),]
    c_corr   <- d[which(d$op == "~~"),]
    
    c_direct$row <- match(c_direct$lhs, dat$lvarnames) - dat$n_lagged
    c_direct$col <- match(c_direct$rhs, dat$lvarnames)
    c_corr$row <- match(c_corr$lhs, dat$lvarnames) - dat$n_lagged
    c_corr$col <- match(c_corr$rhs, dat$lvarnames)
    
    sample_counts <- matrix(0, ncol = (dat$n_vars_total), nrow = (dat$n_vars_total - dat$n_lagged))
    sample_counts[cbind(c_direct$row, c_direct$col)] <- c_direct$xcount
    sample_counts_corr <- matrix(0, ncol = (dat$n_vars_total), nrow = (dat$n_vars_total - dat$n_lagged))
    sample_counts_corr[cbind(c_corr$row, c_corr$col)] <- c_corr$xcount
    colnames(sample_counts) <- dat$varnames
    rownames(sample_counts) <- dat$varnames[(dat$n_lagged+1):(dat$n_vars_total)]
    colnames(sample_counts_corr) <- dat$varnames
    rownames(sample_counts_corr) <- dat$varnames[(dat$n_lagged+1):(dat$n_vars_total)]
    
    if (dat$plot){
      
      sample_colors <- matrix(NA, ncol = (dat$n_vars_total), nrow = (dat$n_vars_total-dat$n_lagged))
      sample_colors[cbind(c_direct$row, c_direct$col)] <- c_direct$color
      sample_colors_corr <- matrix(NA, ncol = (dat$n_vars_total), nrow = (dat$n_vars_total-dat$n_lagged))
      sample_colors_corr[cbind(c_corr$row, c_corr$col)] <- c_corr$color
      
      sample_paths  <- t(sample_counts)/dat$n_subj
      sample_paths_corr <- t(sample_counts_corr)/dat$n_subj
      
      lagged     <- sample_paths[1:(dat$n_lagged), ]
     
      contemp    <- sample_paths[(dat$n_lagged+1):(dat$n_vars_total), ]
      plot_vals  <- rbind(w2e(lagged), w2e(contemp))
      is_lagged  <- c(rep(TRUE, sum(lagged != 0)),
                      rep(FALSE, sum(contemp != 0)))
      
      samp_colors <- t(sample_colors)
      colors      <- c(samp_colors[1:(dat$n_lagged), ],
                       samp_colors[(dat$n_lagged+1):(dat$n_vars_total), ])
      colors      <- colors[!is.na(colors)]
      
      samp_plot <- try(qgraph(plot_vals,
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

     samp_plot$group_plot_paths$graphAttributes$Edges$width <- (plot_vals[,3])*7.137138  
 
      
      samp_colors_corr <- t(sample_colors_corr)
      #commented out by lan 2.10.2020
      # colors_corr      <- c(samp_colors_corr[1:(dat$n_lagged), ],
      #                samp_colors_corr[(dat$n_lagged+1):(dat$n_vars_total), ])
      #commented out by Katie 12.1.2022
      # colors_corr     <- samp_colors_corr[(dat$n_lagged+1):(dat$n_vars_total)]
      colors_corr      <- samp_colors_corr[!is.na(samp_colors_corr)]
     
      if (sum(sample_paths_corr)>0){
        corr   <- sample_paths_corr[(dat$n_lagged+1):(dat$n_vars_total), ]
        plot_vals_cov  <- w2e(corr)
        samp_plot_cov <- try(qgraph(plot_vals_cov,
                                    layout       = "circle",
                                    edge.labels  = FALSE,
                                    edge.color   = colors_corr,
                                    parallelEdge = TRUE,
                                    fade         = FALSE,
                                    arrows       = FALSE,
                                    labels       = 
                                      dat$varnames[(dat$n_lagged+1):(dat$n_vars_total)],
                                    label.cex    = 2,
                                    DoNotPlot    = TRUE))
      
      samp_plot_cov$graphAttributes$Edges$width <- (plot_vals_cov[,3])*7.137138  
      } else {      
        samp_plot_cov <- NULL
}
      
      if (!is.null(dat$out)){
        pdf(file.path(dat$out, "summaryPathsPlot.pdf"))
        plot(samp_plot)
        dev.off()
        if(sum(sample_paths_corr)>0){
        pdf(file.path(dat$out, "summaryCovPlot.pdf"))
        plot(samp_plot_cov)
        dev.off()}
      }
      
    } else {
      samp_plot <- NULL
      samp_plot_cov <- NULL
    }
    
    # 8.13.22 kad: Create df for paths set to 0 by user if applicable
    zero.paths.df <- NULL
    if(!is.null(dat$zero.paths)){
      for(path in dat$zero.paths){
        # Set values for vars that exist in coefs table
        path.split <- strsplit(path, "~")[[1]]
        lhs <- path.split[1]; rhs <- path.split[2]; op <- "~"
        est.std <- se <- ci.lower <- ci.upper <- 0
        z <- pvalue <- id <- level <- color <- NA
        param <- path
        
        # Combine, replicate for each person, set id names
        row <- data.frame(lhs,op,rhs,est.std,se,z,pvalue,ci.lower,ci.upper,id,param,level,color)
        df <- data.frame(lapply(row, rep, length(names(store$coefs))))
        df$id <- names(store$coefs)
        
        # Update df
        zero.paths.df <- rbind(zero.paths.df,df)
      }
    }
    
    # 8.13.22 kad: Combine paths set to 0 with regular coefs for output
    coefs <- rbind(summarize$coefs,zero.paths.df)
    
    indiv_paths     <- coefs[, c("id", "lhs", "op", "rhs", "est.std", 
                                 "se", "z", "pvalue", "level")]
    indiv_paths$lhs <- recode.vars(indiv_paths$lhs, dat$lvarnames, dat$varnames)
    indiv_paths$rhs <- recode.vars(indiv_paths$rhs, dat$lvarnames, dat$varnames)
    indiv_paths     <- indiv_paths[order(indiv_paths$id, indiv_paths$level), ]
    colnames(indiv_paths) <- c("file", "lhs","op", "rhs", "beta", "se", 
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
    
    # Write path counts matrix to output
    if (!is.null(dat$out)){ #& length(coefs[,1]) > 0){ # commented out stl 11.20.17
      write.csv(indiv_paths, file.path(dat$out, "indivPathEstimates.csv"),
                row.names = FALSE)
      write.csv(sample_counts, file.path(dat$out,
                                         "summaryPathCountsMatrix.csv"),
                row.names = FALSE)
      
      ### If hybrid is true or VAR is true, also write output for covariance counts
      if(dat$hybrid|dat$VAR){
      write.csv(sample_counts_corr, file.path(dat$out,
                                         "summaryCovCountsMatrix.csv"),
                row.names = FALSE)
      }
      
      # 6.19.21 kad: if HRF estimates have been calculated from convolved vars,
      # output these as individual files in the individual directory
      if(!is.null(dat$rf_est)){
        for(k in 1:dat$n_subj){
          rf_indiv <- dat$rf_est[[k]]
          write.csv(rf_indiv, file.path(dat$ind_dir, 
                                         paste0(dat$file_order[k,2], 
                                                "EstRF.csv")), row.names = TRUE)
        }
      }
      
      write.csv(fits, file.path(dat$out, "summaryFit.csv"), row.names = FALSE)
      if (dat$subgroup)
      write.table(sub$sim, file.path(dat$out, "similarityMatrix.csv"), sep = ",", col.names = FALSE, row.names = FALSE)
    }
    
  } else {
    # 8.13.22 kad: Create df for paths set to 0 by user if applicable
    zero.paths.df <- NULL
    if(!is.null(dat$zero.paths)){
      for(path in dat$zero.paths){
        # Set values for vars that exist in coefs table
        path.split <- strsplit(path, "~")[[1]]
        lhs <- path.split[1]; rhs <- path.split[2]; op <- "~"
        est.std <- se <- ci.lower <- ci.upper <- 0
        z <- pvalue <- NA
        
        # Combine and update df
        row <- data.frame(lhs,op,rhs,est.std,se,z,pvalue,ci.lower,ci.upper)
        zero.paths.df <- rbind(zero.paths.df,row)
      }
    }
    
    indiv_paths <- store$coefs[[1L]]
    indiv_paths <- rbind(indiv_paths,zero.paths.df) # 8.13.22 kad: Combine paths set to 0 with regular coefs for output
    indiv_paths$file <- "all"
    indiv_paths$lhs  <- recode.vars(indiv_paths$lhs, dat$lvarnames, dat$varnames)
    indiv_paths$rhs  <- recode.vars(indiv_paths$rhs, dat$lvarnames, dat$varnames)
    indiv_paths      <- indiv_paths[order(indiv_paths$file), ]
    indiv_paths      <- indiv_paths[ ,c("file", "lhs","op", "rhs", "est.std", 
                                        "se", "z", "pvalue")]
    colnames(indiv_paths) <- c("file", "lhs", "op", "rhs", "beta", "se", "z", "pval")
    
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
    sample_counts_corr <- NULL
    samp_plot_cov <- NULL
    samp_plot     <- NULL
  }
  
  dx <- list()
  if(diagnos){
    dx[[1]]<- dat
    dx[[2]] <- grp
    dx[[3]] <- store
    names(dx) <- c("dat", "grp", "store")
    }

    
  res <- list(fit           = fits,
              param_est     = indiv_paths,
              samp_plot     = samp_plot,
              samp_plot_cov = samp_plot_cov,
              sub_plots     = summarize$sub_plots,
              sub_plots_cov  = summarize$sub_plots_cov, 
              sample_counts = sample_counts,
              sample_counts_cov =    sample_counts_corr,
              sub_counts    = summarize$sub_counts,
              sub_counts_cov = summarize$sub_counts_cov,
              dx)  
  return(res)
  
}
