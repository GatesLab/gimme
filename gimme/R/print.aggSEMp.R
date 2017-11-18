#'@method plot aggSEMp
#'@export
plot.aggSEMp <- function(x, ...){
    plot(x$plot)
    invisible(x$plot)
}

#'@method print aggSEMp
#'@export
print.aggSEMp <- function(x, estimates = FALSE, fitMeasures = FALSE, ...){
    if (estimates == TRUE){
      cat("Coefficients for final model", "\n")
      print(x$path_se_est, row.names = F)
      invisible(x$path_se_est)
    } else if (estimates == FALSE){
      ind <- x$path_est_mat
      colnames(ind) <- x$path_est_mat
      rownames(ind) <- x$path_est_mat
      ind <- round(ind, digits = 2)
      ind_lag <- ind[(x$n_rois+1):(x$n_rois*2), 1:x$n_rois]
      ind_con <- ind[(x$n_rois+1):(x$n_rois*2), (x$n_rois+1):(x$n_rois*2)]
      cat("\n")
      cat("Lagged Matrix for all", "\n")
      print(ind_lag)
      cat("\n")
      cat("Contemporaneous Matrix for all", "\n")
      print(ind_con)
      invisible(ind)
    }
    if (fitMeasures == TRUE){
      cat("Fit for all", "\n")
      print.data.frame(x$fit, row.names = F)
    }
  }
  