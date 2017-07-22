#'@method plot aggSEMp
#'@export
plot.aggSEMp <- function(x, ...){
    plot(x$g)
    invisible(x$g)
}

#'@method print aggSEMp
#'@export
print.aggSEMp <- function(x, estimates = FALSE, fitMeasures = FALSE, ...){
    if (estimates == TRUE){
      cat("Coefficients for final model", "\n")
      print(x$e, row.names = F)
      invisible(x$e)
    } else if (estimates == FALSE){
      ind <- x$a
      colnames(ind) <- x$b
      rownames(ind) <- x$b
      ind <- round(ind, digits = 2)
      ind_lag <- ind[(x$c+1):(x$c*2), 1:x$c]
      ind_con <- ind[(x$c+1):(x$c*2), (x$c+1):(x$c*2)]
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
      print.data.frame(x$d, row.names = F)
    }
  }
  