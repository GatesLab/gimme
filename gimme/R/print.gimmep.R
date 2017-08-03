#'@method plot gimmep
#'@export
plot.gimmep <- function(x, file = NULL, subgroup = NULL, ...){
  if (is.null(file) & is.null(subgroup)) {
    cat("Please specify a file id for individual plots. Otherwise, summary plot is presented.")
    plot(x$g)
    invisible(x$g)
  } else if (is.null(subgroup)){
    a <- x$f[[file]]
    plot(a)
    invisible(a)
  } else if (is.null(file)){
    ## insert code here to grab subgroup plot
    a <- x$h[[subgroup]]
    if (!is.list(a)){
      cat("Subgroup", subgroup, "contains one individual. No subgroup plot provided.") 
    } else {
      plot(a)
      invisible(a)
    }
  }
}

#'@method print gimmep
#'@export
print.gimmep <- function(x, file = NULL, subgroup = NULL, 
                         mean = FALSE, estimates = FALSE, fitMeasures = FALSE, ...){
  if (!is.null(file)){
    if (estimates == TRUE){
      ind <- x$e[x$e$file %in% file, ]
      ind <- ind[, !names(ind) %in% c("subgroup")]
      ind[ ,4:6] <- round(ind[ ,4:6], digits = 3)
      cat("Coefficients for", file, "\n")
      print(ind, row.names = F)
      invisible(ind)
    } else if (estimates == FALSE){
      ind <- x$a[[file]]
      colnames(ind) <- x$b
      rownames(ind) <- x$b
      ind <- round(ind, digits = 2)
      ind_lag <- ind[(x$c+1):(x$c*2), 1:x$c]
      ind_con <- ind[(x$c+1):(x$c*2), (x$c+1):(x$c*2)]
      cat("\n")
      cat("Lagged Matrix for", file, "\n")
      print(ind_lag)
      cat("\n")
      cat("Contemporaneous Matrix for", file, "\n")
      print(ind_con)
      invisible(ind)
    }
    if (fitMeasures == TRUE){
      allfit <- x$d[, !names(x$d) %in% c("modularity")]
      ind    <- allfit[allfit$file %in% file, ]
      cat("Fit for file", file, "\n")
      print.data.frame(ind, row.names = F)
    }
  }
  if (!is.null(subgroup)){
    if (mean == TRUE){
      subidx   <- x$d[,c("file", "subgroup")]
      files    <- subidx[subidx$subgroup == subgroup, ]$file
      subfiles <- x$a[files]
      if (length(files) == 1){
      cat("Subgroup", subgroup, "contains one individual. No average subgroup matrix provided.") 
      } else {
      s        <- apply(simplify2array(subfiles), 1:2, mean, na.rm = TRUE)
      colnames(s) <- x$b
      rownames(s) <- x$b
      s     <- round(s, digits = 2)
      s     <- s[(x$c+1):(x$c*2), ]
      s_lag <- s[ , 1:x$c]
      s_con <- s[ , (x$c+1):(x$c*2)]
      cat("\n")
      cat("Lagged Average Matrix for Subgroup", subgroup, "\n")
      print(s_lag)
      cat("\n")
      cat("Contemporaneous Average Matrix for Subgroup", subgroup, "\n")
      print(s_con)
      invisible(s)
      }
    } else if (mean == FALSE & estimates == TRUE){
      # INSERT ESTIMATES FOR THAT SUBGROUP
      subest <- x$e
      subest <- subest[subest$subgroup == subgroup, ]
      subest[ ,4:6] <- round(subest[ ,4:6], digits = 3)
      cat("Coefficients for individuals in subgroup", subgroup, "\n")
      print(subest, row.names = F)
      invisible(subest)
    } else if (mean == FALSE & estimates == FALSE){
        sub <- x$k[[subgroup]]
        if (is.null(sub)){
          cat("Subgroup", subgroup, "contains one individual. No subgroup matrix provided.") 
        } else {
          sub_lag <- sub[ , 1:x$c]
          sub_con <- sub[ , (x$c+1):(x$c*2)]
          cat("\n")
          cat("Lagged Count Matrix for Subgroup", subgroup, "\n")
          print(sub_lag)
          cat("\n")
          cat("Contemporaneous Count Matrix for Subgroup", subgroup, "\n")
          print(sub_con)
          invisible(sub)
        }
    }
    if (fitMeasures == TRUE){
      # RETURN FITMEASURES FOR MEMBERS OF THAT SUBGROUP
      allfit <- x$d[, !names(x$d) %in% c("modularity")]
      allfit[allfit$subgroup == subgroup, ]
      cat("Fit for individuals in subgroup", subgroup, "\n")
      print.data.frame(allfit, row.names = F)
    }
  }
  if (is.null(file) & is.null(subgroup)) {
    if (mean == FALSE & estimates == FALSE){
      cat("Please specify a file id for individual coefficient matrix. ", "\n", 
          "Otherwise, a summary count matrix is presented below.", "\n")
      all <- x$j
      all_lag <- all[ , 1:x$c]
      all_con <- all[ , (x$c+1):(x$c*2)]
      cat("\n")
      cat("Lagged Count Matrix for Sample", "\n")
      print(all_lag)
      cat("\n")
      cat("Contemporaneous Count Matrix for Sample", "\n")
      print(all_con)
      invisible(all)
    } else if (mean == FALSE & estimates == TRUE){
      print(x$e) 
      invisible(x$e)
    }
    if (mean == TRUE){
      cat("Please specify a file id for individual coefficient matrix. ", "\n", 
          "Otherwise, a summary average matrix is presented below.", "\n")
      all2 <- apply(simplify2array(x$a), 1:2, mean, na.rm = TRUE)
      colnames(all2) <- x$b
      rownames(all2) <- x$b
      all2 <- round(all2, digits = 2)
      all2_lag <- all2[(x$c+1):(x$c*2), 1:x$c]
      all2_con <- all2[(x$c+1):(x$c*2), (x$c+1):(x$c*2)]
      cat("\n")
      cat("Lagged Average Matrix for Sample", "\n")
      print(all2_lag)
      cat("\n")
      cat("Contemporaneous Average Matrix for Sample", "\n")
      print(all2_con)
    }
    if (fitMeasures == TRUE){
      allfit <- x$d[, !names(x$d) %in% c("modularity")]
      cat("Fit for sample", "\n")
      print.data.frame(allfit, row.names = F)
      invisible(allfit)
    }
  }
}
  # else if (length(file) > 1){
  #   ind <- x$a[file]
  #   for (i in 1:length(file)){
  #     w <- ind[[i]]
  #     colnames(w) <- x$b
  #     rownames(w) <- x$b
  #     w <- round(w, digits = 2)
  #     w_lag <- w[(x$c+1):(x$c*2), 1:x$c]
  #     w_con <- w[(x$c+1):(x$c*2), (x$c+1):(x$c*2)]
  #     cat("\n")
  #     cat("Lagged Coefficient Matrix for", file[i], "\n")
  #     print(w_lag)
  #     cat("\n")
  #     cat("Contemporaneous Coefficient Matrix for", file[i], "\n")
  #     print(w_con)
  #   }
  #   invisible(w)
  # } 




