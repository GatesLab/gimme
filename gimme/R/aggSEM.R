#' @name aggSEM
#' @aliases aggSEM
#' @title Group-level structural equation model search.
#' @description Concatenates all individual-level data files and fits a group model to the data.
#' @usage
#' aggSEM(data   = "",
#'        out    = "",
#'        sep    = "",
#'        header = ,
#'        ar     = TRUE,
#'        plot   = TRUE,
#'        paths  = NULL)
#' @param data The path to the directory where the data files are located, or the name of the
#' list containing each individual's time series. Each file or matrix must contain one matrix 
#' for each individual containing a T (time) by 
#' p (number of variables) matrix where the columns represent variables and
#' the rows represent time. 
#' @param out The path to the directory where the results will be stored (optional). If specified,
#' a copy of output files will be replaced in directory. If directory at specified path does not
#' exist, it will be created.
#' @param sep The spacing of the data files. "" indicates space-delimited, 
#' "/t" indicates tab-delimited, "," indicates comma delimited. Only necessary to specify
#' if reading data in from physical directory.
#' @param header Logical. Indicate TRUE for data files with a header. Only necessary to specify
#' if reading data in from physical directory.
#' @param ar Logical. If TRUE, begins search for group model with autoregressive (AR) paths open. Defaults
#' to TRUE.
#' @param plot Logical. If TRUE, graphs depicting relations among variables of interest will automatically
#' be created. For aggregate-level plot, red paths represent positive weights and blue paths represent negative weights. Defaults to TRUE.
#' @param paths \code{lavaan}-style syntax containing paths with which to begin model estimation. That is, Y~X indicates that Y is regressed on X, or X predicts Y. If no header is used,
#' then variables should be referred to with V followed (with no separation) by the column number. If a
#' header is used, variables should be referred to using variable names. To reference lag variables, "lag"
#' should be added to the end of the variable name with no separation. Defaults to NULL.
#' @details
#'  In main output directory:
#'  \itemize{
#'  \item{\strong{allBetas}} Matrix. Contains estimates for each path in the aggregate-level model. The row variable is the outcome and the column variable is the predictor variable.
#'  \item{\strong{allStdErrors}} Matrix. Contains standard errors for each path in the aggregate-level model. The row variable is the outcome and the column variable is the predictor variable.
#'  \item{\strong{allPathEstimates}} {Contains estimate, standard error, p-value, and z-value for each path for the concatenated data.}
#'  \item{\strong{summaryFit}} {Contains model fit information for the aggregate-level model.}
#'  \item{\strong{summaryPathsPlot}} Contains aggregate-level plot. Red paths represent positive weights and blue paths represent negative weights.
#' }
#' @author Stephanie Lane
#' @examples
#' fit <- aggSEM(data = ts)
#' print(fit, fitMeasures = TRUE)
#' @export
aggSEM <- function(data,
                   out = NULL,
                   sep = NULL,
                   header = NULL,
                   ar    = TRUE,
                   plot  = TRUE,
                   paths = NULL){

  agg.internal <- function(setup.out){

  subjects = setup.out$subjects
  varnames = setup.out$varnames
  syntax   = setup.out$syntax
  plot     = setup.out$plot
  # files    = list.files(setup.out$data,full.names=TRUE)
  header   = setup.out$header
  sep      = setup.out$sep
  data_list  = setup.out$ts_list
  subgroup = setup.out$subgroup
  agg      = setup.out$agg
  plot     = setup.out$plot

  # data.all <- data.frame()
  # for (k in 1:subjects){
  #   data.file <- read.data(files[k],
  #                          header = header,
  #                          sep    = sep)
  #   data.all  <- rbind(data.all,data.file)
  # }
  
  data.all <- do.call(rbind, data_list)
  
  colnames(data.all) <- c(varnames)

  addind.out <- addind(done            = 0,
                       evaluate        = 1,
                       syntax          = syntax,
                       data.file       = data.all,
                       setup.out       = setup.out)

  evalind.out <- evalind(addind.out = addind.out,
                         setup.out  = setup.out,
                         data.file  = data.all)

  fixfitind.out <- fixfitind(setup.out   = setup.out,
                             evalind.out = evalind.out,
                             data.file   = data.all)

  final.fit.out <- final.fit(setup.out     = setup.out,
                             fixfitind.out = fixfitind.out,
                             data.file     = data.all,
                             k             = 1)

  all_ind_paths <- final.fit.out$ind.paths
  all_elem <- final.fit.out$ind.elements
  all_fit  <- as.matrix(final.fit.out$ind.fit)
  all_fit[,1]  <- "all"
  colnames(all_fit) <- c("subject", "chisq", "df", "pval",
                         "rmsea", "srmr", "nnfi", "cfi", "status")
  row.names(all_elem) <- NULL
  
  ind_plot <- final.fit.out$ind_plot
  
  list <- list("all_elem" = all_elem,
               "all_fit"      = all_fit,
               "all_syntax"   = NULL,
               "all_diff_sub" = FALSE,
               "all_ind_paths" = all_ind_paths,
               "ind_plot"      = ind_plot)
  return(list)
}

  setup.out <- setup(data     = data,
                     sep      = sep,
                     header   = header,
                     out      = out,
                     plot     = plot,
                     ar       = ar,
                     paths    = paths,
                     groupcutoff = NULL,
                     subcutoff = NULL,
                     subgroup = FALSE,
                     ind      = FALSE,
                     agg      = TRUE)

  agg.internal.out <- agg.internal(setup.out = setup.out)

  wrapup.out <- wrapup(indsem.internal.out = agg.internal.out,
                       setup.out           = setup.out)

  print.gimme.aggSEM(z=setup.out)
  
  final <- list(a = agg.internal.out$all_ind_paths,
                b = setup.out$varnames,
                c = setup.out$rois,
                d = wrapup.out$fit,
                e = wrapup.out$param_est,
                f = NULL,
                g = agg.internal.out$ind_plot,
                h = NULL,
                i = NULL)

  class(final) <- "aggSEMp"
  
  invisible(final)
  
}

print.gimme.aggSEM <- function(z){
  writeLines("aggSEM finished running normally")
  if (!is.null(z$out)) writeLines(paste("output is stored in", z$out))
}
################################################################################
################################################################################
