#' @name aggSEM
#' @aliases aggSEM
#' @title Group-level structural equation model search.
#' @description Concatenates all individual-level data files and 
#' fits a group model to the data.
#' @usage
#' aggSEM(data   = "",
#'        out    = "",
#'        sep    = "",
#'        header = ,
#'        ar     = TRUE,
#'        plot   = TRUE,
#'        paths  = NULL,
#'        exogenous = NULL)
#' @param data The path to the directory where the data files are located, 
#' or the name of the list containing each individual's time series. 
#' Each file or matrix must contain one matrix 
#' for each individual containing a T (time) by 
#' p (number of variables) matrix where the columns represent variables and
#' the rows represent time. If in list form, each item in the list (i.e., matrix) must be named.
#' @param out The path to the directory where the results will 
#' be stored (optional). If specified, a copy of output files will be replaced
#' in directory. If directory at specified path does not
#' exist, it will be created.
#' @param sep The spacing of the data files when data are in a directory. "" indicates space-delimited, 
#' "/t" indicates tab-delimited, "," indicates comma delimited. 
#' Only necessary to specify if reading data in from physical directory.
#' @param header Logical. Indicate TRUE for data files with a header. 
#' Only necessary to specify if reading data in from physical directory.
#' @param ar Logical. If TRUE, begins search for group model with 
#' autoregressive (AR) paths open. Defaults to TRUE.
#' @param plot Logical. If TRUE, figures depicting relations among variables 
#' of interest will automatically be created. For aggregate-level plot, red 
#' paths represent positive weights and blue paths represent negative weights. 
#' Dashed lines denote lagged relations (lag 1) and solid lines are contemporaneous 
#' (lag 0). Defaults to TRUE.
#' @param paths \code{lavaan}-style syntax containing paths with which to begin
#' model estimation (optional). That is, Y~X indicates that Y is regressed on X, or X
#' predicts Y. If no header is used, then variables should be referred to
#' with V followed (with no separation) by the column number. If a
#' header is used, variables should be referred to using variable names. 
#' To reference lag variables, "lag" should be added to the end of the variable 
#' name with no separation. Defaults to NULL.
#' @param exogenous Vector of variable names to be treated as exogenous (optional).  
#' That is, exogenous variable X can predict Y  but cannot be predicted by Y.  
#' If no header is used, then variables should be referred to with V followed 
#' (with no separation) by the column number. If a header is used, variables 
#' should be referred to using variable names.  Defaults to NULL.
#' @details
#'  In main output directory:
#'  \itemize{
#'  \item{\strong{allBetas}} Matrix. Contains estimates for each path in the
#'  aggregate-level model. The row variable is the outcome and the column 
#'  variable is the predictor variable.
#'  \item{\strong{allStdErrors}} Matrix. Contains standard errors for each path 
#'  in the aggregate-level model. The row variable is the outcome and the column
#'  variable is the predictor variable.
#'  \item{\strong{allPathEstimates}} {Contains estimate, standard error,
#'  p-value, and z-value for each path for the concatenated data.}
#'  \item{\strong{summaryFit}} {Contains model fit information for the 
#'  aggregate-level model.}
#'  \item{\strong{summaryPathsPlot}} Contains aggregate-level plot. Red paths 
#'  represent positive weights and blue paths represent negative weights.
#' }
#' @author Stephanie Lane
#' @examples
#' \dontrun{
#' exFit <- aggSEM(data = ts)
#' }
#' \dontshow{
#' load(system.file("extdata", "sysdata.rda", package = "gimme"))
#' }
#' plot(exFit)
#' @export

aggSEM <- function(data,
                   out    = NULL,
                   sep    = NULL,
                   header = NULL,
                   ar     = TRUE,
                   plot   = TRUE,
                   paths  = NULL,
                   exogenous = NULL){
  
  ind      = NULL # appease CRAN check
  grp      = NULL # appease CRAN check
  sub_spec = NULL # appease CRAN check

  
  dat  <- setup(data        = data,
                sep         = sep,
                header      = header,
                out         = out,
                plot        = plot,
                ar          = ar,
                paths       = paths,
                exogenous   = exogenous,
                groupcutoff = NULL,
                subcutoff   = NULL,
                subgroup    = FALSE,
                ind         = FALSE,
                agg         = TRUE)

  store <- indiv.search(dat, 
                        grp = NULL, 
                        ind)

  final <- final.org(dat, 
                     grp, 
                     ind = store$ind, 
                     sub, 
                     sub_spec, 
                     store)
  
  res <- list(path_est_mat = store$betas, 
              varnames     = dat$varnames, 
              n_rois       = dat$n_rois,
              fit          = final$fit, 
              path_se_est  = final$param_est, 
              plot         = store$plots, 
              vcov         = store$vcov)
  
  print.gimme.aggSEM(z = dat)
  
  class(res) <- "aggSEMp"
  
  invisible(res)
}

print.gimme.aggSEM <- function(z){
  writeLines("aggSEM finished running normally")
  if (!is.null(z$out)) writeLines(paste("output is stored in", z$out))
}
