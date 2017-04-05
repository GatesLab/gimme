#' @name gimmeSEM
#' @aliases gimme gimmeSEM
#' @title Group iterative multiple model estimation.
#' @description This function identifies structural equation models for each 
#' individual that consist of both group-level and individual-level paths.
#' @usage
#' gimmeSEM(data     = "",
#'          out      = "",
#'          sep      = "",
#'          header   = ,
#'          ar       = TRUE,
#'          plot     = TRUE,
#'          subgroup = FALSE,
#'          paths    = NULL,
#'          groupcutoff = .75,
#'          subcutoff   = .5,
#'          diagnos  = FALSE)
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
#' @param ar Logical. If TRUE, begins search for group model with 
#' autoregressive (AR) paths open. Defaults to TRUE.
#' @param paths \code{lavaan}-style syntax containing paths with which
#' to begin model estimation. That is, Y~X indicates that Y 
#' is regressed on X, or X predicts Y. If no header is used,
#' then variables should be referred to with V followed (with no separation) 
#' by the column number. If a
#' header is used, variables should be referred to using variable names. 
#' To reference lag variables, "lag" should be added to the end of the variable 
#' name with no separation. Defaults to NULL.
#' @param plot Logical. If TRUE, graphs depicting relations among variables 
#' of interest will automatically be
#' created. For individual-level plots, red paths represent positive weights 
#' and blue paths represent negative weights. 
#' For the group-level plot, black represents group-level paths, grey represents 
#' individual-level paths, and (if subgroup = TRUE) 
#' green represents subgroup-level paths. For the group-level plot, 
#' the width of the edge corresponds to the count. Defaults to TRUE.
#' @param subgroup Logical. If TRUE, subgroups are generated based on similarities 
#' in model features using the \code{walktrap.community} 
#' function from the \code{igraph} package.
#' @param groupcutoff Cutoff value for group- level paths. Defaults to .75, 
#' indicating that a path must be significant across 75\% of individuals to be 
#' included as a group-level path.
#' @param subcutoff Cutoff value for subgroup- level paths. Defaults to .5, 
#' indicating that a path must be significant across at least 50\% of the
#' individuals in a subgroup to be considered a subgroup-level path.
#' @param diagnos In development. Defaults to FALSE.
#' @details
#'  In main output directory:
#'  \itemize{
#'  \item{\strong{indivPathEstimates}} {Contains estimate, standard error, 
#'  p-value, and z-value for each path for each individual. 
#'  If subgroup = TRUE and subgroups are found, then a column is present containing the 
#'  subgroup membership for each individual. Also contains the level at which each path 
#'  was estimated: group, subgroup, or individual.}
#'  \item{\strong{summaryFit}} {Contains model fit information for individual-level models. 
#'  If subgroups are requested, this file also contains the subgroup membership 
#'  for each individual.}
#'  \item{\strong{summaryPathCountMatrix}} Contains counts of total number of paths, 
#'  both contemporaneous and lagged, estimated for the sample. The row variable is the 
#'  outcome and the column variable is the predictor variable.
#'  \item{\strong{summaryPathCounts}} {Contains summary count information for paths 
#'  identified at the group-, subgroup (if subgroup = TRUE), and individual-level.}
#'  \item{\strong{summaryPathsPlot}} {Produced if plot = TRUE. Contains figure with group, 
#'  subgroup (if subgroup = TRUE), and individual-level paths 
#'  for the sample. Black paths are group-level, green paths are subgroup-level, 
#'  and grey paths are individual-level, where the thickness of the line
#'  represents the count.}
#'  }
#'  In subgroup output directory (if subgroup = TRUE):
#'  \itemize{
#'  \item{\strong{subgroup\emph{k}PathCounts}} Contains counts of relations 
#'  among lagged and contemporaneous variables for the \strong{\emph{k}}th subgroup.
#'  \item{\strong{subgroup\emph{k}Plot}} Contains plot of group, subgroup, 
#'  and individual level paths for the \strong{\emph{k}}th subgroup. 
#'  Black represents group-level paths, grey represents individual-level paths, 
#'  and green represents subgroup-level paths.
#'  }
#'  Note: if a subgroup of size n = 1 is discovered, subgroup-level output is not produced. \cr

#'  In individual output directory (where \strong{\emph{id}} represents the 
#'  original file name for each individual):
#'  \itemize{
#'  \item{\strong{\emph{id}Betas}} Contains individual-level estimates
#'   of each path for each individual.
#'  \item{\strong{\emph{id}StdErrors}} Contains individual-level standard errors 
#'  for each path for each individual.
#'  \item{\strong{\emph{id}Plot}} Contains individual-level plots. Red paths 
#'  represent positive weights and blue paths represent negative weights.
#' }
#' @references Gates, K.M. & Molenaar, P.C.M. (2012). Group search algorithm 
#' recovers effective connectivity maps for individuals 
#' in homogeneous and heterogeneous samples. NeuroImage, 63, 310-319.
#' @author Stephanie Lane
#' @examples
#'  \dontrun{
#' paths <- 'V2 ~ V1
#'           V3 ~ V4lag'
#'
#' fit <- gimmeSEM(data     = simData,
#'                 out      = "C:/simData_out",
#'                 subgroup = TRUE)
#'                 
#' print(fit, mean = TRUE)
#' print(fit, subgroup = 1, mean = TRUE)
#' print(fit, file = "group_1_1", estimates = TRUE)
#' print(fit, subgroup = 2, fitMeasures = TRUE)
#' plot(fit, file = "group_1_1")
#'  }
#' @keywords gimmeSEM
#' @export gimme gimmeSEM
gimmeSEM <- gimme <- function(data           = NULL,
                              out            = NULL,
                              sep            = NULL,
                              header         = NULL,
                              ar             = TRUE,
                              plot           = TRUE,
                              subgroup       = FALSE,
                              paths          = NULL,
                              groupcutoff    = .75,
                              subcutoff      = .5,
                              diagnos        = FALSE){
  
  setup.out         <- setup(data                 = data,
                             sep                  = sep,
                             header               = header,
                             out                  = out,
                             plot                 = plot,
                             ar                   = ar,
                             paths                = paths,
                             subgroup             = subgroup,
                             ind                  = FALSE,
                             agg                  = FALSE,
                             groupcutoff          = groupcutoff,
                             subcutoff            = subcutoff)
  
  if (diagnos == TRUE){
    dir.create(file.path(out, "diagnos"))
    saveRDS(setup.out, file.path(out, "diagnos", "01_setup.RDS"))
  }
  
  ## the miSEM function is responsible for the search procedure
  ## it is used at both the group and subgroup levels
  miSEM.out         <- miSEM(setup.out            = setup.out,
                             previous.out         = setup.out,
                             subgroup.step        = FALSE,
                             subjects             = NULL,
                             ts_list              = NULL,
                             syntax               = NULL,
                             subgroup_paths       = NULL,
                             second_round         = FALSE)
  
  if (diagnos == TRUE){
    saveRDS(miSEM.out, file.path(out, "diagnos", "02_miSEM.RDS"))
  }
  
  ## the evalbetas function is responsible for pruning paths 
  ## which may have become nonsignificant for the majority of paths
  ## it is used at both the group and subgroup levels
  evalbetas.out <- evalbetas(setup.out            = setup.out,
                             previous.out         = miSEM.out,
                             subgroup.step        = FALSE,
                             s                    = NULL,
                             subjects             = NULL,
                             ts_list              = NULL,
                             post_sub_prune       = FALSE,
                             evalbetas.out        = NULL)
  ## begin subgroup steps ##
  if (diagnos == TRUE){
    saveRDS(evalbetas.out, file.path(out, "diagnos", "03_evalbetas.RDS"))
  }
  
  if (subgroup == TRUE){
    ## runs group-level model for each individual, 
    ## creates similarity matrix, obtains subgroup assignments
    subsetup.out         <- subsetup(setup.out    = setup.out,
                                     previous.out = evalbetas.out)
    
    if (diagnos == TRUE){
      saveRDS(subsetup.out, file.path(out, "diagnos", "04_subsetup.RDS"))
    }
    
    ## searches for paths common to each subgroup
    ## using same general procedure as group-level search, 
    ## but using different cutoff for what constitutes "majority"
    miSEMsub.out         <- miSEMsub(subsetup.out = subsetup.out,
                                     setup.out    = setup.out,
                                     previous.out = evalbetas.out,
                                     second_round = FALSE,
                                     evalbetassub.out = NULL)
    if (diagnos == TRUE){
      saveRDS(miSEMsub.out, file.path(out, "diagnos", "05_miSEMsub.RDS"))
    }
    ## prunes paths which may have become nonsignificant 
    ## for the majority of the subgroup
    evalbetassub.out <- evalbetassub(subsetup.out  = subsetup.out,
                                     previous.out  = miSEMsub.out,
                                     setup.out     = setup.out,
                                     evalbetas.out = evalbetas.out)
    
    if (diagnos == TRUE){
      saveRDS(evalbetassub.out, file.path(out, "diagnos", "06_evalbetassub.RDS"))
    }
    
    ## searches for paths to prune which may have become nonsignificant
    ## for the group-level model now that subgroup-level paths
    ## have been obtained
    prune.group.post.out <- evalbetas(setup.out      = setup.out,
                                      previous.out   = evalbetassub.out,
                                      subgroup.step  = FALSE,
                                      s              = NULL,
                                      subjects       = NULL,
                                      ts_list        = NULL,
                                      post_sub_prune = TRUE,
                                      evalbetas.out  = evalbetas.out)
    
    if (diagnos == TRUE){
      saveRDS(prune.group.post.out, file.path(out, "diagnos", "07_prune_group_post.RDS"))
    }
    
    ## only runs if a group-level path was removed in the previous step
    if (evalbetas.out$n_group_paths != prune.group.post.out$n_group_paths){
      ## searches for any final subgroup-level paths following the group-level pruning
      miSEMsub.round2.out <- miSEMsub(subsetup.out = subsetup.out,
                                      setup.out    = setup.out,
                                      previous.out = prune.group.post.out,
                                      second_round = TRUE,
                                      evalbetassub.out = evalbetassub.out)
      
      if (diagnos == TRUE){
        saveRDS(miSEMsub.round2.out, file.path(out, "diagnos", "08_miSEMsub_round2.RDS"))
      }
      
      ## only run another round of subgroup-level pruning if another subgroup-level path got added
      check.again <- any(miSEMsub.round2.out$unique_syntax[,4] != evalbetassub.out$sub.paths[,2])
      if (check.again == TRUE){
        evalbetassub.round2.out <- evalbetassub(subsetup.out  = subsetup.out,
                                                previous.out  = miSEMsub.round2.out,
                                                setup.out     = setup.out,
                                                evalbetas.out = evalbetas.out)
        # replace first round with second round results
        evalbetassub.out <- evalbetassub.round2.out
      } else {
        # if evalbetas isn't run again, grabs most recent syntax from miSEMsub.round2
        # for use in indsem
        combine_syntax                <- merge(evalbetassub.out$syntax_sub,
                                               miSEMsub.round2.out$unique_syntax,
                                               by.x = "membership", by.y = "V1")
        combine_syntax                <- as.matrix(combine_syntax[-c(3,5,6)])
        colnames(combine_syntax)      <- c("membership", "files", "syntax")
        # replace subgroup-level syntax from post-subgrouping-pruning
        # with subgroup-level syntax from round 2 subgroup-level search
        evalbetassub.out$syntax_sub   <- combine_syntax
        evalbetassub.out$final_syntax <- unique(combine_syntax[ ,c(1,3)])
      }
    }
    
    if (diagnos == TRUE){
      saveRDS(evalbetassub.out, file.path(out, "diagnos", "09_evalbetassub_final.RDS"))
    }
    
  }
  
  if (subgroup == FALSE) evalbetassub.out <- NULL
  
  ## this is the individual-level search, adds path for each individual
  ## runs one person at a time  indsem.internal.out <- list
  indsem.internal.out <- indsem.internal(setup.out        = setup.out,
                                         evalbetassub.out = evalbetassub.out,
                                         evalbetas.out    = evalbetas.out)
  
  if (diagnos == TRUE){
    saveRDS(indsem.internal.out, file.path(out, "diagnos", "10_indSEM_internal.RDS"))
  }
  
  ## grabs information from individual-level search (once complete) 
  ## prints summary output and makes tables
  wrapup.out <- wrapup(indsem.internal.out = indsem.internal.out,
                       setup.out           = setup.out)
  
  if (diagnos == TRUE){
    saveRDS(wrapup.out, file.path(out, "diagnos", "11_wrapup.RDS"))
  }
  
  print.gimme(x = subsetup.out,
              y = subgroup,
              z = setup.out)
  
  final <- list(a = indsem.internal.out$all_ind_paths, 
                b = setup.out$varnames, 
                c = setup.out$rois,
                d = wrapup.out$fit, 
                e = wrapup.out$param_est, 
                f = indsem.internal.out$all_plots, 
                g = wrapup.out$all_plot, 
                h = wrapup.out$sub_plots,
                i = setup.out$subgroup,
                j = wrapup.out$all_counts,
                k = wrapup.out$sub_paths,
                l = indsem.internal.out$vcov_params)
  
  class(final) <- "gimmep"
  
  invisible(final)
}

print.gimme <- function(x, y, z){
  writeLines("gimme finished running normally")
  if (!is.null(z$out)) writeLines(paste("output is stored in", z$out))
  if (y == TRUE) {
    writeLines(paste("Number of subgroups =", x$n_subgroups))
    writeLines(paste("Modularity =", round(x$modularity, digits = 5)))
  }
}

## this setup function creates many values that later code refers back to
setup <- function (data,
                   sep,
                   header,
                   out,
                   plot,
                   ar,
                   paths,
                   subgroup,
                   agg,
                   ind,
                   groupcutoff,
                   subcutoff) {
  
  ## check that data argument is specified
  if (is.null(data)){
    stop(paste0("gimme ERROR: neither a data directory nor a data list is specified. ",
                "Please either specify a directory of data or a list of individual data files."))
  }
  
  ## code to create list of individual data files if used from directory
  ## or create lagged version if list already provided by user
  if (!is.list(data)){
    files <- list.files(data, full.names = TRUE)
    ## add error messages if a data directory is specified but sep or header is missing
    if (is.null(sep)){
      stop(paste0("gimme ERROR: a data directory is specified but a sep argument is not. ",
                  "Please specify a sep argument before continuing."))
    }
    if (is.null(header)){
      stop(paste0("gimme ERROR: a data directory is specified but a header argument is not. ",
                  "Please specify a logical value for header before continuing."))
    }
    ## if above checks are passed, read in files and create list
    ts_list <- list()
    for (i in 1:length(files)){
      ts_ind       <- read.table(files[i], sep = sep, header = header)
      ts_list[[i]] <- ts_ind 
    }
    varnames       <- colnames(ts_ind)
    names(ts_list) <- tools::file_path_sans_ext(basename(files))
    rois           <- ncol(ts_ind)
  } else if (is.list(data)){
    ts_list  <- data
    varnames <- colnames(ts_list[[1]])
    rois     <- ncol(ts_list[[1]])
  }
  
  lvarnames <- character()
  # simplify creation of variable names
  varnames  <- c(paste0(varnames[1:rois], "lag"), varnames)
  lvarnames <- c(paste0("VAR", seq(1:rois), "lag"), paste0("VAR", seq(1:rois)))
  
  ## go back through list and create lagged variables
  for (p in 1:length(ts_list)){
    all          <- ts_list[[p]]
    first        <- all[1:(nrow(all)-1), ]
    second       <- all[2:nrow(all), ]
    ts_lc        <- data.frame(first, second)
    colnames(ts_lc) <- varnames
    ts_list[[p]] <- ts_lc
  }
  
  ## code to check whether output directory specified by the user actually exists
  ## if it doesn't exist, the code wouldn't crash until the individual-level search
  ## this checks and throws an error at the start
  
  # flag to add check if argument is null
  if (is.null(out)){
    cat("gimme MESSAGE: No output directory specified. All output should be directed to an object.", "\n")
  } else if (file_test(op = "-d", out) == FALSE) {
    cat("gimme MESSAGE: specified output directory doesn't exist. Attempting to create now.", "\n")
    dir.create(out, recursive = TRUE)
    if (dir.exists(out)){
      cat("gimme MESSAGE: output directory successfully created.", "\n")
    } else {
      stop("gimme ERROR: unable to create output directory. Please specify a different file path")
    }
  }
  
  subjects         <- length(ts_list)
  trackparts       <- matrix(0, nrow = subjects, ncol = 2)
  trackparts[,1]   <- seq(1:subjects)
  trackparts[,2]   <- names(ts_list)
  all              <- ts_list[[1]]
  # check to see if number of columns is equal to 1. 
  # if so, user likely misspecified the sep argument
  if (rois == 1) {
    stop(paste0("gimme ERROR: only one column of data read in. ",
                "Check if sep argument properly specified."))
  }
  ## check to see if all individuals have same number of columns. 
  ## otherwise, code will fail with ambiguous error later
  cols         <- numeric()
  missingCols  <- numeric()
  constantCols <- logical()
  numericCols  <- logical()
  for (k in 1:subjects){
    data.file <- ts_list[[k]]
    cols[k]   <- ncol(data.file)
    missingCols[k] <- sum(colSums(is.na(data.file)) < nrow(data.file))
    constantCols[k] <- any(apply(data.file, 2, sd, na.rm = TRUE) == 0)
    numericCols[k]  <- any(apply(data.file, 2, is.numeric) == FALSE)
  }
  if (subjects != 1) {
    if (sd(cols) != 0) {
      stop(paste0('gimme ERROR: not all data files have the same number of columns. ',
                  'Please fix or remove file before continuing.'))
    }
    if (sd(missingCols) != 0) {
      stop(paste0('gimme ERROR: at least one data file contains a column with all NA. ',
                  'Please fix or remove files listed below before continuing. \n', 
                  paste0(names(ts_list)[missingCols != cols], collapse = "\n")))
    }
    if (any(cols != missingCols)) {
      stop(paste0('gimme ERROR: at least one data file contains a column with all NA. ',
                  'Please fix or remove file before continuing.'))
    }  
    if (any(constantCols == TRUE)){
      stop(paste0('gimme ERROR: at least one data file contains a column with constant values. ',
                  'Please fix or remove files listed below before continuing. \n', 
                  paste0(names(ts_list)[constantCols == TRUE], collapse = "\n")))
    }
    if (any(numericCols == TRUE)){
      stop(paste0('gimme ERROR: at least one data file contains a column with non-numeric values. ',
                  'Please fix or remove files listed below before continuing. \n', 
                  paste0(names(ts_list)[numericCols == TRUE], collapse = "\n")))
    }
  } 
  if (subjects == 1 & ind == FALSE) {
    stop(paste0('gimme ERROR: only one subject detected in data directory. ',
                'Please use indSEM function instead.'))
  }
  vars              <- rois*2
  cutoffind         <- qchisq(.99, 1)
  n_group_paths     <- 0
  
  x                 <- seq(1:vars)
  y                 <- substring(lvarnames, 4)
  individual        <- file.path(out, "individual")
  subgroup_dir      <- file.path(out, "subgroup")
  fitted            <- file.path(out, "fitted")
  
  if (subgroup == TRUE & !is.null(out))  {
    dir.create(subgroup_dir, showWarnings = FALSE)
  }
  if (agg      == FALSE & !is.null(out)) {
    dir.create(individual, showWarnings = FALSE)
  }
  if (plot == TRUE) {
    plot.names <- varnames[(rois+1):(rois*2)]
  } else {
    plot.names <- ""
  }
  
  #------------------------------------------------------------------------------#
  # prepare paths if semigimme is specified
  if (!is.null(paths))
  {
    prep.paths <- function(paths, varnames, lvarnames)
    {
      table    <- lavParTable(paths)
      # only include paths in the syntax which are specified free by the user
      # allows the possibility for the user to fix certain paths to zero
      tableFree    <- table[table$op == "~" & table$free != 0, ]
      dvsFree      <- recoderFunc(tableFree$lhs,varnames,lvarnames)
      ivsFree      <- recoderFunc(tableFree$rhs,varnames,lvarnames)
      if (nrow(tableFree) != 0){
        vsFree       <- paste0(dvsFree, "~", ivsFree)
      } else vsFree = NULL
      # table up the paths which are fixed to a certain value by the user
      tableFixed   <- table[table$op == "~" & table$free == 0,]
      if (nrow(tableFixed) > 0){
        dvsFixed     <- recoderFunc(tableFixed$lhs,varnames,lvarnames)
        ivsFixed     <- recoderFunc(tableFixed$rhs,varnames,lvarnames)
        vsFixed      <- paste0(dvsFixed, "~", ivsFixed)
      } else {
        vsFixed = NULL
      }
      list = list(paths  = vsFree,
                  remove = vsFixed)
      return(list)
    }
    remove <- prep.paths(paths, varnames, lvarnames)$remove  
    paths  <- prep.paths(paths, varnames, lvarnames)$paths
  } else {
    remove <- NULL
  }
  #------------------------------------------------------------------------------#
  
  ## code below creates the starting null syntax file
  ## creation of syntax simplified Oct 2016
  
  # set up single indicator latent variables 
  line1 <- paste0(lvarnames, "=~1*", varnames, collapse = "\n") 
  
  # covary all exogenous lagged variables
  line2 <- paste(capture.output(for (i in 1:rois) {for (j in 1:i){
    cat(lvarnames[i],"~~",lvarnames[j],
        sep="","\n")}}),collapse="\n")
  
  # estimate variance of single-indicator latent variables 
  line3 <- paste0(lvarnames[(rois+1):vars], "~~", lvarnames[(rois+1):vars], collapse = "\n") 
  
  # freely estimate autoregressive relationships if ar = TRUE
  # if ar = FALSE, set up nonsense paths fixed to zero 
  # to open beta matrix for MI calculation in lavaan
  if (ar == TRUE) {
    line4 <- paste0(lvarnames[(rois+1):vars], "~", lvarnames[1:rois], collapse = "\n")
  } else {
    line4 <- paste0(lvarnames[1:rois], "~0*", lvarnames[(rois+1):vars], collapse = "\n")
  }
  
  if (!is.null(paths))
  {
    line5 <- paste0(paths, collapse = "\n") 
  }
  
  if (is.null(paths)){
    syntax            <- paste(line1 ,line2, line3, line4, sep = "\n")
  } else {
    syntax            <- paste(line1, line2, line3, line4, line5, sep = "\n")
  }
  ## end creating syntax
  
  ## create list of paths that make sense to gimme to open
  ## for example, it doesn't include paths where time would predict time-1
  candidate_paths   <- apply(expand.grid(lvarnames[(rois+1):vars], 
                                         lvarnames[1:vars]), 1, paste, collapse = "~")
  
  # if path specified by a user is fixed to a certain value, 
  # we want to remove it from consideration when looking at MIs, 
  # especially if that path was fixed to zero
  candidate_paths <- candidate_paths[!candidate_paths %in% remove]
  
  ## just creates list of AR paths so that later code doesn't kick them out
  ## code is designed so that if user selects ar = TRUE, they stay in
  ## even if they become nonsignificant for the majority
  ar.paths <- paste0(lvarnames[(rois+1):vars], "~", lvarnames[1:rois]) 
  
  # if user specifies paths, add them to the list of ar.paths
  # this ensures that they don't get kicked out
  if (!is.null(paths)){
    ar.paths <- c(ar.paths, paths) 
  }
  
  n_fixed_paths <- length(ar.paths)
  
  list <- list("subjects"          = subjects,
               "rois"              = rois,
               "x"                 = x,
               "y"                 = y,
               "varnames"          = varnames,
               "trackparts"        = trackparts,
               "vars"              = vars,
               "lvarnames"         = lvarnames,
               "cutoffind"         = cutoffind,
               "fitted"            = fitted,
               "individual"        = individual,
               "plot.names"        = plot.names,
               "syntax"            = syntax,
               "candidate_paths"   = candidate_paths,
               "ar.paths"          = ar.paths,
               "ar"                = ar,
               "out"               = out,
               "plot"              = plot,
               "n_group_paths"     = n_group_paths,
               "subgroup_dir"      = subgroup_dir,
               "n_fixed_paths"     = n_fixed_paths,
               "subgroup"          = subgroup,
               "agg"               = agg,
               "groupcutoff"       = groupcutoff,
               "subcutoff"         = subcutoff,
               "ts_list"           = ts_list)
  return(list)
}


## model for lavaan
fit.model <- function (varnames,
                       syntax,
                       data.file) {
  
  fit <- tryCatch(lavaan::lavaan(syntax,
                                 data            = data.file,
                                 model.type      = "sem",
                                 missing         = "fiml",
                                 estimator       = "ml",
                                 int.ov.free     = TRUE,
                                 int.lv.free     = FALSE,
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

## quick function that bypasses the need to read in data and write it back out
## just reads in file and creates the lagged vars
# read.data <- function (file,
#                        sep,
#                        header) {
#   all       <- as.matrix(read.table(file, sep = sep, header = header))
#   first     <- all[1:(nrow(all)-1), ]
#   second    <- all[2:nrow(all), ]
#   data.file <- data.frame(first, second)
#   return(data.file)
# }

count_excellent <- function(indices){
  rmsea     <- indices[4]
  srmr      <- indices[5]
  cfi       <- indices[6]
  nnfi      <- indices[7]
  rmseaE    <- ifelse(rmsea < .05, 1, 0)
  srmrE     <- ifelse(srmr  < .05, 1, 0)
  cfiE      <- ifelse(cfi   > .95, 1, 0)
  nnfiE     <- ifelse(nnfi  > .95, 1, 0)
  excellent <- sum(rmseaE, srmrE, cfiE, nnfiE)
  return(excellent)
}

## this function is responsible for the group-level search procedure
miSEM <- function (setup.out,
                   previous.out,
                   subgroup.step,
                   subjects,
                   ts_list,
                   syntax,
                   subgroup_paths,
                   second_round) {
  
  vars            <- setup.out$vars
  varnames        <- setup.out$varnames
  candidate_paths <- setup.out$candidate_paths
  ar              <- setup.out$ar
  groupcutoff     <- setup.out$groupcutoff
  subcutoff       <- setup.out$subcutoff
  
  n_group_paths <- previous.out$n_group_paths
  
  param = NULL
  Freq  = NULL
  
  mi.index            <- matrix(1:((vars-1)*vars), nrow = ((vars-1)*vars), ncol = 1)
  
  ## the miSEM function is used in several stages of gimme 
  ## it sets up some values depending on what stage it's in.
  ## if subgroup.step == FALSE, it means that it's the first round (group-level search).
  ## subgroup.step == TRUE & second_round == FALSE means it's using the function to find
  ## group-level paths within the subgroups, or in other words, it's the subgroup-level search
  ## subgroup.step == TRUE & second_round == TRUE means that the group-level search and subgroup-level
  ## search have already taken place, but that it's going back through to check for any additional
  ## subgroup-level paths. this only happens if group-level paths were pruned after the subgroup-level
  ## search took place the first go-around.
  if (subgroup.step == FALSE) {
    n_group_paths    <- 0
    n_subgroup_paths <- 0
    subjects         <- setup.out$subjects
    ts_list          <- setup.out$ts_list 
    syntax           <- previous.out$syntax
  }
  
  if (subgroup.step == TRUE & second_round == FALSE){
    n_subgroup_paths <- 0
    subgroup_paths   <- character()
  }
  
  if (subgroup.step == TRUE & second_round == TRUE){
    subgroup_paths   <- unlist(strsplit(subgroup_paths, "[,]"))
    n_subgroup_paths <- length(subgroup_paths)
  }
  
  cutoff               <- qchisq(.95, subjects)
  cutoffgroup          <- qchisq(1-.05/subjects, 1)
  continue <- 1
  
  while (continue == 1) {
    mi.list            <- matrix(0, nrow = vars*(vars-1)*subjects, ncol = 6)
    colnames(mi.list)  <- c("subject", "index", "lhs", "op", "rhs", "mi")
    count_converge     <- 0
    ## section below just iteratively fits model, grabs MIs, selects one significant
    ## for most people, adds to syntax, continues until no path sig for majority (75%)
    for (k in 1:subjects){
      # data.file           <- read.data(file   = files[k],
      #                                  sep    = sep,
      #                                  header = header)
      # colnames(data.file) <- c(varnames)
      fit                 <- fit.model(syntax    = syntax,
                                       data.file = ts_list[[k]])
      
      check_npd            <- any(grepl("error", class(fit)) == TRUE)
      if (check_npd == FALSE) {
        check_zero_se <- sum(lavInspect(fit,"se")$beta, na.rm = TRUE) == 0
      } else check_zero_se <- TRUE
      
      if (ar == FALSE & n_group_paths == 0){check_zero_se <- FALSE}
      
      if (check_npd == FALSE & check_zero_se == FALSE) {
        all_mi         <- tryCatch(modindices(fit), error = function(e) e)
        check_singular <- any(grepl("singular", all_mi) == TRUE)
        check_error    <- any(grepl("error", class(all_mi)) == TRUE)
        converge       <- lavInspect(fit, "converged")
      } else {
        if (subgroup.step == FALSE) {
          writeLines(paste("group-level search, subject", k, "nonconvergence"))
        } else {
          writeLines(paste("subgroup-level search, subject", k, "nonconvergence"))
        }
        check_singular <- TRUE
        converge       <- FALSE
        check_error    <- TRUE
      }
      
      if (check_singular == FALSE & converge == TRUE & check_zero_se == FALSE & 
          check_error == FALSE) {
        if (subgroup.step == FALSE) {
          writeLines(paste("group-level search, subject", k))
        } else {
          writeLines(paste("subgroup-level search, subject", k))
        }
        mi                <- as.matrix(all_mi[all_mi$op == "~",])[ ,c("lhs","op","rhs","mi")]
        padLength         <- ((vars-1)*vars) - nrow(mi)
        pad               <- matrix(NA, nrow = padLength, ncol = 4)
        mi                <- rbind(mi, pad)
        count_converge    <- count_converge + 1
      }
      # if it doesn't converge or is computationally singular or NPD
      if (converge == FALSE | check_singular == TRUE | check_npd == TRUE | check_zero_se == TRUE) {
        mi                <- matrix(NA, nrow = ((vars-1)*vars), ncol = 4)
      }
      #stacking matrices
      mi.subject          <- matrix(k, nrow = ((vars-1)*vars), ncol = 1)
      mi.list[(((nrow(mi)*k)-nrow(mi))+1):(nrow(mi)*k),(1:6)] <-
        as.matrix(cbind(mi.subject, mi.index, mi), 
                  rownames.force = FALSE)
    }
    
    ## this section just sorts the matrix of MIs and grabs the best one
    mi.all                <- as.data.frame(mi.list[complete.cases(mi.list),])
    mi.all$param          <- paste0(mi.all$lhs, mi.all$op, mi.all$rhs)
    mi.all                <- mi.all[-c(3:5)]
    mi.all                <- subset(mi.all, param %in% candidate_paths)
    mi.all[,3]            <- as.numeric(as.character(mi.all[,3]))
    mi.all[,1]            <- as.numeric(as.character(mi.all[,1]))
    mi.all                <- transform(mi.all, sum = ave(mi, param, FUN = sum))
    mi.high               <- subset(mi.all, mi > cutoffgroup)
    # make sure that at least one MI is above cutoff, or it will error
    if (nrow(mi.high) != 0){
      mi.count              <- subset(as.data.frame(table(mi.high$param)), Freq > 0)
      mi.high.count         <- subset(mi.high, !duplicated(param))
      mi.merge              <- merge(x = mi.high.count, y = mi.count,
                                     by.x = "param", by.y = "Var1")
      paramadd              <- mi.merge[order(-mi.merge$Freq, -mi.merge$sum),][1,1]
      ## list of every individual's MIs for the selected path
      y                     <- subset(mi.all, param == paramadd, select = mi)$mi
      drop_element          <- ifelse(max(y) < cutoffgroup, TRUE, FALSE)
      prop                  <- sum(y > cutoffgroup)/count_converge
      halt <- length(y) <= subjects/2
    } else {
      halt <- TRUE
      prop <- 0
      drop_element = TRUE
    }
    
    ## stops the search if less than half of subjects are terminating normally
    
    
    if (subgroup.step == TRUE)  cutoffprop <- subcutoff
    if (subgroup.step == FALSE) cutoffprop <- groupcutoff
    
    if (prop <= cutoffprop | drop_element == TRUE | halt == TRUE) {
      continue <-0
    } else {
      continue            <- 1
      syntax              <- paste(syntax, as.name(paramadd), sep = "\n")
      if (subgroup.step == FALSE) n_group_paths    <- n_group_paths + 1
      if (subgroup.step == TRUE) {
        n_subgroup_paths <- n_subgroup_paths + 1
        subgroup_paths       <- append(subgroup_paths,paramadd)
      }
    }
    
  }  #end of while continue>0
  list <- list("syntax"           = syntax,
               "n_group_paths"    = n_group_paths,
               "n_subgroup_paths" = n_subgroup_paths,
               "subgroup_paths"   = subgroup_paths)
  return(list)
}
## this function is used for pruning, both for group-level and subgroup-level
## it basically runs the models again and looks for paths that are now nonsignificant
## for the majority
evalbetas <- function (setup.out,
                       previous.out,
                       subgroup.step,
                       s,
                       subjects,
                       ts_list,
                       post_sub_prune,
                       evalbetas.out) {
  
  ar.paths          = setup.out$ar.paths
  ar                = setup.out$ar
  varnames          = setup.out$varnames
  rois              = setup.out$rois
  n_fixed_paths     = setup.out$n_fixed_paths
  groupcutoff       = setup.out$groupcutoff
  subcutoff         = setup.out$subcutoff
  
  op    = NULL
  sig   = NULL
  param = NULL
  z     = NULL
  
  ## flags are similar to miSEM
  if (subgroup.step == FALSE & post_sub_prune == FALSE) {
    bad              = ifelse(identical(setup.out$syntax, previous.out$syntax), 0, 1)
    syntax           = previous.out$syntax
    subjects         = setup.out$subjects
    n_subgroup_paths = 0
    ts_list          = setup.out$ts_list
    n_group_paths    = previous.out$n_group_paths
    if (n_group_paths == 0) bad = 0
  }
  
  if (subgroup.step == TRUE){
    bad              = ifelse(previous.out$unique_syntax[s,3] == 0, 0, 1)
    n_subgroup_paths = as.numeric(previous.out$unique_syntax[s,3])
    syntax           = previous.out$unique_syntax[s,2]
    subgroup_paths   = previous.out$unique_syntax[s,4]
    subgroup_paths   = unlist(strsplit(subgroup_paths, "[,]"))
    n_group_paths    = previous.out$n_group_paths
  }
  
  if (post_sub_prune == TRUE){
    ts_list     <- setup.out$ts_list
    bad         <- 1
    group       <- unlist(strsplit(evalbetas.out$syntax, "\n"))
    start       <- unlist(strsplit(setup.out$syntax, "\n"))
    group_paths <- group[-(1:length(start))]
    subjects    <- setup.out$subjects
    syntax.all  <- as.matrix(previous.out$syntax_sub$syntax)
    files       <- as.matrix(previous.out$syntax_sub$files)
    n_group_paths <- length(group_paths)
    if (n_group_paths == 0) bad = 0
  }
  
  cutoffz           = abs(qnorm(.05/subjects))
  
  #getting coefficients for final model
  while (bad == 1) {
    if (post_sub_prune == FALSE){
      list.all  <- matrix(NA, nrow = (n_fixed_paths + n_group_paths + 
                                        n_subgroup_paths)*subjects, ncol = 2)
    } else {
      list.all  <- matrix(NA, nrow = n_group_paths*subjects, ncol = 2)
    }
    
    colnames(list.all) <- c("param", "z")
    count_converge     <- 0
    
    for (k in 1:subjects){
      if (post_sub_prune == TRUE) syntax  <- syntax.all[k, ]
      
      # data.file           <- read.data(file   = files[k],
      #                                  sep    = sep,
      #                                  header = header)
      # colnames(data.file) <- c(varnames)
      fit                 <- fit.model(syntax    = syntax,
                                       data.file = ts_list[[k]])
      
      check_npd             <- any(grepl("error", class(fit)) == TRUE)
      check_zero_se  <- sum(lavInspect(fit, "se")$beta, na.rm = TRUE) == 0
      
      if (check_npd == FALSE & check_zero_se == FALSE) {
        all_mi         <- tryCatch(modindices(fit), error = function(e) e)
        check_singular <- any(grepl("singular", all_mi) == TRUE)
        converge       <- lavInspect(fit, "converged")
        check_error    <- any(grepl("error", class(all_mi)) == TRUE)
      } else {
        check_singular <- TRUE
        converge       <- FALSE
        check_error    <- TRUE
      }
      
      if (check_singular == FALSE & converge == TRUE & 
          check_zero_se == FALSE & check_error == FALSE) {
        ## printing the step
        if (subgroup.step == FALSE) {
          writeLines(paste("group-level pruning, subject", k))
        } else {
          writeLines(paste("subgroup-level pruning, subject", k))
        }
        z.list         <- subset(standardizedSolution(fit), op == "~")
        z.list$param   <- paste0(z.list$lhs, z.list$op, z.list$rhs)
        z.list         <- z.list[, c("param", "z")]
        z.list         <- z.list[complete.cases(z.list), ]
        if (post_sub_prune == TRUE) z.list <- z.list[z.list$param %in% group_paths, ]
        # below line added because a group path could have an NA std. error
        # for a person, which could result in it not being in the z.list object
        # this would then cause the object to be empty and throw an error later
        if (nrow(z.list) == 0 & post_sub_prune == TRUE){
          z.list <- matrix(NA, nrow = n_group_paths, ncol = 2)
        } 
        if (nrow(z.list) == 0 & post_sub_prune == FALSE){
          z.list <- matrix(NA, nrow = (rois + n_group_paths), ncol = 2)
        }
        count_converge <- count_converge + 1
      }
      # if it doesn't converge
      if (converge == FALSE | check_singular == TRUE | check_npd == TRUE | 
          check_zero_se == TRUE | check_error == TRUE) {
        if (post_sub_prune == FALSE){
          z.list <- matrix(NA, nrow = (rois + n_group_paths), ncol = 2)
        } else {
          z.list <- matrix(NA, nrow = n_group_paths, ncol = 2)
        }
        ## printing the step
        if(subgroup.step == FALSE) {
          writeLines(paste("group-level pruning, subject", k, "nonconvergence"))
        } else {
          writeLines(paste("subgroup-level pruning, subject", k, "nonconvergence"))
        }
      }
      list.all[(((nrow(z.list)*k)-nrow(z.list)+1):(nrow(z.list)*k)), ] <- as.matrix(z.list)
    }
    
    list.all     <- as.data.frame(list.all)
    list.all$z   <- as.numeric(as.character(list.all$z))
    list.all$sig <- ifelse(abs(list.all$z) > cutoffz, 1, 0)
    list.all     <- transform(list.all, sum = ave(sig, param, FUN = sum))
    # add line to not include AR as
    if (ar == TRUE){
      list.all <- subset(list.all, !(param %in% ar.paths))
    }
    # only those paths
    # added unique to subgrouping can be removed in evalbetas
    if (subgroup.step == TRUE) list.all <- subset(list.all, param %in% subgroup_paths)
    if ((nrow(list.all) == 0) == TRUE) bad <- 0
    paramdrop <- as.character(list.all[which.min(list.all$sum), 1])
    #create histogram of t/z values
    
    y2   <- abs(subset(list.all, param == paramdrop, select = z)$z)
    
    prop <- sum(y2 >= cutoffz)/count_converge
    
    drop_element <- ifelse(max(y2) < cutoffz, TRUE, FALSE)
    
    if (subgroup.step == TRUE)   cutoffprop <- subcutoff
    if (subgroup.step == FALSE)  cutoffprop <- groupcutoff
    
    if (prop <= cutoffprop | drop_element == TRUE) {
      ## drops offending path from syntax character string
      syntax <- unlist(strsplit(syntax, "[\n]"))
      syntax <- syntax[!syntax %in% paramdrop]
      syntax <- paste0(syntax, collapse = "\n")
      bad    <- 1
      if (subgroup.step == TRUE) {n_subgroup_paths <- n_subgroup_paths - 1
      if ((n_subgroup_paths == 0) == TRUE) bad <- 0
      subgroup_paths <- subgroup_paths[!subgroup_paths %in% paramdrop]
      }
      if (subgroup.step == FALSE) {n_group_paths <- n_group_paths - 1
      if ((n_group_paths == 0) == TRUE) bad <- 0
      subgroup_paths <- NULL
      }
      if (post_sub_prune == TRUE) {
        for (i in 1:nrow(syntax.all)){
          syntax.all.paths <- unlist(strsplit(syntax.all[i,],"\n"))
          syntax.all.paths <- syntax.all.paths[!syntax.all.paths %in% paramdrop]
          syntax.all[i,]   <- paste0(syntax.all.paths, collapse="\n")
        }
      }
    } else {bad <- 0
    if (subgroup.step == FALSE) subgroup_paths <- NULL
    }
    if (length(y2) <= subjects/2) bad <- 0
  }
  
  
  if (subgroup.step == FALSE) subgroup_paths <- NULL
  if (subgroup.step == TRUE)  subgroup_paths <- paste0(subgroup_paths, collapse = ",")
  if (post_sub_prune == TRUE) {
    syntax <- syntax.all
    n_subgroup_paths <- NULL
  }
  
  list <- list("syntax"               = syntax,
               "n_group_paths"    = n_group_paths,
               "n_subgroup_paths" = n_subgroup_paths,
               "subgroup_paths"       = subgroup_paths)
  return(list)
}

#this function adds individual paths until two of four fit indices are excellent for individual
addind <- function (done,
                    evaluate,
                    syntax,
                    data.file,
                    setup.out) {
  
  varnames        = setup.out$varnames
  cutoffind       = setup.out$cutoffind
  candidate_paths = setup.out$candidate_paths
  ar              = setup.out$ar
  
  param = NULL

  count.ind.paths     <- 0
  vec.MI              <- character()
  
  while (done == 0) { # done <- 0

    fit <- fit.model(syntax    = syntax,
                     data.file = data.file)
    
    check_npd            <- any(grepl("error",class(fit)) == TRUE)
    if (check_npd == FALSE){
      check_zero_se <- sum(lavInspect(fit,"se")$beta, na.rm = TRUE) == 0
      converge             <- lavInspect(fit, "converged")
    } else {
      check_zero_se <- FALSE
      converge <- FALSE
    } 
    if (ar == FALSE & count.ind.paths == 0) check_zero_se <- FALSE
    
    if (check_npd == FALSE & check_zero_se == FALSE & converge == TRUE) {
      all_mi          <- tryCatch(modindices(fit, op = "~"), error = function(e) e)  
      check_singular  <- any(grepl("singular", all_mi) == TRUE)
      empty_mi        <- ifelse(nrow(all_mi) == 0, TRUE, FALSE)
      if (check_singular == TRUE) empty_mi <- TRUE
      converge        <- lavInspect(fit, "converged")
      check_error     <- any(grepl("error", class(all_mi))==TRUE)
      if (converge == FALSE){
        check_fit <- FALSE
      } else if (converge == TRUE){
        check_fit <- any(is.na(fitMeasures(fit,c("chisq","df","pvalue","rmsea",
                                                 "srmr","nnfi","cfi"))))
      }
    } else {
      check_singular <- TRUE
      converge       <- FALSE
      check_error    <- TRUE
      check_fit      <- TRUE
      empty_mi       <- TRUE
    }
    if (check_singular == FALSE & converge == TRUE & check_npd == FALSE & 
        check_zero_se == FALSE & check_error == FALSE & 
        check_fit == FALSE & empty_mi == FALSE) {
      indMI          <- as.matrix(all_mi[,c("lhs","op","rhs","mi")])
      indMI          <- as.data.frame(indMI[complete.cases(indMI),])
      indMI$param    <- paste0(indMI$lhs, indMI$op, indMI$rhs)
      indMI          <- subset(indMI, param %in% candidate_paths)
      indMI$mi       <- as.numeric(as.character(indMI$mi))
      indparamadd    <- indMI[which.max(indMI$mi),5]
      indparamaddval <- indMI[which.max(indMI$mi),4]
      indices        <- fitMeasures(fit,c("chisq","df","pvalue",
                                          "rmsea","srmr", "nnfi","cfi"))
      excellent <- count_excellent(indices)
      if (excellent >= 2) {
        done   <- 1
        fixfit <- 0
      } else if (excellent <2) {
        done <- 0
        if (indparamaddval >= cutoffind) {
          syntax          <- paste(syntax, as.name(indparamadd), sep = "\n")
          count.ind.paths <- count.ind.paths + 1
          vec.MI          <- append(vec.MI,indparamadd)
        } else {
          done   <- 1
          fixfit <- 1
        }
      }
    } else {done = 1; fixfit = 0}
    fit <- fit.model(syntax    = syntax,
                     data.file = data.file)
    
    check_npd            <- any(grepl("error",class(fit)) == TRUE)
    if (check_npd == FALSE){
      check_zero_se <- sum(lavInspect(fit,"se")$beta, na.rm = TRUE) == 0
    } else check_zero_se = TRUE
    
    if (check_npd == FALSE & check_zero_se == FALSE) {
      all_mi            <- tryCatch(modindices(fit), error=function(e) e)
      check_singular    <- any(grepl("singular", all_mi) == TRUE)
      empty_mi          <- ifelse(nrow(all_mi) == 0, TRUE, FALSE)
      if (check_singular == TRUE) empty_mi <- TRUE
      converge          <- lavInspect(fit, "converged")
      check_error       <- any(grepl("error", class(all_mi)) == TRUE)
      if (converge == FALSE){
        check_fit <- FALSE
      } else if (converge == TRUE){
        check_fit <- any(is.na(fitMeasures(fit,c("chisq","df","pvalue","rmsea",
                                                 "srmr","nnfi","cfi"))))
      }
    } else {
      check_singular <- TRUE
      converge       <- FALSE
      check_error    <- TRUE
      check_fit      <- TRUE
      empty_mi       <- TRUE
    }
    
    if (converge == TRUE & check_singular == FALSE & check_npd == FALSE &
        check_zero_se == FALSE & check_error == FALSE & 
        check_fit == FALSE & empty_mi == FALSE) {
      indices   <- fitMeasures(fit,c("chisq","df","pvalue","rmsea",
                                     "srmr", "nnfi","cfi"))
      excellent <- count_excellent(indices)
      if (excellent >= 2) {done <- 1; fixfit <- 0}
    } else {done <- 1;  fixfit <- 0; evaluate <- 0}
    if (converge == FALSE)      done <- 1
    if (check_singular == TRUE) done <- 1
    if (check_fit == TRUE)      done <- 1
    if (check_error == TRUE)    done <- 1
    if (check_npd == TRUE)      done <- 1
    if (empty_mi == TRUE)       done <- 1
  }
  if (count.ind.paths==0) evaluate <- 0
  list <- list("evaluate" = evaluate,
               "fixfit"   = fixfit,
               "syntax"   = syntax,
               "vec.MI"   = vec.MI)
  return(list)
}

## this function prunes paths that have become nonsignificant for that individual,
## but group- and subgroup-level paths aren't candidates for pruning
evalind <- function (addind.out,
                     setup.out,
                     data.file) {
  
  evaluate  = addind.out$evaluate
  varnames  = setup.out$varnames
  syntax    = addind.out$syntax
  fixfit    = addind.out$fixfit
  vec.MI    = addind.out$vec.MI
  
  op    = NULL
  param = NULL
  
  
  while (evaluate==1) {
    fit <- fit.model(syntax    = syntax,
                     data.file = data.file)
    
    check_npd            <- any(grepl("error",class(fit))==TRUE)
    check_zero_se <- sum(lavInspect(fit,"se")$beta,na.rm=TRUE)==0
    
    if (check_npd == FALSE & check_zero_se == FALSE) {
      all_mi         <- tryCatch(modindices(fit), error = function(e) e)
      check_singular <- any(grepl("singular", all_mi) == TRUE)
      converge       <- lavInspect(fit, "converged")
      check_error    <- any(grepl("error", class(all_mi)) == TRUE)
    } else {
      check_singular <- TRUE
      converge       <- FALSE
      check_error    <- TRUE
    }
    
    if (check_singular==FALSE & converge==TRUE & check_npd==FALSE & check_error==FALSE) {
      indlist        <- subset(standardizedSolution(fit), op == "~")
      indlist$param  <- paste0(indlist$lhs, indlist$op, indlist$rhs)
      indlist        <- indlist[c("param","z")]
      indlist        <- indlist[complete.cases(indlist),]
      indlist        <- as.data.frame(indlist)
      indlist$z      <- as.numeric(as.character(indlist$z))
      indlist        <- subset(indlist, param %in% vec.MI)
      parampruneposs <- as.character(indlist[which.min(indlist$z),1])
      parampruneval  <- as.numeric(indlist[which.min(indlist$z),2])
      prune          <- ifelse(parampruneval>1.96, 0, 1)
      if (nrow(indlist) == 0) prune <- 0
      if (prune == 1) {
        syntax       <- unlist(strsplit(syntax, "[\n]"))
        syntax       <- syntax[!syntax %in% parampruneposs]
        syntax       <- paste0(syntax, collapse = "\n")
        paramprune   <- parampruneposs
        evaluate     <- 1
      } else {evaluate <- 0; fixfit <- 0}
      fit <- fit.model(syntax=syntax,
                       data.file=data.file)
      check_npd            <- any(grepl("error",class(fit))==TRUE)
      check_zero_se <- sum(lavInspect(fit,"se")$beta,na.rm=TRUE)==0
      
      if (check_npd == FALSE & check_zero_se==FALSE) {
        all_mi            <- tryCatch(modindices(fit),error=function(e) e)
        check_singular      <- any(grepl("singular",all_mi)==TRUE)
        converge            <- lavInspect(fit, "converged")
        check_error         <- any(grepl("error",class(all_mi))==TRUE)
      } else {
        check_singular <- TRUE
        converge       <- FALSE
        check_error    <- TRUE
      }
      
      if (converge==TRUE & check_singular==FALSE & check_npd==FALSE & 
          check_zero_se==FALSE & check_error == FALSE){
        indices      <- fitMeasures(fit,c("chisq","df","pvalue",
                                          "rmsea","srmr",
                                          "nnfi","cfi"))
        excellent <- count_excellent(indices)
        fixfit    <- ifelse(excellent >= 2, 0, 1)
        # if it doesn't converge
      }
      if (converge == FALSE | check_singular == TRUE) {evaluate <- 0; fixfit <- 0}
    }
  }
  list <- list("fixfit" = fixfit,
               "syntax" = syntax)
  return(list)
}


## this adds paths, even if not significant, until excellent model is obtained
fixfitind <- function (setup.out,
                       evalind.out,
                       data.file) {
  
  fixfit          = evalind.out$fixfit
  varnames        = setup.out$varnames
  syntax          = evalind.out$syntax
  candidate_paths = setup.out$candidate_paths
  ar              = setup.out$ar
  
  param = NULL
  
  while (fixfit == 1) {
    fit <- fit.model(varnames  = varnames,
                     syntax    = syntax,
                     data.file = data.file)
    
    check_npd            <- any(grepl("error", class(fit)) == TRUE)
    ## only for fixfit, temporary solution
    if (ar == FALSE){
      check.se.zero   <- sum(lavInspect(fit, "se")$beta, na.rm = TRUE) == 0
      check.beta.zero <- sum(lavInspect(fit, "est")$beta, na.rm = TRUE) == 0
      if (check.beta.zero == FALSE & check.se.zero == TRUE) {
        check_zero_se <- TRUE
      } else check_zero_se <- FALSE
    } else check_zero_se <- sum(lavInspect(fit, "se")$beta, na.rm = TRUE) == 0
    
    if (check_npd == FALSE & check_zero_se == FALSE) {
      all_mi            <- tryCatch(modindices(fit), error = function(e) e)
      check_singular      <- any(grepl("singular", all_mi) == TRUE)
      converge            <- lavInspect(fit, "converged")
      check_error         <- any(grepl("error", class(all_mi)) == TRUE)
    } else {
      check_singular <- TRUE
      converge       <- FALSE
      check_error    <- FALSE
    }
    
    if (check_singular == FALSE & converge == TRUE & check_npd == FALSE & 
        check_zero_se == FALSE & check_error == FALSE) {
      indMI        <- as.matrix(all_mi[all_mi$op == "~", ])[,c("lhs","op","rhs","mi")]
      indMI        <- as.data.frame(indMI[complete.cases(indMI),])
      indMI$param  <- paste0(indMI$lhs, indMI$op, indMI$rhs)
      indMI        <- subset(indMI, param %in% candidate_paths)
      indMI$mi     <- as.numeric(as.character(indMI$mi))
      indparamadd  <- indMI[which.max(indMI$mi),5]
      syntax       <- paste(syntax, as.name(indparamadd), sep = "\n")
      fit          <- fit.model(syntax    = syntax,
                                data.file = data.file)
      converge       <- lavInspect(fit, "converged")
      all_mi       <- tryCatch(modindices(fit), error = function(e) e)
      check_singular <- any(grepl("singular", all_mi) == TRUE)
      check_error    <- any(grepl("error", class(all_mi)) == TRUE)
      if (check_singular == FALSE & converge == TRUE & check_error == FALSE) {
        indices      <- fitMeasures(fit,c("chisq", "df", "pvalue",
                                          "rmsea", "srmr", "nnfi", "cfi"))
        df     <- indices[2]
        if (df == 0) fixfit <- 0
        excellent <- count_excellent(indices)
        fixfit    <- ifelse(excellent >= 2, 0, 1)
      }
      # if it doesn't converge
      if (converge == FALSE | check_singular == TRUE | check_error == TRUE) fixfit <- 0
    } else fitfix <- 0
  }
  return(syntax)
}

## this fits the final model and extracts relevant params that we'll summarize for user later
## also creates/prints individual-level matrices and plots as it runs
final.fit <- function(setup.out,
                      fixfitind.out,
                      data.file,
                      k){
  
  varnames   = setup.out$varnames
  lvarnames  = setup.out$lvarnames
  syntax     = fixfitind.out
  trackparts = setup.out$trackparts
  rois       = setup.out$rois
  plot.names = setup.out$plot.names
  x          = setup.out$x
  y          = setup.out$y
  fitted     = setup.out$fitted
  individual = setup.out$individual
  out        = setup.out$out
  plot       = setup.out$plot
  agg        = setup.out$agg
  
  ind_plot = NA
  op = NULL
  
  fit <- fit.model(syntax    = syntax,
                   data.file = data.file)
  
  check_npd            <- any(grepl("error", class(fit)) == TRUE)
  if (check_npd == FALSE){
    check_zero_se      <- sum(lavInspect(fit, "se")$beta, na.rm = TRUE) == 0
  } else check_zero_se = FALSE
  
  if (check_npd == FALSE & check_zero_se == FALSE) {
    all_mi            <- tryCatch(modindices(fit), error = function(e) e)
    check_singular      <- any(grepl("singular", all_mi) == TRUE)
    converge            <- lavInspect(fit, "converged")
    check_error         <- any(grepl("error", class(all_mi)) == TRUE)
  } else {
    check_singular <- TRUE
    converge       <- FALSE
    check_error    <- TRUE
  }
  if (converge == TRUE) last.converge <- TRUE
  
  ## update code here to remove last element
  if (converge == FALSE | check_error == TRUE) {
    last.converge <- FALSE
    syntax <- unlist(strsplit(syntax, "[\n]"))
    syntax <- syntax[-length(syntax)]
    syntax <- paste0(syntax, collapse = "\n")
    fit    <- fit.model(syntax    = syntax,
                        data.file = data.file)
    
    check_npd            <- any(grepl("error", class(fit)) == TRUE)
    if (check_npd == F) {
      check_zero_se <- sum(lavInspect(fit, "se")$beta, na.rm = TRUE) == 0
    } else check_zero_se <- TRUE
    
    if (check_npd == FALSE & check_zero_se == FALSE) {
      all_mi         <- tryCatch(modindices(fit), error = function(e) e)
      check_singular <- any(grepl("singular", all_mi) == TRUE)
      converge       <- lavInspect(fit, "converged")
      check_error    <- any(grepl("error", class(all_mi)) == TRUE)
    } else {
      check_singular <- TRUE
      converge       <- FALSE
      check_error    <- TRUE
    }
  }
  
  ind.fit <- matrix(NA, nrow = 1, ncol = 10)
  if (converge == TRUE & check_singular == FALSE) {
    vcov_ind <- lavInspect(fit, "vcov.std.all")
    indices <- fitMeasures(fit,c("chisq","df","pvalue","rmsea",
                                 "srmr", "nnfi","cfi", "bic"))
    # insert indfit
    ind.fit[1,2] <- round(indices[1], digits = 4)
    ind.fit[1,3] <- indices[2]
    ind.fit[1,4] <- round(indices[3], digits = 4)
    ind.fit[1,5] <- round(indices[4], digits = 4)
    ind.fit[1,6] <- round(indices[5], digits = 4)
    ind.fit[1,7] <- round(indices[7], digits = 4)
    ind.fit[1,8] <- round(indices[6], digits = 4)
    ind.fit[1,9] <- round(indices[8], digits = 4)
    if (last.converge == FALSE) {
      ind.fit[1,10] <- "last known convergence"
    } else {
      ind.fit[1,10] <- "converged normally"
    }
    
    indlist              <- subset(standardizedSolution(fit), op == "~")
    indlist$param        <- paste0(indlist$lhs, indlist$op, indlist$rhs)
    indlist              <- indlist[complete.cases(indlist), ]
    indlist              <- as.data.frame(indlist)
    indsubject           <- matrix(k, nrow = (nrow(indlist)), ncol = 1)
    colnames(indsubject) <- c("subject")
    indlist              <- cbind(indsubject, indlist)
    #creating individual-level beta matrices
    individual.paths     <- matrix(0, nrow = (rois*2), ncol = (rois*2))
    individual.SEs       <- matrix(0, nrow = (rois*2), ncol = (rois*2))
    indlist$row          <- substring(indlist$lhs, 4)
    indlist$col          <- substring(indlist$rhs, 4)
    indrows              <- indlist$row
    indcols              <- indlist$col
    ## multiple steps to recode column names to appropriate numerics
    
    indcols     <- as.numeric(recoderFunc(indcols, y, x))
    indrows     <- as.numeric(recoderFunc(indrows, y, x))
    
    indbetas    <- as.numeric(as.character(indlist$est.std))
    indSEs      <- as.numeric(as.character(indlist$se))
    indelements <- nrow(indlist)
    
    for (s in 1:indelements){
      individual.paths[indrows[s], indcols[s]] <- indbetas[s]
    }
    
    for (q in 1:indelements){
      individual.SEs[indrows[q], indcols[q]] <- indSEs[q]
    }
    
    individual.paths <- round(individual.paths, digits = 4)
    individual.paths[is.na(individual.paths)] <- 0
    
    individual.SEs <- round(individual.SEs,digits=4)
    individual.SEs[is.na(individual.SEs)] <- 0
    
    individual.paths.all         <- individual.paths[(rois+1):(rois*2),1:(rois*2)]
    rownames(individual.paths.all) <- varnames[(rois+1):(rois*2)]
    colnames(individual.paths.all) <- varnames[1:(rois*2)]
    individual.SEs             <- individual.SEs[(rois+1):(rois*2),1:(rois*2)]
    rownames(individual.SEs)   <- varnames[(rois+1):(rois*2)]
    colnames(individual.SEs)   <- varnames[1:(rois*2)]
    
    if (agg == TRUE & !is.null(out)) {
      write.csv(individual.paths.all, file = file.path(out,"allBetas.csv"),     row.names = TRUE)
      write.csv(individual.SEs,       file = file.path(out,"allStdErrors.csv"), row.names = TRUE)
    } else if (agg == FALSE & !is.null(out)){
      write.csv(individual.paths.all, file = file.path(individual, 
                                                       paste0(trackparts[k,2],"Betas.csv")),    
                row.names = TRUE)
      write.csv(individual.SEs,       file = file.path(individual, 
                                                       paste0(trackparts[k,2],"StdErrors.csv")), 
                row.names = TRUE)
    }
    
    if (plot == TRUE){
      individual.paths.t     <- t(individual.paths)
      Lagged                 <- individual.paths.t[1:(rois),(rois+1):(rois*2)]
      Contemporaneous        <- individual.paths.t[(rois+1):(rois*2),(rois+1):(rois*2)]
      eLagged                <- W2E(Lagged)
      eContemporaneous       <- W2E(Contemporaneous)
      isLagged               <- c(rep(TRUE, sum(Lagged != 0)), rep(FALSE, sum(Contemporaneous != 0)))
      plotind                <- file.path(individual,paste0(trackparts[k,2],"Plot.pdf"))
      if (agg == TRUE) {
        plotind <- file.path(out, "summaryPathsPlot.pdf")
      }
      ind_plot <- tryCatch(qgraph(rbind(eLagged, eContemporaneous),
                                  layout              = "circle",
                                  lty                 = ifelse(isLagged,2, 1),
                                  edge.labels         = F,
                                  curve               = FALSE,
                                  parallelEdge        = TRUE,
                                  fade                = FALSE,
                                  posCol              = "red",
                                  negCol              = "blue",
                                  labels              = plot.names,
                                  label.cex           = 3,
                                  edge.label.cex      = 1.5,
                                  edge.label.position = .3,
                                  DoNotPlot           = TRUE),error=function(e) e)
      donotplot <- ifelse("error" %in% class(ind_plot), TRUE, FALSE)
      if (!is.null(out) & !donotplot){
        pdf(plotind)
        plot(ind_plot)
        dev.off()
      }
    }
  }
  if (converge == FALSE) {
    ind.fit[1,10] <- "nonconvergence"
    indlist      <- data.frame()
    individual.paths <- data.frame()
    vcov_ind <- NA
  }
  if (check_singular == TRUE) {
    ind.fit[1,10] <- "computationally singular"
    indlist      <- data.frame()
    individual.paths <- data.frame()
    vcov_ind <- NA
  }
  if (check_error == TRUE) {
    ind.fit[1,10] <- "error"
    indlist      <- data.frame()
    individual.paths <- data.frame()
    vcov_ind <- NA
  }
  list <- list("ind.fit"      = ind.fit,
               "ind.elements" = indlist,
               "syntax"       = syntax,
               "ind.paths"    = individual.paths,
               "ind_plot"     = ind_plot,
               "vcov_ind"     = vcov_ind)
  return(list)
}

## puts all of them together
indsem.internal <- function(setup.out,
                            evalbetassub.out,
                            evalbetas.out){
  
  subjects         = setup.out$subjects
  varnames         = setup.out$varnames
  trackparts       = setup.out$trackparts
  plot             = setup.out$plot
  modularity       = evalbetassub.out$modularity
  subgroup         = setup.out$subgroup
  ts_list          = setup.out$ts_list
  
  all_elem          <- data.frame()
  all_fit           <- matrix(NA, nrow=subjects, ncol=10)
  all_syntax        <- matrix(NA, nrow=subjects,ncol=4)
  all_ind_paths     <- list()
  all_plots         <- list()
  vcov_params       <- list()
  colnames(all_fit) <- c("subject", "chisq", "df", "pval",
                         "rmsea", "srmr", "nnfi", "cfi", "bic", "status")
  
  #files            <- list.files(data, full.names=TRUE)
  
  all.diff.subgroups <- length(unique(evalbetassub.out$syntax_sub[,1])) == subjects
  
  if (subgroup == TRUE & all.diff.subgroups == FALSE){
    all <- unlist(strsplit(evalbetassub.out$final_syntax[,2][1], "\n"))
    # if subgroup paths are empty, then just leave final syntax as above "all"
    if (is.na(evalbetassub.out$sub.paths[,2][1]) == FALSE){
      sub <- unlist(strsplit(evalbetassub.out$sub.paths[,2][1],","))
      # if syntax has no subgroup paths, use as is
      if (length(sub) == 0){
        all <- paste0(all, collapse = "\n")
        # if syntax has subgroup paths, remove them to get group syntax
      } else {
        all <- head(all, -length(sub))
        all <- paste0(all, collapse = "\n")
      }
    }
    sub.info <- as.data.frame(evalbetassub.out$syntax_sub, stringsAsFactors = FALSE)
    sorted   <- sub.info[order(sub.info$files), ]
  }
  
  if (subgroup == TRUE & all.diff.subgroups == TRUE){
    all      <- evalbetassub.out$final_syntax[,2][1]
    sub.info <- as.data.frame(evalbetassub.out$syntax_sub, stringsAsFactors = FALSE)
    sorted   <- sub.info[order(sub.info$files), ]
  }
  
  if (subgroup == FALSE){
    all <- evalbetas.out$syntax
  }
  
  if (is.null(evalbetas.out)){
    all <- setup.out$syntax
  }
  
  for (k in 1:subjects) {
    
    writeLines(paste("individual-level search, subject", k))
    
    if (subgroup == TRUE) {
      syntax    <- sorted[k,3]
      
      # data.file <- read.data(file   = sorted[k,2],
      #                        sep    = sep,
      #                        header = header)
      
      data.file <- ts_list[[sorted[k,2]]]
      
      if (is.na(sorted[k,1]) == T) syntax <- all # subgroup assignment check
      
    } else {
      
      # data.file <- read.data(file   = files[k],
      #                        sep    = sep,
      #                        header = header)
      
      data.file <- ts_list[[k]]
      
      syntax    <- all
    }
    
    addind.out <- addind(done            = 0,
                         evaluate        = 1,
                         syntax          = syntax,
                         data.file       = data.file,
                         setup.out       = setup.out)
    
    evalind.out <- evalind(addind.out = addind.out,
                           setup.out  = setup.out,
                           data.file  = data.file)
    
    fixfitind.out <- fixfitind(setup.out   = setup.out,
                               evalind.out = evalind.out,
                               data.file   = data.file)
    
    final.fit.out <- final.fit(setup.out     = setup.out,
                               fixfitind.out = fixfitind.out,
                               data.file     = data.file,
                               k             = k)
    
    all_elem           <- rbind(all_elem, final.fit.out$ind.elements)
    all_fit[k,]        <- as.matrix(final.fit.out$ind.fit)
    all_syntax[k,3]    <- as.matrix(final.fit.out$syntax)
    all_syntax[k,1]    <- k
    all_ind_paths[[k]] <- final.fit.out$ind.paths
    all_plots[[k]]     <- final.fit.out$ind_plot
    vcov_params[[k]]   <- final.fit.out$vcov_ind
    # names(all_ind_paths[k]) <- trackparts[k,2]
  }
  
  all_fit[,1]    <- trackparts[,2]
  #all_syntax[,2] <- files
  all_syntax[,2] <- trackparts[,2]
  all_syntax[,4] <- all
  colnames(all_syntax) <- c("subject", "files", "syntax_ind", "syntax_group")
  names(all_ind_paths) <- trackparts[ ,2]
  names(all_plots)     <- trackparts[ ,2]
  names(vcov_params)   <- trackparts[ ,2]
  
  if (subgroup == TRUE){
    colnames(sorted)[3] <- c("syntax_sub")
    #add code here to arrange group, subgroup, and individual-level paths
    all.syntax_sub     <- merge(all_syntax, sorted, by = "files")
    all_fit            <- cbind(all_fit, as.numeric(as.character(sorted[,1])))
    modularity         <- append(modularity, rep("",nrow(all_fit)-1))
    all_fit            <- cbind(all_fit,modularity)
    colnames(all_fit)  <- c("subject", "chisq", "df", "pval",
                            "rmsea", "srmr", "nnfi", "cfi", "bic",
                            "status", "subgroup", "modularity")
    
  } else all.syntax_sub <- as.data.frame(all_syntax)
  
  list <- list("all_elem" = all_elem,
               "all_fit"      = all_fit,
               "all_syntax"   = all.syntax_sub,
               "all_diff_sub" = all.diff.subgroups,
               "all_ind_paths" = all_ind_paths, 
               "all_plots"     = all_plots,
               "vcov_params"   = vcov_params)
  return(list)
}

## sorts and organizes information, prints summary files
wrapup <- function(indsem.internal.out,
                   setup.out){
  
  all_elem     = indsem.internal.out$all_elem
  all_fit      = indsem.internal.out$all_fit
  all_syntax   = indsem.internal.out$all_syntax
  all_diff_sub = indsem.internal.out$all_diff_sub
  rois         = setup.out$rois
  out          = setup.out$out
  plot.names   = setup.out$plot.names
  subjects     = setup.out$subjects
  x            = setup.out$x
  y            = setup.out$y
  varnames     = setup.out$varnames
  lvarnames    = setup.out$lvarnames
  subgroup_dir = setup.out$subgroup_dir
  plot         = setup.out$plot
  agg          = setup.out$agg
  subgroup     = setup.out$subgroup
  all_ind_paths = indsem.internal.out$all_ind_paths

  present = NULL
  sig     = NULL
  est.std = NULL
  param   = NULL
  membership = NULL
  
  if (subgroup == TRUE & all_diff_sub == FALSE){
    all_merged_sub_final <- data.frame()
    sub_final            <- character()
    sub_ind_paths        <- character()
    sub_fullsub_paths    <- character()
    all_merged           <- merge(all_elem, all_syntax, by = "subject")
    n_subgroups          <- length(unique(all_merged$membership))
    
    for (p in 1:n_subgroups){
      all_merged_sub         <- subset(all_merged, membership == p)
      if (nrow(all_merged_sub) != 0){
        ## add code to add column for "group","subgroup","individual"
        n_sub_subjects         <- length(unique(all_merged_sub$files))
        all_merged_sub$est.std <- as.numeric(as.character(all_merged_sub$est.std))
        all_merged_sub         <- transform(all_merged_sub, 
                                            mean.beta = (ave(est.std, param, FUN = sum))/n_sub_subjects)
        all_merged_sub$sig     <- ifelse(abs(all_merged_sub$z) > 1.96, 1, 0)
        all_merged_sub$present <- ifelse(all_merged_sub$sig < 2, 1, 0)
        all_merged_sub         <- transform(all_merged_sub, sum.sig = ave(sig, param, FUN = sum))
        all_merged_sub         <- transform(all_merged_sub, count   = ave(present, param, FUN = sum))
        
        if (n_sub_subjects != 1){
          ## using syntax to identify group, subgroup, individual level paths
          group    <- unlist(strsplit(as.vector(unique(all_merged_sub$syntax_group)), "\n"))
          sub      <- unlist(strsplit(as.vector(unique(all_merged_sub$syntax_sub)), "\n"))
          ind      <- unlist(strsplit(as.vector(unique(all_merged_sub$syntax_ind)), "\n"))
          ind      <- unique(subset(ind, !ind %in% sub))
          sub      <- subset(sub, !sub %in% group)
          sub_fullsub_paths <- append(sub_fullsub_paths,sub)
        } else {
          
          # if there's only one person in the subgroup, then the ind paths are the sub paths
          # only used to determine whether or not a subgroup path is promoted to a group path
          group    <- unlist(strsplit(as.vector(unique(all_merged_sub$syntax_group)),"\n"))
          ind      <- unlist(strsplit(as.vector(unique(all_merged_sub$syntax_ind)),"\n"))
          ind      <- subset(ind, !ind %in% group)
          sub      <- ind
          sub_ind_paths <- append(sub_ind_paths, sub)
        }
        
        label  <- c(rep("group",    length(group)), 
                    rep("subgroup", length(sub)), 
                    rep("ind",      length(ind)))
        color  <- c(rep("black",    length(group)), 
                    rep("green3",   length(sub)),   
                    rep("gray50",   length(ind)))
        param  <- c(group, sub, ind)
        levels <- data.frame(param, label, color)
        
        all_merged_sub       <- merge(all_merged_sub, levels, by = "param")
        all_merged_sub_final <- rbind(all_merged_sub_final, all_merged_sub)
        sub_final            <- append(sub_final, sub)
      }
    }
    
    ## this portion checks to see if there is a subgroup path that exists for each subgroup
    ## if it does, it bumps it up to the group level
    a            <- table(sub_final)
    sub_to_group <- names(a[a == n_subgroups])
    sub_stay_sub <- names(a[a != n_subgroups])
    all_merged_sub_final$label[all_merged_sub_final$param %in% sub_to_group] <- "group"
    all_merged_sub_final$color[all_merged_sub_final$param %in% sub_to_group] <- "black"
    
    
    sub_plots <- list()
    sub_paths <- list()
    ## code to create subgroup-level plots and matrices
    for (p in 1:n_subgroups){
      all_merged_sub   <- subset(all_merged_sub_final, membership == p)
      if (nrow(all_merged_sub) !=0){
        n_sub_subjects   <- length(unique(all_merged_sub$files))
        sub.paths        <- matrix(0, nrow = (rois*2), ncol = (rois*2))
        sub.colors       <- matrix(NA, nrow = (rois*2), ncol = (rois*2))
        elements         <- nrow(all_merged_sub)
        rows             <- substring(all_merged_sub$lhs, 4)
        cols             <- substring(all_merged_sub$rhs, 4)
        cols             <- as.numeric(recoderFunc(cols, y, x))
        rows             <- as.numeric(recoderFunc(rows, y, x))
        counts           <- as.numeric(all_merged_sub$count)
        beta.weights     <- as.numeric(all_merged_sub$mean.beta)
        colors           <- as.character(all_merged_sub$color)
        for (m in 1:elements){
          sub.paths [rows[m], cols[m]] <- counts[m]
          sub.colors[rows[m], cols[m]] <- colors[m]
        }
        sub.paths <- round(sub.paths,digits = 4)
        sub.paths[is.na(sub.paths)] <- 0
        sub.paths.all <- sub.paths[(rois+1):(rois*2),(1:(rois*2))]
        rownames(sub.paths.all) <- varnames[(rois+1):(rois*2)]
        colnames(sub.paths.all) <- varnames[1:(rois*2)]
        
        if (n_sub_subjects != 1 & !is.null(out))  {
          write.csv(sub.paths.all, file = file.path(subgroup_dir, 
                                                    paste0("subgroup", p, 
                                                           "PathCountsMatrix.csv")), 
                    row.names = TRUE)
        }
        
        if (plot == TRUE & n_sub_subjects != 1){
          sub.paths.t      <- t(sub.paths)
          sub.paths.t      <- sub.paths.t/n_sub_subjects
          Lagged           <- sub.paths.t[1:(rois),(rois+1):(rois*2)]
          Contemporaneous  <- sub.paths.t[(rois+1):(rois*2),(rois+1):(rois*2)]
          eLagged          <- W2E(Lagged)
          eContemporaneous <- W2E(Contemporaneous)
          sub.colors.t     <- t(sub.colors)
          colorLagged           <- sub.colors.t[1:(rois),(rois+1):(rois*2)]
          colorContemporaneous  <- sub.colors.t[(rois+1):(rois*2),(rois+1):(rois*2)]
          color.list       <- c(colorLagged,colorContemporaneous)
          color.list       <- color.list[!is.na(color.list)]
          isLagged         <- c(rep(TRUE,nrow(eLagged)), rep(FALSE,nrow(eContemporaneous)))
          # if an element is repeated, change its curve value
          curve            <- rep(1, length(isLagged))
          curve[which(duplicated(rbind(eLagged, eContemporaneous)[,1:2]))] <- .5
          plotsub          <- file.path(subgroup_dir,paste0("subgroup",p,"Plot.pdf"))
          sub_plot <- tryCatch(qgraph(rbind(eLagged,eContemporaneous),
                                      layout              = "circle",
                                      lty                 = ifelse(isLagged,2, 1),
                                      edge.labels         = FALSE,
                                      # curve               = curve,
                                      curve               = FALSE,
                                      parallelEdge        = TRUE,
                                      labels              = plot.names,
                                      edge.color          = color.list,
                                      fade                = FALSE,
                                      label.cex           = 3,
                                      edge.label.cex      = 1.5,
                                      edge.label.position = .3,
                                      DoNotPlot           = TRUE), error=function(e) e)
          if (!is.null(out)){
            pdf(plotsub)
            plot(sub_plot)
            dev.off()
          }
        } else {
          sub_plot  <- NULL
          sub.paths.all <- NULL
        } 
      } else {
        sub_plots <- NULL
        sub.paths.all <- NULL
      }
      sub_plots[[p]] <- sub_plot
      sub_paths[[p]] <- sub.paths.all
    }
    
    ## quick fix to make sure subgroup paths are represented as subgroup,
    ## not individual. this could happen if it's subgroup for some and individual for others
    all_elem_final       <- all_merged_sub_final
    all_elem_final$label <- as.character(all_elem_final$label)
    
    ## adjust count of paths, prepare for summary output
    a1 <- aggregate(present ~ label + param + membership, data = all_elem_final, sum)
    a1$label <- ifelse(a1$label == 'subgroup',
                       paste0(a1$label, a1$membership),as.character(a1$label))
    a2 <- aggregate(present ~ param + label, data = a1, sum)
    a2 <- a2[order(-a2$present, a2$label),]
    ## now transpose object to have column per label
    a2 <- reshape(a2, timevar = "label", idvar = "param", direction = "wide")
    a2[is.na(a2)] <- 0
    
    all_elem_final <- transform(all_elem_final, present = ave(present, param, FUN = sum))
    all_elem_final <- transform(all_elem_final, sum.sig = ave(sig, param, FUN = sum))
    # this line causes issues with membership, maybe come back, today = 09.19.2016
    # all_elem_final <- subset(all_elem_final, !duplicated(param)) ## this line problem
    
    # obtains shortened list of each unique path + label + subgroup combination for later merging
    all_elem_final <- all_elem_final[row.names(unique(all_elem_final[c("param", "label", "membership")])),]
    
    # grab relevant info from all_elem_final for later use
    # when merging subgroup and level information into indivPathEstimates
    # must take place here. otherwise, below code that labels a subgroup level path
    # as subgroup across ALL individuals will introduce incorrect information
    label_info   <- all_elem_final[,c("lhs", "rhs", "label", "membership")]
    
    ## removes individual-level paths from this list for subgroup size 1
    ## this way, no paths are marked as subgroup-level that only existed for one person
    ## in the final summary graphic
    sub_fullsub_paths <- unique(subset(sub_fullsub_paths, !sub_fullsub_paths %in% sub_to_group))
    # prioritizes any subgroup path for summary graphic; shows as green instead of gray
    all_elem_final$label[all_elem_final$param %in% sub_fullsub_paths] <- "subgroup"
    all_elem_final$color[all_elem_final$param %in% sub_fullsub_paths] <- "green3"
  } else {
    ## insert code for group and ind if all_diff_sub = TRUE
    sub_plots <- NULL
    sub_paths <- NULL
    all_elem_final         <- as.data.frame(all_elem)
    all_elem_final$est.std <- as.numeric(as.character(all_elem_final$est.std))
    if (agg == FALSE) {
      group    <- unlist(strsplit(as.vector(unique(all_syntax$syntax_group)), "\n"))
      ind      <- unlist(strsplit(as.vector(unique(all_syntax$syntax_ind)), "\n"))
      ind      <- unique(subset(ind, !ind %in% group))
      label  <- c(rep("group", length(group)), rep("ind", length(ind)))
      color  <- c(rep("black", length(group)), rep("gray50", length(ind)))
      param  <- c(group, ind)
      levels <- data.frame(param, label, color)
      all_elem_final <- transform(all_elem,
                                  mean.beta = (ave(est.std, param, FUN = sum))/subjects)
      all_elem_final$sig     <- ifelse(abs(all_elem_final$z) > 1.96, 1, 0)
      all_elem_final$present <- ifelse(all_elem_final$sig < 2, 1, 0)
      all_elem_final         <- transform(all_elem_final, sum.sig = ave(sig, param, FUN = sum))
      all_elem_final         <- transform(all_elem_final, present = ave(present, param, FUN = sum))
      
      all_elem_final         <- subset(all_elem_final, !duplicated(param))
      all_elem_final         <- merge(all_elem_final, levels, by = "param")
      a1 <- aggregate(present ~ label + param ,data = all_elem_final,sum)
      a1 <- a1[order(-a1$present,a1$label),]
      ## now transpose object to have column per label
      a1 <- reshape(a1, timevar = "label", idvar = "param", direction = "wide")
      a1[is.na(a1)] <- 0
    }
  }
  if (agg == FALSE){
    final.paths      <- matrix(0, nrow = (rois*2), ncol = (rois*2))
    final.colors     <- matrix(NA, nrow = (rois*2), ncol = (rois*2))
    elements         <- nrow(all_elem_final)
    rows             <- substring(all_elem_final$lhs, 4)
    cols             <- substring(all_elem_final$rhs, 4)
    cols             <- as.numeric(recoderFunc(cols, y, x))
    rows             <- as.numeric(recoderFunc(rows, y, x))
    counts           <- as.numeric(all_elem_final$present)
    beta.weights     <- as.numeric(all_elem_final$mean.beta)
    colors           <- as.character(all_elem_final$color)
    ## if all subgroups are different, just use black for group and gray for ind
    if (subgroup == TRUE & all_diff_sub == TRUE) { # begin
      colors <- character()
      colors <- replace(colors, counts == subjects, "black")
      colors <- replace(colors, counts != subjects, "grey50")
    } # end
    for (m in 1:elements){
      final.paths[rows[m], cols[m]] <- counts[m]
      final.colors[rows[m], cols[m]] <- colors[m]
    }
    final.paths[is.na(final.paths)] <- 0
    final.paths.all <- final.paths[(rois+1):(rois*2),1:(rois*2)]
    
    rownames(final.paths.all) <- varnames[(rois+1):(rois*2)]
    colnames(final.paths.all) <- varnames[1:(rois*2)]
    
    if (!is.null(out)){
      write.csv(final.paths.all, file = file.path(out, "summaryPathCountMatrix.csv"), row.names = TRUE)
    }
    
    if (subgroup == FALSE | all_diff_sub == TRUE){
      all_elem_summary  <- a1
    } else {
      all_elem_summary  <- a2
    }
    
    if (plot == TRUE){
      final.paths.t    <- t(final.paths)
      final.paths.t    <- final.paths.t/subjects
      Lagged           <- final.paths.t[1:(rois),(rois+1):(rois*2)]
      Contemporaneous  <- final.paths.t[(rois+1):(rois*2),(rois+1):(rois*2)]
      eLagged          <- W2E(Lagged)
      eContemporaneous <- W2E(Contemporaneous)
      final.colors.t        <- t(final.colors)
      colorLagged           <- final.colors.t[1:(rois),(rois+1):(rois*2)]
      colorContemporaneous  <- final.colors.t[(rois+1):(rois*2),(rois+1):(rois*2)]
      color.list       <- c(colorLagged,colorContemporaneous)
      color.list       <- color.list[!is.na(color.list)]
      isLagged         <- c(rep(TRUE, nrow(eLagged)), rep(FALSE, nrow(eContemporaneous)))
      plotfinal        <- file.path(out,"summaryPathsPlot.pdf")
      # if an element is repeated, change its curve value
      curve                  <- rep(1, length(isLagged))
      curve[which(duplicated(rbind(eLagged, eContemporaneous)[,1:2]))] <- .5
      #if (subgroup==TRUE){edge.color=color.list}else{edge.color=NULL}
      edge.color=color.list
      all_plot <- tryCatch(qgraph(rbind(eLagged,eContemporaneous),
                                  layout              = "circle",
                                  lty                 = ifelse(isLagged,2, 1),
                                  edge.labels         = FALSE,
                                  # curve               = curve,
                                  curve               = FALSE,
                                  parallelEdge        = TRUE,
                                  labels              = plot.names,
                                  edge.color          = edge.color,
                                  fade                = FALSE,
                                  label.cex           = 3,
                                  edge.label.cex      = 1.5,
                                  edge.label.position = .3,
                                  DoNotPlot           = TRUE),error=function(e) e)
      if (!is.null(out)){
        pdf(plotfinal)
        plot(all_plot)
        dev.off()
      }
    } else all_plot <- NULL
  } else {
    final.paths.all <- NULL
    all_plot        <- NULL
  }
  
  if (nrow(all_elem) != 0) {
    all_elem        <- all_elem[c("subject", "param", "est.std", "se", "pvalue", "z")]
    all_elem[ ,3:6] <- round(all_elem[ ,3:6], digits = 4)
    colnames(setup.out$trackparts) <- c("subject", "file")
    all_elem   <- merge(setup.out$trackparts, all_elem, 
                        by.x = "subject", by.y = "subject")[,-1]
    all_elem   <- all_elem[order(all_elem$file),]
    
    if (agg == TRUE) {
      all_elem[,1] <- "all"
      all_elem_summary <- data.frame()
    }
    
    dv           <- sapply(strsplit(all_elem$param,"~"),"[",1)
    iv           <- sapply(strsplit(all_elem$param,"~"),"[",2)
    dvs.ivs      <- as.data.frame(cbind(dv, iv))
    all_elem <- cbind(dvs.ivs, all_elem)
    all_elem <- all_elem[c("file","dv","iv","est.std","se","z","pvalue")]
    
    if (agg == FALSE){
      dv        <- sapply(strsplit(all_elem_summary$param,"~"),"[",1)
      iv        <- sapply(strsplit(all_elem_summary$param,"~"),"[",2)
      dvs.ivs   <- as.data.frame(cbind(dv, iv))
      all_elem_summary <- cbind(dvs.ivs, all_elem_summary)
      all_elem_summary$param <- NULL
    }
    
    # sub back variable names into all_elem summary and all_elem
    for (v in 1:length(varnames)){
      all_elem$dv <- recoderFunc(all_elem$dv, lvarnames, varnames)
      all_elem$iv <- recoderFunc(all_elem$iv, lvarnames, varnames)
      if (agg == FALSE) {
        all_elem_summary$dv <- recoderFunc(all_elem_summary$dv, lvarnames, varnames)
        all_elem_summary$iv <- recoderFunc(all_elem_summary$iv, lvarnames, varnames)
      }
    }
    colnames(all_elem)[4]       <- "beta"
    colnames(all_elem)[7]       <- "pval"
    row.names(all_elem_summary) <- NULL
    row.names(all_elem)         <- NULL
    all_fit                         <- as.data.frame(all_fit)
    colnames(all_fit)[1]            <- "file"
    if (!is.null(out)){
      write.csv(all_fit, file.path(out, "summaryFit.csv"), row.names = FALSE)
    }
    # ------- merge level and subgroup information into indivPathEstimates----------
    ###############################################################################
    if (subgroup == TRUE & all_diff_sub == FALSE){
      # merge subgroup assignment into indivPathEstimates file
      all_elem <- merge(all_elem, all_fit[, c("file", "subgroup")], by = "file")
      # sub in variable names to be consistent with all_elem
      for (v in 1:length(varnames)){
        label_info$lhs   <- recoderFunc(label_info$lhs, lvarnames, varnames)
        label_info$rhs   <- recoderFunc(label_info$rhs, lvarnames, varnames)
      }
      sub_ind_to_merge <- label_info[!label_info$label %in% "group", ]
      group_to_merge   <- label_info[label_info$label %in% "group", ][,-4]
      # merge all elements with subgroup and ind paths
      all_elem.a <- merge(all_elem, sub_ind_to_merge, 
                          by.x = c("dv",  "iv",  "subgroup"),
                          by.y = c("lhs", "rhs", "membership"),
                          all.y = TRUE)
      # merge all elements with group paths
      all_elem.b <- merge(all_elem, group_to_merge,
                          by.x = c("dv", "iv"),
                          by.y = c("lhs", "rhs"), 
                          all.y = TRUE)
      # put back together
      all_elem <- rbind(all_elem.a, all_elem.b)
      all_elem <- all_elem[,c("file", "dv", "iv", "beta",
                              "se", "z", "pval", "subgroup", "label")]
      colnames(all_elem)[9] <- c("level")
      all_elem   <- all_elem[order(all_elem$file, all_elem$level),]
      # new line below added 9.19.16, removes identical rows
      # identical rows may happen because we overly preserve information above
      all_elem   <- unique(all_elem)
    }
    
    # merge level into indivPathEstimates for subgroup = FALSE or all_diff_sub = TRUE
    if (agg == FALSE){
      if (subgroup == FALSE | all_diff_sub == TRUE){
        label_info         <- all_elem_final[,c("lhs", "rhs", "label")]
        for (v in 1:length(varnames)){
          label_info$lhs   <- recoderFunc(label_info$lhs, lvarnames, varnames)
          label_info$rhs   <- recoderFunc(label_info$rhs, lvarnames, varnames)
        }
        all_elem.b <- merge(all_elem, label_info,
                            by.x = c("dv", "iv"),
                            by.y = c("lhs", "rhs"), 
                            all = TRUE)
        all_elem   <- all_elem.b
        all_elem   <- all_elem[,c("file", "dv", "iv", "beta",
                                  "se", "z", "pval", "label")]
        colnames(all_elem)[8] <- c("level")
        # if any paths labeled as subgroup accidentally remain in the all_diff_sub case, 
        # just label them as individual
        all_elem$level[all_elem$level == "subgroup"] <- "ind"
        all_elem <- all_elem[order(all_elem$file, all_elem$level),]
      }
    }
    ###############################################################################
    ## end merging information into indivPathEstimates--------------------------
    
    if (agg == FALSE & !is.null(out)) {
      write.csv(all_elem, file.path(out, "indivPathEstimates.csv"), row.names = FALSE)
    } else if (agg == TRUE & !is.null(out)) {
      write.csv(all_elem, file.path(out, "allPathEstimates.csv"), row.names = FALSE)
    }
    if (agg == FALSE & !is.null(out)){
      write.csv(all_elem_summary, file.path(out, "summaryPathCounts.csv"), row.names=FALSE)
    }
    
  } else {
    if (agg == T) message(paste0("aggSEM: Model could not be estimated ",
                                 "due to nonpositive definite matrix"))
    if (agg == F) message(paste0("gimmeSEM: No models could be estimated ",
                                 "due to nonpositive definite matrices for all individuals"))
  }
  # just added on 11/3/16
  est <- list(fit       = all_fit,
              param_est = all_elem,
              all_plot  = all_plot,
              sub_plots = sub_plots,
              all_counts = final.paths.all,
              sub_paths  = sub_paths)
  return(est)
  # final <- list(all_ind = all_ind_paths)
  # class(final) <- "gimmep"
  # return(final)
}

################################################################################
# misem_evalbetas_sub_funcs.R
################################################################################
## these two functions use the miSEM and evalbetas functions
## at the subgroup level

miSEMsub <- function (subsetup.out,
                      setup.out,
                      previous.out,
                      second_round,
                      evalbetassub.out) {
  
  ts_list = setup.out$ts_list 
  membership = NULL
  
  if (second_round == FALSE) miSEMsub.out = NULL
  n_subgroups         = subsetup.out$n_subgroups
  subgroup_membership = subsetup.out$subgroup_membership
  varnames            = setup.out$varnames
  syntax              = previous.out$syntax
  n_group_paths   = previous.out$n_group_paths
  
  unique_syntax       <- matrix(seq(1:n_subgroups), ncol = 4, nrow = n_subgroups)
  subgroup_membership <- as.data.frame(subgroup_membership)
  if (second_round == TRUE) {
    syntax.all <- cbind(subsetup.out$subgroup_membership$membership, previous.out$syntax)
    syntax.all <- syntax.all[!duplicated(syntax.all), ]
  }
  
  for (s in 1:n_subgroups){
    if (second_round == TRUE){
      indx   <- syntax.all[,1] == s
      syntax <- syntax.all[indx,2]
      indx2  <- evalbetassub.out$sub.paths[,1] == s
      subgroup_paths <- evalbetassub.out$sub.paths[indx2,2]
      names(subgroup_paths) <- NULL
      
    }
    sub.subjids     <- as.character(subset(subgroup_membership, membership == s)[ ,1])
    subjects        <- length(sub.subjids)
    ts_sub          <- ts_list[sub.subjids]
    
    if (subjects > 1){
      subgroup.miSEM <- miSEM(setup.out         = setup.out,
                              previous.out      = previous.out,
                              subgroup.step     = TRUE,
                              subjects          = subjects,
                              ts_list           = ts_sub,
                              syntax            = syntax,
                              subgroup_paths    = subgroup_paths,
                              second_round      = second_round)
      
      unique_syntax[s,2] <- subgroup.miSEM$syntax
      unique_syntax[s,3] <- subgroup.miSEM$n_subgroup_paths
      unique_syntax[s,4] <- paste0(subgroup.miSEM$subgroup_paths, collapse = ",")
      
    } else {
      unique_syntax[s,2] <- syntax
      unique_syntax[s,3] <- 0
      unique_syntax[s,4] <- ""}
  }
  list <- list("unique_syntax"     = unique_syntax,
               "n_group_paths" = n_group_paths)
  return(list)
}

evalbetassub <- function (subsetup.out,
                          previous.out,
                          setup.out,
                          evalbetas.out) {
  
  n_subgroups         = subsetup.out$n_subgroups
  subgroup_membership = subsetup.out$subgroup_membership
  modularity          = subsetup.out$modularity
  ar.paths            = setup.out$ar.paths
  ar                  = setup.out$ar
  out                 = setup.out$out
  ts_list             = setup.out$ts_list
  
  membership = NULL
  
  subgroup_membership <- as.data.frame(subgroup_membership)
  final_syntax        <- matrix(seq(1:n_subgroups), ncol=2, nrow = n_subgroups)
  sub.paths           <- matrix(seq(1:n_subgroups), ncol=2, nrow = n_subgroups)
  
  for (s in 1:n_subgroups){
    sub.subjids     <- as.character(subset(subgroup_membership, membership == s)[ ,1])
    subjects        <- length(sub.subjids)
    ts_sub          <- ts_list[sub.subjids]
    
    if (subjects > 1){
      evalbetas.out.sub <- evalbetas(setup.out      = setup.out,
                                     previous.out   = previous.out,
                                     subgroup.step  = TRUE,
                                     s              = s,
                                     subjects       = subjects,
                                     ts_list        = ts_sub,
                                     post_sub_prune = FALSE)
      
      final_syntax[s,2] <- evalbetas.out.sub$syntax
      if (length(evalbetas.out.sub$subgroup_paths) != 0) {
        sub.paths[s,2]  <- evalbetas.out.sub$subgroup_paths
      } else {
        sub.paths[s,2] <- NA
      }
    } else { 
      final_syntax[s,2] <- evalbetas.out$syntax
      sub.paths[s,2]   <- NA
    }
  }
  
  colnames(final_syntax) <- c("membership", "syntax")
  syntax_sub             <- merge(subgroup_membership, 
                                  final_syntax, by = "membership", all = TRUE)
  syntax_sub             <- syntax_sub[order(syntax_sub$subjid), ]
  syntax_sub$files       <- as.character(syntax_sub$subjid)
  syntax_sub$syntax      <- as.character(syntax_sub$syntax)
  list <- list("syntax_sub"   = syntax_sub,
               "final_syntax" = final_syntax,
               "sub.paths"    = sub.paths,
               "modularity"   = modularity)
  return(list)
}

## runs each individual using group-level model 
subsetup <- function (setup.out,
                      previous.out) {
  
  op         = NULL
  modularity = NULL
  param      = NULL
  
  subjects        = setup.out$subjects
  varnames        = setup.out$varnames
  vars            = setup.out$vars
  candidate_paths = setup.out$candidate_paths
  n_fixed_paths   = setup.out$n_fixed_paths
  ts_list         = setup.out$ts_list
  syntax          = previous.out$syntax
  syntax_paths    = unlist(strsplit(syntax, "\n"))
  n_group_paths   = previous.out$n_group_paths
  out             <- setup.out$out
  cutoffind       <- setup.out$cutoffind
  cutoffgroup     <- qchisq(1-.05/subjects, 1)
  cutoffmi        <- qchisq(1-(.05/((vars*(vars-1)/2)*subjects)), 1)
  
  mi_matrix    <- matrix(0, ncol = subjects, nrow = ((vars*(vars-1)/2)))
  epc_matrix   <- matrix(0, ncol = subjects, nrow = ((vars*(vars-1)/2)))
  p_matrix     <- matrix(0, ncol = subjects, nrow = (n_group_paths + n_fixed_paths))
  beta_matrix  <- matrix(0, ncol = subjects, nrow = (n_group_paths + n_fixed_paths))
  ## create matrix of MIs individual by individual after evalbetas
  
  for (k in 1:subjects)
  {
    # data.file           <- read.data(file   = files[k],
    #                                  sep    = sep,
    #                                  header = header)
    # colnames(data.file) <- c(varnames)
    fit <- fit.model(syntax    = syntax,
                     data.file = ts_list[[k]])
    
    check_npd <- any(grepl("error", class(fit)) == TRUE)
    
    if (check_npd == FALSE) {
      all_mi         <- tryCatch(modindices(fit), error = function(e) e)
      check_singular <- any(grepl("singular", all_mi) == TRUE)
      empty_mi       <- ifelse(nrow(all_mi) == 0, TRUE, FALSE)
      if (check_singular == TRUE) empty_mi <- TRUE
      converge       <- lavInspect(fit, "converged")
    } else {
      check_singular <- TRUE
      converge       <- FALSE
      empty_mi       <- TRUE
    }
    
    # printing replication
    writeLines(paste("subgroup search, subject", k))
    ## grab MIs and grab beta values; put into preallocated matrices
    if (check_singular == FALSE & converge == TRUE & empty_mi == FALSE) {
      mi           <- all_mi[all_mi$op == "~",][,c("lhs","op","rhs","mi","epc")]
      mi$param     <- paste0(mi$lhs, mi$op, mi$rhs)
      mi           <- mi[-c(1:3)]
      mi           <- subset(mi, param %in% candidate_paths)
      
      mi$mi        <- as.numeric(as.character(mi$mi))
      mi$epc       <- as.numeric(as.character(mi$epc))
      
      mi$mi[mi$param %in% syntax_paths]  <- 0
      mi$epc[mi$param %in% syntax_paths] <- 0
      
      mi$thresh    <- ifelse(mi$mi > cutoffmi, 1, 0)
      beta         <- as.data.frame(subset(standardizedSolution(fit), op == "~"))
      beta$thresh  <- ifelse(beta$pvalue < (.05/subjects), 1, 0)
    }
    # if it doesn't converge or is computationally singular
    if (converge == FALSE | check_singular == TRUE | check_npd == TRUE | empty_mi == TRUE) {
      mi_matrix[,k]     <- NA
      epc_matrix[,k]    <- NA
      p_matrix[,k]      <- NA
      beta_matrix[,k]   <- NA
    } else {
      mi_matrix[1:(nrow(mi)),k]  <- mi$thresh
      epc_matrix[1:(nrow(mi)),k] <- mi$epc
      p_matrix[,k]               <- beta$thresh
      beta_matrix[,k]            <- beta$est.std
    }
  }
  
  all_val_matrix           <- rbind(epc_matrix, beta_matrix)
  all_p_matrix             <- rbind(mi_matrix, p_matrix)
  
  colnames(all_val_matrix) <- names(ts_list)
  colnames(all_p_matrix)   <- names(ts_list)
  
  all_val_matrix        <- all_val_matrix[,colSums(is.na(all_val_matrix)) != nrow(all_val_matrix)]
  all_p_matrix          <- all_p_matrix[,colSums(is.na(all_p_matrix)) != nrow(all_p_matrix)]
  names                 <- colnames(all_val_matrix)
  
  sim <- matrix(0, ncol = ncol(all_val_matrix), nrow = ncol(all_val_matrix))
  
  ## creates similarity matrix to pass to walktrap
  ## counts if sign and significance are the same for EPC or beta weight
  for (i in 1:ncol(all_p_matrix)){
    for (j in 1:ncol(all_p_matrix)){
      ind1 <- as.matrix(all_p_matrix[ ,i])
      ind2 <- as.matrix(all_p_matrix[ ,j])
      val1 <- as.matrix(all_val_matrix[ ,i])
      val2 <- as.matrix(all_val_matrix[ ,j])
      sim[i,j] <- sum(ind1 == 1 & ind2 == 1 & sign(val1) == sign(val2), na.rm = TRUE)
    }
  }
  
  colnames(sim) <- names
  sim           <- sim - min(sim, na.rm = TRUE)
  diag(sim)     <- 0
  ## here is where walktrap runs on the similarity matrix and gives you subgroups
  g                   <- graph.adjacency(sim, mode = "undirected")
  membership          <- walktrap.community(g, steps = 4)$membership
  modularity_value    <- modularity(walktrap.community(g, steps = 4))
  
  subgroup_membership <- cbind(names, membership)
  n_subgroups         <- length(table(subgroup_membership[ ,2]))
  filenames           <- as.matrix(names(ts_list))
  subgroup_membership <- merge(filenames, subgroup_membership,
                               by.x = "V1", by.y = "names", all.x = TRUE)
  
  colnames(subgroup_membership) <- c("subjid", "membership")
  list <- list("epc_beta_matrix"     = all_val_matrix,
               "epc_beta_p_matrix"   = all_p_matrix,
               "sim"                 = sim,
               "subgroup_membership" = subgroup_membership,
               "n_subgroups"         = n_subgroups,
               "modularity"          = modularity_value)
  return(list)
}

################################################################################
# helper functions
################################################################################
W2E <-function(x) cbind(which(x != 0, arr.ind = TRUE), x[x != 0])

recoderFunc <- function(data,
                        oldvalue,
                        newvalue){
  if (is.factor(data))     data     <- as.character(data)
  if (is.factor(oldvalue)) oldvalue <- as.character(oldvalue)
  if (is.factor(newvalue)) newvalue <- as.character(newvalue)
  newvec <- data
  for (i in unique(oldvalue)) newvec[data == i] <- newvalue[oldvalue == i]
  newvec
}
################################################################################
