#' @name gimmeSEM
#' @aliases gimme gimmeSEM
#' @title Group iterative multiple model estimation.
#' @description This function identifies structural equation models for each
#' individual that consist of both group-level and individual-level paths.
#' @usage
#' gimmeSEM(data        = "",
#'          out         = "",
#'          sep         = "",
#'          header      = ,
#'          ar          = TRUE,
#'          plot        = TRUE,
#'          subgroup    = FALSE,
#'          sub_feat    = "lag & contemp",
#'          confirm_subgroup = NULL,
#'          paths       = NULL,
#'          exogenous   = NULL,
#'          ex_lag      = FALSE,
#'          mult_vars   = NULL,
#'          mean_center_mult = FALSE,
#'          standardize = FALSE,
#'          groupcutoff = .75,
#'          subcutoff   = .5,
#'          diagnos     = FALSE)
#' @param data The path to the directory where the data files are located,
#' or the name of the list containing each individual's time series. Each file
#' or matrix must contain one matrix for each individual containing a T (time)
#' by p (number of variables) matrix where the columns represent variables and
#' the rows represent time. Individuals must have the same variables (p)
#' but can have different lengths of observations (T).
#' @param out The path to the directory where the results will be stored
#' (optional). If specified,
#' a copy of output files will be replaced in directory. If directory at
#' specified path does not exist, it will be created.
#' @param sep The spacing of the data files. 
#' "" indicates space-delimited,
#' "/t" indicates tab-delimited, "," indicates comma delimited. Only necessary
#' to specify if reading data in from physical directory.
#' @param header Logical. Indicate TRUE for data files with a header. Only
#' necessary to specify if reading data in from physical directory.
#' @param ar Logical. If TRUE, begins search for group model with
#' autoregressive (AR) paths freed for estimation. Defaults to TRUE.
#' @param paths \code{lavaan}-style syntax containing paths with which
#' to begin model estimation (optional). That is, Y~X indicates that Y
#' is regressed on X, or X predicts Y. If no header is used,
#' then variables should be referred to with V followed (with no separation)
#' by the column number. If a
#' header is used, variables should be referred to using variable names.
#' To reference lag variables, "lag" should be added to the end of the variable
#' name with no separation. Defaults to NULL.
#' @param exogenous Vector of variable names to be treated as exogenous (optional).
#' That is, exogenous variable X can predict Y but cannot be predicted by Y.
#' If no header is used, then variables should be referred to with V followed
#' (with no separation) by the column number. If a header is used, variables
#' should be referred to using variable names. Defaults to NULL.
#' @param ex_lag Logical.  If true, lagged variables are created for exogenous variables.  
#' Defaults to FALSE.
#' @param mult_vars Vector of variable names to be multiplied to explore bilinear/modulatory
#' effects (optional). All multiplied variables will be treated as exogenous (X can predict
#' Y but cannot be predicted by Y). Within the vector, multiplication of two variables should be
#' indicated with an asterik (e.g. V1*V2). If no header is used, variables should be referred to with 
#' V followed by the column number (with no separation). If a header is used, each variable should be
#' referred to using variable names. If multiplication with the lag 1 of a variable is desired, the 
#' variable name should be followed by "lag" with no separation (e.g. V1*V2lag). Note that if
#' multiplied variables are desired, at least one variable in the dataset must be specified as exogenous.
#' Defaults to NULL.
#' @param mean_center_mult Logical. If TRUE, the variables indicated in mult_vars will be mean-centered
#' before being multiplied together. Defaults to FALSE. 
#' @param standardize Logical. If TRUE, all variables will be standardized to have a mean of zero and a
#' standard deviation of one. Defaults to FALSE. 
#' @param plot Logical. If TRUE, graphs depicting relations among variables
#' of interest will automatically be
#' created. Solid lines represent contemporaneous relations (lag 0) and dashed lines reflect 
#' lagged relations (lag 1). For individual-level plots, red paths represent positive weights
#' and blue paths represent negative weights. Width of paths corresponds to estimated path weight.
#' For the group-level plot, black represents group-level paths, grey represents
#' individual-level paths, and (if subgroup = TRUE)
#' green represents subgroup-level paths. For the group-level plot,
#' the width of the edge corresponds to the count. Defaults to TRUE.
#' @param subgroup Logical. If TRUE, subgroups are generated based on
#' similarities in model features using the \code{walktrap.community}
#' function from the \code{igraph} package. Defaults to FALSE. 
#' @param confirm_subgroup Dataframe. Option only available when subgroup = TRUE. Dataframe should contain two columns. The first
#' column should specify file labels (the name of the data files without file extension), 
#' and the second should contain integer values (beginning at 1) 
#' specifying the subgroup membership for each individual.
#' function from the \code{igraph} package. Defaults to TRUE. 
#' @param sub_feature Option to indicate feature(s) used to subgroup individuals. Defaults to
#' "lag & contemp" for lagged and contemporaneous, which is the original method. Can use 
#' "lagged" or "contemp" to subgroup solely on features related to lagged and contemporaneous 
#' relations, respectively.
#' @param confirm_subgroup Dataframe. If subgroup is also TRUE, option to provide
#' subgroup labels contained in the dataframe. Dataframe has 2 columns,
#' the first referring to file labels (without extensions), and the second an integer variable referring to subgroup label.
#' @param groupcutoff Cutoff value for group-level paths. Defaults to .75,
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
#'  If subgroup = TRUE and subgroups are found, then a column is present
#'  containing the subgroup membership for each individual. Also contains the
#'  level at which each path was estimated: group, subgroup, or individual.}
#'  \item{\strong{summaryFit}} {Contains model fit information for individual-
#'  level models. If subgroups are requested, this file also contains the
#'  subgroup membership for each individual.}
#'  \item{\strong{summaryPathCountMatrix}} Contains counts of total number of
#'  paths, both contemporaneous and lagged, estimated for the sample. The row
#'  variable is the outcome and the column variable is the predictor variable.
#'  \item{\strong{summaryPathCounts}} {Contains summary count information for
#'  paths identified at the group-, subgroup (if subgroup = TRUE), and
#'  individual-level.}
#'  \item{\strong{summaryPathsPlot}} {Produced if plot = TRUE. Contains figure
#'  with group, subgroup (if subgroup = TRUE), and individual-level paths
#'  for the sample. Black paths are group-level, green paths are subgroup-level,
#'  and grey paths are individual-level, where the thickness of the line
#'  represents the count.}
#'  }
#'  In subgroup output directory (if subgroup = TRUE):
#'  \itemize{
#'  \item{\strong{subgroup\emph{k}PathCounts}} Contains counts of relations
#'  among lagged and contemporaneous variables for the \strong{\emph{k}}th
#'  subgroup.
#'  \item{\strong{subgroup\emph{k}Plot}} Contains plot of group, subgroup,
#'  and individual level paths for the \strong{\emph{k}}th subgroup.
#'  Black represents group-level paths, grey represents individual-level paths,
#'  and green represents subgroup-level paths.
#'  }
#'  Note: if a subgroup of size n = 1 is discovered, subgroup-level output is
#'  not produced. \cr

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
#' @references Lane, S.T. & Gates, K.M. (2017). Automated selection of robust
#' individual-level structural equation models for time series data.
#' Structural Equation Modeling.
#' @author Stephanie Lane
#' @examples
#'  \dontrun{
#' paths <- 'V2 ~ V1
#'           V3 ~ V4lag'
#'
#' fit <- gimmeSEM(data     = simData,
#'                 out      = "C:/simData_out",
#'                 subgroup = TRUE, 
#'                 paths    = paths)
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
                              sub_feature    = "lag & contemp",
                              confirm_subgroup = NULL,
                              paths          = NULL,
                              exogenous      = NULL,
                              ex_lag         = FALSE,
                              mult_vars      = NULL,
                              mean_center_mult = FALSE,
                              standardize    = FALSE,
                              groupcutoff    = .75,
                              subcutoff      = .5,
                              diagnos        = FALSE){

  sub_membership = NULL

  dat         <- setup(data                 = data,
                       sep                  = sep,
                       header               = header,
                       out                  = out,
                       plot                 = plot,
                       ar                   = ar,
                       paths                = paths,
                       exogenous            = exogenous,
                       ex_lag               = ex_lag,
                       mult_vars            = mult_vars,
                       mean_center_mult     = mean_center_mult,
                       standardize          = standardize,
                       subgroup             = subgroup,
                       ind                  = FALSE,
                       agg                  = FALSE,
                       groupcutoff          = groupcutoff,
                       subcutoff            = subcutoff)
  
  #Error Check for Confirm Subgroup Labels
  if(subgroup & !is.null(confirm_subgroup)){
    if(dim(confirm_subgroup)[[2]] != 2){
      stop(paste0("gimme ERROR: confirmatory subgroup dataframe is not a two column dataframe.",
                  " Please ensure that the confirmatory subgroup dataframe consists of a column of filenames and a column of community assignments."))
    }
    if(length(match(confirm_subgroup[,1], (dat$file_order)$names)) != dim(confirm_subgroup)[[1]]){
      stop(paste0("gimme ERROR: confirmatory subgroup dataframe contains mismatched filenames.",
                  " Please ensure that the confirmatory subgroup filenames match the data filenames, sans extensions (Example: sub_1000, not sub_1000.csv)"))
    }
    if(!is.numeric(confirm_subgroup[,2])){
      stop(paste0("gimme ERROR: confirmatory subgroup assignments are non-numeric.",
                  " Please ensure that the confirmatory subgroup assignments are integer valued, beginning from 1. (Example: 1, 2, 3, 4)"))
    }
  }
  
  
  grp <- list("n_group_paths" = 0,
              "n_fixed_paths" = length(dat$fixed_paths),
              "group_paths"   = c())

  s1  <- search.paths(base_syntax  = dat$syntax,
                     fixed_syntax = NULL,
                     add_syntax   = grp$group_paths,
                     n_paths      = grp$n_group_paths,
                     data_list    = dat$ts_list,
                     elig_paths   = dat$candidate_paths,
                     prop_cutoff  = dat$group_cutoff,
                     n_subj       = dat$n_subj,
                     chisq_cutoff = qchisq(1-.05/dat$n_subj, 1),
                     subgroup_stage = FALSE)

  grp[c("n_group_paths", "group_paths")] <- s1

  prune <- ifelse(grp$n_group_paths != 0, TRUE, FALSE)

  if (prune){
    s2 <- prune.paths(base_syntax  = dat$syntax,
                      fixed_syntax = NULL,
                      add_syntax   = grp$group_paths,
                      data_list    = dat$ts_list,
                      n_paths      = grp$n_group_paths,
                      n_subj       = dat$n_subj,
                      prop_cutoff  = dat$group_cutoff,
                      elig_paths   = grp$group_paths,
                      subgroup_stage = FALSE)
    grp[c("n_group_paths", "group_paths")] <- s2
  }

  # determine subgroup assignments, if requested

  if (subgroup){
    sub <- determine.subgroups(base_syntax  = c(dat$syntax, grp$group_paths),
                               data_list    = dat$ts_list,
                               n_subj       = dat$n_subj,
                               chisq_cutoff = dat$chisq_cutoff_mi_epc,
                               file_order   = dat$file_order,
                               elig_paths   = dat$candidate_paths,
                               confirm_subgroup = confirm_subgroup,
                               out_path     = dat$out, 
                               sub_feature  = sub_feature)

  # begin subgroup-level search for paths ------------------------------------ #

  sub_spec <- vector("list", sub$n_subgroups)

  for (s in 1:sub$n_subgroups){

    sub_s <- list(sub_paths     = character(),
                  n_sub_paths   = 0,
                  sub_s_subjids = subset(sub$sub_mem,
                                         sub_membership == s)[ ,"names"],
                  n_sub_subj    = sum(sub$sub_mem$sub_membership == s,
                                      na.rm = TRUE),
                  sub_membership    = s)

    if (sub_s$n_sub_subj > 1){
      s4 <- search.paths(base_syntax  = dat$syntax,
                         fixed_syntax = grp$group_paths,
                         add_syntax   = character(),
                         n_paths      = 0,
                         data_list    = dat$ts_list[sub_s$sub_s_subjids],
                         elig_paths   = dat$candidate_paths,
                         prop_cutoff  = dat$sub_cutoff,
                         n_subj       = sub_s$n_sub_subj,
                         chisq_cutoff = qchisq(1-.05/sub_s$n_sub_subj, 1),
                         subgroup_stage = TRUE)
      sub_s[c("n_sub_paths", "sub_paths")] <- s4
    }
    sub_spec[[s]] <- sub_s
  }
  # end subgroup-level search for paths -------------------------------------- #

  # begin subgroup-level pruning --------------------------------------------- #
  for (s in 1:sub$n_subgroups){
    prune <- sub_spec[[s]]$n_sub_paths != 0 & sub_spec[[s]]$n_sub_subj != 1
    if(prune){
      s5 <- prune.paths(base_syntax  = dat$syntax,
                        fixed_syntax = grp$group_paths,
                        add_syntax   = sub_spec[[s]]$sub_paths,
                        data_list    = dat$ts_list[sub_spec[[s]]$sub_s_subjids],
                        n_paths      = sub_spec[[s]]$n_sub_paths,
                        n_subj       = sub_spec[[s]]$n_sub_subj,
                        prop_cutoff  = dat$sub_cutoff,
                        elig_paths   = sub_spec[[s]]$sub_paths,
                        subgroup_stage = TRUE)
      sub_spec[[s]][c("n_sub_paths", "sub_paths")] <- s5
    }
  }

  # begin second-round group-level pruning ----------------------------------- #
  prune <- any(lapply(sub_spec, FUN = function(x) x$n_sub_paths != 0) == TRUE)

  sub_spec_comb <- do.call(rbind, sub_spec)
  ind           <- merge(sub$sub_mem, sub_spec_comb, "sub_membership", all.x = TRUE)
  ind           <- ind[order(ind$index),]
  ind$sub_paths[is.na(ind$sub_paths)] <- ""
  temp_count    <- grp$n_group_paths

  if (prune){
    s6 <- prune.paths(base_syntax  = dat$syntax,
                      fixed_syntax = ind$sub_paths,
                      add_syntax   = grp$group_paths,
                      data_list    = dat$ts_list,
                      n_paths      = grp$n_group_paths,
                      n_subj       = dat$n_subj,
                      prop_cutoff  = dat$group_cutoff,
                      elig_paths   = grp$group_paths,
                      subgroup_stage = FALSE)

    grp[c("n_group_paths", "group_paths")] <- s6
  }

  if (temp_count != grp$n_group_paths){
    temp_sub_spec <- sub_spec
    for (s in 1:sub$n_subgroups){
      if (sub_spec[[s]]$n_sub_subj > 1){
        s7 <- search.paths(base_syntax  = dat$syntax,
                           fixed_syntax = grp$group_paths,
                           add_syntax   = sub_spec[[s]]$sub_paths,
                           n_paths      = sub_spec[[s]]$n_sub_paths,
                           data_list    =
                             dat$ts_list[sub_spec[[s]]$sub_s_subjids],
                           elig_paths   = dat$candidate_paths,
                           prop_cutoff  = dat$sub_cutoff,
                           n_subj       = sub_spec[[s]]$n_sub_subj,
                           chisq_cutoff =
                             qchisq(1-.05/sub_spec[[s]]$n_sub_subj, 1),
                           subgroup_stage = FALSE)
        sub_spec[[s]][c("n_sub_paths", "sub_paths")] <- s7
      }
    }

    if (!identical(temp_sub_spec, sub_spec)){
      for (s in 1:sub$n_subgroups){
        prune <- temp_sub_spec[[s]]$n_sub_paths != sub_spec[[s]]$n_sub_paths
        if(prune){
          s8 <- prune.paths(base_syntax  = dat$syntax,
                            fixed_syntax = grp$group_paths,
                            add_syntax   = sub_spec[[s]]$sub_paths,
                            data_list    =
                              dat$ts_list[sub_spec[[s]]$sub_s_subjids],
                            n_paths      = sub_spec[[s]]$n_sub_paths,
                            n_subj       = sub_spec[[s]]$n_sub_subj,
                            prop_cutoff  = dat$sub_cutoff,
                            elig_paths   = sub_spec[[s]]$sub_paths,
                            subgroup_stage = FALSE)
          sub_spec[[s]][c("n_sub_paths", "sub_paths")] <- s8
        }
      }
    }
  }

  sub_spec_comb <- do.call(rbind, sub_spec)
  ind           <- merge(sub$sub_mem, sub_spec_comb, "sub_membership", all.x = TRUE)
  ind$sub_paths[is.na(ind$sub_paths)] <- ""
  ind           <- ind[order(ind$index),]

  } else {
    # create ind object here if no subgrouping takes place
    sub      <- NULL
    sub_spec <- NULL
    ind      <- dat$file_order
  }

  # individual-level search
  store <- indiv.search(dat, grp, ind)

  print.gimme(x = sub,
              y = subgroup,
              z = dat)

  # wrap-up and create output
  final <- final.org(dat,
                     grp,
                     ind = store$ind,
                     sub,
                     sub_spec,
                     store)

  # these objects are used in print.gimmep
  # if you change an object name here, 
  # you need to change it in the print.gimmep.R
  res <- list(path_est_mats   = store$betas,
              varnames        = dat$varnames,
              n_rois          = dat$n_rois,
              fit             = final$fit,
              path_se_est     = final$param_est,
              plots           = store$plots,
              group_plot      = final$samp_plot,
              sub_plots       = final$sub_plots,
              subTF           = subgroup,
              path_counts     = final$sample_counts,
              path_counts_sub = final$sub_counts,
              vcov            = store$vcov,
              sim_matrix      = sub$sim, 
              syntax          = dat$syntax
              )
  class(res) <- "gimmep"

  invisible(res)
}

print.gimme <- function(x, y, z){
  writeLines("gimme finished running normally")
  if (!is.null(z$out)) writeLines(paste("output is stored in", z$out))
  if (y == TRUE) {
    writeLines(paste("Number of subgroups =", x$n_subgroups))
    writeLines(paste("Modularity =", round(x$modularity, digits = 5)))
  }
}
