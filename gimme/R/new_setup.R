## this setup function creates many values that later code refers back to
setup <- function (data,
                   sep,
                   header,
                   out,
                   plot,
                   ar,
                   paths,
                   exogenous,
                   mult_vars,
                   mean_center_mult,
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
    files <- list.files(data, full.names = TRUE) #Creates a list of all files in the dir pointed at
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
      ts_list[[i]] <- ts_ind #creates a list of the TS of all files, with each sub as a slice
    }
    varnames       <- colnames(ts_ind) 
    names(ts_list) <- tools::file_path_sans_ext(basename(files)) 
    rois           <- ncol(ts_ind) 
    n_endog        <- ncol(ts_ind) ### ADDED TO DEAL WITH "ROIS" ISSUE LATER IN THE SCRIPT
  } else if (is.list(data)){
    ts_list  <- list()
    ts_list  <- data
    rois     <- ncol(ts_list[[1]])
    n_endog  <- ncol(ts_list[[1]])
    varnames <- colnames(ts_list[[1]])
    if (is.null(varnames)){
      varnames <- c(paste0("V", seq(1,rois))) ###HERE THE VARS IN A LIST ARE NAMED X1 INSTEAD OF V1-- THIS IS MY PROBLEM; CHANGED TO V TO TRY A BANDAID FIX FOR RIGHT NOW
      ts_list <- lapply(ts_list, function(x) { 
        colnames(x)<-varnames 
        x 
      })
      }
  }
  
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
  
  ### This loop reads in mulitplied variables if there are any, and splits them and determines which variables are being multiplied.
  ### For each subject it finds these variables in their data, mean centers it if desired, and then multiplies the two together,
  ### Renames it based on the original inputted name, and then binds it to the subject's data. Once this is done for each sub,
  ### It will do the same for next multiplied variable, if there are any. 
  ### TO ACCOUNT FOR LIST PROBLEM, WILL SOMEONE HAVE TO ADD IN A LOOP/RECODE/CHECK HERE THAT CHANGES THEM TO V INSTEAD OF X
  if (!is.null(mult_vars)){
    for(i in 1:length(mult_vars)){ 
      mult_pairs <- mult_vars[[i]]
      vars_to_mult <- strsplit(mult_pairs, "*", fixed = TRUE)
      vars_to_mult_mat <- unlist(vars_to_mult)
      factor_1 <- vars_to_mult_mat[1]
      factor_2 <- vars_to_mult_mat[2]
      for(p in 1:length(ts_list)){
        all <- ts_list[[p]]
        if (mean_center_mult == TRUE){
          var_1 <- all[,factor_1]
          var_1_center <- scale(var_1, scale = FALSE)
          var_2 <- all[,factor_2]
          var_2_center <- scale(var_2, scale = FALSE)
          multiplied_center <- var_1_center*var_2_center
          df_tobind <- data.frame(multiplied_center)
          colnames(df_tobind) <- paste0(factor_1,"by",factor_2)
          all_appended <- cbind(all,df_tobind)
          ts_list[[p]] <- all_appended
        } else{
          var_1 <- all[,factor_1]
          var_2 <- all[,factor_2]
          multiplied <- var_1*var_2
          df_tobind <- data.frame(multiplied)
          colnames(df_tobind) <- paste0(factor_1,"by",factor_2)
          all_appended <- cbind(all,df_tobind)
          ts_list[[p]] <- all_appended
        }
      }
    }
  }
 ###Problem line is now currently that it is unable to index the column values based on names when I run 
  ##from gimmesem(). I get this error:  Error in `[.data.frame`(all, , factor_1) : undefined columns selected .
  ##However it works perfectly when I just run this code....not sure what to do! 
  
  ### Added this loop to rename the multiplied variables to their latent lvarnames, since this is what is used for candidate paths
  lmult_pairs <- NULL
  if (!is.null(mult_vars)){
    for(i in 1:length(mult_vars)){ 
    mult_pairs <- mult_vars[[i]]
    vars_to_mult <- strsplit(mult_pairs, "*", fixed = TRUE)
    vars_to_mult_mat <- unlist(vars_to_mult)
    lvars_to_mult <- recode.vars(vars_to_mult_mat, varnames, lvarnames)
    lfactor_1 <- lvars_to_mult[1]
    lfactor_2 <- lvars_to_mult[2]
    lmult_name <- paste0(lfactor_1,"by",lfactor_2)
    lmult_pairs[[i]] <- lmult_name
    }
  }   
  n_bilinear <- length(lmult_pairs) ###Added to count the number of bilinear/multiplied variables
  
  lexog<- NULL
  if (!is.null(exogenous)){
    lexog <- recode.vars(exogenous, varnames, lvarnames) 
  }
  n_exog <- length(lexog) ###Added to count the number of exogenous variables
  
  lexogenous <- c(lexog,lmult_pairs)
  n_exog_total <- length(lexogenous) ###Added to count the combined total number of exogenous variables
  
  all <- ts_list[[1]] ###I only named it this because it's used below (at L162) for a reason that's unclear to me so decided to keep it. Not sure if I can just move that up here?
  varnames <- colnames(all) ###renamed to account for any new mult vars
  lvarnames <- c(lvarnames,lmult_pairs) ###renamed to account for any new mult vars. Combines original latent names plus latent names for multiplied vars
  n_vars_total <- length(varnames)
  n_contemporaneous <- n_endog + n_bilinear
  
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
  
  subjects         <- length(ts_list) #Num of subs is just the num of subs in ts_list
  file_order       <- data.frame(index = c(seq(1:subjects)), 
                                 names = c(names(ts_list)), 
                                 stringsAsFactors = FALSE) #indexes the subjects from 1:num sub based on order in ts_list
  
  all              <- ts_list[[1]] #Not sure yet why it's doing this, this is not called again. Have no idea what it's doing here. 

  if (rois == 1) {
    stop(paste0("gimme ERROR: only one column of data read in. ",
                "Check if sep argument properly specified."))
  }
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
  if (subjects == 1 & !ind) {
    stop(paste0('gimme ERROR: only one subject detected in data directory. ',
                'Please use indSEM function instead.'))
  }
  
  vars              <- length(varnames) ###Changed from rois*2 to length of ALL varnames, to account for any new multiplied vars
  cutoffind         <- qchisq(.99, 1)
  n_group_paths     <- 0
  
  x                 <- seq(1:vars) #makes the x variables the var NUMBERS only (eg 1, 2, 3,...,10)
  y                 <- substring(lvarnames, 4) #this makes the y variables the lvarname SUFFIXES only, by taking the lvarnames begining at char. 4 (so 1lag, 2lag,...,4,5)
  ###This will have to be changed if we need this, but do we need this???? What is this doing????
  individual        <- file.path(out, "individual")
  subgroup_dir      <- file.path(out, "subgroup")
  
  if (subgroup & !is.null(out))  {
    dir.create(subgroup_dir, showWarnings = FALSE)
  }
  if (!agg & !is.null(out)) {
    dir.create(individual, showWarnings = FALSE)
  }
  if (plot) {
    plot.names <- varnames[(rois+1):(vars)] 
    ###THIS MIGHT HAVE TO BE CHANGED AGAIN
  } else {
    plot.names <- ""
  }
  
  # check to make sure variables in exogenous argument exist in data
  if(!is.null(exogenous)){
    for(exog in exogenous)
    if (!exog %in% varnames){
      stop(paste0('gimme ERROR: Exogenous variable name not in data column names
                  Please fix.'))
    }
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
      dvsFree      <- recode.vars(tableFree$lhs, varnames, lvarnames)
      ivsFree      <- recode.vars(tableFree$rhs, varnames, lvarnames)
      
      # check if any exogenous variables have been incorrectly specified
      # for free paths
      if(!is.null(exogenous)){
      for (exog in lexogenous){
        if (exog %in% dvsFree){
          stop(paste0('gimme ERROR: an exogenous variable was treated as endogenous in 
                      specified paths.  Please remove variable from exogenous list or 
                      correct path specification'))
          #This is because exog can't be treated as DV. We will have to add this/make sure other things we're adding (e.g. multiplied paths) are also considered in this
        }
      }
      }
      
      if (nrow(tableFree) != 0){
        vsFree       <- paste0(dvsFree, "~", ivsFree)
      } else vsFree <- NULL
      # table up the paths which are fixed to a certain value by the user
      tableFixed   <- table[table$op == "~" & table$free == 0,]
      if (nrow(tableFixed) > 0){
        dvsFixed     <- recode.vars(tableFixed$lhs, varnames, lvarnames)
        
        # check if any exogenous variables have been incorrectly specified
        # for fixed paths
        if (!is.null(lexogenous)){
        for (exog in lexogenous){
          if (exog %in% dvsFixed){
            stop(paste0('gimme ERROR: an exogenous variable was treated as endogenous in 
                        specified paths.  Please remove variable from exogenous list or 
                        correct path specification'))
          }
        }
        }
        
        ivsFixed     <- recode.vars(tableFixed$rhs, varnames, lvarnames)
        vsFixed      <- paste0(dvsFixed, "~", ivsFixed)
      } else {
        vsFixed <- NULL
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
  line1 <- paste0(lvarnames, "=~1*", varnames) 
  
  # covary all exogenous lagged variables
  line2 <- apply(expand.grid.unique(lvarnames[1:rois], lvarnames[1:rois], incl.eq = TRUE), 
                 1, paste, collapse = "~~")
  
  # estimate variance of single-indicator latent variables 
  line3 <- paste0(lvarnames[(rois+1):vars], "~~", lvarnames[(rois+1):vars]) 
  
  # freely estimate autoregressive relationships if ar = TRUE
  # if ar = FALSE, set up nonsense paths fixed to zero 
  if (ar == TRUE) {
    line4 <- paste0(lvarnames[(rois+1):vars], "~", lvarnames[1:rois])
    ## creates list of AR paths so that later code doesn't kick them out
    fixed_paths <- paste0(lvarnames[(rois+1):vars], "~", lvarnames[1:rois]) 
  } else {
    line4 <- paste0(lvarnames[1:rois], "~0*", lvarnames[(rois+1):vars])
    fixed_paths <- NULL
  }
  
  
  syntax <- c(line1, line2, line3, line4)
  
  if (!ar & is.null(paths)) {
    covzero <- NULL
    for (i in (rois+2):vars) {
      for (j in (rois+1):(i-1)){
        covzero <- c(covzero, paste0(lvarnames[i],"~~0*", lvarnames[j]))
      }
    }
    syntax <- c(syntax, covzero)
  }
  
  
  if (!is.null(paths)) syntax <- c(syntax, paths)
  
  ## create list of paths that make sense to gimme to open
  candidate_paths   <- apply(expand.grid(lvarnames[(rois+1):vars], 
                                         lvarnames[1:vars]), 1, paste, collapse = "~")
  ###So what this does now is uses lvarnames[(rois+1):vars] to list VAR1..VAR5 (no lags) PLUS the two multiplied,
  ###Then uses lvarnames[1:vars] to list ALL the variables and regresses them each on each from the first list.
  ###Is this okay? It's not problematic that the mult vars are listed twice, because the next part takes care
  ###of removing where they're DVs, right?
 
   # if path specified by a user is fixed to a certain value, 
  # remove it from consideration when looking at MIs
  candidate_paths <- candidate_paths[!candidate_paths %in% remove]
  
  ## create list of impossible exogenous paths
  exog_paths<-NULL
  if(!is.null(lexogenous)){
  exog_paths <- apply(expand.grid(lexogenous[1:length(lexogenous)],
                  lvarnames[1:length(lvarnames)]), 1, paste, collapse = "~")
  }
  
   # remove impossible exogenous paths from candidate paths, fixed paths, and syntax
  if(!is.null(exog_paths)){
  candidate_paths <- candidate_paths[!candidate_paths %in% exog_paths]
  }
  
  if(!is.null(exog_paths)){
    fixed_paths <- fixed_paths[!fixed_paths %in% exog_paths]
  }
  
  if(!is.null(exog_paths)){
    syntax <- syntax[!syntax %in% exog_paths]
  }
  
  if(!is.null(exog_paths)){
    syntax <- c(syntax,paste0(lexogenous[1:n_exog_total], "~0*", lvarnames[1]))
  }
  # if user specifies paths, add them to the list of fixed paths
  if (!is.null(paths)) fixed_paths <- c(fixed_paths, paths) 
  
  dat <- list("ar" = ar, 
              "out"= out,
              "plot" = plot,
              "subgroup" = subgroup,
              "agg" = agg,
              "n_subj" = subjects,
              "n_rois" = n_endog, ### Changed to n_endog because n_rois is used in functions to refer to only endogenous vars!
              "n_exog" = n_exog,
              "n_bilinear" = n_bilinear,
              "n_exog_total" = n_exog_total,
              "n_vars_total" = n_vars_total,
              "n_contemporaneous" = n_contemporaneous,
              "varnames" = varnames,
              "lvarnames" = lvarnames,
              "cutoffind" = cutoffind,
              "subgroup_dir" = subgroup_dir,
              "ind_dir"   = individual,
              "syntax"     = syntax,
              "candidate_paths" = candidate_paths,
              "fixed_paths"  = fixed_paths,
              "group_cutoff" = groupcutoff,
              "sub_cutoff" = subcutoff,
              "ts_list" = ts_list,
              "file_order" = file_order,
              "chisq_cutoff_mi_epc" = qchisq(1-(.05/((vars*(vars-1)/2)*subjects)), 1))
  return(dat)
}
  
