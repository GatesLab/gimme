#' Prepare data containing augmented variables
#' @param data A list of individual data sets with augmented columns ending
#' in "_AUG".
#' @keywords internal
setupAugmentedData <- function(data) {
  if (!is.list(data)) {
    stop("gimme ERROR: augmented data must be a list with one data frame or matrix per individual.")
  }
  if (length(data) == 0) {
    stop("gimme ERROR: augmented data must contain at least one individual data set.")
  }

  if (is.null(names(data))) {
    names(data) <- paste0("subj", seq_along(data))
  }

  augmented_vars <- lapply(data, function(x) {
    varnames <- colnames(x)
    if (is.null(varnames)) {
      stop("gimme ERROR: dataAUG = TRUE requires column names so augmented variables can be identified.")
    }
    sub("_AUG$", "", varnames[grepl("_AUG$", varnames)])
  })

  data <- lapply(data, function(x) {
    colnames(x) <- sub("_AUG$", "", colnames(x))
    if (anyDuplicated(colnames(x))) {
      stop("gimme ERROR: removing '_AUG' creates duplicate variable names within an individual.")
    }
    x
  })

  reference_names <- colnames(data[[1]])
  data <- lapply(data, function(x) {
    if (!setequal(colnames(x), reference_names)) {
      stop("gimme ERROR: all data entries must contain the same variables after removing '_AUG'.")
    }
    x[, reference_names, drop = FALSE]
  })

  names(augmented_vars) <- names(data)
  list(data = data, augmented_vars = augmented_vars)
}

#' Identify relations prohibited for an augmented variable
#' @param lhs Left-hand-side variable names.
#' @param op Lavaan operators.
#' @param rhs Right-hand-side variable names.
#' @param augmented_vars Base names of augmented variables for one individual.
#' @keywords internal
isAugmentedPath <- function(lhs, op, rhs, augmented_vars) {
  if (length(augmented_vars) == 0) {
    return(rep(FALSE, length(lhs)))
  }

  augmented_names <- c(augmented_vars, paste0(augmented_vars, "lag"))
  involves_augmented <- lhs %in% augmented_names | rhs %in% augmented_names
  allowed_ar <- op == "~" & lhs %in% augmented_vars & rhs == paste0(lhs, "lag")
  allowed_variance <- op == "~~" & lhs == rhs
  allowed_intercept <- op == "~" & rhs == "1"

  involves_augmented & !allowed_ar & !allowed_variance & !allowed_intercept
}

#' Fix prohibited augmented-variable paths to zero for one individual
#' @param syntax Lavaan model syntax.
#' @param augmented_vars Base names of augmented variables for one individual.
#' @keywords internal
constrainAugmentedSyntax <- function(syntax, augmented_vars) {
  if (length(augmented_vars) == 0 || length(syntax) == 0) {
    return(syntax)
  }

  vapply(syntax, function(line) {
    compact <- gsub("\\s+", "", line)
    match <- regexec("^([^~=]+)(~~|~)(?:[^*]+\\*)?([^~]+)$", compact)
    parts <- regmatches(compact, match)[[1]]

    if (length(parts) == 0) {
      return(line)
    }

    lhs <- parts[2]
    op <- parts[3]
    rhs <- parts[4]
    if (isAugmentedPath(lhs, op, rhs, augmented_vars)) {
      paste0(lhs, op, "0*", rhs)
    } else {
      line
    }
  }, character(1), USE.NAMES = FALSE)
}

#' Remove MI/EPC evidence for prohibited augmented-variable paths
#' @param mis A modification-index data frame.
#' @param augmented_vars Base names of augmented variables for one individual.
#' @keywords internal
maskAugmentedMIs <- function(mis, augmented_vars) {
  if (length(augmented_vars) == 0 || length(mis) == 1 && is.na(mis)) {
    return(mis)
  }

  prohibited <- isAugmentedPath(mis$lhs, mis$op, mis$rhs, augmented_vars)
  mis$augmented_prohibited <- prohibited
  zero_columns <- intersect(c("mi", "epc", "sepc.lv", "sepc.all", "sepc.nox"), names(mis))
  mis[prohibited, zero_columns] <- 0
  p_columns <- intersect(c("pvalue", "p"), names(mis))
  mis[prohibited, p_columns] <- 1
  mis
}

#' Remove significance evidence for prohibited augmented-variable paths
#' @param zs A standardized-solution data frame.
#' @param augmented_vars Base names of augmented variables for one individual.
#' @keywords internal
maskAugmentedZs <- function(zs, augmented_vars) {
  if (length(augmented_vars) == 0 || length(zs) == 1 && is.na(zs)) {
    return(zs)
  }

  prohibited <- isAugmentedPath(zs$lhs, zs$op, zs$rhs, augmented_vars)
  zero_columns <- intersect(c("est.std", "se", "z"), names(zs))
  zs[prohibited, zero_columns] <- 0
  p_columns <- intersect(c("pvalue", "p"), names(zs))
  zs[prohibited, p_columns] <- 1
  zs
}
