## 4 Different runs are tested:
#   (1) Typical
#   (2) HRF
#   (3) Subgrouping
#   (4) Latent Variables

## Each run has 2 tests:
#   (1) The "expect_equal" line tests that the beta estimates for the current 
# version are equal to the path estimates from the source of truth version, with
# a tolerance of 1e-5 (i.e. differences smaller than 1e-5 are ignored)
#   (2) The "expect_identical" line tests that the paths recovered by the
# current version are the same as the paths recovered by the source of truth
# version.
# 
# ## The "source of truth" version is the version of GIMME on CRAN as of Oct 2023.
# 
# 
# test_that("Run 1 gives expected results", {
#   run1_sot <- readRDS("./rds/run1_path_matrix.rds")
#   run1_paths_sot <-readRDS("./rds/run1_path_counts.rds")
#   run1 <- gimme(data = gimme::ts)
#   expect_equal(run1[["path_est_mats"]], run1_sot, tolerance = 1e-5)
#   expect_identical(run1[["path_counts"]], run1_paths_sot)
# })


# 
# test_that("Run 2 gives expected results", {
#   run2_sot <- readRDS("rds/run2_path_matrix.rds")
#   run2_paths_sot <-readRDS("rds/run2_path_counts.rds")
#   run2 <- gimme(data = gimme::HRFsim,
#                 ar = TRUE,
#                 exogenous = "V5",
#                 conv_vars = "V5",
#                 conv_length = 16,
#                 conv_interval = 1,
#                 mult_vars = "V4*V5",
#                 mean_center_mult = TRUE
#   )
#   expect_equal(run2[["path_est_mats"]], run2_sot, tolerance = 1e-5)
#   expect_identical(run2[["path_counts"]], run2_paths_sot)
# })
# 
# 
# # 
# test_that("Run 3 gives expected results", {
#   run3_sot <- readRDS("./rds/run3_path_matrix.rds")
#   run3_paths_sot <-readRDS("./rds/run3_path_counts.rds")
#   run3 <- gimme(data = gimme::simData,
#                 subgroup = TRUE)
#   expect_equal(run3[["path_est_mats"]], run3_sot, tolerance = 1e-5)
#   expect_identical(run3[["path_counts"]], run3_paths_sot)
# })
# 
# 
# test_that("Run 4 gives expected results", {
#   run4_sot <- readRDS("rds/run4_path_matrix.rds")
#   run4_paths_sot <-readRDS("rds/run4_path_counts.rds")
#   lv_model_all <- "
#   L1 =~ V1 + V2 + V3 
#   L2 =~ V4 + V5 + V6
#   L3 =~ V7 + V8 + V9
#   "
#   run4 <- gimme(data = gimme::simDataLV,
#                 subgroup = TRUE,
#                 lv_model = lv_model_all,
#                 lv_estimator = "miiv",
#                 lv_score = "regression",
#                 lv_final_estimator = "miiv")
#   expect_equal(run4[["path_est_mats"]], run4_sot, tolerance = 1e-3)
#   expect_identical(run4[["path_counts"]], run4_paths_sot)
# })

test_that("highest.mi standard stop_crit stops once fit is adequate", {
  mi_list <- list(data.frame(lhs = "y", op = "~", rhs = "x", mi = 5, epc = 0.2))
  indices <- c(chisq = 1, df = 1, pvalue = .1, rmsea = .01, srmr = .01, nnfi = .99, cfi = .99)
  
  res <- gimme:::highest.mi(
    mi_list = mi_list,
    indices = indices,
    elig_paths = "y~x",
    prop_cutoff = NULL,
    n_subj = 1,
    chisq_cutoff = 3.84,
    stop_crit = "standard",
    allow.mult = FALSE,
    ms_tol = 1e-5,
    hybrid = FALSE,
    dir_prop_cutoff = 0
  )
  
  expect_true(isTRUE(res$goodfit))
  expect_true(is.na(res$add_param))
})

test_that("highest.mi model fit stop_crit can add a nonsignificant path", {
  mi_list <- list(data.frame(lhs = "y", op = "~", rhs = "x", mi = 2, epc = 0.2))
  indices <- c(chisq = 10, df = 1, pvalue = .001, rmsea = .20, srmr = .20, nnfi = .50, cfi = .50)
  
  res <- gimme:::highest.mi(
    mi_list = mi_list,
    indices = indices,
    elig_paths = "y~x",
    prop_cutoff = NULL,
    n_subj = 1,
    chisq_cutoff = 3.84,
    stop_crit = "model fit",
    allow.mult = FALSE,
    ms_tol = 1e-5,
    hybrid = FALSE,
    dir_prop_cutoff = 0
  )
  
  expect_false(isTRUE(res$goodfit))
  expect_identical(res$add_param, "y~x")
})

test_that("highest.mi significance stop_crit keeps adding significant paths", {
  mi_list <- list(data.frame(lhs = "y", op = "~", rhs = "x", mi = 5, epc = 0.2))
  indices <- c(chisq = 1, df = 1, pvalue = .1, rmsea = .01, srmr = .01, nnfi = .99, cfi = .99)
  
  res <- gimme:::highest.mi(
    mi_list = mi_list,
    indices = indices,
    elig_paths = "y~x",
    prop_cutoff = NULL,
    n_subj = 1,
    chisq_cutoff = 3.84,
    stop_crit = "significance",
    allow.mult = FALSE,
    ms_tol = 1e-5,
    hybrid = FALSE,
    dir_prop_cutoff = 0
  )
  
  expect_true(isTRUE(res$goodfit))
  expect_identical(res$add_param, "y~x")
})

test_that("dataAUG processing records augmented variables and restores common names", {
  data_aug <- list(
    person1 = data.frame(x_AUG = 1:4, y = 5:8),
    person2 = data.frame(y_AUG = 1:4, x = 5:8)
  )

  res <- gimme:::setupAugmentedData(data_aug)

  expect_identical(res$augmented_vars$person1, "x")
  expect_identical(res$augmented_vars$person2, "y")
  expect_identical(colnames(res$data$person1), c("x", "y"))
  expect_identical(colnames(res$data$person2), c("x", "y"))
})

test_that("setup uses data when dataAUG is TRUE", {
  data_aug <- list(
    person1 = data.frame(x_AUG = 1:4, y = 5:8),
    person2 = data.frame(y_AUG = 1:4, x = 5:8)
  )

  ctrl <- list(dataAUG = TRUE, sep = NULL, header = NULL)
  raw <- gimme:::setupDataLists(data_aug, ctrlOpts = ctrl)
  res <- gimme:::setupAugmentedData(raw)

  expect_identical(res$augmented_vars$person1, "x")
  expect_identical(res$augmented_vars$person2, "y")
  expect_identical(colnames(res$data$person2), c("x", "y"))
  expect_identical(formals(gimme)$dataAUG, FALSE)
  expect_identical(formals(indSEM)$dataAUG, FALSE)
})

test_that("augmented variables only retain their autoregressive path", {
  lhs <- c("x", "y", "x", "x", "x")
  op <- c("~", "~", "~", "~~", "~~")
  rhs <- c("xlag", "xlag", "y", "y", "x")

  prohibited <- gimme:::isAugmentedPath(lhs, op, rhs, "x")

  expect_identical(prohibited, c(FALSE, TRUE, TRUE, TRUE, FALSE))
})

test_that("augmented paths are zeroed in syntax and MI evidence", {
  syntax <- c("x~xlag", "y~xlag", "x~y", "x~~y", "x~~x", "x~1")
  constrained <- gimme:::constrainAugmentedSyntax(syntax, "x")

  expect_identical(
    constrained,
    c("x~xlag", "y~0*xlag", "x~0*y", "x~~0*y", "x~~x", "x~1")
  )

  mis <- data.frame(
    lhs = c("x", "y", "x"),
    op = "~",
    rhs = c("xlag", "xlag", "y"),
    mi = c(3, 8, 9),
    epc = c(.1, .4, .5),
    pvalue = c(.08, .001, .001)
  )
  masked <- gimme:::maskAugmentedMIs(mis, "x")

  expect_identical(masked$mi, c(3, 0, 0))
  expect_identical(masked$epc, c(.1, 0, 0))
  expect_identical(masked$pvalue, c(.08, 1, 1))
})

test_that("highest.mi does not select a zeroed augmented path", {
  mi_list <- list(gimme:::maskAugmentedMIs(data.frame(
    lhs = "y", op = "~", rhs = "x", mi = 0, epc = 0
  ), "x"))
  indices <- c(
    chisq = 10, df = 1, pvalue = .001, rmsea = .20,
    srmr = .20, nnfi = .50, cfi = .50
  )

  res <- gimme:::highest.mi(
    mi_list = mi_list,
    indices = indices,
    elig_paths = "y~x",
    prop_cutoff = NULL,
    n_subj = 1,
    chisq_cutoff = 3.84,
    stop_crit = "model fit",
    allow.mult = FALSE,
    ms_tol = 1e-5,
    hybrid = FALSE,
    dir_prop_cutoff = 0
  )

  expect_true(is.na(res$add_param))
})
