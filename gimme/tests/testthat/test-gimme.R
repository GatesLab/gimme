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

## The "source of truth" version is the version of GIMME on CRAN as of Oct 2023.


test_that("Run 1 gives expected results", {
  run1_sot <- readRDS("rds/run1_path_matrix.rds")
  run1_paths_sot <-readRDS("rds/run1_path_counts.rds")
  run1 <- gimme(data = gimme::ts)
  expect_equal(run1[["path_est_mats"]], run1_sot, tolerance = 1e-5)
  expect_identical(run1[["path_counts"]], run1_paths_sot)
})



test_that("Run 2 gives expected results", {
  run2_sot <- readRDS("rds/run2_path_matrix.rds")
  run2_paths_sot <-readRDS("rds/run2_path_counts.rds")
  run2 <- gimme(data = gimme::HRFsim,
                ar = TRUE,
                exogenous = "V5",
                conv_vars = "V5",
                conv_length = 16,
                conv_interval = 1,
                mult_vars = "V4*V5",
                mean_center_mult = TRUE
  )
  expect_equal(run2[["path_est_mats"]], run2_sot, tolerance = 1e-5)
  expect_identical(run2[["path_counts"]], run2_paths_sot)
})



test_that("Run 3 gives expected results", {
  run3_sot <- readRDS("rds/run3_path_matrix.rds")
  run3_paths_sot <-readRDS("rds/run3_path_counts.rds")
  run3 <- gimme(data = gimme::simData,
                subgroup = TRUE)
  expect_equal(run3[["path_est_mats"]], run3_sot, tolerance = 1e-5)
  expect_identical(run3[["path_counts"]], run3_paths_sot)
})


test_that("Run 4 gives expected results", {
  run4_sot <- readRDS("rds/run4_path_matrix.rds")
  run4_paths_sot <-readRDS("rds/run4_path_counts.rds")
  lv_model_all <- "
  L1 =~ V1 + V2 + V3 
  L2 =~ V4 + V5 + V6
  L3 =~ V7 + V8 + V9
  "
  run4 <- gimme(data = gimme::simDataLV,
                subgroup = TRUE,
                lv_model = lv_model_all,
                lv_estimator = "miiv",
                lv_score = "regression",
                lv_final_estimator = "miiv")
  expect_equal(run4[["path_est_mats"]], run4_sot, tolerance = 1e-3)
  expect_identical(run4[["path_counts"]], run4_paths_sot)
})

