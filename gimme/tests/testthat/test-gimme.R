
# Initial Test Passed
test_that("Run 1 gives expected results", {
  run1_sot <- readRDS("rds/run1_path_matrix.rds")
  run1_paths_sot <-readRDS("rds/run1_path_counts.rds")
  run1 <- gimme(data = gimme::ts)
  expect_equal(run1[["path_est_mats"]], run1_sot, tolerance = 1e-5)
  expect_identical(run1[["path_counts"]], run1_paths_sot)
})


# Error - test not run
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


# Initial Test failed - differences smaller than 2 decimal places
# Change code to allow this unless one of the elements is zero
test_that("Run 3 gives expected results", {
  run3_sot <- readRDS("rds/run3_path_matrix.rds")
  run3_paths_sot <-readRDS("rds/run3_path_counts.rds")
  run3 <- gimme(data = gimme::simData,
                subgroup = TRUE)
  expect_equal(run3[["path_est_mats"]], run3_sot, tolerance = 1e-5)
  expect_identical(run3[["path_counts"]], run3_paths_sot)
})

# Initial Test Passed
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
  expect_equal(run4[["path_est_mats"]], run4_sot, tolerance = 1e-5)
  expect_identical(run4[["path_counts"]], run4_paths_sot)
})

