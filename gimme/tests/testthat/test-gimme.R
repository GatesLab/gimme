context("gimme")
skip_on_cran()

test_that("gimme output is the same", {
  
  dir <- tempdir()
  
  sol <- gimme(data="C:/Users/hpike/Desktop/practicedata/input", 
               out=dir, 
               sep="", 
               header=FALSE, 
               ar=TRUE, 
               plot=FALSE, 
               subgroup=TRUE, 
               paths=NULL)

  ##get from data output after making code changes
  modularity <- read.csv(file.path(dir, "summaryFit.csv", fsep=.Platform$file.sep))
  modularity <- modularity[!colnames(modularity) %in% "subgroup"]
    
  ##get when know code is working
  modularity2 <- structure(list(file = structure(1:5, .Label = c("ts1", "ts2", "ts3", "ts4", "ts5"), 
                class = "factor"), chisq = c(5.9521, 1.2039, 1.2499, 10.038, 7.1295), df = c(6L, 6L, 5L, 7L, 6L), 
                pval = c(0.4286, 0.9767, 0.94, 0.1864, 0.309), rmsea = c(0, 0, 0, 0.0941, 0.062), srmr = 
                c(0.0493, 0.0136, 0.0158, 0.0349, 0.0436), nnfi = c(1, 1, 1, 0.9817, 0.993), cfi = 
                c(1.0007, 1.0784, 1.0755, 0.9607, 0.9824), status = structure(c(1L, 1L, 1L, 1L, 1L), 
                .Label = "converged normally", class = "factor"), subgroup = 1:5, modularity = c(0L, NA, NA, NA, NA)), 
                .Names = c("file", "chisq", "df", "pval", "rmsea", "srmr", "nnfi", "cfi", "status", "subgroup", 
                "modularity"), class = "data.frame", row.names = c(NA, -5L))
  modularity2 <- modularity2[!colnames(modularity2) %in% "subgroup"]
  
  expect_that(modularity, equals(modularity2))
})





























