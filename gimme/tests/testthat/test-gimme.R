context("gimme")
skip_on_cran()

test_that("gimme output is the same", {
  
  sol <- gimme(data="C:/Users/hpike/Desktop/practicedata/input", 
               out="C:/Users/hpike/Desktop/practicedataout/output", 
               sep="", 
               header=FALSE, 
               ar=TRUE, 
               plot=FALSE, 
               subgroup=TRUE, 
               paths=NULL)
  
  out="C:/Users/hpike/Desktop/practicedataout/output"

  ##get from data output after making code changes
  modularity <- read.csv(file.path(out, "summaryFit.csv", fsep=.Platform$file.sep))
  modularity <- modularity[!colnames(modularity) %in% "subgroup"]
    
  ##get when know code is working
  modularity2 <- read.csv("C:/Users/hpike/Desktop/gimmeTroubleshooting/practiceoutput.csv")
  modularity2 <- modularity2[!colnames(modularity2) %in% "subgroup"]
  
  expect_that(modularity, equals(modularity2))
})





























