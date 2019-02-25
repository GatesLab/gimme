## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, eval = F)

## ----1-------------------------------------------------------------------
#  install.packages("tools") #to use the function to get file names sans extension
#  #get filenames in the folder without extension
#  filename <- file_path_sans_ext(list.files(path = "t_120_n_25_v_5", full.names = FALSE))
#  subgroup <- c(1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2) #create the subgroup
#  confirm_dataframe <- data.frame(filename, subgroup)              #create the dataframe

## ----2-------------------------------------------------------------------
#  gimme(                  # can use "gimme" or "gimmeSEM"
#      data = 't_120_n_25_v_5',          # source directory where your data are
#      out = 'SampleOutput',            # output directory where you'd like your output to go
#      sep = ",",           # how data are separated. "" for space; "," for comma, "/t" for tab-delimited
#      header = FALSE,          # TRUE or FALSE, is there a header
#      ar = TRUE,          # TRUE (default) or FALSE, start with autoregressive paths open
#      plot = TRUE,        # TRUE (default) or FALSE, generate plots
#      subgroup = TRUE,    # Must be TRUE to perform confirmatory subgrouping
#      confirm_subgroup = confirm_dataframe, # confirm_dataframe is the dataframe constructed previously
#      paths = NULL,       # option to list paths that will be group-level (semi-confirmatory)
#      groupcutoff = .75,  # the proportion that is considered the majority at the group level
#      subcutoff = .75      # the proportion that is considered the majority at the subgroup level
#      )

