library(gimme)

#datalocation <- '/Users/gateskm/Dropbox/gateslab/GG_Workshop Materials/gimme/Data example/sub_2_prop_0.5_n_25_rep_1' #Mac
datalocation <- 'C:\\Users\\gateskm\\Dropbox\\gateslab\\GG_Workshop Materials\\gimme\\Data example\\t_120_n_25_v_5' #windows


library(tools) # necessary to get file names without extension
list_files <- as.matrix(file_path_sans_ext(list.files(datalocation)))

##### RUNNING GIMME WITH DEFAULTS #######
# AR = TRUE; SUBGROUP = FALSE; PLOT = TRUE; group cutoff = 0.75
output <- gimmeSEM(data        = datalocation,  #notice that we are saving the output to an object; optional
         #out         = "/Users/gateslab/Desktop/t_120_n_25_v_5 2 OUT", # optional output folder
         sep         = ",", # options: "/t" for tab-delimited, "" for space-delimited, "," for comma
         header      = FALSE # no default; user must specify
         )

##### MAKING A LIST #####
datalist <- list()
for (p in 1:length(list_files))
  datalist[[p]] <- read.table(paste0(datalocation, "/", list_files[p], ".txt"), header = FALSE, sep = ",")

names(datalist) <- list_files # you must name the slices in the list, or else gimme won't work.

# running gimme with list and defaults #
output <- gimmeSEM(data        = datalist,  #notice that we are saving the output to an object; optional
                   out         = "/Users/gateslab/Desktop/t_120_n_25_v_5 2 OUT" # optional output folder
                   )

##### SETTING AR TO FALSE #####
# Note:  it is recommended to keep AR = TRUE
output <- gimmeSEM(data        = datalist,  #notice that we are saving the output to an object; optional
                   out         = "/Users/gateslab/Desktop/t_120_n_25_v_5 2 OUT", # optional output folder
                   AR          = FALSE
                   )

##### FINDING DATA-DRIVEN SUBGROUPS #####
# The "subgroups" argument specifies if you would like unsupervised classification of individuals 
# based on their model features.
output <- gimmeSEM(data        = datalist,  #notice that we are saving the output to an object; optional
                 #  out         = "/Users/gateslab/Desktop/t_120_n_25_v_5 2 OUT", # optional output folder
                   subgroup   = TRUE
                   )

##### SPECIFYING CONFIRMATORY SUBGROUPS #####

subgroup_ids <- matrix(, length(datalist), 2) # make a matrix that has 2 columns
subgroup_ids[,1] <- list_files
subgroup_ids[1:12,2] <- 1
subgroup_ids[13:25,2] <- 2

View(subgroup_ids)

output_subconf <- gimmeSEM(data        = datalist,
         confirm_subgroup = subgroup_ids
        )

##### SPECIFYING EXOGENEOUS PATHS #####
# These are paths that can predict but cannot be predicted. 
exog_vars <- c("V1", "V5")

output_exog <- gimmeSEM(data        = datalist,
         exogenous   = exog_vars
         )


##### SPECIFYING CONFIRMATORY PATHS ####
# These are paths that are present at the beginning of the search. They will be estimated for everyone in the sample. 
# lavaan-style syntax is used. 
library(lavaan)
conf_paths <- 'V1 ~ V2 
               V3 ~ V4lag'

output_confpaths <- gimmeSEM(data        = datalist,
                    paths       = conf_paths
                    )
