## ---- eval = FALSE-------------------------------------------------------
#  gimme(                  # can use "gimme" or "gimmeSEM"
#        data = '',          # list object or source directory where your data are 
#        out = '',           # output directory where you'd like your output to go
#        sep = "",           # how data are separated. "" for space; "," for comma, "/t" for tab-delimited
#        header = ,          # TRUE or FALSE, is there a header
#        ar = TRUE,          # TRUE (default) or FALSE, start with autoregressive paths open
#        plot = TRUE,        # TRUE (default) or FALSE, generate plots
#        subgroup = FALSE,   # TRUE or FALSE (default), cluster individuals based on similarities in effects
#        paths = NULL,       # option to list paths that will be group-level (semi-confirmatory)
#        groupcutoff = .75,  # the proportion that is considered the majority at the group level
#        subcutoff = .75,     # the proportion that is considered the majority at the subgroup level
#        VAR       = FALSE,  # TRUE or FALSE (default), option to use VAR model instead of uSEM
#        hybrid    = FALSE   # TRUE or FALSE (default), option to use hybrid-VAR model instead of uSEM
#)  
