pkgname <- "gimme"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('gimme')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("aggSEM")
### * aggSEM

flush(stderr()); flush(stdout())

### Name: aggSEM
### Title: Group-level structural equation model search.
### Aliases: aggSEM

### ** Examples

data(ts1,ts2,ts3,ts4,ts5)
input.path  <- file.path(tempdir(),"input")
output.path <- file.path(tempdir(),"output")
dir.create(input.path)
dir.create(output.path)
write.table(ts1,file.path(input.path,"ts1.txt"),col.names=FALSE,row.names=FALSE)
write.table(ts2,file.path(input.path,"ts2.txt"),col.names=FALSE,row.names=FALSE)
write.table(ts3,file.path(input.path,"ts3.txt"),col.names=FALSE,row.names=FALSE)
write.table(ts4,file.path(input.path,"ts4.txt"),col.names=FALSE,row.names=FALSE)
write.table(ts5,file.path(input.path,"ts5.txt"),col.names=FALSE,row.names=FALSE)
aggSEM(data   = input.path,
       sep    = "",
       header = FALSE,
       out    = output.path,
       ar     = TRUE,
       plot   = TRUE,
       paths  = NULL)



cleanEx()
nameEx("gimme")
### * gimme

flush(stderr()); flush(stdout())

### Name: gimme
### Title: Group iterative multiple model estimation.
### Aliases: gimme
### Keywords: gimme

### ** Examples

## Not run: 
##D paths <- 'V2 ~ V1
##D           V3 ~ V4lag'
##D 
##D gimme(data     = "C:/data100",
##D       sep      = ",",
##D       header   = FALSE,
##D       out      = "C:/data100_gimme_out",
##D       ar       = TRUE,
##D       plot     = TRUE,
##D       paths    = paths,
##D       subgroup = FALSE)
##D  
## End(Not run)



cleanEx()
nameEx("indSEM")
### * indSEM

flush(stderr()); flush(stdout())

### Name: indSEM
### Title: Individual-level structural equation model search.
### Aliases: indSEM
### Keywords: indSEM

### ** Examples

## Not run: 
##D indSEM(data   = "C:/data100",
##D        sep    = ",",
##D        header = FALSE,
##D        out    = "C:/data100_indSEM_out",
##D        ar     = TRUE,
##D        plot   = TRUE,
##D        paths  = NULL)
##D  
## End(Not run)



### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
