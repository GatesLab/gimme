#' @name gimmeInteractive
#' @aliases gimmeInteractive gimmeGUI
#' @title Graphical User Interface designed for use with the gimme package
#' @description This function launches a graphical user interface which calls the main functions
#' for the gimme package. gimmeInteractive() requires the package gWidgetsRGtk2, which requires the GTK+ library. 
#' Please install this package before using gimmeInteractive().
#' @usage gimmeInteractive()
#' @details Returns a simple text file, guiSyntax.txt, containing a copy of the syntax submitted by the GUI. Placed in output folder specified by the user.
#' @author Hallie Pike, Stephanie Lane
#' @export gimmeInteractive gimmeGUI
gimmeInteractive <- gimmeGUI <- function() {
# added code to require the package RGtk2 and GTK+
 if (!requireNamespace("gWidgetsRGtk2", quietly = TRUE)) {
    stop("gWidgetsRGtk2 needed for this function to work. Please install it, along with the GTK+ library, before running gimmeInteractive().",
         call. = FALSE)
 }
  
  ssrun = NULL
  gui.env          <- new.env()
  gui.env$Header   <- TRUE
  gui.env$Plot     <- TRUE
  gui.env$Ar       <- TRUE
  gui.env$Subgroup <- FALSE
  gui.env$Paths    <- NULL
  gui.env$groupcutoff <- .75
  gui.env$subcutoff   <- .5

  introwin <- gwindow(title="gimmeInteractive", visible=TRUE, width = 600, height = 200)
  gwin     <- gwindow(title="gimmeSEM", visible=FALSE)
  iwin     <- gwindow(title="indSEM", visible=FALSE)
  awin     <- gwindow(title="aggSEM", visible=FALSE)

  introgroup <- ggroup(horizontal=FALSE, container=introwin)
  agroup     <- ggroup()
  bgroup     <- ggroup(container=agroup)
  lyt        <- glayout(container=agroup)

  gimmebutton  <- gbutton("                                                gimmeSEM \n (Default choice for robust model selection with heterogeneous data)",
                         container=introgroup)
  indSEMbutton <- gbutton("                                                                indSEM \n (Conduct individual-level model search using no shared information from sample)",
                          container=introgroup)
  aggbutton    <- gbutton("                                                                         aggSEM \n (Concatenate data across all individuals and arrive at one model; assumes strict homogeneity)",
                       container=introgroup)

  size(gimmebutton)  <- c(15,60)
  size(indSEMbutton) <- c(15,60)
  size(aggbutton)    <- c(15,60)

  code1  <- lyt[21,1, anchor=c(1,0),  expand=TRUE] <- glabel("data = ", container=lyt)
  code2  <- lyt[21,2, anchor=c(-1,0), expand=TRUE] <- glabel(NULL, container=lyt)
  code7  <- lyt[22,1, anchor=c(1,0),  expand=TRUE] <- glabel("out = ", container=lyt)
  code8  <- lyt[22,2, anchor=c(-1,0), expand=TRUE] <- glabel(NULL, container=lyt)
  code3  <- lyt[23,1, anchor=c(1,0),  expand=TRUE] <- glabel("sep = ", container=lyt)
  code4  <- lyt[23,2, anchor=c(-1,0), expand=TRUE] <- glabel(NULL, container=lyt)
  code5  <- lyt[24,1, anchor=c(1,0),  expand=TRUE] <- glabel("header = ", container=lyt)
  code6  <- lyt[24,2, anchor=c(-1,0), expand=TRUE] <- glabel(TRUE, container=lyt)
  code9  <- lyt[25,1, anchor=c(1,0),  expand=TRUE] <- glabel("ar = ", container=lyt)
  code10 <- lyt[25,2, anchor=c(-1,0), expand=TRUE] <- glabel(TRUE, container=lyt)
  code11 <- lyt[26,1, anchor=c(1,0),  expand=TRUE] <- glabel("plot = ", container=lyt)
  code12 <- lyt[26,2, anchor=c(-1,0), expand=TRUE] <- glabel(TRUE, container=lyt)
  code13 <- lyt[27,1, anchor=c(1,0),  expand=TRUE] <- glabel("subgroup = ", container=lyt, visible=FALSE)
  code14 <- lyt[27,2, anchor=c(-1,0), expand=TRUE] <- glabel(FALSE, container=lyt, visible=FALSE)
  code15 <- lyt[28,1, anchor=c(1,0),  expand=TRUE] <- glabel("groupcutoff = ", container=lyt, visible=FALSE)
  code16 <- lyt[28,2, anchor=c(-1,0), expand=TRUE] <- glabel(.75, container=lyt, visible=FALSE)
  code17 <- lyt[29,1, anchor=c(1,0),  expand=TRUE] <- glabel("subcutoff = ", container=lyt, visible = F)
  code18 <- lyt[29,2, anchor=c(-1,0), expand=TRUE] <- glabel(.5, container=lyt, visible = F)
  
  back <- gimage(stock.id="go-back",
                 dirname="stock",
                 container=bgroup,
                 size="large_toolbar",
                 anchor=c(0,1),
                 handler=function(h,...) {
                   #visible(introwin) <- TRUE
                   if (visible(gwin)==TRUE) delete(gwin, agroup)
                   if (visible(awin)==TRUE) delete(awin, agroup)
                   if (visible(iwin)==TRUE) delete(iwin, agroup)
                   visible(separ6) <- TRUE
                   visible(subgroupLabel) <- TRUE
                   visible(subgroup) <- TRUE
                   if (visible(iwin) == TRUE | visible(awin) == TRUE) {
                     visible(code13) <- FALSE
                     visible(code14) <- FALSE
                   }
                   visible(gwin) <- FALSE
                   visible(iwin) <- FALSE
                   visible(awin) <- FALSE
                 })

  datalabel <- lyt[2,1] <- glabel("Where are the data files located?",
                                  container=lyt)
  data <- lyt[2,2] <- gfilebrowse(text="Select a folder...",
                                  type="selectdir",
                                  quote=TRUE,
                                  container=lyt,
                                  action=code2,
                                  handler=function(h,...) {
                                    oldVal <- svalue(h$obj)
                                    str <- gsub("'", replacement="", oldVal, fixed=TRUE)
                                    svalue(h$action) <- dQuote(str)
                                    assign("Data", str, envir=gui.env)
                                  })
  separ1 <- lyt[3,1:2] <- gseparator(container=lyt)
  
  outLabel <- lyt[4,1] <- glabel("Where would you like results to be saved?", container=lyt)
  out <- lyt[4,2] <- gfilebrowse(text="Select a folder...",
                                 type="selectdir",
                                 quote=TRUE,
                                 container=lyt,
                                 action=code8,
                                 handler=function(h,...) {
                                   oldVal <- svalue(h$obj)
                                   str <- gsub("'", replacement="", oldVal, fixed=TRUE)
                                   svalue(h$action) <- dQuote(str)
                                   assign("Out", str, envir=gui.env)
                                 })
  
  separ3 <- lyt[5,1:2] <- gseparator(container=lyt)

  sepLabel <- lyt[6,1] <- glabel("What is the format of the data files?",
                                 container=lyt)
  sep <- lyt[6,2] <- gcombobox(c("", "Space Delimited", "Tab Delimited", "Comma Delimited"),
                               container=lyt,
                               action=code4,
                               handler=function(h,...) {
                                 oldVal <- svalue(h$obj)
                                 if (oldVal=="Space Delimited") str <- ""
                                 else if (oldVal=="Tab Delimited") str <- "\t"
                                 else if (oldVal=="Comma Delimited") str <- ","
                                 svalue(h$action) <- dQuote(str)
                                 if (oldVal == "Tab Delimited"){
                                   svalue(h$action) <- paste0("'", "\\", "t", "'")
                                 }
                                 assign("Sep", str, envir=gui.env)
                               })
  separ2 <- lyt[7,1:2] <- gseparator(container=lyt)

  headerLabel <- lyt[8,1] <- glabel("Do the data files have a header?",
                                    container=lyt)
  header <- lyt[8,2] <- gradio(c("Yes", "No"),
                               container=lyt,
                               action=code6,
                               horizontal=TRUE,
                               handler=function(h,...) {
                                 oldVal <- svalue(h$obj)
                                 if (oldVal=="Yes") bool <- TRUE
                                 else if (oldVal=="No") bool <- FALSE
                                 svalue(h$action) <- bool
                                 assign("Header", bool, envir=gui.env)
                               })

  separ5 <- lyt[9,1:2] <- gseparator(container=lyt)

  arLabel <- lyt[10,1] <- glabel("Should gimme start with autoregressive effects?",
                                 container=lyt)
  ar <- lyt[10,2] <- gradio(c("Yes", "No"),
                            container=lyt,
                            action=code10,
                            selected=1,
                            horizontal=TRUE,
                            handler=function(h,...) {
                              oldVal <- svalue(h$obj)
                              if (oldVal=="Yes") bool <- TRUE
                              else if (oldVal=="No") bool <- FALSE
                              svalue(h$action) <- bool
                              assign("Ar", bool, envir=gui.env)
                            })
  separ4 <- lyt[11,1:2] <- gseparator(container=lyt)

  plotLabel <- lyt[12,1] <- glabel("Would you like plots?", container=lyt)
  plot <- lyt[12,2] <- gradio(c("Yes", "No"),
                              container=lyt,
                              action=code12,
                              selected=1,
                              horizontal=TRUE,
                              handler=function(h,...) {
                                oldVal <- svalue(h$obj)
                                if (oldVal=="Yes") bool <- TRUE
                                else if (oldVal=="No") bool <- FALSE
                                svalue(h$action) <- bool
                                assign("Plot", bool, envir=gui.env)
                              })
  separ6 <- lyt[13,1:2] <- gseparator(container=lyt)

  subgroupLabel <- lyt[14,1] <- glabel("Should gimme subgroup individuals?",
                                       container=lyt)
  subgroup <- lyt[14,2] <- gradio(c("Yes", "No"),
                                  container=lyt,
                                  action=code14,
                                  selected=2,
                                  horizontal=TRUE,
                                  handler=function(h,...) {
                                    oldVal <- svalue(h$obj)
                                    if (oldVal=="Yes") bool <- TRUE
                                    else if (oldVal=="No") bool <- FALSE
                                    svalue(h$action) <- bool
                                    assign("Subgroup", bool, envir=gui.env)
                                  })
  separ7 <- lyt[15,1:2] <- gseparator(container=lyt)

  pathwayLabel <- lyt[16,1] <- glabel("Would you like to define pathways?",
                                      container=lyt)
  pathway <- lyt[16,2] <- gradio(c("Yes", "No"),
                                 container=lyt,
                                 selected = 2,
                                 horizontal=TRUE,
                                 handler=function(h,...){
                                   oldVal <- svalue(h$obj)
                                   if (oldVal=="Yes") bool <- TRUE
                                   else if (oldVal=="No") bool <- FALSE
                                   assign("Pathway", bool, envir=gui.env)
                                 })
  separ8 <- lyt[17,1:2] <- gseparator(container=lyt)
## inserted by stl
  groupcutLabel <- lyt[18,1] <- glabel("Group path cutoff value:",
                                        container = lyt)
  groupcut <- lyt[18,2] <- gslider(from = 0.5, to = 1, by = 0.01, value = 0.75,
                                    container = lyt,
                                    action = code16,
                                    handler = function(h,...) {
                                      oldVal <- svalue(h$obj)
                                      svalue(h$action) <- oldVal
                                      assign("groupcutoff", oldVal, envir = gui.env)
                                    })
  separ9 <- lyt[19,1:2] <- gseparator(container = lyt)
  subcutLabel <- lyt[20,1] <- glabel("Subgroup path cutoff value:",
                                      container = lyt)
  subcut <- lyt[20,2] <- gslider(from = 0.5, to = 1, by = 0.01, value = 0.5,
                                  container = lyt,
                                  action = code18,
                                  handler = function(h,...) {
                                    oldVal <- svalue(h$obj)
                                    svalue(h$action) <- oldVal
                                    assign("subcutoff", oldVal, envir = gui.env)
                                  })
## end inserted by stl  
  
  srun <- lyt[30,1:2] <- gbutton("Define Pathways",
                                 container=lyt,
                                 handler=function(h,...) {
                                   sswin <- gwindow("Lagged Effects", visible=TRUE)
                                   fgroup <- ggroup(horizontal=FALSE,
                                                    container=sswin)
                                   cgroup <- ggroup(horizontal=FALSE,
                                                    container=fgroup)
                                   dgroup <- ggroup(horizontal=FALSE,
                                                    container=fgroup)
                                   egroup <- ggroup(horizontal=FALSE,
                                                    container=fgroup,
                                                    visible=FALSE)
                                   backb <- gimage(stock.id="go-back",
                                                   dirname="stock",
                                                   container=cgroup,
                                                   size="large_toolbar",
                                                   anchor=c(-1,1),
                                                   handler=function(h,...) {
                                                     if (svalue(sswin)=="Lagged Effects") {
                                                       visible(sswin) <- FALSE
                                                     }
                                                     else if(svalue(sswin)=="Contemporaneous Effects") {
                                                       visible(egroup) <- FALSE
                                                       visible(dgroup) <- TRUE
                                                       svalue(sswin) <- "Lagged Effects"
                                                       delete(egroup, cg)
                                                       delete(egroup, ssrun)
                                                     }
                                                   })
                                   path <- paste(gui.env$Data, list.files(gui.env$Data)[1], sep="/")
                                   dir <- gsub("//", path, replacement="/", fixed=TRUE)
                                   dir <- gsub("\\", path, replacement="/", fixed=TRUE)
                                   assign("Dir", dir, envir=gui.env)
                                   if (gui.env$Sep == "" | gui.env$Sep == ",") {
                                     values <- c(colnames(read.csv(gui.env$Dir, header=gui.env$Header, sep=gui.env$Sep)))
                                   }
                                   else if (gui.env$Sep=="\t") {
                                     values <- c(colnames(read.table(gui.env$Dir, header=gui.env$Header, sep=gui.env$Sep)))
                                   }
                                   ag <- glayout(container=dgroup, spacing=5)
                                   ag[1,1, anchor=c(1,0)] <- "Predicted:   "
                                   ag[1,2, anchor=c(-1,0)] <- "Predictor"
                                   for (i in 1:length(values)) {
                                     ag[i+1, 1, anchor=c(1,0)] <- paste(values[i], ":    ", sep="")
                                     ag[i+1, 2, anchor=c(-1,0)] <- gcheckboxgroup(values, horizontal=TRUE, container=ag)
                                   }
                                   cg <- glayout(spacing=5)
                                   cg[1,1, anchor=c(1,0)] <- "Predicted:   "
                                   cg[1,2, anchor=c(-1,0)] <- "Predictor"
                                   for (i in 1:length(values)) {
                                     cg[i+1, 1, anchor=c(1,0)] <- paste(values[i], ":    ", sep="")
                                     cg[i+1, 2, anchor=c(-1,0)] <- gcheckboxgroup(values[-i], horizontal=TRUE, container=cg)
                                   }
                                   crun <- gbutton("Next", container=dgroup)
                                   ssrun <- gbutton("Done")

                                   pathhand <- addHandlerClicked(ssrun, handler=function(h,...){
                                     patha <- NULL
                                     pathc <- NULL
                                     for (i in 1:length(values)) {
                                       oldValA <- svalue(ag[i+1,2])
                                       oldValC <- svalue(cg[i+1,2])
                                       if (length(oldValA)!=0) {
                                         AVal <- paste(values[i], svalue(ag[i+1,2]), sep="~")
                                         patha <- c(patha, AVal)
                                         pathaa <- paste0(patha, "lag")
                                         assign("patha", pathaa, envir=gui.env)
                                       }
                                       if (length(oldValC)!=0) {
                                         CVal <- paste(values[i], svalue(cg[i+1,2]), sep="~")
                                         pathc <- c(pathc, CVal)
                                         assign("pathc", pathc, envir=gui.env)
                                       }
                                     }
                                     # grab paths from gui environment
                                     oldVal     <- c(gui.env$patha, gui.env$pathc)
                                     assign("Paths", c(oldVal), envir=gui.env)
                                     enabled(ssrun) <- FALSE
                                   })

                                   addHandlerClicked(crun, action=ssrun, handler=function(h,...){
                                     svalue(sswin) <- "Contemporaneous Effects"
                                     visible(dgroup) <- FALSE
                                     visible(egroup) <- TRUE
                                     add(egroup, cg)
                                     add(egroup, ssrun)
                                     blockHandler(ssrun, pathhand)
                                     if (visible(iwin)==TRUE) svalue(h$action) <- "Run indSEM"
                                     else if (visible(awin)==TRUE) svalue(h$action) <- "Run aggSEM"
                                     unblockHandler(ssrun, pathhand)
                                     unblockHandler(grun)
                                   })
                                   addHandlerClicked(ssrun, handler=function(h,...){
                                   dispose(sswin)  
                                   unblockHandler(grun)
                                   enabled(grun) <- TRUE
                                   })
                                   
                                 })
  grun <- lyt[31,1:2] <- gbutton(NULL,
                                 container=lyt)

  enabled(srun) <- FALSE
  enabled(grun) <- TRUE

  grunhand <- addHandlerClicked(grun, handler=function(h,...){
    if (visible(gwin)==TRUE) {
      dataPrint <- gsub("\\\\", "/", gui.env$Data)
      outPrint  <- gsub("\\\\", "/", gui.env$Out)
      if (gui.env$Sep == "\t"){
        sepPrint <- paste0("'", "\\", "t", "'")
      } else {
        sepPrint <- paste0("'", gui.env$Sep, "'")
      }
      fileConn  <- file(file.path(gui.env$Out, "guiSyntax.txt"), "w")
      write(paste0("gimmeSEM(data        = ", "'", dataPrint, "'", ","), fileConn, append = TRUE)
      write(paste0("         ", "out         = ", "'", outPrint, "'", ","), fileConn, append = TRUE)
      write(paste0("         ", "sep         = ", sepPrint, ","), fileConn, append = TRUE)
      write(paste0("         ", "header      = ", ifelse((gui.env$Header==TRUE), "TRUE", "FALSE"), ","), fileConn, append = TRUE)
      write(paste0("         ", "ar          = ", ifelse((gui.env$Ar==TRUE), "TRUE", "FALSE"), ","), fileConn, append = TRUE)
      write(paste0("         ", "plot        = ", ifelse((gui.env$Plot==TRUE), "TRUE", "FALSE"), ","), fileConn, append = TRUE)
      write(paste0("         ", "subgroup    = ", ifelse((gui.env$Subgroup==TRUE), "TRUE", "FALSE"), ","), fileConn, append = TRUE)
      write(paste0("         ", "paths       = ", ifelse((is.null(gui.env$Paths)), "NULL", paste0("c(", paste0(gui.env$Paths, collapse = ","), ")")), ","), fileConn, append = TRUE)
      write(paste0("         ", "groupcutoff = ", gui.env$groupcutoff, ","), fileConn, append = TRUE)
      write(paste0("         ", "subcutoff   = ", ifelse((gui.env$Subgroup==TRUE), gui.env$subcutoff, "NULL"), ")"), fileConn, append = TRUE)
      close(fileConn)
      gimme(data=gui.env$Data, out=gui.env$Out, sep=gui.env$Sep,header=gui.env$Header,
            ar=gui.env$Ar,plot=gui.env$Plot,subgroup=gui.env$Subgroup,paths=gui.env$Paths,
            groupcutoff = gui.env$groupcutoff, subcutoff = gui.env$subcutoff)
    }
    else if (visible(iwin)==TRUE) {
      dataPrint <- gsub("\\\\", "/", gui.env$Data)
      outPrint  <- gsub("\\\\", "/", gui.env$Out)
      if (gui.env$Sep == "\t"){
        sepPrint <- "\"/t\""
      } else {
        sepPrint <- paste0("'", gui.env$Sep, "'")
      }
      fileConn  <- file(file.path(gui.env$Out, "guiSyntax.txt"), "w")
      write(paste0("indSEM(data        = ", "'", dataPrint, "'", ","), fileConn, append = TRUE)
      write(paste0("       ", "out         = ", "'", outPrint, "'", ","), fileConn, append = TRUE)
      write(paste0("       ", "sep         = ", sepPrint, ","), fileConn, append = TRUE)
      write(paste0("       ", "header      = ", ifelse((gui.env$Header==TRUE), "TRUE", "FALSE"), ","), fileConn, append = TRUE)
      write(paste0("       ", "ar          = ", ifelse((gui.env$Ar==TRUE), "TRUE", "FALSE"), ","), fileConn, append = TRUE)
      write(paste0("       ", "plot        = ", ifelse((gui.env$Plot==TRUE), "TRUE", "FALSE"), ","), fileConn, append = TRUE)
      write(paste0("       ", "paths       = ", ifelse((is.null(gui.env$Paths)), "NULL", paste0("c(", paste0(gui.env$Paths, collapse = ","), ")")), ")"), fileConn, append = TRUE)
      close(fileConn)
      indSEM(data=gui.env$Data,out=gui.env$Out, sep=gui.env$Sep,header=gui.env$Header,
             ar=gui.env$Ar,plot=gui.env$Plot, paths = gui.env$Paths)
    }
    else if (visible(awin)==TRUE) {
      dataPrint <- gsub("\\\\", "/", gui.env$Data)
      outPrint  <- gsub("\\\\", "/", gui.env$Out)
      if (gui.env$Sep == "\t"){
        sepPrint <- "\"/t\""
      } else {
        sepPrint <- paste0("'", gui.env$Sep, "'")
      }
      fileConn  <- file(file.path(gui.env$Out, "guiSyntax.txt"), "w")
      write(paste0("aggSEM(data        = ", "'", dataPrint, "'", ","), fileConn, append = TRUE)
      write(paste0("       ", "out         = ", "'", outPrint, "'", ","), fileConn, append = TRUE)
      write(paste0("       ", "sep         = ", sepPrint, ","), fileConn, append = TRUE)
      write(paste0("       ", "header      = ", ifelse((gui.env$Header==TRUE), "TRUE", "FALSE"), ","), fileConn, append = TRUE)
      write(paste0("       ", "ar          = ", ifelse((gui.env$Ar==TRUE), "TRUE", "FALSE"), ","), fileConn, append = TRUE)
      write(paste0("       ", "plot        = ", ifelse((gui.env$Plot==TRUE), "TRUE", "FALSE"), ","), fileConn, append = TRUE)
      write(paste0("       ", "paths       = ", ifelse((is.null(gui.env$Paths)), "NULL", paste0("c(", paste0(gui.env$Paths, collapse = ","), ")")), ")"), fileConn, append = TRUE)
      close(fileConn)
      aggSEM(data=gui.env$Data, out=gui.env$Out, sep=gui.env$Sep, header=gui.env$Header,
             ar=gui.env$Ar, plot=gui.env$Plot, paths = gui.env$Paths)
    }
  })

  addHandlerChanged(pathway, handler=function(h,...) {
    if (svalue(pathway)=="Yes" & !is.null(gui.env$Data) & !is.null(gui.env$Out)) {
      enabled(srun) <- TRUE
      enabled(grun) <- FALSE
    }
    else if (svalue(pathway)=="No") {
      enabled(srun) <- FALSE
      enabled(grun) <- TRUE
    }
  })
  
  # if the user hasn't yet specified directories, 
  # don't allow them to specify paths, as it will fail
  # without being able to read in the first data file
  # below code ensures the above works properly
  addHandlerChanged(data, handler=function(h,...) {
    if (!is.null(gui.env$Out) & !is.null(gui.env$Sep) & svalue(pathway)=="Yes"){
      enabled(srun) <- TRUE
      enabled(grun) <- FALSE
    }
    else if (!is.null(gui.env$Out) | !is.null(gui.env$Sep) | svalue(pathway)=="No") {
      enabled(srun) <- FALSE
      enabled(grun) <- TRUE
    }
  })
  
  addHandlerChanged(out, handler=function(h,...) {
    if (!is.null(gui.env$Data) & !is.null(gui.env$Sep) & svalue(pathway)=="Yes"){
      enabled(srun) <- TRUE
      enabled(grun) <- FALSE
    }
    else if (!is.null(gui.env$Data) | !is.null(gui.env$Sep) | svalue(pathway)=="No") {
      enabled(srun) <- FALSE
      enabled(grun) <- TRUE
    }
  })
  
  addHandlerChanged(sep, handler=function(h,...) {
    if (!is.null(gui.env$Data) & !is.null(gui.env$Out) & svalue(pathway)=="Yes"){
      enabled(srun) <- TRUE
      enabled(grun) <- FALSE
    }
    else if (!is.null(gui.env$Data) | !is.null(gui.env$Out) | svalue(pathway)=="No") {
      enabled(srun) <- FALSE
      enabled(grun) <- TRUE
    }
  })
  
  if (svalue(pathway) == "Yes"){
    addHandlerClicked(ssrun, handler=function(h,...){
      enabled(grun) <- TRUE
      enabled(srun) <- FALSE
    })
  }

  addHandlerClicked(gimmebutton, handler=function(h,...){
    visible(introwin) <- FALSE
    visible(gwin)     <- TRUE
    add(gwin, agroup)
    visible(code13) <- TRUE
    visible(code14) <- TRUE
    visible(code15) <- TRUE
    visible(code16) <- TRUE
    visible(code17) <- FALSE
    visible(code18) <- FALSE
    enabled(subcut) <- FALSE
    blockHandler(grun, grunhand)
    svalue(h$action) <- "Run gimmeSEM"
    unblockHandler(grun, grunhand)
  }, action=grun)

  addHandlerChanged(subgroup, handler=function(h,...) {
    if (svalue(subgroup) == "Yes") {
      enabled(subcut) <- TRUE
      enabled(grun)   <- TRUE
      visible(code17) <- TRUE
      visible(code18) <- TRUE
    }
    else if (svalue(subgroup) == "No") {
      enabled(subcut) <- FALSE
      enabled(grun)   <- TRUE
      visible(code17) <- FALSE
      visible(code18) <- FALSE
    }
  })
  
  
  addHandlerClicked(indSEMbutton, handler=function(h,...){
    visible(introwin) <- FALSE
    visible(iwin)     <- TRUE
    visible(separ6)   <- FALSE
    visible(subgroupLabel) <- FALSE
    visible(subgroup) <- FALSE
    visible(code13)   <- FALSE
    visible(code14)   <- FALSE
    visible(code15)   <- FALSE
    visible(code16)   <- FALSE
    visible(code17)   <- FALSE
    visible(code18)   <- FALSE
    visible(separ9)   <- FALSE
    visible(groupcut) <- FALSE
    visible(subcut)   <- FALSE
    visible(groupcutLabel) <- FALSE
    visible(subcutLabel)   <- FALSE
    add(iwin, agroup)
    blockHandler(grun, grunhand)
    svalue(h$action) <- "Run indSEM"
    unblockHandler(grun, grunhand)
  }, action=grun)

  addHandlerClicked(aggbutton, handler=function(h,...){
    visible(introwin) <- FALSE
    visible(awin)     <- TRUE
    visible(separ6)   <- FALSE
    visible(subgroupLabel) <- FALSE
    visible(subgroup) <- FALSE
    visible(groupcut) <- FALSE
    visible(subcut)   <- FALSE
    visible(separ9)   <- FALSE
    visible(groupcutLabel) <- FALSE
    visible(subcutLabel)   <- FALSE
    visible(code13) <- FALSE
    visible(code14) <- FALSE
    visible(code15) <- FALSE
    visible(code16) <- FALSE
    visible(code17) <- FALSE
    visible(code18) <- FALSE
    add(awin, agroup)
    blockHandler(grun, grunhand)
    svalue(h$action) <- "Run aggSEM"
    unblockHandler(grun, grunhand)
  }, action=grun)

addHandlerClicked(grun, handler=function(h,...){
  visible(gwin) <- FALSE
  visible(awin) <- FALSE
  visible(iwin) <- FALSE
})

}

