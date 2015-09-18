setup <- function (data,
                   sep,
                   header,
                   out,
                   plot,
                   ar) {
  
  files            <- list.files(data, full.names=TRUE)
  subjects         <- length(files)
  trackparts       <- matrix(0,nrow=subjects,ncol=2)
  trackparts[,1]   <- seq(1:subjects)
  trackparts[,2]   <- sapply(strsplit(basename(files),"\\."), 
                             function(x) paste(x[1:(length(x)-1)], collapse=".")) 
  
  all               <- as.matrix(read.table(files[1],sep=sep,header=header))
  rois              <- ncol(all)
  vars              <- rois*2
  varnames          <- colnames(all)
  varnames          <- rep(varnames,2)
  cutoffind         <- qchisq(.99,1)
  lvarnames         <- character()
  count.group.paths <- 0
  
  for (j in 1:(rois)) varnames[j]           <- paste(varnames[j],"lag",sep="")
    
  for (j in 1:rois) lvarnames[j]            <- paste("VAR",j,"lag",sep="")
  
  for (j in (rois+1):(rois*2)) lvarnames[j] <- paste("VAR",j-rois,sep="")
  
  x                 <- seq(1:vars)
  y                 <- substring(lvarnames,4)
  individual        <- file.path(out,"individual")
  subgroup.dir          <- file.path(out,"subgroup")
  betas             <- file.path(individual,"betas")
  SEs               <- file.path(individual,"SEs")
  fitted            <- file.path(out,"fitted")

  dir.create(subgroup.dir,      showWarnings=FALSE)
  dir.create(individual,    showWarnings=FALSE)
  dir.create(betas,         showWarnings=FALSE) 
  dir.create(SEs,           showWarnings=FALSE)
  dir.create(fitted,        showWarnings=FALSE)
  
  if (plot==TRUE) {plots      <- file.path(individual,"plots")
                   dir.create(plots,showWarnings=FALSE)
                   plot.names <-  varnames[(rois+1):(rois*2)]
  } else {plots <- ""; plot.names <-""}
  
  line1             <- paste(capture.output(for (i in 1:(rois*2)){ 
    cat(lvarnames[i],"=~","1*",varnames[i],
        sep="","\n")}),collapse="\n")
  line2             <- paste(capture.output(for (i in 1:rois) {for (j in 1:i){ 
    cat(lvarnames[i],"~~",lvarnames[j],
        sep="","\n")}}),collapse="\n")
  line3             <- paste(capture.output(for (i in (rois+1):(rois*2)) {
    cat (lvarnames[i],"~~",lvarnames[i],
         sep="","\n")}),collapse="\n")
  
  if (ar==TRUE) {line4 <- paste(capture.output(for (j in 1:rois) {
    cat (lvarnames[j+rois],"~",lvarnames[j],
         sep="","\n")}),collapse="\n")}
  
  if (ar==FALSE) {line4 <- paste(capture.output(for (j in 1:rois) {
    cat (lvarnames[j],"~0*",lvarnames[j+rois],
         sep="","\n")}),collapse="\n")}
  
  syntax            <- paste(line1,line2,line3,line4,sep="\n")
  candidate.paths   <- capture.output(
    for (i in (rois+1):(rois*2)) 
    {for (j in 1:(rois*2)) 
    {cat(lvarnames[i],"~",lvarnames[j],
         sep="","\n")}})

  ar.paths          <- capture.output(for (i in 1:rois){ 
    cat(lvarnames[i+rois],"~",lvarnames[i],
        sep="","\n")})
  
  list <- list("subjects"          = subjects,
               "rois"              = rois, 
               "x"                 = x,
               "y"                 = y,
               "varnames"          = varnames,
               "trackparts"        = trackparts,
               "vars"              = vars,
               "lvarnames"         = lvarnames, 
               "cutoffind"         = cutoffind,
               "betas"             = betas,
               "SEs"               = SEs,
               "fitted"            = fitted,
               "plots"             = plots,
               "plot.names"        = plot.names,
               "syntax"            = syntax,
               "candidate.paths"   = candidate.paths,
               "ar.paths"          = ar.paths,
               "data"              = data,
               "ar"                = ar,
               "sep"               = sep,
               "header"            = header,
               "out"               = out,
               "plot"              = plot,
               "count.group.paths" = count.group.paths,
               "subgroup.dir"      = subgroup.dir)
  return(list)
}

fit.model <- function (varnames, 
                       syntax, 
                       data.file) {
  
  fit <- tryCatch(lavaan(syntax, 
                    data            = data.file, 
                    model.type      = "sem", 
                    sample.nobs     = time,
                    missing         = "fiml", 
                    estimator       = "ml",
                    int.ov.free     = TRUE,
                    int.lv.free     = FALSE,
                    auto.fix.first  = TRUE,
                    auto.var        = TRUE,
                    auto.cov.lv.x   = TRUE, 
                    auto.th         = TRUE,
                    auto.delta      = TRUE,
                    auto.cov.y      = FALSE,
                    auto.fix.single = TRUE,
                    warn            = FALSE),
                  error=function(e) e)
  return(fit)
}

read.data <- function (file,
                       sep,
                       header) {
  all    <- as.matrix(read.table(file,sep=sep,header=header))
  first  <- all[1:(nrow(all)-1),]
  second <- all[2:nrow(all),]
  full   <- cbind(first,second)
  data.file <- as.data.frame(full)
  return(data.file)
}


miSEM <- function (setup.out,
                   previous.out,
                   subgroup.step,
                   subjects,
                   files,
                   syntax,
                   subgroup.paths,
                   second.round) {
  
  vars            <- setup.out$vars
  varnames        <- setup.out$varnames
  candidate.paths <- setup.out$candidate.paths

  ar              <- setup.out$ar
  sep             <- setup.out$sep
  header          <- setup.out$header
  data            <- setup.out$data
    
  count.group.paths <- previous.out$count.group.paths
  
  mi.index             <- matrix(1:((vars-1)*vars),nrow=((vars-1)*vars),ncol=1)
  
  if (subgroup.step == FALSE) {
    count.group.paths <- 0
    count.subgroup.paths <- 0
    subjects          <- setup.out$subjects
    files  <- list.files(data, full.names=TRUE)
    syntax            <- previous.out$syntax
  }
  
  if (subgroup.step == TRUE & second.round==FALSE){
  count.subgroup.paths <- 0
  subgroup.paths       <- character()
  }
  
  if (subgroup.step == TRUE & second.round==TRUE){
    subgroup.paths       <- unlist(strsplit(subgroup.paths, "[,]"))
    count.subgroup.paths <- length(subgroup.paths)
  }
  
  cutoff               <- qchisq(.95,subjects)
  cutoffgroup          <- qchisq(1-.05/subjects,1)
  
  continue <- 1

  while (continue == 1) {
    mi.list            <- matrix(0, nrow=vars*(vars-1)*subjects, ncol=6)
    colnames(mi.list)  <- c("subject", "index", "lhs", "op", "rhs", "mi")
    count.converge   <- 0
    for (k in 1:subjects)
    {
      data.file           <- read.data(file   = files[k],
                                       sep    = sep,
                                       header = header)
      time                <- nrow(data.file)
      colnames(data.file) <- c(varnames)
      fit                 <- fit.model(varnames  = varnames,  
                                       syntax    = syntax,
                                       data.file = data.file)
      
      check.npd           <- any(grepl("error",class(fit))==TRUE)
      check.not.identified <- sum(lavInspect(fit,"se")$beta,na.rm=TRUE)==0
    
      if (check.npd == FALSE & check.not.identified==FALSE) { 
        
      singular            <- tryCatch(modindices(fit),error=function(e) e)    
      check.singular      <- any(grepl("singular",singular)==TRUE)

      converge            <- lavInspect(fit, "converged")
      } else {
        if (subgroup.step==FALSE) {writeLines(paste("group-level search, subject", k, "nonconvergence"))
        } else {writeLines(paste("subgroup-level search, subject", k, "nonconvergence"))
        }
        check.singular <- TRUE
        converge       <- FALSE
      }

      if (check.singular==FALSE & converge==TRUE & check.not.identified==FALSE) {
        if (subgroup.step==FALSE) {writeLines(paste("group-level search, subject", k))
        } else {writeLines(paste("subgroup-level search, subject", k))
        }
                                                    
        mi                <- as.matrix(singular[singular$op == "~",])[,1:4]
        count.converge    <- count.converge + 1
      }
      # if it doesn't converge or is computationally singular or NPD 
      if (converge==FALSE | check.singular==TRUE | check.npd == TRUE | check.not.identified==TRUE) {
        mi                <- matrix(NA, nrow=((vars-1)*vars), ncol=4)
      }
      #stacking matrices
      mi.subject          <- matrix(k,nrow=((vars-1)*vars),ncol=1)
      mi.list[(((nrow(mi)*k)-nrow(mi))+1):(nrow(mi)*k),(1:6)] <- 
        as.matrix(cbind(mi.subject,mi.index,mi),rownames.force=FALSE) 
    }
    
    mi.all                <- as.data.frame(mi.list[complete.cases(mi.list),])
    mi.all$param          <- paste(mi.all$lhs, mi.all$op, mi.all$rhs, sep="")
    mi.all                <- mi.all[-c(3:5)]
    mi.all                <- subset(mi.all,param %in% candidate.paths) 
    mi.all[,3]            <- as.numeric(as.character(mi.all[,3]))
    mi.all[,1]            <- as.numeric(as.character(mi.all[,1]))
    mi.all                <- transform(mi.all, sum = ave(mi, param, FUN=sum))
    mi.high               <- subset(mi.all, mi>cutoffgroup)
    mi.count              <- subset(as.data.frame(table(mi.high$param)),Freq>0)
    mi.high.count         <- subset(mi.high, !duplicated(param))
    mi.merge              <- merge(x=mi.high.count, y=mi.count, 
                                   by.x="param", by.y="Var1")
    paramadd              <- mi.merge[order(-mi.merge$Freq, -mi.merge$sum),][1,1]
    y                     <- subset(mi.all, param==paramadd, select=mi)$mi
    drop.element          <- ifelse(max(y)<cutoffgroup,TRUE,FALSE)
    prop                  <- sum(y > cutoffgroup)/count.converge
    
    # create code here to check number of people in subgroup
    
    if (subgroup.step==TRUE){
      n.paths               <- length(y)
      if (n.paths < 5) {
        prop                <- sum(y>cutoffgroup)/n.paths
        if (prop>.75){integratelow   = 0
                      integratehigh  = 1}
        if (prop<=.75){integratelow  = 1
                       integratehigh = 0}
      }
    }
    
    if (subgroup.step==TRUE)  cutoffprop <- .5
    if (subgroup.step==FALSE) cutoffprop <- .75
    
    halt <- length(y)<=subjects/2

    if (prop <= cutoffprop | drop.element==TRUE | halt==TRUE) {continue <-0 
    } else {
      continue            <- 1
      syntax              <- paste(syntax,as.name(paramadd),sep="\n")
      if (subgroup.step==FALSE) count.group.paths    <- count.group.paths + 1
      if (subgroup.step==TRUE) {count.subgroup.paths <- count.subgroup.paths + 1
                                subgroup.paths       <- append(subgroup.paths,paramadd)}               
    }

  }  #end of while continue>0
  list <- list("syntax"               = syntax,
               "count.group.paths"    = count.group.paths,
               "count.subgroup.paths" = count.subgroup.paths,
               "subgroup.paths"       = subgroup.paths)
  return(list)
}

evalbetas <- function (setup.out,
                       previous.out,
                       subgroup.step,
                       s,
                       subjects,
                       files,
                       post.sub.prune,
                       evalbetas.out) {
  
  data              = setup.out$data
  ar.paths          = setup.out$ar.paths
  ar                = setup.out$ar
  varnames          = setup.out$varnames
  rois              = setup.out$rois
  sep               = setup.out$sep
  header            = setup.out$header
  


  if (subgroup.step==FALSE & post.sub.prune == FALSE) {
    bad                  = ifelse(identical(setup.out$syntax,previous.out$syntax),0,1)
    syntax               = previous.out$syntax
    subjects             = setup.out$subjects
    count.subgroup.paths = 0
    files                <- list.files(data, full.names=TRUE)
    count.group.paths = previous.out$count.group.paths
  }
  
  if (subgroup.step==TRUE){
    bad = ifelse(previous.out$unique.syntax[s,3]==0,0,1)  
    count.subgroup.paths = as.numeric(previous.out$unique.syntax[s,3])
    syntax               =  previous.out$unique.syntax[s,2]
    subgroup.paths       <- previous.out$unique.syntax[s,4]
    subgroup.paths       <- unlist(strsplit(subgroup.paths, "[,]"))
    count.group.paths    = previous.out$count.group.paths
  }

  if (post.sub.prune == TRUE){
    bad         = 1
    group       <- unlist(strsplit(evalbetas.out$syntax,"\n"))
    start       <- unlist(strsplit(setup.out$syntax,"\n"))
    group.paths <- group[-(1:length(start))]
    subjects    <- setup.out$subjects
    syntax.all  <- as.matrix(previous.out$syntax.sub$syntax)
    files       <- as.matrix(previous.out$syntax.sub$files)
    count.group.paths <- length(group.paths)  
  }
  
  cutoffz           = abs(qnorm(.05/subjects))
  if (count.group.paths==0) bad = 0
  #getting coefficients for final model
  while (bad == 1) {
    if (post.sub.prune == FALSE) list.all <- matrix(NA, nrow=(rois+count.group.paths+count.subgroup.paths)*subjects, ncol=2)
    if (post.sub.prune == TRUE) list.all <- matrix(NA,nrow=count.group.paths*subjects,ncol=2)
    colnames(list.all) <- c("param", "z")
    count.converge     <- 0
    
    for (k in 1:subjects)
    {
      if (post.sub.prune==TRUE) syntax  <- syntax.all[k,]
      
      data.file           <- read.data(file   = files[k],
                                       sep    = sep,
                                       header = header)
      time                <- nrow(data.file)
      colnames(data.file) <- c(varnames)
      fit                 <- fit.model(varnames  = varnames, 
                                       syntax    = syntax,
                                       data.file = data.file)
      
      check.npd             <- any(grepl("error",class(fit))==TRUE)
      check.not.identified  <- sum(lavInspect(fit,"se")$beta,na.rm=TRUE)==0
      
      if (check.npd == FALSE & check.not.identified==FALSE) { 
        singular            <- tryCatch(modindices(fit),error=function(e) e)
        check.singular      <- any(grepl("singular",singular)==TRUE)
        converge            <- lavInspect(fit, "converged")
      } else {
        check.singular <- TRUE
        converge       <- FALSE
      }
      
      if (check.singular==FALSE & converge==TRUE & check.not.identified==FALSE) {
        ## printing the step
        if (subgroup.step==FALSE) {writeLines(paste("group-level search, subject", k))
        } else {writeLines(paste("subgroup-level search, subject", k))
        }
       
        z.list         <- subset(standardizedSolution(fit),op=="~")
        z.list$param   <- paste(z.list$lhs, z.list$op, z.list$rhs, sep="")
        z.list         <- z.list[c("param","z")]
        z.list         <- z.list[complete.cases(z.list),]
        if (post.sub.prune == TRUE) z.list <- z.list[z.list$param %in% group.paths,]   
        count.converge <- count.converge + 1
      }
      # if it doesn't converge
      if (converge == FALSE | check.singular == TRUE | check.npd == TRUE | check.not.identified==TRUE) {
        z.list <- matrix(NA, nrow=(rois+count.group.paths), ncol=2) 
        if (post.sub.prune == TRUE) z.list <- matrix(NA, nrow=length(group.paths), ncol=2) 
        ## printing the step
        if(subgroup.step==FALSE) {writeLines(paste("group-level search, subject", k, "nonconvergence"))
        } else {writeLines(paste("subgroup-level search, subject", k, "nonconvergence"))
        }
      }  
      list.all[(((nrow(z.list)*k)-nrow(z.list)+1):(nrow(z.list)*k)),] <- as.matrix(z.list)
    }
    
    list.all     <- as.data.frame(list.all)
    list.all$z   <- as.numeric(as.character(list.all$z))
    list.all$sig <- ifelse(abs(list.all$z)>cutoffz, 1, 0)
    list.all     <- transform(list.all, sum=ave(sig, param, FUN=sum))
    # add line to not include AR as 
    if (ar==TRUE){
      list.all <- subset(list.all,!(param %in% ar.paths))
    }
    # insert part here for subgroup.step where only those paths
    # added unique to subgrouping can be removed in evalbetas
    if (subgroup.step==TRUE) list.all <- subset(list.all,param %in% subgroup.paths)
    if ((nrow(list.all)==0)==TRUE) bad <- 0
    paramdrop <- as.character(list.all[which.min(list.all$sum),1])
    #create histogram of t/z values

    y2   <- abs(subset(list.all, param==paramdrop, select=z)$z) 
    
    prop <- sum(y2>=cutoffz)/count.converge
    
    drop.element <- ifelse(max(y2)<cutoffz,TRUE,FALSE)
    
    if (subgroup.step == TRUE) cutoffprop <- .5
    if (subgroup.step == FALSE)  cutoffprop <- .75
    
    if (prop <= cutoffprop | drop.element==TRUE) {
      syntax <- unlist(strsplit(syntax, "[\n]"))
      syntax <- syntax[!syntax %in% paramdrop]
      syntax <- paste(syntax,sep="",collapse="\n")
      bad    <- 1
      if (subgroup.step==TRUE) {count.subgroup.paths <- count.subgroup.paths - 1
                                if ((count.subgroup.paths==0)==TRUE) bad <- 0
                                subgroup.paths <- subgroup.paths[!subgroup.paths %in% paramdrop]
      }
      if (subgroup.step==FALSE) {count.group.paths <- count.group.paths - 1
                                 if ((count.group.paths==0)==TRUE) bad <- 0
                                 subgroup.paths <- NULL
      }
      if (post.sub.prune==TRUE) {     
        for (i in 1:nrow(syntax.all)){
          syntax.all.paths <- unlist(strsplit(syntax.all[i,],"\n"))
          syntax.all.paths <- syntax.all.paths[!syntax.all.paths %in% paramdrop]
          syntax.all[i,]   <- paste(syntax.all.paths,sep="",collapse="\n")
        }
      }
    } else {bad <- 0
            if (subgroup.step==FALSE) subgroup.paths <- NULL
    }
    if (length(y2)<=subjects/2) bad <- 0
  }

  
  if (subgroup.step==FALSE) subgroup.paths <- NULL
  if (subgroup.step==TRUE)  subgroup.paths <- paste(subgroup.paths,sep="",collapse=",")
  if (post.sub.prune==TRUE) {syntax <- syntax.all
                             count.subgroup.paths <- NULL
  }
  
  list <- list("syntax"               =syntax,
               "count.group.paths"    =count.group.paths,
               "count.subgroup.paths" =count.subgroup.paths,
               "subgroup.paths"       =subgroup.paths)
  return(list)
}

