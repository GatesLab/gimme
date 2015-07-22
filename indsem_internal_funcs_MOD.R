addind <- function (done,
                    evaluate,
                    syntax,
                    data.file,
                    setup.out) {
  
  varnames        = setup.out$varnames
  cutoffind       = setup.out$cutoffind
  candidate.paths = setup.out$candidate.paths
  
  while (done==0) {
    count.ind.paths     <- 0
    vec.MI              <- character()
    time                <- nrow(data.file)
    colnames(data.file) <- c(varnames)
    
    fit <- fit.model(varnames  = varnames,
                     syntax    = syntax,
                     data.file = data.file)
    
    check.npd            <- any(grepl("error",class(fit))==TRUE)
    check.not.identified <- sum(lavInspect(fit,"se")$beta,na.rm=TRUE)==0
    
    if (check.npd == FALSE & check.not.identified==FALSE) { 
      
      singular            <- tryCatch(modindices(fit),error=function(e) e)    
      check.singular      <- any(grepl("singular",singular)==TRUE)
      converge            <- lavInspect(fit, "converged")
    } else {
      check.singular <- TRUE
      converge       <- FALSE
    }
    
    if (check.singular==FALSE & converge==TRUE & check.npd == FALSE & check.not.identified==FALSE) {
      indMI          <- as.matrix(singular[singular$op == "~",])[,1:4]
      indMI          <- as.data.frame(indMI[complete.cases(indMI),])
      indMI$param    <- paste(indMI$lhs, indMI$op, indMI$rhs, sep="")
      indMI          <- subset(indMI,param %in% candidate.paths)
      indMI$mi       <- as.numeric(as.character(indMI$mi))
      indparamadd    <- indMI[which.max(indMI$mi),5]
      indparamaddval <- indMI[which.max(indMI$mi),4]
      indices        <- fitMeasures(fit,c("chisq","df","pvalue","rmsea","srmr",
                                          "nnfi","cfi"))  
      rmsea   <- indices[4]
      srmr    <- indices[5]
      cfi     <- indices[6]
      nnfi    <- indices[7]
      rmseaE  <- ifelse(rmsea<.05, 1, 0)
      srmrE   <- ifelse(srmr<.05, 1, 0)
      cfiE    <- ifelse(cfi>.95, 1, 0)
      nnfiE   <- ifelse(nnfi>.95, 1, 0)
      excellent <- sum(rmseaE, srmrE, cfiE, nnfiE)
      if (excellent >= 2) {
        done   <- 1 
        fixfit <- 0
      } else if (excellent <2) {
        done <- 0
        if (indparamaddval >= cutoffind) {
          syntax          <- paste(syntax,as.name(indparamadd),sep="\n")
          count.ind.paths <- count.ind.paths + 1
          vec.MI          <- append(vec.MI,indparamadd)
        } else {
          done   <- 1 
          fixfit <- 1
        }
      }
    } else {done=1;fixfit=0}
    fit <- fit.model(varnames  = varnames,
                     syntax    = syntax,
                     data.file = data.file)
    
    check.npd           <- any(grepl("error",class(fit))==TRUE)
    check.not.identified <- sum(lavInspect(fit,"se")$beta,na.rm=TRUE)==0
    
    if (check.npd == FALSE & check.not.identified==FALSE) { 
      singular            <- tryCatch(modindices(fit),error=function(e) e)
      check.singular      <- any(grepl("singular",singular)==TRUE)
      converge            <- lavInspect(fit, "converged")
    } else {
      check.singular <- TRUE
      converge       <- FALSE
    }
    
    if (converge==TRUE & check.singular==FALSE & check.npd ==FALSE & check.not.identified==FALSE){
      indices   <- fitMeasures(fit,c("chisq","df","pvalue","rmsea","srmr",
                                     "nnfi","cfi"))  
      rmsea     <- indices[4]
      srmr      <- indices[5]
      cfi       <- indices[6]
      nnfi      <- indices[7]
      rmseaE    <- ifelse(rmsea<.05, 1, 0)
      srmrE     <- ifelse(srmr<.05, 1, 0)
      cfiE      <- ifelse(cfi>.95, 1, 0)
      nnfiE     <- ifelse(nnfi>.95, 1, 0)
      excellent <- sum(rmseaE, srmrE, cfiE, nnfiE)
      if (excellent >= 2) {done <- 1; fixfit <- 0}
    } else {done <- 1;  fixfit <- 0;evaluate <- 0}  
    if (converge==FALSE) done <- 1
    if (check.singular==TRUE) done <-1
    if (check.npd == TRUE) done <- 1
  }
  if (count.ind.paths==0) evaluate <- 0
  list <- list("evaluate" = evaluate,
               "fixfit"   = fixfit,
               "syntax"   = syntax,
               "vec.MI"   = vec.MI)
  return(list)
}  


evalind <- function (addind.out,
                     setup.out,
                     data.file) {
  
  evaluate  = addind.out$evaluate
  varnames  = setup.out$varnames
  syntax    = addind.out$syntax
  fixfit    = addind.out$fixfit
  vec.MI    = addind.out$vec.MI
  
  while (evaluate==1) {
    fit <- fit.model(varnames  = varnames, 
                     syntax    = syntax,
                     data.file = data.file)  
    
    check.npd            <- any(grepl("error",class(fit))==TRUE)
    check.not.identified <- sum(lavInspect(fit,"se")$beta,na.rm=TRUE)==0
    
    if (check.npd == FALSE & check.not.identified==FALSE) { 
      singular            <- tryCatch(modindices(fit),error=function(e) e)
      check.singular      <- any(grepl("singular",singular)==TRUE)
      converge            <- lavInspect(fit, "converged")
    } else {
      check.singular <- TRUE
      converge       <- FALSE
    }
    
    if (check.singular==FALSE & converge==TRUE & check.npd==FALSE) {
      indlist        <- subset(standardizedSolution(fit),op=="~")
      indlist$param  <- paste(indlist$lhs, indlist$op, indlist$rhs, sep="")
      indlist        <- indlist[c("param","z")]
      indlist        <- indlist[complete.cases(indlist),]
      indlist        <- as.data.frame(indlist)
      indlist$z      <- as.numeric(as.character(indlist$z))
      indlist        <- subset(indlist,param %in% vec.MI)
      parampruneposs <- as.character(indlist[which.min(indlist$z),1])
      parampruneval  <- as.numeric(indlist[which.min(indlist$z),2])
      prune          <- ifelse(parampruneval>1.96, 0, 1)
      if (nrow(indlist)==0) prune <- 0
      if (prune == 1) {
        syntax       <- unlist(strsplit(syntax, "[\n]"))
        syntax       <- syntax[!syntax %in% parampruneposs]
        syntax       <- paste(syntax,sep="",collapse="\n")
        paramprune   <- parampruneposs
        evaluate     <- 1
      } else {evaluate <- 0;fixfit <- 0}
      fit <- fit.model(varnames=varnames,  
                       syntax=syntax,
                       data.file=data.file) 
      check.npd            <- any(grepl("error",class(fit))==TRUE)
      check.not.identified <- sum(lavInspect(fit,"se")$beta,na.rm=TRUE)==0
      
      if (check.npd == FALSE & check.not.identified==FALSE) { 
        
        singular            <- tryCatch(modindices(fit),error=function(e) e)    
        check.singular      <- any(grepl("singular",singular)==TRUE)
        converge            <- lavInspect(fit, "converged")
      } else {
        check.singular <- TRUE
        converge       <- FALSE
      }
      
      if (converge==TRUE & check.singular==FALSE & check.npd==FALSE & check.not.identified==FALSE){
        indices      <- fitMeasures(fit,c("chisq","df","pvalue","rmsea","srmr",
                                          "nnfi","cfi"))  
        rmsea     <- indices[4]
        srmr      <- indices[5]
        cfi       <- indices[6]
        nnfi      <- indices[7]
        rmseaE    <- ifelse(rmsea<.05, 1, 0)
        srmrE     <- ifelse(srmr<.05, 1, 0)
        cfiE      <- ifelse(cfi>.95, 1, 0)
        nnfiE     <- ifelse(nnfi>.95, 1, 0)
        excellent <- sum(rmseaE, srmrE, cfiE, nnfiE)
        if (excellent >= 2) fixfit <-0 else fixfit <-1
        # if it doesn't converge
      }
      if (converge==FALSE | check.singular==TRUE) {evaluate <- 0; fixfit <- 0}
    }
  }
  list <- list("fixfit" = fixfit,
               "syntax" = syntax)
  return(list)
}

fixfitind <- function (setup.out,
                       evalind.out,
                       data.file) {
  
  fixfit          = evalind.out$fixfit
  varnames        = setup.out$varnames
  syntax          = evalind.out$syntax
  candidate.paths = setup.out$candidate.paths
  
  while (fixfit==1) {
    fit <- fit.model(varnames  = varnames,  
                     syntax    = syntax,
                     data.file = data.file)
    
    check.npd            <- any(grepl("error",class(fit))==TRUE)
    check.not.identified <- sum(lavInspect(fit,"se")$beta,na.rm=TRUE)==0
    
    if (check.npd == FALSE & check.not.identified==FALSE) { 
      singular            <- tryCatch(modindices(fit),error=function(e) e)
      check.singular      <- any(grepl("singular",singular)==TRUE)
      converge            <- lavInspect(fit, "converged")
    } else {
      check.singular <- TRUE
      converge       <- FALSE
    }
    
    if (check.singular == FALSE & converge == TRUE & check.npd==FALSE & check.not.identified==FALSE) {
      indMI        <- as.matrix(singular[singular$op == "~",])[,1:4]
      indMI        <- as.data.frame(indMI[complete.cases(indMI),])
      indMI$param  <- paste(indMI$lhs, indMI$op, indMI$rhs, sep="")
      indMI        <- subset(indMI,param %in% candidate.paths)
      indMI$mi     <- as.numeric(as.character(indMI$mi))
      indparamadd  <- indMI[which.max(indMI$mi),5]
      syntax       <- paste(syntax,as.name(indparamadd),sep="\n")     
      fit          <- fit.model(varnames = varnames,  
                                syntax = syntax,
                                data.file = data.file)
      converge       <- lavInspect(fit, "converged")
      singular       <- tryCatch(modindices(fit),error=function(e) e)
      check.singular <- any(grepl("singular",singular)==TRUE)
      if (check.singular==FALSE & converge==TRUE) {
        indices      <- fitMeasures(fit,c("chisq","df","pvalue","rmsea","srmr",
                                          "nnfi","cfi"))  
        rmsea  <- indices[4]
        srmr   <- indices[5]
        cfi    <- indices[6]
        nnfi   <- indices[7]
        df     <- indices[2]
        rmseaE <- ifelse(rmsea<.05, 1, 0)
        srmrE  <- ifelse(srmr<.05, 1, 0)
        cfiE   <- ifelse(cfi>.95, 1, 0)
        nnfiE  <- ifelse(nnfi>.95, 1, 0)
        if (df == 0) fixfit <- 0
        excellent <- sum(rmseaE, srmrE, cfiE, nnfiE)
        if (excellent >= 2) fixfit<-0 else fixfit<-1    
      }  
      # if it doesn't converge
      if (converge==FALSE | check.singular == TRUE) fixfit      <- 0
    }
  }
  return(syntax)
}

final.fit <- function(setup.out,
                      fixfitind.out,
                      data.file,
                      plot,
                      agg,
                      k){
  
  varnames   = setup.out$varnames
  lvarnames  = setup.out$lvarnames
  syntax     = fixfitind.out
  trackparts = setup.out$trackparts
  rois       = setup.out$rois
  plot.names = setup.out$plot.names
  x          = setup.out$x
  y          = setup.out$y
  betas      = setup.out$betas
  SEs        = setup.out$SEs
  plots      = setup.out$plots
  fitted     = setup.out$fitted
  
  fit <- fit.model(varnames  = varnames,
                   syntax    = syntax,
                   data.file = data.file)
  
  check.npd            <- any(grepl("error",class(fit))==TRUE)
  check.not.identified <- sum(lavInspect(fit,"se")$beta,na.rm=TRUE)==0
  
  if (check.npd == FALSE & check.not.identified==FALSE) { 
    singular            <- tryCatch(modindices(fit),error=function(e) e)
    check.singular      <- any(grepl("singular",singular)==TRUE)
    converge            <- lavInspect(fit, "converged")
  } else {
    check.singular <- TRUE
    converge       <- FALSE
  }
  if (converge==TRUE) last.converge <- TRUE 
  
  
  ## update code here to remove last element
  if (converge==FALSE) {
    last.converge <- FALSE
    syntax <- unlist(strsplit(syntax, "[\n]"))
    syntax <- syntax[-length(syntax)]
    syntax <- paste(syntax,sep="",collapse="\n")
    fit    <- fit.model(varnames  = varnames,  
                        syntax    = syntax,
                        data.file = data.file)
    
    check.npd            <- any(grepl("error",class(fit))==TRUE)
    check.not.identified <- sum(lavInspect(fit,"se")$beta,na.rm=TRUE)==0
    
    if (check.npd == FALSE & check.not.identified==FALSE) { 
      singular            <- tryCatch(modindices(fit),error=function(e) e)
      check.singular      <- any(grepl("singular",singular)==TRUE)
      converge            <- lavInspect(fit, "converged")
    } else {
      check.singular <- TRUE
      converge       <- FALSE
    }
  }
  
  ind.fit <- matrix(NA,nrow=1,ncol=9)
  if (converge==TRUE & check.singular==FALSE) {
    indices <- fitMeasures(fit,c("chisq","df","pvalue","rmsea","srmr",
                                 "nnfi","cfi"))  
    rmsea  <- indices[4]
    srmr   <- indices[5]
    cfi    <- indices[6]
    nnfi   <- indices[7]
    chisq  <- indices[1]
    df     <- indices[2]
    pval   <- indices[3]
    rmseaE <- ifelse(rmsea<.05, 1, 0)
    srmrE  <- ifelse(srmr <.05, 1, 0)
    cfiE   <- ifelse(cfi  >.95, 1, 0)
    nnfiE  <- ifelse(nnfi >.95, 1, 0)
    # insert indfit
    ind.fit[1,2] <- round(chisq,digits=4)
    ind.fit[1,3] <- df
    ind.fit[1,4] <- round(pval, digits=4)
    ind.fit[1,5] <- round(rmsea,digits=4)
    ind.fit[1,6] <- round(srmr, digits=4)
    ind.fit[1,7] <- round(nnfi, digits=4)
    ind.fit[1,8] <- round(cfi,  digits=4)
    if (last.converge==FALSE) {ind.fit[1,9] <- "last known convergence"
    } else {ind.fit[1,9] <- "converged normally"}
    
    indlist              <- subset(standardizedSolution(fit),op=="~")
    indlist$param        <- paste(indlist$lhs, indlist$op, indlist$rhs, sep="")
    indlist              <- indlist[complete.cases(indlist),]
    indlist              <- as.data.frame(indlist)
    indsubject           <- matrix(k,nrow=(nrow(indlist)), ncol=1)
    colnames(indsubject) <- c("subject")
    indlist <- cbind(indsubject, indlist)
    #creating individual-level beta matrices
    individual.paths     <- matrix(0,nrow=(rois*2), ncol=(rois*2))
    individual.SEs       <- matrix(0,nrow=(rois*2), ncol=(rois*2))
    indlist$row          <- substring(indlist$lhs,4)
    indlist$col          <- substring(indlist$rhs,4)
    indrows              <- indlist$row
    indcols              <- indlist$col
    ## multiple steps to recode column names to appropriate numerics
    
    indcols     <- as.numeric(recoderFunc(indcols,y,x))
    indrows     <- as.numeric(recoderFunc(indrows,y,x))
    
    indbetas    <- as.numeric(as.character(indlist$est.std))
    indSEs      <- as.numeric(as.character(indlist$se))
    indelements <- nrow(indlist)
    
    for (s in 1:indelements){
      individual.paths[indrows[s], indcols[s]] <- indbetas[s]
    }
    
    for (q in 1:indelements){
      individual.SEs[indrows[q], indcols[q]] <- indSEs[q]
    }
    
    
    if (agg==TRUE) {
      write.csv(individual.paths.lag, file=file.path(betas,"all_lagged.csv"), row.names=TRUE)
      write.csv(individual.paths.con, file=file.path(betas,"all_contemp.csv"), row.names=TRUE)
      write.csv(individual.SEs.lag, file=file.path(SEs,"all_lagged.csv"), row.names=TRUE)
      write.csv(individual.SEs.con, file=file.path(SEs,"all_contemp.csv"), row.names=TRUE) 
    }
    
    individual.paths[is.na(individual.paths)] <- 0 
    individual.paths.lag <- individual.paths[(rois+1):(rois*2),(1:rois)]
    individual.paths.con <- individual.paths[(rois+1):(rois*2),(rois+1):(rois*2)]
    
    rownames(individual.paths.lag) <- lvarnames[1:rois]
    colnames(individual.paths.lag) <- lvarnames[1:rois]
    rownames(individual.paths.con) <- lvarnames[(rois+1):(rois*2)]
    colnames(individual.paths.con) <- lvarnames[(rois+1):(rois*2)]
    
    individual.SEs[is.na(individual.SEs)] <- 0 
    individual.SEs.lag <- individual.SEs[(rois+1):(rois*2),(1:rois)]
    individual.SEs.con <- individual.SEs[(rois+1):(rois*2),(rois+1):(rois*2)]
    rownames(individual.SEs.lag) <- lvarnames[1:rois]
    colnames(individual.SEs.lag) <- lvarnames[1:rois]
    rownames(individual.SEs.con) <- lvarnames[(rois+1):(rois*2)]
    colnames(individual.SEs.con) <- lvarnames[(rois+1):(rois*2)]
    
    if (agg==TRUE) {
      write.csv(individual.paths.lag, file=file.path(betas,"all_lagged.csv"), row.names=TRUE)
      write.csv(individual.paths.con, file=file.path(betas,"all_contemp.csv"), row.names=TRUE)
      write.csv(individual.SEs.lag, file=file.path(SEs,"all_lagged.csv"), row.names=TRUE)
      write.csv(individual.SEs.con, file=file.path(SEs,"all_contemp.csv"), row.names=TRUE) 
    } else {
      
      write.csv(individual.paths.lag, file=file.path(betas,paste(trackparts[k,2],"lagged.csv",sep="_")), row.names=TRUE)
      write.csv(individual.paths.con, file=file.path(betas,paste(trackparts[k,2],"contemp.csv",sep="_")), row.names=TRUE)   
      write.csv(individual.SEs.lag, file=file.path(SEs,paste(trackparts[k,2],"lagged.csv",sep="_")), row.names=TRUE)
      write.csv(individual.SEs.con, file=file.path(SEs,paste(trackparts[k,2],"contemp.csv",sep="_")), row.names=TRUE)
      
    }
    
    
    individual.paths.dich <- matrix(0,nrow=(rois*2), ncol=(rois*2))
    indbetas.dich         <- ifelse(indbetas!=0,1,0)
    for (r in 1:indelements){
      individual.paths.dich[indrows[r], indcols[r]] <- indbetas.dich[r]
    }
    inddichoutput                <- paste (fitted,"/",trackparts[k,2],".txt",sep="")
    if (agg==TRUE) inddichoutput <- paste (fitted,"/all.csv",sep="")
    individual.paths.dich[is.na(individual.paths.dich)] <- 0 
    write.table(individual.paths.dich, 
                file      = inddichoutput, 
                sep       = ",", 
                row.names = FALSE, 
                col.names = FALSE)
    
    if (plot==TRUE){
      individual.paths.t     <- t(individual.paths)
      Lagged                 <- individual.paths.t[1:(rois),(rois+1):(rois*2)]
      Contemporaneous        <- individual.paths.t[(rois+1):(rois*2),(rois+1):(rois*2)]
      eLagged                <- W2E(Lagged)
      eContemporaneous       <- W2E(Contemporaneous)
      isLagged               <- c(rep(TRUE,nrow(eLagged)), rep(FALSE,nrow(eContemporaneous)))
      plotind                <- paste(plots,"/",trackparts[k,2],".pdf",sep="")
      if (agg==TRUE) plotind <- file.path(plots,"all.pdf")
      pdf(plotind)
      tryCatch(qgraph(rbind(eLagged,eContemporaneous),
                      layout              = "circle", 
                      lty                 = ifelse(isLagged,2, 1),
                      edge.labels         = F, 
                      curve               = TRUE, 
                      fade                = FALSE,
                      posCol              = "red",
                      negCol              = "blue", 
                      labels              = plot.names,
                      label.cex           = 3, 
                      edge.label.cex      = 1.5,
                      edge.label.position = .3),error=function(e) e)
      dev.off()
    } 
  }
  if (converge==FALSE) {
    ind.fit[1,9] <- "nonconvergence"
    indlist      <- data.frame()
  }
  if (check.singular==TRUE) {
    ind.fit[1,9] <- "computationally singular"
    indlist      <- data.frame()
  }
  list <- list("ind.fit"      = ind.fit,
               "ind.elements" = indlist,
               "syntax"       = syntax)
  return(list)
}


indsem.internal <- function(setup.out,
                            subgroup,
                            evalbetassub.out,
                            evalbetas.out){
  
  subjects         = setup.out$subjects
  varnames         = setup.out$varnames
  trackparts       = setup.out$trackparts
  header           = setup.out$header
  data             = setup.out$data
  sep              = setup.out$sep
  plot             = setup.out$plot
  
  all.elements      <- data.frame()
  all.fit           <- matrix(NA, nrow=subjects, ncol=9)
  all.syntax        <- matrix(NA, nrow=subjects,ncol=4)
  colnames(all.fit) <- c("subject", "chisq", "df", "pval", 
                         "rmsea", "srmr", "nnfi", "cfi", "status")
  
  files            <- list.files(data, full.names=TRUE)
  
  all.diff.subgroups <- length(unique(evalbetassub.out$syntax.sub[,1]))==subjects
  
  if (subgroup==TRUE & all.diff.subgroups==FALSE){
    all <- unlist(strsplit(evalbetassub.out$final.syntax[,2][1],"\n"))
    sub <- unlist(strsplit(evalbetassub.out$sub.paths[,2][1],","))
    all <- head(all,-length(sub))
    all <- paste(all,sep="",collapse="\n")
    
    sub.info <- as.data.frame(evalbetassub.out$syntax.sub,stringsAsFactors=FALSE)
    sorted   <- sub.info[order(sub.info$files),]
  }
  
  if (subgroup==TRUE & all.diff.subgroups==TRUE){
    all      <- evalbetassub.out$final.syntax[,2][1]
    sub.info <- as.data.frame(evalbetassub.out$syntax.sub,stringsAsFactors=FALSE)
    sorted   <- sub.info[order(sub.info$files),]
  }
  
  if (subgroup==FALSE){
    all <- evalbetas.out$syntax
  }
  
  for (k in 1:subjects) {
    
    writeLines(paste("individual-level search, subject", k))
    
    if (subgroup==TRUE) {
      syntax    <- sorted[k,3]
      
      data.file <- read.data(file   = sorted[k,2],
                             sep    = sep,
                             header = header)
      
      if (is.na(sorted[k,1])==T) syntax <- all
      
    } else {
      
      data.file <- read.data(file   = files[k],
                             sep    = sep,
                             header = header)
      syntax    <- all 
    }
    
    colnames(data.file) <- c(varnames)
    
    addind.out <- addind(done            = 0,
                         evaluate        = 1,
                         syntax          = syntax,
                         data.file       = data.file,
                         setup.out       = setup.out)
    
    evalind.out <- evalind(addind.out = addind.out,
                           setup.out  = setup.out,
                           data.file  = data.file)
    
    fixfitind.out <- fixfitind(setup.out   = setup.out,
                               evalind.out = evalind.out,
                               data.file   = data.file)
    
    final.fit.out <- final.fit(setup.out     = setup.out,
                               fixfitind.out = fixfitind.out,
                               data.file     = data.file,
                               k             = k,
                               agg           = FALSE,
                               plot          = plot)
    
    all.elements    <- rbind(all.elements,final.fit.out$ind.elements)
    all.fit[k,]     <- as.matrix(final.fit.out$ind.fit)
    all.syntax[k,3] <- as.matrix(final.fit.out$syntax) 
    all.syntax[k,1] <- k
  }
  
  all.fit[,1]    <- trackparts[,2]
  all.syntax[,2] <- files
  all.syntax[,4] <- all
  colnames(all.syntax) <- c("subject","files","syntax.ind","syntax.group")
  
  if (subgroup==TRUE){
    colnames(sorted)[3] <- c("syntax.sub")    
    #add code here to arrange group, subgroup, and individual-level paths
    all.syntax.sub <- merge(all.syntax,sorted,by="files")
    all.fit           <- cbind(all.fit, as.numeric(as.character(sorted[,1])))
    colnames(all.fit) <- c("subject", "chisq", "df", "pval", 
                           "rmsea", "srmr", "nnfi", "cfi", "status","subgroup")
    
  } else all.syntax.sub <- all.syntax
  
  list <- list("all.elements" = all.elements,
               "all.fit"      = all.fit,
               "all.syntax"   = all.syntax.sub,
               "all.diff.sub" = all.diff.subgroups)
  return(list)
}

wrapup <- function(indsem.internal.out,
                   setup.out,
                   agg,
                   subgroup){
  
  all.elements = indsem.internal.out$all.elements
  all.fit      = indsem.internal.out$all.fit
  all.syntax   = indsem.internal.out$all.syntax
  all.diff.sub = indsem.internal.out$all.diff.sub
  rois         = setup.out$rois
  out          = setup.out$out
  plot.names   = setup.out$plot.names
  subjects     = setup.out$subjects
  plots        = setup.out$plots
  x            = setup.out$x
  y            = setup.out$y
  header       = setup.out$header
  varnames     = setup.out$varnames
  lvarnames    = setup.out$lvarnames
  subgroup.dir = setup.out$subgroup.dir
  plot         = setup.out$plot
  
  present = NULL
  sig     = NULL
  est.std = NULL
  param   = NULL
  
  if (subgroup == TRUE & all.diff.sub == FALSE){
    all.merged.sub.final <- data.frame() 
    sub.final            <- character()
    sub.ind.paths        <- character()
    all.merged           <- merge(all.elements,all.syntax,by="subject")
    n.subgroups          <- length(unique(all.merged$membership))
    
    for (p in 1:n.subgroups){
      all.merged.sub         <- subset(all.merged,membership==p)
      ## add code to add column for "group","subgroup","individual"
      sub.subjects           <- length(unique(all.merged.sub$files))
      all.merged.sub$est.std <- as.numeric(as.character(all.merged.sub$est.std))
      all.merged.sub         <- transform(all.merged.sub,mean.beta = (ave(est.std, param, FUN=sum))/sub.subjects)
      all.merged.sub$sig     <- ifelse(abs(all.merged.sub$z)>1.96, 1, 0)
      all.merged.sub$present <- ifelse(all.merged.sub$sig<2,1,0)
      all.merged.sub         <- transform(all.merged.sub, sum.sig = ave(sig, param, FUN=sum))
      all.merged.sub         <- transform(all.merged.sub, count   = ave(present, param, FUN=sum))
      
      if (sub.subjects != 1){
        ## using syntax to identify group, subgroup, individual level paths
        group    <- unlist(strsplit(as.vector(unique(all.merged.sub$syntax.group)),"\n"))
        sub      <- unlist(strsplit(as.vector(unique(all.merged.sub$syntax.sub)),"\n"))
        ind      <- unlist(strsplit(as.vector(unique(all.merged.sub$syntax.ind)),"\n"))
        ind      <- unique(subset(ind, !ind %in% sub))
        sub      <- subset(sub, !sub %in% group)  
      } else {
        # if there's only one person in the subgroup, then the ind paths are the sub paths
        # only used to determine whether or not a subgroup path is promoted to a group path
        group    <- unlist(strsplit(as.vector(unique(all.merged.sub$syntax.group)),"\n"))
        ind      <- unlist(strsplit(as.vector(unique(all.merged.sub$syntax.ind)),"\n"))
        ind      <- subset(ind, !ind %in% group)
        sub      <- ind
        sub.ind.paths <- append(sub.ind.paths,sub)
      }
      
      label  <- c(rep("group",length(group)),rep("subgroup",length(sub)),rep("ind",length(ind)))
      color  <- c(rep("black",length(group)),rep("green3",length(sub)),rep("gray50",length(ind)))
      param  <- c(group,sub,ind)
      levels <- data.frame(param,label,color)
      
      all.merged.sub       <- merge(all.merged.sub,levels,by="param")
      all.merged.sub.final <- rbind(all.merged.sub.final,all.merged.sub)
      sub.final            <- append(sub.final,sub)
    }
    
    ## this portion checks to see if there is a subgroup path that exists for each subgroup
    ## if it does, it bumps it up to the group level
    a            <- table(sub.final)
    sub.to.group <- names(a[a==n.subgroups])
    sub.stay.sub <- names(a[a!=n.subgroups])
    all.merged.sub.final$label[all.merged.sub.final$param==sub.to.group] <- "group"
    all.merged.sub.final$color[all.merged.sub.final$param==sub.to.group] <- "black"
    
    ## code to create subgroup-level plots and matrices
    for (p in 1:n.subgroups){
      all.merged.sub   <- subset(all.merged.sub.final,membership==p) 
      sub.subjects     <- length(unique(all.merged.sub$files))
      sub.paths        <- matrix(0,nrow=(rois*2), ncol=(rois*2))
      sub.colors       <- matrix(NA,nrow=(rois*2), ncol=(rois*2))
      elements         <- nrow(all.merged.sub)
      rows             <- substring(all.merged.sub$lhs, 4)
      cols             <- substring(all.merged.sub$rhs, 4)
      cols             <- as.numeric(recoderFunc(cols,y,x))
      rows             <- as.numeric(recoderFunc(rows,y,x))
      counts           <- as.numeric(all.merged.sub$count)
      beta.weights     <- as.numeric(all.merged.sub$mean.beta)
      colors           <- as.character(all.merged.sub$color)
      for (m in 1:elements){
        sub.paths [rows[m], cols[m]] <- counts[m]
        sub.colors[rows[m], cols[m]] <- colors[m]
      }
      sub.paths[is.na(sub.paths)] <- 0
      sub.paths.lag <- sub.paths[(rois+1):(rois*2),(1:rois)]
      sub.paths.con <- sub.paths[(rois+1):(rois*2),(rois+1):(rois*2)]
      rownames(sub.paths.lag) <- lvarnames[1:rois]
      colnames(sub.paths.lag) <- lvarnames[1:rois]
      rownames(sub.paths.con) <- lvarnames[(rois+1):(rois*2)]
      colnames(sub.paths.con) <- lvarnames[(rois+1):(rois*2)]
      
      if (sub.subjects != 1)  {
        write.csv(sub.paths.lag,file.path(subgroup.dir,paste("subgroup",p,"_lagged.csv",sep="")))
        write.csv(sub.paths.con,file.path(subgroup.dir,paste("subgroup",p,"_contemp.csv",sep="")))
      }
      
      if (plot==TRUE & sub.subjects != 1){
        sub.paths.t      <- t(sub.paths)
        sub.paths.t      <- sub.paths.t/sub.subjects
        Lagged           <- sub.paths.t[1:(rois),(rois+1):(rois*2)]
        Contemporaneous  <- sub.paths.t[(rois+1):(rois*2),(rois+1):(rois*2)]
        eLagged          <- W2E(Lagged)
        eContemporaneous <- W2E(Contemporaneous)
        sub.colors.t     <- t(sub.colors)
        colorLagged           <- sub.colors.t[1:(rois),(rois+1):(rois*2)]
        colorContemporaneous  <- sub.colors.t[(rois+1):(rois*2),(rois+1):(rois*2)]  
        color.list       <- c(colorLagged,colorContemporaneous)
        color.list       <- color.list[!is.na(color.list)]
        isLagged         <- c(rep(TRUE,nrow(eLagged)), rep(FALSE,nrow(eContemporaneous)))
        plotsub          <- file.path(subgroup.dir,paste("subgroup",p,"_plot.pdf",sep=""))
        pdf(plotsub)
        tryCatch(qgraph(rbind(eLagged,eContemporaneous),
                        layout              = "circle", 
                        lty                 = ifelse(isLagged,2, 1),
                        edge.labels         = FALSE, 
                        curve               = TRUE,
                        labels              = plot.names,
                        edge.color          = color.list,
                        fade                = FALSE,
                        label.cex           = 3,  
                        edge.label.cex      = 1.5, 
                        edge.label.position = .3),error=function(e) e)
        dev.off()
      }
    }
    
    ## quick fix to make sure subgroup paths are represented as subgroup,
    ## not individual. this could happen if it's subgroup for some and individual for others
    all.elements.final       <- all.merged.sub.final
    all.elements.final$label <- as.character(all.elements.final$label)
    ## removes individual-level paths from this list for subgroup size 1
    ## this way, no paths are marked as subgroup-level that only existed for one person
    ## in the final summary graphic
    sub.stay.sub <- subset(sub.stay.sub,!sub.stay.sub %in% sub.ind.paths)
    #
    all.elements.final$label[all.elements.final$param %in% sub.stay.sub] <- "sub"
    all.elements.final$color[all.elements.final$param %in% sub.stay.sub] <- "green3"
    ## adjust count of paths 
    all.elements.final <- transform(all.elements.final, count   = ave(present, param, FUN=sum))
    all.elements.final <- transform(all.elements.final, sum.sig = ave(sig, param, FUN=sum))
    all.elements.final <- subset(all.elements.final, !duplicated(param))

    a1 <- aggregate(present ~ label + param + membership,data=all.merged.sub.final,sum)
    a1$label <- ifelse(a1$label=='subgroup',
                          paste(a1$label,a1$membership,sep=""),as.character(a1$label))
    a2 <- aggregate(present~param+label,data=a1,sum)
    a2 <- a2[order(-a2$present,a2$label),]
    ## now transpose object to have column per label
    a2 <- reshape(a2,timevar="label",idvar="param",direction="wide")
    a2[is.na(a2)] <- 0
    
    
  } else {
    all.elements.final         <- as.data.frame(all.elements)
    all.elements.final$est.std <- as.numeric(as.character(all.elements.final$est.std))
    if (agg==FALSE) {all.elements.final <- transform(all.elements, 
                                                     mean.beta = (ave(est.std, param, FUN=sum))/subjects)}
    all.elements.final$sig     <- ifelse(abs(all.elements.final$z)>1.96, 1, 0)
    all.elements.final$present <- ifelse(all.elements.final$sig<2,1,0)
    all.elements.final         <- transform(all.elements.final, sum.sig = ave(sig, param, FUN=sum))
    all.elements.final         <- transform(all.elements.final, count = ave(present, param, FUN=sum))
    all.elements.final         <- subset(all.elements.final, !duplicated(param))
  }
  

  
  if (agg==FALSE){
    final.paths      <- matrix(0,nrow=(rois*2), ncol=(rois*2))
    final.betas      <- matrix(0,nrow=(rois*2), ncol=(rois*2))
    final.paths.dich <- matrix(0,nrow=(rois*2), ncol=(rois*2))
    final.colors     <- matrix(NA,nrow=(rois*2), ncol=(rois*2))
    elements         <- nrow(all.elements.final)
    rows             <- substring(all.elements.final$lhs, 4)
    cols             <- substring(all.elements.final$rhs, 4)
    cols             <- as.numeric(recoderFunc(cols,y,x))
    rows             <- as.numeric(recoderFunc(rows,y,x))
    counts           <- as.numeric(all.elements.final$count)
    beta.weights     <- as.numeric(all.elements.final$mean.beta)
    if (subgroup==TRUE) colors <- as.character(all.elements.final$color)
    for (m in 1:elements){
      final.paths[rows[m], cols[m]] <- counts[m]
      if (subgroup==TRUE) final.colors[rows[m], cols[m]] <- colors[m]
    }
    final.paths[is.na(final.paths)] <- 0
    final.paths.lag <- final.paths[(rois+1):(rois*2),(1:rois)]
    final.paths.con <- final.paths[(rois+1):(rois*2),(rois+1):(rois*2)]
    rownames(final.paths.lag) <- lvarnames[1:rois]
    colnames(final.paths.lag) <- lvarnames[1:rois]
    rownames(final.paths.con) <- lvarnames[(rois+1):(rois*2)]
    colnames(final.paths.con) <- lvarnames[(rois+1):(rois*2)]
    write.csv(final.paths.lag, file=file.path(out,"finalpaths_lagged.csv"), row.names=TRUE)
    write.csv(final.paths.con, file=file.path(out,"finalpaths_contemp.csv"), row.names=TRUE)
    
    if (subgroup==FALSE | all.diff.sub==TRUE){
      all.elements.summary <- all.elements.final[c("param","mean.beta","sum.sig","count")]
    } else {
#       all.elements.summary <- all.elements.final[c("param","mean.beta","sum.sig","count","label")]
#       all.elements.summary  <- all.elements.summary[order(-all.elements.summary$count),]
      all.elements.summary  <- a2
    }
    
    if (plot==TRUE){
      final.paths.t    <- t(final.paths)
      final.paths.t    <- final.paths.t/subjects
      Lagged           <- final.paths.t[1:(rois),(rois+1):(rois*2)]
      Contemporaneous  <- final.paths.t[(rois+1):(rois*2),(rois+1):(rois*2)]
      eLagged          <- W2E(Lagged)
      eContemporaneous <- W2E(Contemporaneous)
      if (subgroup==TRUE){
        final.colors.t     <- t(final.colors)
        colorLagged           <- final.colors.t[1:(rois),(rois+1):(rois*2)]
        colorContemporaneous  <- final.colors.t[(rois+1):(rois*2),(rois+1):(rois*2)]  
        color.list       <- c(colorLagged,colorContemporaneous)
        color.list       <- color.list[!is.na(color.list)]
      }
      isLagged         <- c(rep(TRUE,nrow(eLagged)), rep(FALSE,nrow(eContemporaneous)))
      plotfinal        <- file.path(out,"/finalpaths.pdf")
      if (subgroup==TRUE){edge.color=color.list}else{edge.color=NULL}
      pdf(plotfinal)
      ## insert error wrapper here
      tryCatch(qgraph(rbind(eLagged,eContemporaneous),
                      layout              = "circle", 
                      lty                 = ifelse(isLagged,2, 1),
                      edge.labels         = FALSE, 
                      curve               = TRUE,
                      labels              = plot.names, 
                      edge.color          = edge.color,
                      fade                = FALSE,
                      label.cex           = 3,  
                      edge.label.cex      = 1.5, 
                      edge.label.position = .3),error=function(e) e)
      dev.off()
    }
  }
  



  all.elements   <- all.elements[c("subject","param","est.std","se","pvalue","z")]
  colnames(setup.out$trackparts) <- c("subject","ID")
  all.elements   <- merge(setup.out$trackparts,all.elements,by.x="subject",by.y="subject")[,-1]
  all.elements   <- all.elements[order(all.elements$ID),]


  if (agg==TRUE) {all.elements[,1] <- "all"; all.elements.summary <- data.frame()}
  # insert code here to sub back variable names into all.elements summary and all.elements
  if (header==TRUE){
    for (v in 1:length(varnames)){
      all.elements$param <- gsub(lvarnames[v],varnames[v],all.elements$param)
      if (agg==FALSE) all.elements.summary$param <- gsub(lvarnames[v],varnames[v],all.elements.summary$param)
    }
  }
  row.names(all.elements.summary) <- NULL
  row.names(all.elements)         <- NULL
  all.fit                         <- as.data.frame(all.fit)
  write.csv(all.fit,file.path(out,"allfit.csv"),row.names=FALSE)
  write.csv(all.elements,file.path(out,"all.elements.csv"),row.names=FALSE)
  write.csv(all.elements.summary,file.path(out,"all.elements.summary.csv"),row.names=FALSE)
}

