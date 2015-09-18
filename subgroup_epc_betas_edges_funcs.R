subsetup <- function (setup.out,
                      previous.out) {
  
  subjects        = setup.out$subjects
  varnames        = setup.out$varnames
  vars            = setup.out$vars
  candidate.paths = setup.out$candidate.paths

  
  syntax          = previous.out$syntax
  syntax.paths    = unlist(strsplit(syntax,"\n"))
  group.paths     = previous.out$count.group.paths
  sep             <- setup.out$sep
  header          <- setup.out$header
  data            <- setup.out$data
  out             <- setup.out$out
  cutoffind       <- setup.out$cutoffind
  cutoffgroup     <- qchisq(1-.05/subjects,1)
  cutoffmi        <- qchisq(1-(.05/((vars*(vars-1)/2)*subjects)),1)
  
  files           <- list.files(data, full.names=TRUE)
  trackparts      <- setup.out$trackparts
  mi.matrix    <- matrix(0,ncol=subjects,nrow=((vars*(vars-1)/2)))
  epc.matrix   <- matrix(0,ncol=subjects,nrow=((vars*(vars-1)/2)))
  p.matrix     <- matrix(0,ncol=subjects,nrow=(group.paths+(vars/2)))
  beta.matrix  <- matrix(0,ncol=subjects,nrow=(group.paths+(vars/2)))
  ## create matrix of MIs individual by individual after evalbetas
 
  for (k in 1:subjects)
  {
    data.file           <- read.data(file   = files[k],
                                     sep    = sep,
                                     header = header)
    time                <- nrow(data.file)
    colnames(data.file) <- c(varnames)
    fit <- fit.model(varnames  = varnames, 
                     syntax    = syntax,
                     data.file = data.file)
    
    check.npd           <- any(grepl("error",class(fit))==TRUE)
    
    if (check.npd == FALSE) { 
      singular            <- tryCatch(modindices(fit),error=function(e) e)
      check.singular      <- any(grepl("singular",singular)==TRUE)
      converge            <- lavInspect(fit, "converged")
    } else {
      check.singular <- TRUE
      converge       <- FALSE
    }
    
    # printing replication
    writeLines(paste("subgroup search, subject", k))
    
    if (check.singular==FALSE & converge==TRUE) {
      mi           <- as.matrix(singular[singular$op == "~",])[,c(1:5)]
      mi           <- as.data.frame(mi)
      mi$param     <- paste(mi$lhs, mi$op, mi$rhs, sep="")
      mi           <- mi[-c(1:3)]
      mi           <- subset(mi,param %in% candidate.paths)

      mi$mi        <- as.numeric(as.character(mi$mi))
      mi$epc       <- as.numeric(as.character(mi$epc))
      
      mi$mi[mi$param %in% syntax.paths]  <- 0
      mi$epc[mi$param %in% syntax.paths] <- 0
      
      mi$thresh    <- ifelse(mi$mi>cutoffmi,1,0)
      beta         <- as.data.frame(subset(standardizedSolution(fit),op=="~"))
      beta$thresh  <- ifelse(beta$pvalue<(.05/subjects),1,0)
    }
    # if it doesn't converge or is computationally singular
    if (converge==FALSE | check.singular==TRUE | check.npd==TRUE) {
      mi.matrix[,k]     <- NA
      epc.matrix[,k]    <- NA
      p.matrix[,k]      <- NA
      beta.matrix[,k]   <- NA
    } else {
      mi.matrix[,k]     <- mi$thresh
      epc.matrix[,k]    <- mi$epc
      p.matrix[,k]      <- beta$thresh
      beta.matrix[,k]   <- beta$est.std
    }
  }
  
  all.val.matrix         <- rbind(epc.matrix,beta.matrix)
  all.p.matrix           <- rbind(mi.matrix,p.matrix)
  
  colnames(all.val.matrix) <- files
  colnames(all.p.matrix)   <- files
  
  all.val.matrix        <- all.val.matrix[,colSums(is.na(all.val.matrix)) != nrow(all.val.matrix)]
  all.p.matrix          <- all.p.matrix[,colSums(is.na(all.p.matrix)) != nrow(all.p.matrix)]
  names                 <- colnames(all.val.matrix)
  
  sim <- matrix(0,ncol=ncol(all.val.matrix),nrow=ncol(all.val.matrix))
  
  for (i in 1:ncol(all.p.matrix)){
    for (j in 1:ncol(all.p.matrix)){
      ind1 <- as.matrix(all.p.matrix[,i])
      ind2 <- as.matrix(all.p.matrix[,j])
      val1 <- as.matrix(all.val.matrix[,i])
      val2 <- as.matrix(all.val.matrix[,j])
      sim[i,j] <- sum(ind1 == 1 & ind2 == 1 & sign(val1)==sign(val2))
    }
  }
  
  colnames(sim) <- names
  sim           <- sim - min(sim)
  diag(sim)     <- 0
  
  g                   <- graph.adjacency(sim,mode="undirected")
  membership          <- walktrap.community(g,steps=4)$membership
  
  modstem           <- file.path(out,"modularity.txt")
  modularity.value  <- modularity(walktrap.community(g,steps=4))
  write.table(modularity.value,modstem,row.names=FALSE,col.names=FALSE)
  
  subgroup.membership <- cbind(names,membership)
  n.subgroups         <- length(table(subgroup.membership[,2]))
  files               <- as.matrix(files)
  subgroup.membership <- merge(files,subgroup.membership,
                               by.x="V1",by.y="names",all.x=TRUE)
  colnames(subgroup.membership) <- c("files","membership")
  list <- list("mi.matrix"           = mi.matrix,
               "cor.mi"              = sim,
               "subgroup.membership" = subgroup.membership,
               "n.subgroups"         = n.subgroups,
               "modularity"          = modularity.value)
  return(list)
}
