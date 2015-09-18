miSEMsub <- function (subsetup.out,
                      setup.out,
                      previous.out,
                      second.round,
                      evalbetassub.out) {
  
  if (second.round==FALSE) miSEMsub.out = NULL
  n.subgroups         = subsetup.out$n.subgroups
  subgroup.membership = subsetup.out$subgroup.membership
  varnames            = setup.out$varnames
  syntax              = previous.out$syntax
  count.group.paths   = previous.out$count.group.paths
  
  unique.syntax       <- matrix(seq(1:n.subgroups),ncol=4,nrow=n.subgroups)
  subgroup.membership <- as.data.frame(subgroup.membership)
  if (second.round==TRUE) {syntax.all <- cbind(subsetup.out$subgroup.membership$membership,
                                           previous.out$syntax)
                           syntax.all <- syntax.all[!duplicated(syntax.all),]
                           
  }
  
  for (s in 1:n.subgroups){
    if (second.round==TRUE){indx   <- syntax.all[,1]==s
                            syntax <- syntax.all[indx,2]
                            indx2  <- evalbetassub.out$sub.paths[,1]==s
                            subgroup.paths <- evalbetassub.out$sub.paths[indx2,2]
                            names(subgroup.paths) <- NULL
                            
    }
    sub.out.files     <- as.character(subset(subgroup.membership,membership==s)[,1])
    subjects          <- length(sub.out.files)
    if (subjects > 1){
      subgroup.miSEM <- miSEM(setup.out         = setup.out,
                              previous.out      = previous.out,
                              subgroup.step     = TRUE,
                              subjects          = subjects,
                              files             = sub.out.files,
                              syntax            = syntax,
                              subgroup.paths    = subgroup.paths,
                              second.round      = second.round)
      
      unique.syntax[s,2] <- subgroup.miSEM$syntax
      unique.syntax[s,3] <- subgroup.miSEM$count.subgroup.paths
      unique.syntax[s,4] <- paste(subgroup.miSEM$subgroup.paths,sep="",collapse=",")
      
    } else {unique.syntax[s,2] <- syntax
            unique.syntax[s,3] <- 0
            unique.syntax[s,4] <- ""}
  }
  list <- list("unique.syntax"     = unique.syntax,
               "count.group.paths" = count.group.paths)
  return(list)
}

evalbetassub <- function (subsetup.out,
                          previous.out,
                          setup.out,
                          evalbetas.out) {
  
  n.subgroups         = subsetup.out$n.subgroups
  subgroup.membership = subsetup.out$subgroup.membership
  ar.paths            = setup.out$ar.paths
  ar                  = setup.out$ar
  out                 = setup.out$out
  
  subgroup.membership <- as.data.frame(subgroup.membership)
  final.syntax        <- matrix(seq(1:n.subgroups),ncol=2,nrow=n.subgroups)
  sub.paths           <- matrix(seq(1:n.subgroups),ncol=2,nrow=n.subgroups)
  
  for (s in 1:n.subgroups){
    sub.out.files     <- as.character(subset(subgroup.membership,membership==s)[,1])
    subjects          <- length(sub.out.files)
    
    if (subjects > 1){
      evalbetas.out.sub <- evalbetas(setup.out      = setup.out,
                                     previous.out   = previous.out,
                                     subgroup.step  = TRUE,
                                     s              = s,
                                     subjects       = subjects,
                                     files          = sub.out.files,
                                     post.sub.prune = FALSE)
      
      final.syntax[s,2]      <- evalbetas.out.sub$syntax
      if (length(evalbetas.out.sub$subgroup.paths) !=0) {
        sub.paths[s,2]         <- evalbetas.out.sub$subgroup.paths
      } else sub.paths[s,2] <- NA
    } else { final.syntax[s,2] <- evalbetas.out$syntax
             sub.paths[s,2] <- NA
    }
  }
  
  colnames(final.syntax) <- c("membership","syntax")
  syntax.sub             <- merge(subgroup.membership,final.syntax,by="membership",all=TRUE)
  syntax.sub             <- syntax.sub[order(syntax.sub$files),]
  syntax.sub$files       <- as.character(syntax.sub$files)
  syntax.sub$syntax      <- as.character(syntax.sub$syntax)
  list <- list("syntax.sub"   = syntax.sub,
               "final.syntax" = final.syntax,
               "sub.paths"    = sub.paths)
  return(list)
}
