gimme <- function(data,
                  sep,
                  header,
                  out,
                  plot = TRUE,
                  ar   = FALSE,
                  subgroup,
                  deconvolve_hrf = FALSE,
                  control=list(deconvolve_method="bush")){
  
  setup.out         <- setup(data   = data, 
                             sep    = sep,
                             header = header, 
                             out    = out,
                             plot   = plot,
                             ar     = ar,
                             deconvolve_hrf = deconvolve_hrf,
                             control = control)
  
  miSEM.out         <- miSEM(setup.out            = setup.out,
                             previous.out         = setup.out,
                             subgroup.step        = FALSE,
                             subjects             = NULL,
                             files                = NULL,
                             syntax               = NULL,
                             subgroup.paths       = NULL,
                             second.round         = FALSE)
  
  evalbetas.out <- evalbetas(setup.out            = setup.out,
                             previous.out         = miSEM.out,
                             subgroup.step        = FALSE,
                             s                    = NULL,
                             subjects             = NULL,
                             files                = NULL,
                             post.sub.prune       = FALSE,
                             evalbetas.out        = NULL)
  
  ## begin subgroup steps ##
  if (subgroup==TRUE){
    
    subsetup.out         <- subsetup(setup.out    = setup.out,
                                     previous.out = evalbetas.out)
    
    miSEMsub.out         <- miSEMsub(subsetup.out = subsetup.out,
                                     setup.out    = setup.out,
                                     previous.out = evalbetas.out,
                                     second.round = FALSE,
                                     evalbetassub.out = NULL)
    
    evalbetassub.out <- evalbetassub(subsetup.out  = subsetup.out,
                                     previous.out  = miSEMsub.out,
                                     setup.out     = setup.out,
                                     evalbetas.out = evalbetas.out)
    
    prune.group.post.out <- evalbetas(setup.out      = setup.out,
                                      previous.out   = evalbetassub.out,
                                      subgroup.step  = FALSE,
                                      s              = NULL,
                                      subjects       = NULL,
                                      files          = NULL,
                                      post.sub.prune = TRUE,
                                      evalbetas.out  = evalbetas.out)
    
    ## this only runs if a group path was removed in the post-subgrouping pruning
    if (evalbetas.out$count.group.paths != prune.group.post.out$count.group.paths){
      
      miSEMsub.round2.out <- miSEMsub(subsetup.out = subsetup.out,
                                      setup.out    = setup.out,
                                      previous.out = prune.group.post.out,
                                      second.round = TRUE,
                                      evalbetassub.out = evalbetassub.out)
      
      check.again <- any(miSEMsub.round2.out$unique.syntax[,4] != evalbetassub.out$sub.paths[,2])
        
      if (check.again == TRUE){    
        evalbetassub.round2.out <- evalbetassub(subsetup.out  = subsetup.out,
                                                previous.out  = miSEMsub.round2.out,
                                                setup.out     = setup.out,
                                                evalbetas.out = evalbetas.out)
        
        evalbetassub.out <- evalbetassub.round2.out
      } else {
        # if evalbetas isn't run again, grabs most recent syntax from miSEMsub.round2
        # for use in indsem
        test <- merge(evalbetassub.out$syntax.sub,miSEMsub.round2.out$unique.syntax,by.x="membership",by.y="V1")
        test <- as.matrix(test[-c(3,5,6)])
        colnames(test) <- c("membership","files","syntax")
        evalbetassub.out$syntax.sub   <- test
        evalbetassub.out$final.syntax <- unique(test[,c(1,3)])
      }
    }
    
    sub.paths <- evalbetassub.out$sub.paths
    colnames(sub.paths) <- c("subgroup","subgroup-level paths")
    write.csv(sub.paths,file.path(setup.out$out,"subgroup_paths.csv"),row.names=FALSE)
    ## end subgroup steps ##
  
  }
  
  if (subgroup==FALSE) evalbetassub.out <- NULL
  
  indsem.internal.out <- indsem.internal(setup.out        = setup.out,
                                         subgroup         = subgroup,
                                         evalbetassub.out = evalbetassub.out,
                                         evalbetas.out    = evalbetas.out)

  wrapup.out <- wrapup(indsem.internal.out = indsem.internal.out,
                       setup.out           = setup.out,
                       agg                 = FALSE,
                       subgroup            = subgroup)
  
  print.gimme(x=subsetup.out,
              y=subgroup,
              z=setup.out)
}

print.gimme <- function(x,y,z){
  writeLines("gimme finished running normally")
  writeLines(paste("output is stored in ", z$out))
  if (y == TRUE) {
    writeLines(paste("Number of subgroups =", x$n.subgroups))
    writeLines(paste("Modularity =", round(x$modularity,digits=5)))
  }
}

