#' Perturb networks and evaluate subgroup structures.
#' @param mat A symmetric, weighted matrix or weighted graph
#' 
#' 

perturb <- function(mat, method){
  #library(mcclust) 
  #library(dils)
  #library(igraph)
  #library(ggplot2)
  
  if (is.igraph(mat)) 
    { g <- mat
    
    mat <- as.matrix( get.adjacency(g, attr = "weight"))
    
    } else
    g                   <- graph.adjacency(mat, mode = "undirected", weighted = TRUE)
  
  diag(mat)        <- 0
  
  truemembership          <- walktrap.community(g, steps = 4)$membership
 
  # now randomly perturb 
  n.elements <- length(mat[,1])*(length(mat[,1])-1)/2
  percent <-seq(from=0, to = n.elements, by=round(0.01*(n.elements))) #disrupt 1% at a time
  VI<-matrix(,nrow = 100, ncol = length(percent))
  ARI<-matrix(,nrow = 100,ncol = length(percent))
  modularity_value <- matrix(,nrow = 100, ncol = length(percent))
  dump <- matrix(1, length(mat[,1]), length(mat[,1]))
  diag(dump) <- 0
  eligable <- which(lower.tri(dump)!=0,arr.ind = T)#ensure that diagonal isn't considered, and all potential edges are

  VI[,1] <- 0 #when 0 edges are disrupted no variation of information
  ARI[,1]<- 1 #when 0 edges are disrupted ARI = 1
  
  for(k in 2:length(percent)) { 
    for(p in 1:100){ # for each degree of perturbation run 100 times
      toalter <- sample(1:length(eligable[,1]), percent[k])
      if (length(toalter)>1)
      {
      randomized <- sample(toalter)
      } else {
        randomized <- toalter
        }
      new.v <- mat
      for (l in 1:length(randomized))
      new.v[eligable[toalter[l],1],eligable[toalter[l],2]]<- new.v[[eligable[randomized[l],1],eligable[randomized[l],2]]] 
      # maintain symmetry:
      new.v[eligable[toalter[l],2],eligable[toalter[l],1]]<- new.v[[eligable[randomized[l],2],eligable[randomized[l],1]]] 
      new.g                  <- graph.adjacency(new.v, mode = "undirected", weighted = TRUE)
      membership       <- walktrap.community(new.g, steps = 4)$membership
      modularity_value[p,k] <- modularity(walktrap.community(new.g, steps = 4))
      ## ARI and VI for comparison of this.g and membership in original
      VI[p,k] <- vi.dist(membership, truemembership)
      ARI[p,k] <-arandi(membership, truemembership)
    } 
  }
  
  # plots of ARI and VI compared to original
  percentlab <- percent/n.elements
  meanVI <- colMeans(VI)
  meanARI <- colMeans(ARI)
  meanMod <- colMeans(modularity_value)
  
  # explicitly change 10 and 20 percent of community affiliations to add lines to graph
  comms <- unique(truemembership)
  lengths <- NULL
  for (p in 1:length(comms))
    lengths[p] <- length(which(truemembership == comms[p]))
  max <- which(lengths == max(lengths))
  perc10 <- round(.10*length(truemembership))
  perc20 <- round(.20*length(truemembership))
  tochange <- which(truemembership == comms[max])
  changed10 <- truemembership
  for (p in 1:perc10)
    changed10[tochange[p]] <- comms[sample(which(lengths == min(lengths)[1]), 1)]
  changed20 <- truemembership
  for (p in 1:perc20)
    changed20[tochange[p]] <- comms[sample(which(lengths == min(lengths)[1]), 1)]
  
  plotVI <- qplot(percentlab, meanVI, main = "Comparison of original result against perturbed graphs: VI", xlab = "Proportion Perturbed", ylab = "VI")
  plotVI <- plotVI + geom_hline(mapping = NULL, data = NULL, yintercept = c(vi.dist(changed10, truemembership), vi.dist(changed20, truemembership)))
  
  plotARI <-  qplot(percentlab, meanARI, main = "Comparison of original result against perturbed graphs: ARI", xlab = "Proportion Perturbed", ylab = "Mean ARI",
                    ymin = 0, ymax = 1) 
  plotARI <- plotARI + geom_hline(mapping = NULL, data = NULL, yintercept = c(arandi(changed10, truemembership), arandi(changed20, truemembership)))
  
  #perhaps plot violin plots at each X value rather than mean? need to scale differently (margins)
  
  res <- list(
    VI = VI, 
    ARI = ARI, 
    modularity_value = modularity_value,
    percent = percentlab,
    plotVI = plotVI,
    plotARI = plotARI
    )
  return(res)
}
