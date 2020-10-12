#' @keywords internal 
modmax <- function(x, m, sub_method){
      # sparsify m based on x
      diag(m) <- 0
      m[which(m <= stats::quantile(m[upper.tri(m, diag = FALSE)], x))] <- 0
      olm = m
      m = graph_from_adjacency_matrix(m)
      weights      <- E(m)$weight
        if (sub_method == "Walktrap")
          p = cluster_walktrap(m, weights = weights, steps = 4)
        if (sub_method == "Infomap")
          p = cluster_infomap(m, e.weights = weights)
        if (sub_method == "Edge Betweenness")
          p = cluster_edge_betweenness(m, weights = weights)
        if (sub_method == "Fast Greedy")
          p = cluster_fast_greedy(m, weights = weights)
        if (sub_method == "Label Prop")
          p = cluster_label_prop(m, weights = weights)
        if (sub_method == "Leading Eigen")
          p = cluster_leading_eigen(m, weights = weights)
        if (sub_method == "Louvain")
          p = cluster_louvain(m, weights = weights)
        if (sub_method == "Spinglass")
          p = cluster_spinglass(m, weights = weights)
      modval <- igraph::modularity(p)
      sk = ask = list()
      mem = p$membership
      dim = length(mem)

      for(i in 1:dim){
        sk[i] = 0
        ask[i] = 0
        for(j in 1:dim){
          if(olm[i,j] == 0){NULL}else{
            if(mem[i] == mem[j]){ask[[i]] = ask[[i]] + olm[i,j]}else{
              sk[[i]] = sk[[i]] + olm[i,j]
            }
          }
        }
        sk[[i]] = sk[[i]]/2
        ask[[i]] = sk[[i]] + ask[[i]]/2
      }
      sk = unlist(sk); ask = unlist(ask)

      temp.1 = cbind(mem, sk, ask)
      
      clist = matrix(NA, 1, max(mem))
      for(i in 1:max(temp.1[,1])){
        if(temp.1[mem == i,2] == 0 && temp.1[mem == i,3] == 0){
          clist[1,i] = 1
        }else{
          clist[1,i] = sum(temp.1[mem == i,2]) / 
            sum(temp.1[mem == i,3])  
        }
      }
      simplified = round(clist, 3)
      (conductance = 1 - sum(simplified)/length(simplified))
      return(-1*conductance)
}







