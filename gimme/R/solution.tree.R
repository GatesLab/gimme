#' solution.tree
#' @export solution.tree
solution.tree <- function(x, level = c("group", "individual"), cols = NULL, ids = "all", plot.tree = FALSE){
  
  tree <- gimme:::batch.create.tree(
    x$grp_hist, 
    x$ind_hist, 
    x$ind_fit, 
    x$dat$subgroup, 
    names(x$dat$ts_list),
    x$sub)
  
  if(is.null(cols)){
    cols <- c("stage")
  } else {
    cols <- c("stage", cols)
  }
  
  if(level == "group"){
    
    grp_tree <- tree[[1]]
    
    if(plot.tree){
      final.plot <- plot(grp_tree)
      return(final.plot)
    } else {
      do.call(print, c(grp_tree, cols))
    }
    
  } else {

    ind_tree <-  tree[[2]]
    
    
    if(length(ids == 1) && ids == "all"){
      
      to_print <- seq(1:length(x$dat$ts_list))
      
    } else {
      
      to_print <- which(names(x$dat$ts_list) %in% ids)
      
    }
    
    if(length(ids == 1) && ids != "all"){
      
      if(plot.tree){
        plot(ind_tree[[to_print]])
      } else {
        do.call(print, c(ind_tree, cols))
      }
      
    } else {
      
      for(i in to_print){
        cat(paste0(names(x$dat$ts_list)[i]),"\n\n")
        if(plot.tree){
          plot(ind_tree[[i]])
        } else {
          do.call(print, c(ind_tree[[i]], cols))
        }
      }
      
    }
    
    
  }
    
}
