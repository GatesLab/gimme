#' Write MS-GIMME results to data.frame.
#' @keywords internal 
gimmems.write <- function(x){
  
  out.dir <- x$dat$ctrlOpts$out
  
  #------------------------------------------#
  # Summary fit
  #------------------------------------------#
  
  summaryFit <- as.data.frame(do.call("rbind", 
    
    # grpsol_i                       
    lapply(seq_along(x$ind_fit), function(i){                
      
      as.data.frame(do.call("rbind",
                            
        # grpsol_i_ind_j
        lapply(seq_along(x$ind_fit[[i]]), function(j){         
          
          as.data.frame(do.call("rbind",
                                        
            # grpsol_i_ind_j_indsol_k
            lapply(seq_along(x$ind_fit[[i]][[j]]), function(k){  
            
              grp_sol_i_ind_j_sol_k <-  x$ind_fit[[i]][[j]][[k]]
              
              # we want to add some identifying information
              id_cols <- data.frame(
                "subj" = grp_sol_i_ind_j_sol_k$subj,
                "grp_sol" = i,
                "ind_sol" = k
              )
              
              df <- cbind(id_cols,  t(data.frame(grp_sol_i_ind_j_sol_k$fits)))
              
              # add model fit information
              #df <- cbind(df, t(data.frame(grp_sol_i_ind_j_sol_k$fits)))
              
              df$status <- grp_sol_i_ind_j_sol_k$status
              
              df
            
            })
            
          ))
          
        })
      
      ))
      
    })
    
  ))
  
  # now let's write this df to the ms directory
  write.csv(summaryFit, file.path(out.dir, "summaryFit.csv"), row.names=FALSE)
  
  
  indivPathEstimates <- as.data.frame(do.call("rbind", 
    
    # grpsol_i                       
    lapply(seq_along(x$ind_fit), function(i){                
      
      as.data.frame(do.call("rbind",
                            
        # grpsol_i_ind_j
        lapply(seq_along(x$ind_fit[[i]]), function(j){         
          
          as.data.frame(do.call("rbind",
                                        
            # grpsol_i_ind_j_indsol_k
            lapply(seq_along(x$ind_fit[[i]][[j]]), function(k){  
            
              grp_sol_i_ind_j_sol_k <-  x$ind_fit[[i]][[j]][[k]]
              
              df <- grp_sol_i_ind_j_sol_k$coefs
              
              # we can clean up this df a little. for example,
              # we don't want to include any impossible paths.
              # such as those where lag is on the LHS. 
              
              df <- df[!grepl("lag", df$lhs) & grepl("~", df$op),]
              
              # we want to add some identifying information
              id_cols <- data.frame(
                "subj" = grp_sol_i_ind_j_sol_k$subj,
                "grp_sol" = i,
                "ind_sol" = k
              )
              
              df <- cbind(id_cols, df)
              
              # add model fit information
              #df <- cbind(df, t(data.frame(grp_sol_i_ind_j_sol_k$fits)))
              
              df$op <- df$ci.lower <- df$ci.upper <- NULL
              
              
              colnames(df) <-c("subj", "grp_sol", "ind_sol", 
                               "dv", "iv", "beta", "se", 
                               "z", "pval")
              
              df
            
            })
            
          ))
          
        })
      
      ))
      
    })
    
  ))
  
  # now let's write this df to the ms directory
  write.csv(indivPathEstimates, file.path(out.dir, "indivPathEstimates.csv"), row.names=FALSE)
  
  invisible(indivPathEstimates)

}