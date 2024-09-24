#' @name simulateVAR
#' @title Simulate data from Vector AutoRegression (VAR) models.
#' @description This function simulates data. It allows for structural VAR and VAR data generating models. 
#' @usage
#' simulateVAR(A   = NULL, 
#'             Phi       = NULL, 
#'             Psi       = NULL, 
#'             subAssign = NULL, 
#'             N         = NULL, 
#'             ASign     = "random",  
#'             PhiSign   = "random",  
#'             Obs       = NULL,
#'             indA      = 0.01, 
#'             indPhi    = 0.01,
#'             indPsi    = 0.00)
#' @param A A matrix (for no subgroups) or list of A matrices, with slice # = # of subgroups. 
#' @param Phi Phi matrix (for no subgroups) or list of Phi matrices, with slice # = # of subgroups.  
#' @param Psi matrix (for no subgroups) or list of Psi matrices, with slice # = # of subgroups. 
#' @param subAssign Optional vector of length N that indicates which subgroup each individual is in. 
#' @param N Number of indvidiuals.
#' @param Obs Number of observations (T) per individual. Burn in of 400 is used to generate then discarded.  
#' @param indA Sparsity of individual-level A paths. 0 indicates no individual-level. Use decimals. Default is 
#' 0.01, meaning that each path that is not in the group-level A matrix has a 0.01 chance of being added. 
#' @param indPhi Sparsity of individual-level Phi paths. 0 indicates no individual-level. Use decimals. Default is 
#' 0.01, meaning that each path that is not in the group-level Phi matrix has a 0.01 chance of being added.
#' @param indPsi Sparsity of individual-level Psi paths. 0 indicates no individual-level. Use decimals. Default is 
#' 0, meaning that each path that is not in the group-level Psi matrix has a 0 chance of being added at the ind. level.
#' Individual- level paths added at this rate per individual. 
#' @param ASign Defaults to "random" for ind level paths, with 50 percent chance of positive and 50 percent negative, other option is either "neg" or "pos" which provides all negative or all positive relations, respectively.
#' @param PhiSign Defaults to "random" for ind level paths, with 50 percent chance of positive and 50 percent negative, other option is either "neg" or "pos" which provides all negative or all positive relations, respectively. 
#' @author KM Gates, Ai Ye, Ethan McCormick, & Zachary Fisher 
#' @export simulateVAR

simulateVAR <- function(A         = NULL, 
                        Phi       = NULL, 
                        Psi       = NULL, 
                        subAssign = NULL, 
                        N         = NULL, 
                        ASign     = "random",  
                        PhiSign   = "random",  
                        Obs       = NULL,
                        indA      = 0.01, 
                        indPhi    = 0.01,
                        indPsi    = 0.00) {
  
  AList   <- list()
  PhiList <- list()
  PsiList <- list()
  dataList  <- list()
  subgroups <- unique(subAssign)
  
  ###### Checks and Rules #####
  
  if(!is.list(Phi)){
    PhiList[[1]] <- Phi
  }
  
  if(is.null(Psi)){
    Psi <- matrix(0, nrow = dim(PhiList[[1]])[1], ncol = dim(PhiList[[1]])[2])
    diag(Psi) <- 1
  }
  
  if(is.null(A)){
    stop(paste0("ERROR: No data generating matrix provided for A",
                " Please provide a matrix or list of matrices. They can be zero."))  }
  
  if(is.list(A)) {
    AList   <- A
    PhiList <- Phi
    PsiList <- Psi
    if(is.null(subAssign)){
      stop(paste0("ERROR: Multiple patterns provided with no subgroup assignments",
                  " Please provide subgroup assignments in subAssign argument."))
    }
    if(length(is.list(A)) != length(is.list(Phi))) 
      stop(paste0("ERROR: Different numbers of matrices provided for A and Phi.",
                  " Please ensure there is one A and Phi matrix for each subgroup."))
    if(length(is.list(A)) != length(is.list(Psi))) 
      stop(paste0("ERROR: Different numbers of matrices provided for A and Psi.",
                  " Please ensure there is one A and Psi matrix for each subgroup."))
  } else {
    writeLines("NOTICE: One A matrix provided. No subgroups generated.")
    subAssign <- matrix(1, N, 1)
    AList[[1]]   <- A
    PsiList[[1]] <- Psi
  }
  
  if(is.null(N))
    stop(paste0("ERROR: Please provide N for the number of individuals to generate"))
  if(is.null(Obs))
    stop(paste0("ERROR: Please provide Obs for the number of time points per person"))
  
  ##### Data Generation ####
  
  vars <- length(AList[[1]][1,])
  
  Ind.nonstation <- c()
  
  for (ind in 1:N) {
    go <- 1
    counter <- 0                                   
    while (go > 0 & counter < 100){
      ATemp <- AList[[subAssign[ind]]]
      PhiTemp <- PhiList[[subAssign[ind]]]
      PsiTemp <- PsiList[[subAssign[ind]]]
      
      AMean <- mean(ATemp)
      PhiMean <- mean(PhiTemp)
      PsiMean <- mean(PsiTemp)
      
      if(indA>0){
        ATemp[which(ATemp == 0)] <- stats::rbinom(1,1,indA)
        if(ASign == "random"){
          random.sign <- sample(c(0,1), size=1)
          if(random.sign==0) {ATemp[which(ATemp == 1)] <- -stats::rnorm(1, AMean, 0.3)}
          if(random.sign==1) {ATemp[which(ATemp == 1)] <- stats::rnorm(1, AMean, 0.3)}
        }
        if(ASign == "neg")
          ATemp[which(ATemp == 1)] <- -stats::rnorm(1, AMean, 0.3)
        if(ASign == "pos")
          ATemp[which(ATemp == 1)] <- stats::rnorm(1, AMean, 0.3)
        diag(ATemp) <- 0
      }
      
      if(indPhi>0){
        PhiTemp[which(PhiTemp == 0)] <- stats::rbinom(1,1,indPhi)
        if(PhiSign == "random"){
          random.sign <- sample(c(0,1), size=1)
          if(random.sign==0) {PhiTemp[which(PhiTemp == 1)] <- -stats::rnorm(1, PhiMean, 0.3)}
          if(random.sign==1) {PhiTemp[which(PhiTemp == 1)] <- stats::rnorm(1,PhiMean, 0.3)}
        }
        if(PhiSign == "neg")
          PhiTemp[which(PhiTemp == 1)] <- -stats::rnorm(1, , 0.3)
        if(PhiSign == "pos")
          PhiTemp[which(PhiTemp == 1)] <- stats::rnorm(1, , 0.3)
      }
      
      if(indPsi>0){
        PsiTemp[which(PsiTemp == 0)] <- stats::rbinom(1,1,indPsi)
        PsiTemp[which(PsiTemp == 1)] <- stats::rnorm(1, PsiMean, 0.3)
      }
      
      negA <- solve(diag(vars)-ATemp) 
      
      time  <- matrix(0, nrow = vars, ncol = Obs+400) # 400 burn-in observations
      
      time1 <- matrix(0, nrow = vars, ncol = Obs+400)
      
      noise <- negA %*% t(MASS::mvrnorm(n = (Obs+400),rep(0,vars),PsiTemp, empirical = TRUE))
      
      time[,1] <- noise[,1]
      
      time1[,1] <- negA %*% PhiTemp %*% time[,1] + noise[,1]
      
      time[,2]  <- time1[,1]
      
      for (t in 2:(Obs+400)){
        time1[,t]  <- negA %*% PhiTemp %*% time[,(t)] + noise[,t]
        if (t<(Obs+400))
          time[,(t+1)] <- time1[,t]
      }
      go <- 0
      for (c in 1:length(time[,1])){
        # aTSA no longer an option 1/22/2024
        # if(aTSA::adf.test(time[c,], out = FALSE)$type3[1,3]>0.05)
        #   go <- go + 1
        if(suppressWarnings(tseries::adf.test(time[c,])$p.value)>0.05){
        go <- go + 1
        counter <- sum(counter, 1)
        }
      }
    }
    if(counter == 100){
      Ind.nonstation <- append(Ind.nonstation, ind)
      writeLines(paste0('WARNING: No Stationary Time Series Data Generated for Individual:', ind))
    } else {
      dataList[[ind]] <- t(time[,401:(400+Obs)]) 
      names(dataList)[ind]<- paste0('ind', ind)
    }
  }
  
  
  return(list(A=AList,
              Phi = PhiList, 
              Psi = PsiList, 
              dataList = dataList, 
              subAssign = subAssign))
  
}
