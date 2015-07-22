W2E <-function(x) cbind(which(x!=0,arr.ind=TRUE),x[x!=0])

recoderFunc <- function(data, 
                        oldvalue, 
                        newvalue) {
  if (is.factor(data))     data     <- as.character(data)
  if (is.factor(oldvalue)) oldvalue <- as.character(oldvalue)
  if (is.factor(newvalue)) newvalue <- as.character(newvalue)
  newvec <- data
  for (i in unique(oldvalue)) newvec[data == i] <- newvalue[oldvalue == i]
  newvec
}
