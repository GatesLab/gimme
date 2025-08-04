#' Counts number of excellent fit indices
#' @param indices A vector of fit indices from lavaan.
#' @param rmsea_cutoff Cutoff for RMSEA for an individual model(default is .05; must be between 0.0 and 1.0).
#' @param srmr_cutoff Cutoff for SRMR for an individual model (default is .05; must be between 0.0 and 1.0).
#' @param nnfi_cutoff Cutoff for NNFI for an individual model (default is .95; must be between 0.0 and 1.0).
#' @param cfi_cutoff Cutoff for CFI for an individual model (default is .95; must be between 0.0 and 1.0).
#' @return The number of fit indices that are excellent.
#' @keywords internal 
count.excellent <- function(indices,
                            rmsea_cutoff = .05,
                            srmr_cutoff = .05,
                            nnfi_cutoff = .95,
                            cfi_cutoff = .95) {
  rmseaE    <- ifelse(indices[4] <= rmsea_cutoff, 1, 0)
  srmrE     <- ifelse(indices[5] <= srmr_cutoff, 1, 0)
  nnfiE     <- ifelse(indices[6] >= nnfi_cutoff, 1, 0)
  cfiE      <- ifelse(indices[7] >= cfi_cutoff, 1, 0)
  excellent <- sum(rmseaE, srmrE, nnfiE, cfiE, na.rm = TRUE)
  return(excellent)
}