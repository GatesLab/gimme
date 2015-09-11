##R port of Bush and Cisler 2013, Magnetic Resonance Imaging
##Adapted from the original provided by Keith Bush

## Author:      Keith Bush, PhD
## Institution: University of Arkansas at Little Rock
## Date:        Aug. 9, 2013

deconvolve_nlreg <- function(BOLDobs, kernel, nev_lr=.01, epsilon=.005) {
    ## Description:
    ## This function deconvolves the BOLD signal using Bush 2011 method
    ##
    ## Inputs:
    ##    BOLDobs - observed BOLD timeseries
    ##    kernel  - assumed kernel of the BOLD signal
    ##    nev_lr  - learning rate for the assignment of neural events
    ##    epsilon - relative error change (termination condition)
    ##
    ## Outputs:
    ##    encoding - reconstructed neural events

    ## Determine time series length
    N = length(BOLDobs)

    ##Calc simulation steps related to simulation time
    K = length(kernel)
    A = K - 1 + N
    
    ##Termination Params
    preverror = 1e9 #previous error
    currerror = 0   #current error

    ##Construct random activation vector (fluctuate slightly around zero between -2e-9 and 2e-9)
    activation = rep(2e-9, A)*runif(A) - 1e-9
    
    ##Presolve activations to fit target_adjust as encoding
    max_hrf_id_adjust = which.max(kernel) - 1 #element of kernel 1 before max
    BOLDobs_adjust = BOLDobs[max_hrf_id_adjust:N]
    pre_encoding = BOLDobs_adjust - min(BOLDobs_adjust)
    pre_encoding = pre_encoding/max(pre_encoding) #unit normalize
    encoding = pre_encoding
    activation[K:(K-1 + length(BOLDobs_adjust))] = log(pre_encoding/(1-pre_encoding))

    while (abs(preverror-currerror) > epsilon) {

        ##Compute encoding vector
        encoding = sigmoid(activation)
        
        ##Construct feature space
        feature = generate_feature(encoding,K)

        ##Generate virtual bold response by multiplying feature (N x K) by kernel (K x 1) to get N x 1 estimated response
        ytilde = feature[K:nrow(feature),] %*% kernel

        ##Convert to percent signal change
        meanCurrent = mean(ytilde)
        brf = ytilde - meanCurrent
        brf = brf/meanCurrent
        
        ##Compute dEdbrf
        dEdbrf = brf - BOLDobs

        ##Assume normalization does not impact deriv much.
        dEdy = dEdbrf

        ##Precompute derivative components
        dEde = diag(K) %*% kernel
        back_error = c(rep(0, K-1), dEdy, rep(0, K-1))
        
        ##Backpropagate Errors
        delta = c()
        for (i in 1:A) {
            active = activation[i]
            deda = dsigmoid(active);
            dEda = dEde * deda
            this_error = back_error[i:(i-1+K)]
            delta = c(delta, sum(dEda * this_error))
        }

        ##Update estimate
        activation = activation - nev_lr * delta

        ##Iterate Learning
        preverror = currerror
        currerror = sum(dEdbrf^2)
    }
    return(encoding)
}
  

## Support functions
sigmoid <- function(x) {
    y <- 1/(1+exp(-x))
    return(y)
}

dsigmoid <- function(x) {                        
    y=(1-sigmoid(x))*sigmoid(x)
    return(y)
}

generate_feature <- function(encoding, K) {
    fmatrix = matrix(0, length(encoding), K)
    fmatrix[,1] = encoding

    for (i in 2:K) {
        fmatrix[,i] = c(rep(0, i-1), encoding[1:(length(encoding) - (i-1))])
    }
    return(fmatrix)
}


spm_hrf <- function(TR, P) {
    ## spm_hrf - Returns a hemodynamic response function as a difference of gammas (canonical)
    ##
    ## USAGE:
    ##   hrf = spm_hrf(TR,[p])
    ##
    ##INPUT:
    ## TR   - scan repetition time (seconds)
    ## P    - 1x7 vector of parameters for the response function (two gamma functions)
    ##                                                 (default value, in seconds)
    ## P[1] - delay of response (relative to onset)    (6)
    ## P[2] - delay of undershoot (relative to onset)  (16)
    ## P[3] - dispersion of response                   (1)
    ## P[4] - dispersion of undershoot                 (1)
    ## P[5] - ratio of response to undershoot          (6)
    ## P[6] - onset (seconds)                          (0)
    ## P[7] - length of kernel (seconds)               (32)
    ##
    ## OUTPUT:
    ##  $kernel - hemodynamic response function
    ##  $p      - parameters of the response function
    ##_______________________________________________________________________
    ## Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

    ## Karl Friston
    ## $Id: spm_hrf.m 387 2005-12-17 18:31:23Z klaas $

    ## default parameters
    ##-----------------------------------------------------------------------
    fMRI_T = 16 #microtime resolution is 1/16th of TR

    p = c(6, 16, 1, 1, 6, 0, 32)

    if (!missing(P)) {
        p[1:length(P)] = P
    }

    ## modelled hemodynamic response function - {mixture of Gammas}
    ##-----------------------------------------------------------------------
    dt    = TR/fMRI_T
    u     = 0:(p[7]/dt) - p[6]/dt #sampling grid of HRF in microtime units (e.g., 0:256 for default params and 2s TR)

    ## MH: eliminated use of spm Gpdf in favor of built-in gamma PDF in R. Checked that this yields identical
    ## results, although the third parameter of spm_Gpdf is really rate, not shape.
    hrf   = dgamma(u, shape=p[1]/p[3], rate=dt/p[3]) - dgamma(u, shape=p[2]/p[4], rate=dt/p[4])/p[5]

    hrf   = hrf[0:(p[7]/TR)*fMRI_T + 1] #subsample hrf back onto TR grid
    hrf   = hrf/sum(hrf); #Unit normalize hrf

    return(list(kernel=hrf, p=p))
}
