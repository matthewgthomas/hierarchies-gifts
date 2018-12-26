#######################################################################
## Model comparison and diagnostics
## implement a way to compute loo using sparse model design matrices 
## (because we have a *lot* of observations and a *lot* of groups so using standard matrices eats all the memory)
##
# this function calculates the pointwise log-likelihoods
# here, data and draws will be sparse matrices, so need to convert them 
# see the 'loo' package documentation for futher info about this kind of function
llfun = function(i, data, draws) {
  #draws = as.matrix(draws)  # maybe better to do this once, before calling the function - rather than on every iteration within the function...
  
  # convert the sparse matrix containing the data into a normal matrix and chuck out the 'y' and 'trials' columns
  d2 = data[, -which(colnames(data) %in% c("y", "trials"))]
  d2 = as.matrix(d2)
  
  # calculate the probabilities for the current observation
  eta = draws %*% d2
  eta = plogis(eta)
  # calculate log-likelihood
  dbinom(data[,"y"], size = data[,"trials"], prob = eta, log = TRUE)
}

#' This function takes a model fit with rstanarm and creates the 'args' parameter for loo.function()
#'
#' @param m stanreg, The model
#' @return nlist of arguments (data, draws, S, N)
get_args = function(m) {
  library(Matrix)
  
  # make a sparse data matrix containing the fixed and random effects
  x = get_x(m)
  data = Matrix(data = c(y=d.long$GiftGiven, trials=rep(1, nrow(d.long)), x=x), ncol=2+ncol(x), dimnames=list(NULL, c("y", "trials", dimnames(x)[[2]])))
  z <- get_z(m)
  data = Matrix::cBind(data, z)
  
  # make a sparse matrix of draws from the posterior distribution
  post = Matrix(as.matrix(m))
  draws = post[, seq_len(ncol(x)), drop = FALSE]
  b <- post[, rstanarm:::b_names(colnames(post)), drop = FALSE]
  draws <- Matrix::cBind(draws, b)
  draws = as.matrix(draws)
  
  args = nlist(data = data, draws = draws, S = NROW(draws), N = nrow(data))
  args
}

#' Calculates loo using sparse matrices
#'
#' @param m stanreg, The model
#' @return loo object
#' 
loo_sparse = function(m, cores=1) {
  args = get_args(m)
  # calculate loo for this model
  loo_out = loo(llfun, args = args, cores=cores)
  attr(loo_out, "model")  = deparse(substitute(m))
  return(loo_out)
}

#' Calculate log-likelihood matrix using sparse model design matrices
#'
#' @seealso get_args, llfun
#' @param object stanreg, The model, fit with rstanarm
#'
#' @return matrix of log-likelihoods
#' 
log_lik.sparse <- function(object) {
  args <- get_args(object)

  out <- vapply(
    seq_len(args$N),
    FUN = function(i) {
      as.vector(llfun(
        i = i,
        data = args$data[i,, drop = FALSE],
        draws = args$draws
      ))
    },
    FUN.VALUE = numeric(length = args$S)
  )
  colnames(out) <- rownames(model.frame(object))
  return(out)
}


#' A tweak of loo::model_weights() for cases where `loo` objects have already been calculated for models
#' 
#' All params apart from `loo_list` are the same as in loo::model_weights()
#'
#' @param loo_list list, A bunch of `loo` objects
#' @param method 
#' @param BB 
#' @param BB_n 
#' @param alpha 
#' @param seed 
#' @param optim_method 
#' 
#' @inheritParams loo::model_weights
#' @seealso loo::model_weights
#'
#' @return A vector of optimal model weights.
#' 
model_weights.loo <-function(loo_list, method="stacking",BB=T,BB_n=1000, alpha=1, seed=NULL, optim_method="BFGS")
{
  if (!method %in%c("stacking","pseudobma") )
    stop("Must specify a method in stacking or pseudobma .")
  K<-length(loo_list)                #number of models
  if (K==1)
    stop("Only one model is found.")
  
  N <- nrow(loos[[1]]$pointwise)       #number of data points
  lpd_point<-matrix(NA,N,K)            #point wise log likelihood
  elpd_loo<-rep(NA,K)
  for( k in 1:K){
    #log_likelihood<- log_lik_list[[k]]
    L <- loo_list[[k]]
    lpd_point[,k] <- L$pointwise[,1]    #calculate log(p_k (y_i | y_-i))
    elpd_loo[k]<-L$elpd_loo
  }
  ## 1) stacking on log score
  if (method =="stacking"){
    w_stacking <- loo::stacking_weight(lpd_point, optim_method=optim_method)
    cat("The stacking weights are:\n")
    print(rbind(paste("Model"  ,c(1:K) ), round(w_stacking*100 )/100))
    return(w_stacking)
  }
  else
    if (method =="pseudobma"){
      uwts <- exp( elpd_loo - max( elpd_loo))
      w_loo1 <- uwts / sum(uwts)
      if(BB==F){
        cat("The Pseudo-BMA weights are:\n")
        print(rbind(paste("Model"  ,c(1:K) ),  round(w_loo1*100 )/100))
        return(w_loo1)
      }
      if(BB==T){
        w_loo2  <- loo::pseudobma_weight(lpd_point, BB_n,alpha, seed)
        cat("The Pseudo-BMA+ weights using Bayesian Bootstrap  are:\n ")
        print(rbind(paste("Model",c(1:K) ),  round(w_loo2*100 )/100))
        return (w_loo2 )
      }
    }
}


#' Model comparison with delta_elpd and delta_se
#'
#' @param ... At least two objects returned by \code{\link{loo}} or
#'   \code{\link{waic}}.
#'
#' @return A data.frame of model comparison info.
#' 
#' @details See \code{compare_models} in \pkg{rstanarm} and \code{compare} in \pkg{loo}.
#' 
compare_models_delta <- function(...)
{
  L = list(...)
  #L = list(loo.herd.null, loo.herd.control, loo.herd.deg, loo.herd.deg_int, loo.herd.btwn, loo.herd.eigen)
  
  if (!all(sapply(L, function(x) inherits(x, "loo"))))
    stop("All inputs should have class 'loo'.")
  if (length(L) <= 1L)
    stop("'compare' requires at least two models.")
  
  # retrieve model names from loo objects
  mnames = sapply(L, function(x) attr(x, "name"))
  # mnames <- as.character(match.call(expand.dots = TRUE))[-1L]

  if (length(mnames)==0)
    mnames <- paste0("model", seq_along(L))

  # compute SE of differences between adjacent models from top to bottom in ranking
  # this part was shamelessly nabbed and tweaked from compare() in Richard McElreath's rethinking package
  dSE.matrix <- matrix( NA , nrow=length(L) , ncol=length(L) )
  
  colnames(dSE.matrix) <- mnames
  rownames(dSE.matrix) <- mnames
  
  for ( i in 1:(length(L)-1) ) {
    for ( j in (i+1):length(L) ) {
      # get pointwise looic/waic from each loo object
      loo_ptw1 <- L[[i]]$pointwise[, grep("^elpd", dimnames(L[[i]]$pointwise)[[2]])]
      loo_ptw2 <- L[[j]]$pointwise[, grep("^elpd", dimnames(L[[j]]$pointwise)[[2]])]
      
      dSE.matrix[i,j] <- as.numeric( sqrt( length(loo_ptw1)*var( loo_ptw1 - loo_ptw2 ) ) )
      dSE.matrix[j,i] <- dSE.matrix[i,j]
    }#j
  }#i
  
  IC.list <- abs( unlist( sapply(L, function(x) x[ grep("^elpd", names(x)) ]) ) )
  dIC <- IC.list - min( IC.list )
  
  topm <- which( dIC==0 )
  dSEcol <- dSE.matrix[,topm]
  
  # w.IC <- rethinking::ICweights( IC.list )  # use rethinking package to calculate model weights
  
  result <- data.frame( delta=dIC, d_se=dSEcol )
  rownames(result) <- mnames
  result <- result[ order( result[["delta"]] ) , ]
  
  if (length(L)==2)
  {
    # two loo objects: produce a full model comparison table
    # (this part shamelessly lifted from loo::compare())
    sel <- grep("pointwise|pareto_k", names(L[[1L]]), invert = TRUE)
    x <- sapply(L, function(x) unlist(x[sel]))
    
    colnames(x) <- mnames
    rnms <- rownames(x)
    comp <- x
    patts <- c("^waic$|^looic$", "^se_waic$|^se_looic$", "elpd", "p_")
    row_ord <- unlist(sapply(patts, function(p) grep(p, rownames(comp))),
                      use.names = FALSE)
    col_ord <- order(x[grep("^elpd", rnms), ], decreasing = TRUE)
    comp <- as.data.frame( t(comp[row_ord, col_ord]) )
  } else {
    # three or more loo objects, so use rstanarm's comparison function as is
    #comp <- as.data.frame( rstanarm::compare_models( ... ) )
    comp <- as.data.frame( loo::compare( x=L ) )
  }
  
  # stick the deltas and weights on the end (in the correct order)
  #mnames_sorted = row.names(comp)
  comp <- dplyr::bind_cols(comp, result)
  rownames(comp) <- rownames(result)
  comp
  
  # rm(L, mnames, dSE.matrix, loo_ptw1, loo_ptw2, IC.list, dIC, topm, dSEcol, result, comp, i,j)
}
