#' Wrapper of different combining algorithms
#'
#' @param logpath a list of names of files containing BEAST log outputs
#' @param treepath a list of names of files containing BEAST trees outputs
#' @param algo c("debias", "gaussian", "weiszfeld")
#' @param time.offset a vector of length K, to be added to the times of each subsets
#' @param num.samples only the last num.samples samples are used
#' @param recon.len number of grids for the reconstructed population size
#' @param skip number of lines to skip when reading the BEAST log outputs
#' @param tree.prior c("skyline", "skyride")
#' @param maxiter maximum number of iterations of Weiszfeld algorithm
#' @param fct a factor to be multiplied to the times and population sizes
#'
#' @return reconstructed population size
#' @export
#'
combine = function(logpath, treepath, algo = "debias", time.offset = NULL, num.samples = 1000,
                   recon.len = 100, skip = 3, tree.prior = "skyline", maxiter = 1000, fct = 1){
  K = length(logpath)
  if(algo == "debias"){
    # if(!hetero)
    #   res = combine.debiased(logpath, treepath, skip = skip, tree.prior = tree.prior, recon.len = recon.len)
    # else
      res.debias = combine.debiased.hetero(logpath, treepath, time.offset = time.offset, num.samples = num.samples, recon.len = recon.len, skip = skip)
      time.grid.recon = res.debias[[1]]
      N.recon = res.debias[[2]]
      N.quant.debiased = data.frame(med = apply(N.recon, 2, quantile, prob = 0.5),
                                    lo = apply(N.recon, 2, quantile, prob = 0.025),
                                    hi = apply(N.recon, 2, quantile, prob = 0.975, na.rm = TRUE))
    return(list(time.grid.recon, N.quant.debiased))
  }

  res = process.beast.logs(logpath, treepath, time.offset = time.offset, fct = fct, tree.prior = tree.prior, skip = skip, recon.len = recon.len, M = M, coal.return = TRUE)
  f.all = res[[1]]
  t.all = res[[2]]
  cutoff = res[[3]]
  time.grid.recon = seq(0, cutoff, length = recon.len)
  # Gaussian approx
  if(algo == "gaussian"){
    res.gauss = Gaussian.approx(f.all)
    mu.bar = res.gauss[[1]]
    Sigma.bar = res.gauss[[2]]
    # time.grid.recon = seq(0, cutoff, length = recon.len)
    sigma.diag = sqrt(diag(Sigma.bar))
    N.quant = data.frame(med = mu.bar, lo = mu.bar - sigma.diag * qnorm(0.975), hi = mu.bar + sigma.diag * qnorm(0.975))
    return(list(time.grid.recon, exp(N.quant)))
  }
  if(algo == "weiszfeld"){
    omega.star = Weiszfeld(f.all, maxiter = maxiter)
    draw.samples = NULL
    for(i in 1:K){
      idx = sample(1000, floor(omega.star[i] * 1000), replace = TRUE)
      draw.samples = rbind(draw.samples, t(f.all[i, , idx] ))
    }
    N.quant = data.frame(med = apply(draw.samples, 2, quantile, prob = 0.5),
                         lo = apply(draw.samples, 2, quantile, prob = 0.025),
                         hi = apply(draw.samples, 2, quantile, prob = 0.975))
    return(list(time.grid.recon, exp(N.quant)))
  }
}
