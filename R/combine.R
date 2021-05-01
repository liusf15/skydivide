#' Gaussian approximation
#'
#' @param f.all an array with shape (K, recon.len, num.samples), where K is the number of subsets, recon.len is the number of grid points, num.samples is the number of samples used
#'
#' @return reconstructed population size
#' @export
#'
Gaussian.approx = function(f.all){
  K = dim(f.all)[1]
  recon.len = dim(f.all)[2]
  num.samples = dim(f.all)[3]
  mu.mean = matrix(0, K, recon.len)
  Sigma.mean = array(0, dim = c(K, recon.len, recon.len))
  # compute sample mean and covariance for each subsets
  for(j in 1:K){
    mu.mean[j, ] = apply(f.all[j, , ], 1, mean)
    Sigma.mean[j, , ] = cov(t(f.all[j, , ]))
  }
  # compute weighted average
  mu.bar = matrix(0, recon.len, 1)
  Sigma.bar = matrix(0, recon.len, recon.len)
  for(j in 1:K){
    Sigma.bar = Sigma.bar + solve(Sigma.mean[j, , ] + 0.001 * diag(recon.len))
    mu.bar = mu.bar + solve(Sigma.mean[j, , ] + 0.001 * diag(recon.len)) %*% mu.mean[j, ]
  }
  Sigma.bar = solve(Sigma.bar)
  mu.bar = Sigma.bar %*% mu.bar
  return(list(mu.bar, Sigma.bar))
}

#' Weiszfeld algorithm
#'
#' @param f.all an array with shape (K, recon.len, num.samples), where K is the number of subsets, recon.len is the number of grid points, num.samples is the number of samples used
#' @param maxiter maximum number of iterations for Weiszfeld algorithm
#'
#' @return reconstructed population size
#' @export
#'
Weiszfeld = function(f.all, maxiter = 1000){
  K = dim(f.all)[1]
  recon.len = dim(f.all)[2]
  num.samples = dim(f.all)[3]
  mu.mean = apply(f.all, c(1, 2), mean)
  omega = rep(1 / K, K)
  for(t in 1:maxiter){
    dd = rep(1, K)
    mean.center = apply(omega * mu.mean, 2, sum)
    for(j in 1:K){
      dd[j] = norm(as.matrix(mean.center - mu.mean[j, ]), "2")
    }
    omega.old = omega
    omega = (1 / dd) / (sum(1 / dd))
    if(norm(as.matrix(omega - omega.old)) < 1e-6){
      cat("Weiszfeld algorithm onverged in", t, "steps", "\n")
      break
    }
  }
  return(omega)
}


#' combine by the debiasing procedure
#'
#' @param logpath a list of paths of the log files
#' @param treepath a list of paths of the tree files
#' @param fct a scaling factor, the population sizes and times are multiplied by fct #TODO: add fct
#' @param time.offset a vector of length K, to be added to the times of each subsets
#' @param tree.prior c("skyline", "skyride") #TODO: add skygrid
#' @param num.samples number of samples (from the last sample) to use
#' @param recon.len number of grid points for the reconstructed population size
#' @param skip number of rows to skip when reading the log files
#'
#' @return a time grid of length recon.len and the reconstructed population sizes at corresponding time grids
#' @export
#'
combine.debiased.hetero = function(logpath, treepath, fct = 1, time.offset = NULL, tree.prior = "skyline", num.samples = 1000, recon.len = 100, skip = 3){
  K = length(logpath)
  if(is.null(time.offset))
    time.offset = rep(0, K)
  f.log = list()
  t.all = list()
  group.size = list()
  treefiles = list()
  num.samples = 1000
  n = 0
  # subset.list = c(1, 2, 3, 4, 5, 7, 10)
  # K = length(subset.list)
  for(k in 1:K){ # read in data
    cat("Reading subset", k, "\n")
    tmp = read.table(file=logpath[[k]], skip = skip, header = TRUE)
    total.len = dim(tmp)[1]
    if(num.samples > total.len)
      num.samples = total.len
    # print(total.len)
    tmp = tmp[(total.len - num.samples + 1) : total.len, ]
    if (tree.prior == "skyline")
      pop.idx = grep("*popSize*", names(tmp))
    if (tree.prior == "skyride")
      pop.idx = grep("*logPopSize*", names(tmp))
    f.log[[k]] = tmp[1:num.samples, pop.idx]
    if (tree.prior == "skyline")
      f.log[[k]] = log(f.log[[k]])

    # pop.idx = grep("*logPopSize*", names(tmp))
    # f.log[[k]] = tmp[, pop.idx]
    group.idx = grep("*groupSize*", names(tmp))
    group.size[[k]] = tmp[, group.idx]
    tmp = ape::read.nexus(treepath[[k]])
    treefiles[[k]] = tmp[(length(tmp) - num.samples + 1) : length(tmp)]
    t.all[[k]] = matrix(0, treefiles[[k]][[1]]$Nnode, num.samples)
    n = n + treefiles[[k]][[1]]$Nnode + 1
  }
  # # align times
  # tip.times = c()
  # for(i in 1:K){
  #   all_dates = matrix(unlist(strsplit(treefiles[[i]][[1]]$tip.label, split = "\\|")), ncol=3, byrow=TRUE)[, 3]
  #   all_dates = decimal_date(as.Date(all_dates))
  #   tip.times = c(tip.times, max(all_dates))
  # }
  # last.time = max(tip.times)
  # time.offset = last.time - tip.times

  f.all = array(0, dim = c(K, n - K, num.samples))

  tmp = c()
  recon.len = 100
  N.debiased = matrix(0, num.samples, n - K - 1)
  coaltimes.comb = matrix(0, num.samples, n - K)  # put all the coalescent times together
  for(i in 1:num.samples){
    if(i %% 200 == 0 )
      cat("Processing ", i / num.samples * 100, "%\n", sep = "")
    coal.tmp = c()
    subset.idx = c()
    lineages.sub = c()
    coaltimes.sub = c()
    for(j in 1:K){
      t.sub = phylodyn::summarize_phylo(treefiles[[j]][[i]])
      coal.sub = t.sub$coal_times + time.offset[j]
      coaltimes.sub[[j]] = coal.sub
      samp.sub = t.sub$samp_times + time.offset[j]
      coal_and_samp = c(coal.sub, samp.sub)
      coal.idx = c(rep(1, length(coal.sub)), rep(0, length(samp.sub)))
      sort_time = sort(coal_and_samp, index.return = TRUE)
      coal_and_samp = sort_time$x
      coal.idx = coal.idx[sort_time$ix]
      lineages.tmp = rep(0, length(coal_and_samp))
      for(q in 1:length(samp.sub)){
        lineages.tmp[coal_and_samp > samp.sub[q]] = lineages.tmp[coal_and_samp > samp.sub[q]] + t.sub$n_sampled[q]
      }
      for(q in 1:length(coal.sub)){
        lineages.tmp[coal_and_samp > coal.sub[q]] = lineages.tmp[coal_and_samp > coal.sub[q]] - 1
      }
      lineages.sub[[j]] = lineages.tmp[coal.idx == 1]  # only need lineages at coal times
      # t.all[[j]][, i] = t.sub
      coal.tmp = c(coal.tmp, coal.sub)
      subset.idx = c(subset.idx, rep(j, length(coal.sub)))
    }
    dat = data.frame(coal.tmp, subset.idx)
    dat.sort = sort(dat$coal.tmp, index.return=TRUE)
    coal.tmp = dat.sort$x  # all coal times sorted
    subset.idx = subset.idx[dat.sort$ix]  # subset idx after sorting the coalescent times

    coaltimes.comb[i, ] = coal.tmp

    N.all = matrix(0, n - K, K)
    lineage.final = matrix(0, n - K, K)  # number of lineages in sub trees
    c.sub = matrix(0, n - K, K)
    for(j in 1:K){
      is.merge = which(subset.idx == j)
      n.i = length(is.merge)
      lineages.tmp = lineages.sub[[j]]  # lineages at sub coal times
      lineages.tmp = rep(lineages.tmp, c(diff(c(0, is.merge))))  # lineages at all coal times
      # addition.scale = rev(diff(c(0, is.merge)))
      # lineage.tmp = rep(seq(n.i, 1) * addition.scale, c(diff(c(0, is.merge))))
      # diff.tmp = diff(is.merge)
      if(length(lineages.tmp) < n - K){
        lineages.tmp = c(lineages.tmp, rep(0, n - K - length(lineages.tmp)))
      }
      lineage.final[, j] = lineages.tmp

      c.sub[, j] = lineages.tmp * (lineages.tmp - 1) / 2
      # length(rep(c(diff(c(0, is.merge))), c(diff(c(0, is.merge)))))
      # lineage.sub[is.merge, j] = seq(n.i, 1)

      t.sub = coaltimes.sub[[j]]
      N.sub = exp(rep(unlist(f.log[[j]][i, ]), group.size[[j]][i, ]))
      N.all[, j] = N.standardize(t.sub, N.sub, coal.tmp)
    }
    lineage.tot = apply(lineage.final, 1, sum)
    c.tot = lineage.tot * (lineage.tot - 1) / 2
    N.final = rep(0, n - K)
    num.contri = rep(0, n - K)
    for(j in 1:K){
      msk = c.sub[, j] > 0
      num.contri[msk] = num.contri[msk] + 1
      N.final[msk] = N.final[msk] + N.all[msk, j] * c.tot[msk] / c.sub[msk, j]
    }
    N.debiased[i, ] = N.final[1:(n - K - 1)] / num.contri[1:(n - K - 1)]

  }
  cutoff = quantile(coaltimes.comb[, n - K], 0.5)  # median of all the TMRA
  time.grid.recon = seq(0, cutoff, length = recon.len)

  N.recon = matrix(0, num.samples, recon.len)
  end.idx = n - K - 1
  for(i in 1:num.samples){
    N.recon[i, ] = N.standardize(coaltimes.comb[i, 1:(end.idx - 1)], N.debiased[i, 1:(end.idx - 1)], time.grid.recon)
  }
  return(list(time.grid.recon * fct, N.recon * fct))
}
