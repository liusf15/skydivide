#' combine by debiasing
#'
#' @param logpath a list of paths of the log files
#' @param treepath a list of paths of the tree files
#' @param fct a scaling factor, the population sizes and times are multiplied by fct #TODO: add fct
#' @param skip number of rows to skip when reading the log files
#' @param tree.prior c("skyline", "skyride") #TODO: add skygrid
#' @param num.samples number of samples (from the last sample) to use
#' @param recon.len number of grid points for the reconstructed population size
#' @return time points and the corresponding population sizes
#' @export
#'
combine.debiased = function(logpath, treepath, skip = 3, tree.prior = "skyline", num.samples = 1000, recon.len = 100){
  K = length(logpath)
  f.log = list()
  t.all = list()
  group.size = list()
  treefiles = list()
  n = 0
  for(i in 1:K){ # read in data
    tmp = read.table(file=logpath[[i]], skip = skip, header = TRUE)
    total.len = dim(tmp)[1]
    # print(total.len)
    tmp = tmp[(total.len - num.samples + 1) : total.len, ]
    if (tree.prior == "skyline")
      pop.idx = grep("*popSize*", names(tmp))
    if (tree.prior == "skyride")
      pop.idx = grep("*logPopSize*", names(tmp))
    f.log[[i]] = tmp[1:num.samples, pop.idx]
    if (tree.prior == "skyline")
      f.log[[i]] = log(f.log[[i]])
    group.idx = grep("*groupSize*", names(tmp))
    group.size[[i]] = tmp[, group.idx]
    tmp = ape::read.nexus(treepath[[i]])
    treefiles[[i]] = tmp[(length(tmp) - num.samples + 1) : length(tmp)]
    t.all[[i]] = matrix(0, treefiles[[i]][[1]]$Nnode, num.samples)
    n = n + treefiles[[i]][[1]]$Nnode + 1
  }

  coaltimes.comb = matrix(0, num.samples, n - K)  # put all the coalescent times together
  N.debiased = matrix(0, num.samples, n - K - 1)  # debiased Ne at coaltimes.comb
  for(i in 1:num.samples){
    coal.tmp = c()
    subset.idx = c()

    for(j in 1:K){
      t.sub = cumsum(ape::coalescent.intervals(treefiles[[j]][[i]])$interval.length)
      # t.all[[j]][, i] = t.sub
      coal.tmp = c(coal.tmp, t.sub)
      subset.idx = c(subset.idx, rep(j, length(t.sub)))
    }
    dat = data.frame(coal.tmp, subset.idx)
    dat.sort = sort(dat$coal.tmp, index.return=TRUE)
    coal.tmp = dat.sort$x
    subset.idx = subset.idx[dat.sort$ix]  # subset idx after sorting the coalescent times

    coaltimes.comb[i, ] = coal.tmp

    N.all = matrix(0, n - K, K)  # subset Ne at coaltimes.comb
    lineage.sub = matrix(0, n - K, K)  # record the num of lineages for each subtree
    c.sub = matrix(0, n - K, K)  # n choose 2
    for(j in 1:K){
      is.merge = which(subset.idx == j)
      n.i = length(is.merge)
      lineage.tmp = rep(seq(n.i, 1), c(diff(c(0, is.merge))))
      if(length(lineage.tmp) < n - K){  # complete with zeros
        lineage.tmp = c(lineage.tmp, rep(0, n - K - length(lineage.tmp)))
      }
      lineage.sub[, j] = lineage.tmp
      c.sub[, j] = lineage.tmp * (lineage.tmp - 1) / 2

      t.sub = cumsum(ape::coalescent.intervals(treefiles[[j]][[i]])$interval.length)
      N.sub = exp(rep(unlist(f.log[[j]][i, ]), group.size[[j]][i, ]))  # log to normal scale
      N.all[, j] = N.standardize(t.sub, N.sub, coal.tmp)
    }
    lineage.tot = apply(lineage.sub, 1, sum)  # lineages of whole tree
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
  # end.idx = which(is.na(N.debiased[1,]))[1] - 1
  end.idx = min(unlist(apply(apply(N.debiased, 2, is.na), 1, which))) - 1
  for(i in 1:num.samples){
    N.recon[i, ] = N.standardize(coaltimes.comb[i, 1:end.idx], N.debiased[i, 1:end.idx], time.grid.recon)
  }
  return(list(time.grid.recon, N.recon))
}
