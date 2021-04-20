#' Project a vector of population size to a specified grid
#'
#' @param N.i vector of population sizes at time t.i
#' @param t.i vector of time points at which the population sizes are N.i
#' @param time.grid grid of time
#'
#' @return a vector of population sizes at time.grid
#' @export
#'
N.standardize = function(t.i, N.i, time.grid){
  N.recon = rep(0, length(time.grid))
  grid.width = time.grid[2] - time.grid[1]
  n.cur = 1
  cut = max(which(time.grid <= tail(t.i, 1)))
  for(j in 1:cut){
    # change: find the smallest n.cur such that time.grid[j] <= t.i[n.cur]
    idx = min(which(t.i >= time.grid[j]))
    N.recon[j] = N.i[min(idx, length(N.i))]
    # length(N.recon)
  }
  if(cut < length(N.recon))
    N.recon[(cut + 1) : length(N.recon)] = N.recon[cut]
  return(N.recon)
}

#' test
#'
#' @param logpath a list of paths of the log files
#' @param treepath a list of paths of the tree files
#' @param M number of groups in skyline
#' @param fct a scaling factor, the population sizes and times are multiplied by fct
#' @param skip number of rows to skip when reading the log files
#' @param tree.prior c("skyline", "skyride") #TODO: add skygrid
#' @param num.samples number of samples (from the last sample) to use
#' @param recon.len number of grid points for the reconstructed population size
#' @param coal.return return coalescent times or not
#'
#' @return f.all
#' @export
#'
process.beast.logs = function(logpath, treepath = NULL, M = 15, fct = 1, skip = 3, tree.prior = "skyline", num.samples = 1000, recon.len = 30, coal.return = FALSE){
  K = length(logpath)
  if (tree.prior %in% c("skyline", "skyride")){
    f.log = list()
    t.all = list()
    group.size = list()
    treefiles = list()
    # read in the log and trees
    for(i in 1:K){
      tmp = read.delim(file=logpath[[i]], skip = skip, header = TRUE)
      if (tree.prior == "skyride"){
        pop.idx = grep("*logPopSize*", names(tmp))
      }

      if (tree.prior == "skyline"){
        pop.idx = grep("*popSize*", names(tmp))
      }
      tot.len = dim(tmp)[1]
      tmp = tmp[(tot.len - num.samples + 1) : tot.len, ]
      if (tree.prior == "skyride"){
        f.log[[i]] = tmp[, pop.idx] + log(fct)
      }
      if (tree.prior == "skyline"){
        f.log[[i]] = tmp[, pop.idx] * fct
        f.log[[i]] = log(f.log[[i]])
      }
      group.idx = grep("*groupSize*", names(tmp))
      group.size[[i]] = tmp[, group.idx]
      tmp = read.nexus(treepath[[i]])
      treefiles[[i]] = tmp[(length(tmp) - num.samples + 1) : length(tmp)]
      t.all[[i]] = matrix(0, treefiles[[i]][[1]]$Nnode, num.samples)
    }

    f.all = array(0, dim = c(K, recon.len, num.samples))
    tmp = c()
    for(j in 1:K){
      for(i in 1:num.samples){
        t.sub = cumsum(coalescent.intervals(treefiles[[j]][[i]])$interval.length) * fct

        t.all[[j]][, i] = t.sub
      }
      tmp = c(tmp, t.all[[j]][treefiles[[j]][[1]]$Nnode, ])
    }
    cutoff = quantile(tmp, 0.5)
    time.grid.recon = seq(0, cutoff, length = recon.len)
    for(j in 1:K){
      for(i in 1:num.samples){
        f.all[j, , i] = N.standardize(t.sub, rep(unlist(f.log[[j]][i, ]), group.size[[j]][i, ]), time.grid.recon)
      }
    }
  }

  if (tree.prior == "skygrid"){
    f.all = array(0, dim = c(K, M, num.samples))
    # read in the log and trees
    for(i in 1:K){
      tmp = read.delim(file=logpath[[i]], skip = 4, header = TRUE, nrows = num.samples)
      f.all[i, , ] = t(tmp[, 8:(7 + M)] + log(fct))
    }
  }
  if(!coal.return){
    return(f.all)
  }
  else
    return(list(f.all, t.all, cutoff))
}
