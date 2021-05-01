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
process.beast.logs = function(logpath, treepath = NULL, M = 15, time.offset = NULL, fct = 1, skip = 3, tree.prior = "skyline", num.samples = 1000, recon.len = 100, coal.return = FALSE){
  K = length(logpath)
  if(is.null(time.offset))
    time.offset = rep(0, K)
  if (tree.prior %in% c("skyline", "skyride")){
    f.log = list()
    t.all = list()
    group.size = list()
    treefiles = list()
    # read in the log and trees
    for(i in 1:K){
      tmp = read.table(file=logpath[[i]], skip = skip, header = TRUE)
      if (tree.prior == "skyride"){
        pop.idx = grep("*logPopSize*", names(tmp))
      }

      if (tree.prior == "skyline"){
        pop.idx = grep("*popSize*", names(tmp))
      }
      tot.len = dim(tmp)[1]
      if (num.samples > tot.len)
        num.samples = tot.len
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
      tmp = ape::read.nexus(treepath[[i]])
      treefiles[[i]] = tmp[(length(tmp) - num.samples + 1) : length(tmp)]
      t.all[[i]] = matrix(0, treefiles[[i]][[1]]$Nnode, num.samples)
    }

    f.all = array(0, dim = c(K, recon.len, num.samples))
    tmp = c()
    for(j in 1:K){
      for(i in 1:num.samples){
        t.tmp = phylodyn::summarize_phylo(treefiles[[j]][[i]])
        t.sub = t.tmp$coal_times + time.offset[j]
        # t.sub = cumsum(ape::coalescent.intervals(treefiles[[j]][[i]])$interval.length) * fct

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

#' plot function
#'
#' @param res.skydive output of combine function
#' @param fig.title title of the figure
#' @return a ggplot2 figure
#' @export
#'
plot_skydivide = function(res.skydive, fig.title = "", ylim = c(0.1, 1e5)){
  time.grid.recon = res.skydive[[1]]
  N.quant.debiased = res.skydive[[2]]
  # N.quant.debiased = data.frame(med = apply(N.recon, 2, quantile, prob = 0.5),
  #                               lo = apply(N.recon, 2, quantile, prob = 0.025),
  #                               hi = apply(N.recon, 2, quantile, prob = 0.975, na.rm = TRUE))

  fig = ggplot2::ggplot(N.quant.debiased, ggplot2::aes(x = time.grid.recon, y = med)) +
    ggplot2::geom_line(size = 2) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), alpha = 0.2) +
    ggplot2::scale_y_log10() +
    ggplot2::scale_x_reverse() +
    ggplot2::ylab("N(t)") +
    ggplot2::xlab("Time to present") +
    ggplot2::ggtitle(fig.title) +
    ggplot2::theme(legend.position="none",
          axis.text=ggplot2::element_text(size = 20),
          axis.title=ggplot2::element_text(size = 20),
          plot.title = ggplot2::element_text(size = 20)) +
    ggplot2::coord_cartesian(ylim = ylim)
  return(fig)
}


