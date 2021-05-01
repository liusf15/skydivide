test_that("algorithm works", {
  skip("skip test")
  n = 50
  mu = 10
  K = 2
  traj_list = c("const")#, "bottleneck", "exp")
  seed_list = seq(1, 1)
  for(traj in traj_list){
    for(seed in seed_list){
      cat("Analyzing:", traj, seed, "\n")
      logpath = list()
      treepath = list()
      traj = "const"
      for(k in 1:K){
        logpath[[k]] = sprintf("/Users/sifanliu/Dropbox/Projects/Scalable Bayesian/scalable_mcmc/Rcode/simulations/skyline/samp_res/samp_n_%d_mu_%d_K_%d_seed_%d_%ssub%d.log", n, mu, K, seed, traj, k)
        treepath[[k]] = sprintf("/Users/sifanliu/Dropbox/Projects/Scalable Bayesian/scalable_mcmc/Rcode/simulations/skyline/samp_res/samp_n_%d_mu_%d_K_%d_seed_%d_%ssub%d.trees", n, mu, K, seed, traj, k)
      }
      res = combine(logpath, treepath, "debias",  time.offset = rep(0, 2))
      fig = plot_skydivide(res, ylim = c(0.01, 10))
    }
  }
})
