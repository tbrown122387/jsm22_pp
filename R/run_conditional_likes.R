# /**
#   * RUN MODES
# * 2. Run the Liu-West1 (original auxiliary style) filter for conditional likelihoods (sampling parameters from prior)
# * 5. Run the Liu-West1 (original auxiliary style) filter for conditional likelihoods (sampling parameters from csv)
# * 8. Run the Liu-West2 filter for conditional likelihoods (sampling parameters from prior)
# * 11. Run the Liu-West2 filter for conditional likelihoods (sampling parameters from csv)
# * 14. Run the Particle Swarm (bootstrap filters) algorithm for conditional likelihoods (sampling parameters from prior)
# * 17. Run the Particle Swarm (bootstrap filters) algorithm for conditional likelihoods (sampling parameters from csv)
# * 20. Run standard auxiliary particle filter for conditional likelihoods (with given parameter estimates).
# */


# have to change directory because some filepaths are hardcoded relative style
setwd("~/jsm22_pp/")
source("R/create_configs.R")
setwd("~/jsm22_pp/")

prog <- "./cpp/cmake-build-release/jsmpp_v2"
runModes <- c(2,5,8,11,14,17,20) # see above
dataPath <- "./data/SPY_returns.csv"
outFiles <- c("lw_aux_prior.txt", "lw_aux_csv.txt", "lw2_prior.txt",
              "lw2_csv.txt", "swarm_prior.txt", "swarm_csv.txt", "pf_est.txt")
outFiles <- paste("data/cond_likes/",outFiles, sep ="")

for(idx in seq_along(runModes)){
  cmd <- paste(prog, runModes[idx], dataPath, '>', outFiles[idx])
  cat(idx, ": ", cmd, "\n")
  system(cmd)
}
