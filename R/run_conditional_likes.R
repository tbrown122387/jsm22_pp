# /**
#   * RUN MODES
# * 2. Run the Liu-West1 (original auxiliary style) filter for conditional likelihoods (sampling parameters from prior)
# * 5. Run the Liu-West1 (original auxiliary style) filter for conditional likelihoods (sampling parameters from csv)
# * 8. Run the Liu-West2 filter for conditional likelihoods (sampling parameters from prior)
# * 11. Run the Liu-West2 filter for conditional likelihoods (sampling parameters from csv)
# * 14. Run the Particle Swarm (bootstrap filters) algorithm for conditional likelihoods (sampling parameters from prior)
# * 17. Run the Particle Swarm (bootstrap filters) algorithm for conditional likelihoods (sampling parameters from csv)
# * 20. Run standard bootstrap particle filter for conditional likelihoods (with given parameter estimates).
# */


# have to change directory because some filepaths are hardcoded relative style
setwd("~/jsm22_pp/")

prog <- "./cpp/cmake-build-release/jsmpp_v2"
runModes <- c(2,5,8,11,14,17,20) # see above
dataPath <- "./data/SPY_returns.csv"
outFiles <- c("lw_aux_prior", "lw_aux_csv", "lw2_prior",
              "lw2_csv", "swarm_prior", "swarm_csv", "pf_est")


# run for output
clOutFiles <- paste("data/cond_likes/",outFiles, ".txt", sep ="")
for(idx in seq_along(runModes)){
  cmd <- paste(prog, runModes[idx], dataPath, '>', clOutFiles[idx])
  cat(idx, ": ", cmd, "\n")
  system(cmd)
}

# run for timing
timeOutFiles <- paste("data/cond_likes/", outFiles, "_timing.txt", sep ="")
for(idx in seq_along(runModes)){
  cmd <- paste('{ time', prog, runModes[idx], dataPath, '> /dev/null; } 2>', timeOutFiles[idx])
  cat(idx, ": ", cmd, "\n")
  system(cmd)
}
