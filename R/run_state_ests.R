## Run Modes
# * 1. Run the Liu-West1 (original auxiliary style) filter for state output (sampling parameters from prior)
# * 4. Run the Liu-West1 (original auxiliary style) filter for state output (sampling parameters from csv)
#   * 7. Run the Liu-West2 filter for state output (sampling parameters from prior)
# * 10. Run the Liu-West2 filter for state output (sampling parameters from csv)
#   * 13. Run the Particle Swarm (bootstrap filters) algorithm for state output (sampling parameters from prior)
# * 16. Run the Particle Swarm (bootstrap filters) algorithm for state output (sampling parameters from csv)
#   * 19. Run standard auxiliary particle filter for state output (with given parameter estimates).
# *


# TODO make this happen in Makefile
# have to change directory because some filepaths are hardcoded relative style
setwd("~/jsm22_pp/")

# construct command line input strings
prog <- "./cpp/cmake-build-release/jsmpp_v2"
runModes <- c(1,4,7,10,13,16,19) # see above
dataPath <- "./data/SPY_returns.csv"
outFiles <- c("lw_aux_prior.txt", "lw_aux_csv.txt", "lw2_prior.txt",
              "lw2_csv.txt", "swarm_prior.txt", "swarm_csv.txt", "pf_est.txt")
outFiles <- paste("data/state_estimates/",outFiles, sep ="")

# run all the commands
for(idx in seq_along(runModes)){
  cmd <- paste(prog, runModes[idx], dataPath, '>', outFiles[idx])
  cat(idx, ": ", cmd, "\n")
  system(cmd)
}
