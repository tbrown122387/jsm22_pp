# * 23. Store Liu-West1 (original auxiliary style) parameter posterior samples (sampling from noninformative priors at time 1) Run 100 times and stack results vertically
# * 24. Store Liu-West2 parameter posterior samples (sampling from noninformative prior at time 1) Run 100 times and stack results vertically
 
# TODO make this happen in Makefile
# have to change directory because some filepaths are hardcoded relative style
setwd("~/jsm22_pp/")

prog <- "./cpp/cmake-build-release/jsmpp_v2"
runModes <- c(23, 24) # see above
dataPath <- "./data/SPY_returns.csv" #doesn't matter because hardcoded
outFiles <- c("lw_aux_posterior.txt", "lw2_prior_posterior.txt")
outFiles <- paste("data/posterior_samps/",outFiles, sep ="")


for(idx in seq_along(runModes)){
  cmd <- paste(prog, runModes[idx], dataPath, '>', outFiles[idx])
  cat(idx, ": ", cmd, "\n")
  system(cmd)
}
