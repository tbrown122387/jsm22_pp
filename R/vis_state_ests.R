# /**
#   * RUN MODES
# * 1. Run the Liu-West1 (original auxiliary style) filter for state output (sampling parameters from prior)
# * 4. Run the Liu-West1 (original auxiliary style) filter for state output (sampling parameters from csv)
#   * 7. Run the Liu-West2 filter for state output (sampling parameters from prior)
# * 10. Run the Liu-West2 filter for state output (sampling parameters from csv)
#   * 13. Run the Particle Swarm (bootstrap filters) algorithm for state output (sampling parameters from prior)
# * 16. Run the Particle Swarm (bootstrap filters) algorithm for state output (sampling parameters from csv)
#   * 19. Run standard auxiliary particle filter for state output (with given parameter estimates).
# *
library(ggplot2)
library(reshape2)

# have to change directory because some filepaths are hardcoded relative style
setwd("~/jsm22_pp/")
cNames <- c("lw_aux_prior", "lw_aux_csv", "lw2_prior","lw2_csv", "swarm_prior", "swarm_csv", "pf_est")
outFiles <- paste0(cNames, ".txt")
outFiles <- paste("data/state_estimates/",outFiles, sep ="")
allOutput <- as.data.frame(lapply(outFiles, read.csv, header=F))
#allOutput <- cumsum(allOutput)
colnames(allOutput) <- cNames
allOutput$day <- seq_along(allOutput[,1])
allOutput <- allOutput[1:50,]

meltdf <- melt(allOutput,id.var="day")
ggplot(meltdf,aes(x=day,y=value,colour=variable)) + 
  geom_line() 
# TODO add vertical bar for out of sample
