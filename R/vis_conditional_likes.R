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
library(ggplot2)
library(reshape2)

# have to change directory because some filepaths are hardcoded relative style
setwd("~/jsm22_pp/")
cNames <- c("lw_aux_prior", "lw_aux_csv", "lw2_prior","lw2_csv", "swarm_prior", "swarm_csv", "pf_est")
outFiles <- paste0(cNames, ".txt")
outFiles <- paste("data/cond_likes/",outFiles, sep ="")
allOutput <- as.data.frame(lapply(outFiles, read.csv, header=F))
colnames(allOutput) <- cNames
allOutput$time <- seq_along(allOutput[,1])
#allOutput <- allOutput[51:100,] # temporary line to help visualize stuff by pulling out random chunk


pdf("plots/cond_likes_vis/clike_vis.pdf")
mdf <- melt(allOutput, id.vars = "time")
ggplot(mdf,                            # Draw ggplot2 time series plot
       aes(x = time,
           y = value,
           col = variable)) +
  geom_line()
dev.off()
