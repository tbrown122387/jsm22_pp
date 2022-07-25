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
library(zoo)
library(ggtern)

## have to change directory because some filepaths are hardcoded relative style
setwd("~/jsm22_pp/")
cNames <- c("lw_aux_prior", "lw_aux_csv", "lw2_prior","lw2_csv", "swarm_prior", "swarm_csv", "pf_est")
outFiles <- paste0(cNames, ".txt")
outFiles <- paste("data/cond_likes/",outFiles, sep ="")
allOutput <- as.data.frame(lapply(outFiles, read.csv, header=F))

## remove filtering output on data on and before training data
# warning: the following variable is based off of something 
# that's hardcoded in prices_to_returns.R e.g. trainStop <- "2010-12-31" 
testPriceStart <- "2011-01-01"
numDaysToDisregard <- sum(index(read.csv.zoo("data/SPY.csv", header=T)) < testPriceStart)
numCondLikesToDisregard <- numDaysToDisregard - 2 
# ^ one for conversion from prices to returns, and one from the lagged predictor
allOutput <- allOutput[-(1:numCondLikesToDisregard),]
colnames(allOutput) <- cNames
allOutput$time <- seq_along(allOutput[,1])

## create a cumulative sum data set
allCumsumOutput <- cumsum(allOutput[,colnames(allOutput) != 'time'])
allCumsumOutput$time <- seq_along(allCumsumOutput[,1])


pdf("plots/cond_likes_vis/clike_vis.pdf")
mdf <- melt(allOutput, id.vars = "time")
ggplot(mdf,                            # Draw ggplot2 time series plot
       aes(x = time,
           y = value,
           col = variable)) +
  geom_line()
dev.off()

pdf("plots/cond_likes_vis/cumsum_clike_vis.pdf")
cmdf <- melt(allCumsumOutput, id.vars = "time")
ggplot(cmdf,                            # Draw ggplot2 time series plot
       aes(x = time,
           y = value,
           col = variable)) +
  geom_line()
dev.off()



pdf("plots/cond_likes_vis/clike_ternary_vis_uniforms.pdf")
uniformPriorLLs <- data.frame(allOutput$lw_aux_prior, allOutput$lw2_prior, allOutput$swarm_prior)
uniformPriorLLs <- exp(uniformPriorLLs)/rowSums(exp(uniformPriorLLs))
colnames(uniformPriorLLs) <- c("lw_aux", "lw_v2", "p_swarm")
time <- 1:nrow(allOutput)
ggtern(data=uniformPriorLLs,aes(lw_v2, lw_aux, p_swarm)) + 
  geom_mask() +
  geom_path(aes(colour = time, alpha = .01)) +
  theme_bw() +
  theme_showarrows()
dev.off()


pdf("plots/cond_likes_vis/clike_ternary_vis_csv.pdf")
csvSampsLLs <- data.frame(allOutput$lw_aux_csv, allOutput$lw2_csv, allOutput$swarm_csv)
csvSampsLLs <- exp(csvSampsLLs)/rowSums(exp(csvSampsLLs))
colnames(csvSampsLLs) <- c("lw_aux", "lw_v2", "p_swarm")
time <- 1:nrow(allOutput)
ggtern(data=csvSampsLLs,aes(lw_v2, lw_aux, p_swarm)) + 
  geom_mask() +
  geom_path(aes(colour = time, alpha = .01)) +
  theme_bw() +
  theme_showarrows()
dev.off()