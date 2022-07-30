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
library(zoo)

# have to change directory because some filepaths are hardcoded relative style
setwd("~/jsm22_pp/")
cNames <- c("lw_aux_prior", "lw_aux_csv", "lw2_prior","lw2_csv", "swarm_prior", "swarm_csv", "pf_est")
outFiles <- paste0(cNames, ".txt")
outFiles <- paste("data/state_estimates/",outFiles, sep ="")
allOutput <- as.data.frame(lapply(outFiles, read.csv, header=F))
returns <- read.csv("data/SPY_returns.csv", header=F)[,1]

## remove filtering output on data on and before training data
# warning: the following variable is based off of something 
# that's hardcoded in another file e.g. trainStop <- "2010-12-31" 
testPriceStart <- "2011-01-01"
numDaysToDisregard <- sum(index(read.csv.zoo("data/SPY.csv", header=T)) < testPriceStart)
numCondLikesToDisregard <- numDaysToDisregard - 2 
# ^ one for conversion from prices to returns, and one from the lagged predictor
allOutput <- allOutput[-(1:numCondLikesToDisregard),]
colnames(allOutput) <- cNames
allOutput$day <- seq_along(allOutput[,1])

# make the same adjustment on return data
numExtraReturns <- length(returns) - nrow(allOutput)
returns <- returns[-(1:numExtraReturns)]
returns <- data.frame('returns' = returns, 'day' = seq_along(returns))


# make vis
pdf("plots/state_vis/filter_vis.pdf")
meltdf <- melt(allOutput,id.var="day")
ggplot(meltdf,aes(x=day,y=value,colour=variable)) + 
  geom_line() +
  scale_y_continuous("log-volatility") 
dev.off()


pdf("plots/state_vis/test_returns.pdf")
ggplot(returns,aes(x=day,y=returns)) + 
  geom_line() 
  #geom_line(aes(y = a + returns*b), color = "red") +
  #scale_y_continuous("log-volatility", sec.axis = sec_axis(~ (. - a)/b, name = "returns")) +
  # ggtitle("returns in test data")
dev.off()
