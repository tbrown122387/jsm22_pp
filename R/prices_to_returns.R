
setwd("~/jsm22_pp/")
library(xts)

# write full data for (fast) forecasting algos
SPY <- read.csv.zoo("data/SPY.csv")
returns <- diff(log(coredata(SPY$SPY.Adjusted)))*100
write.table(returns, file = "data/SPY_returns.csv", quote = F, col.names = F, row.names = F)

# write smaller data set for (slow) mcmc estimation
start <- "2010-01-01"
stop <- "2010-12-31"
smallSPY <- SPY[(index(SPY) >= start) & (index(SPY) <= stop)]
smallReturns <- diff(log(coredata(smallSPY$SPY.Adjusted)))*100
write.table(smallReturns, 
            file = "data/SPY_returns_estimation.csv", 
            quote = F, col.names = F, row.names = F)

