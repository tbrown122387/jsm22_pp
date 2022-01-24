
setwd("~/jsm22_pp/")


library(xts)
SPY <- read.csv.zoo("data/SPY.csv")
returns <- diff(log(coredata(SPY$SPY.Adjusted)))*100
write.table(returns, file = "data/SPY_returns.csv", quote = F, col.names = F, row.names = F)
