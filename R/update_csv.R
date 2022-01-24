library(quantmod)
library(xts)

# get old and new data
getSymbols("SPY")

oldDataExists <- file.exists("data/SPY.csv")
if( oldDataExists ){
  cat("old data exists...reading it in\n")
  oldSPY <- read.csv.zoo("data/SPY.csv")
  oldSPY <- as.xts(oldSPY)
}else{
  cat("no old data exists...starting fresh\n")
}

# if old file exists 
# create backup and rewrite
# otherwise write the only data you got 
if( oldDataExists ){
  cat("backing up old data and writing out extended data\n")
  write.zoo(oldSPY, file = "data/SPY.csv.bck", quote = F, sep = ",")
  biggerData <- rbind(oldSPY, SPY)
  biggerData <- biggerData[!duplicated(time(biggerData))]
  write.zoo(biggerData, file = "data/SPY.csv", quote = F, sep = ",")
}else{
  cat("writing new data to file\n")
  write.zoo(SPY, file = "data/SPY.csv", quote = F, sep = ",")
}

