## STEPS
# 1. backup param_samples.csv
# 2. run mcmc
# 3. replace param_samples with most recent samples_* file

# 1. backup param_samples.csv
setwd("~/jsm22_pp/")
system("cp data/param_samples.csv data/param_samples.csv.bck")

# 2. run mcmc
# (be sure you're happy with all the hardcoded stuff in that file, and that it's all compiled)
startTime <- Sys.time()
myCmd <- './cpp/cmake-build-release/jsmpp_v2 22 "data/SPY_returns.csv"'
system(myCmd)
endTime <- Sys.time()

# 3.replace param_samples with most recent samples_* file
allSamplesFiles <- list.files()[grep(pattern = "samples", list.files())]
asDateTimes <- as.POSIXct(substr(allSamplesFiles, 9, 100), format="%Y-%m-%d.%H-%M-%OS")
mostRecent <- allSamplesFiles[which.max(asDateTimes)]
mySecondCmd <- paste("cp", mostRecent, "data/param_samples.csv")
system(mySecondCmd)

# 4. write out time taken
timeTaken <- endTime - startTime
timeMessage <- paste("the amount of time required to generate the mcmc samples was", timeTaken, "seconds")
fileConn<-file("mcmc_time_taken.txt")
writeLines(timeMessage, fileConn)
close(fileConn)
