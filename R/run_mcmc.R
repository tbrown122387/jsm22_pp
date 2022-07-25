## STEPS
# 1. backup param_samples.csv
# 2. run mcmc
# 3. replace param_samples with most recent samples_* file


setwd("~/jsm22_pp/")

# 1. backup file if you need to
if( file.exists("data/posterior_samps/param_samples.csv") ){
    system("cp data/posterior_samps/param_samples.csv data/posterior_samps/param_samples.csv.bck")
}

# 2. run mcmc
# (be sure you're happy with all the hardcoded stuff in that 
# file, that it's all compiled, and that you're cool with config6)
# last argument in myCmd is ignored...it's really found in the config file
startTime <- Sys.time()
myCmd <- './cpp/cmake-build-release/jsmpp_v2 22 data/SPY_returns_estimation.csv'
system(myCmd)
endTime <- Sys.time()

# 3.replace param_samples with most recent samples_* file
allSamplesFiles <- list.files()[grep(pattern = "samples", list.files())]
asDateTimes <- as.POSIXct(substr(allSamplesFiles, 9, 100), format="%Y-%m-%d.%H-%M-%OS")
mostRecent <- allSamplesFiles[which.max(asDateTimes)]
mySecondCmd <- paste("cp", mostRecent, "data/posterior_samps/param_samples.csv")
system(mySecondCmd)

# 4. write out time taken
timeTaken <- endTime - startTime
timeMessage <- paste("the amount of time required to generate the mcmc samples was", timeTaken, "seconds")
fileConn<-file("data/posterior_samps/mcmc_time_taken.txt")
writeLines(timeMessage, fileConn)
close(fileConn)
