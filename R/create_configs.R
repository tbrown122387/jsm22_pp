# changes made here will be propogated to many files
phiLow <- 0.864248
phiHigh <- 0.972859
muLow <- -0.390150  
muHigh <- 0.551337
sigmaLow <- 0.2021094
sigmaHigh <- 0.4718427 
rhoLow <- -0.955950
rhoHigh <- -0.522117
dte <- 5
delta <- .99
paramSamplesFile <- "data/posterior_samps/param_samples.csv"
estDataFile <- "data/SPY_returns_estimation.csv"
numMCMCIters <- 100000
numPFs <- 7
burn <- 100

setwd("~/jsm22_pp/data/")

## check your highs and lows are consistent with the data 
d <- read.csv("posterior_samps/param_samples.csv", header=F)
d[,3] <- sqrt(d[,3])
parEstimates <- colMeans(d[(burn+1):nrow(d),])
phiGood <- (phiLow < parEstimates[1]) & (parEstimates[1] < phiHigh)
muGood <- (muLow < parEstimates[2]) & (parEstimates[2] < muHigh)
sigmaGood <- (sigmaLow < parEstimates[3]) & (parEstimates[3] < sigmaHigh)
rhoGood <- (rhoLow < parEstimates[4]) & (parEstimates[4] < rhoHigh)
stopifnot(phiGood & muGood & sigmaGood & rhoGood)

# config1.csv
# for run-modes 1-3, 7-9
# delta, phi_l, phi_u, mu_l, mu_u, sig_l, sig_u, rho_l, rho_u, dte
options("scipen"=10)    # set high penalty for scientific display
myStr <- paste(delta, phiLow, phiHigh, muLow, muHigh, sigmaLow, sigmaHigh, rhoLow, rhoHigh, dte, sep = ", ")
write.table(myStr, "../configs/config1.csv", append = F, quote = F, col.names = F, row.names = F)

# config 2
# delta, param_samples_filename, dte
myStr <- paste(delta, paramSamplesFile, dte, sep = ",")
write.table(myStr, "../configs/config2.csv", append = F, quote = F, col.names = F, row.names = F)

# config 3
# phi_l, phi_u, mu_l, mu_u, sig_l, sig_u, rho_l, rho_u, dte;
myStr <- paste(phiLow, phiHigh, muLow, muHigh, sigmaLow, sigmaHigh, rhoLow, rhoHigh, dte, sep = ", ")
write.table(myStr, "../configs/config3.csv", append = F, quote = F, col.names = F, row.names = F)

# config 4
# param_samples_filename, dte
myStr <- paste(paramSamplesFile, dte, sep = ",")
write.table(myStr, "../configs/config4.csv", append = F, quote = F, col.names = F, row.names = F)

# config 5
# param_samples_filename, dte
myStr <- paste(c(parEstimates, dte), collapse = ", ")
write.table(myStr, "../configs/config5.csv", append = F, quote = F, col.names = F, row.names = F)

# config 6
# num_mcmc_iters, num_particle_filters
myStr <- paste(numMCMCIters, numPFs, estDataFile, sep = ",")
write.table(myStr, "../configs/config6.csv", append = F, quote = F, col.names = F, row.names = F)
