# visualize MCMC samples in data/param_samples.csv
library(ggplot2)
library(GGally)



#' A Gelman-Rubin convergence diagnostic function.
#'
#' This function calculates the Gelman-Rubin convergence diagnostic (Rhat) on some samples.
#' It assumes each chain is in a separate csv file, and each file contains the samples for
#' all parameters. For each chain, this function splits each column into two, which means
#' "m" comes out to be twice the number of csv files. Finally, this function assumes each
#' csv has a one-line header.
#' @param chain_files vector of string paths for each file of samples
#' @param burn the number of rows you want to discard as burn in. Defaults to 0.
#' @param ... arguments to be passed in to lapply
#' @keywords rhat Gelman-Rubin
#' @export
#' @examples
#' files <- c('/fake/path/samps1.csv', '/fake/path/samps2.csv')
#' rhats(files, burn = 300)
rhats <- function(chain_files, burn = 0, ...){
  # TODO: check all files are hte same shape
  # TODO: print the first couple of rows to make sure they look ok

  if(burn > 0){
    dfs <- lapply(chain_files,
                  function(path) read.csv(path, header=T)[-1:-burn,],
                  ...)
  }else{
    dfs <- lapply(chain_files,
                  function(path) read.csv(path, header=T),
                  ...)
  }
  only_one_column <- is.null(nrow(dfs[[1]]))
  m <- 2*length(chain_files)
  if(only_one_column){
    n <- length(dfs[[1]])
  }else {
    n <- nrow(dfs[[1]])
  }
  n <- (n %/% 2)*2 # in case n is odd
  if(only_one_column){
    dfs_halved <- lapply(dfs, function(df) df[((n/2+1):n)])
    dfs_halved <- c(dfs_halved, lapply(dfs, function(df) df[(1:(n/2))]))
  }else {
    dfs_halved <- lapply(dfs, function(df) df[((n/2+1):n),])
    dfs_halved <- c(dfs_halved, lapply(dfs, function(df) df[(1:(n/2)),]))
  }
  n <- n/2
  if(only_one_column){
    num_params <- 1
  }else {
    num_params <- ncol(dfs_halved[[1]])
  }
  rm(dfs)

  # calculate all rhat statistics for each untransformed univariate parameter
  rhats <- vector(mode="numeric", length = num_params)
  for(param_num in 1:num_params){
    if(only_one_column)
      univariate_chains <- dfs_halved
    else
        univariate_chains <- sapply(dfs_halved, '[', param_num) # pick out the same col from each df
    psi_bar_dot_js <- sapply(univariate_chains, mean)
    psi_bar_dot_dot <- mean(psi_bar_dot_js)
    ss_j <- sapply(univariate_chains, var)
    w <- mean(ss_j)
    b <- var(psi_bar_dot_js)*n
    varhat_plus <- w*(n-1)/n + b/n
    rhats[param_num] <- sqrt(varhat_plus/w)
  }
  return(rhats)
}

#' A acceptance rate function.
#'
#' This function calculates the acceptance rate for samples stored in a data frame.
#' @param df A data frame of samples. Rows correspond with iterations, and columns with parameters.
#' @param burn Number of rows you want to discard as burn in. Defaults to 0.
#' @keywords accept acceptance rate
#' @export
#' @examples
#' fake_samples <- data.frame(rnorm(100))
#' acceptRate(fake_samples)
acceptRate <- function(df, burn=0){
  start <- burn + 1
  mean(abs(diff(df[start:nrow(df),1])) > .000000001)
}


setwd("~/jsm22_pp/")
burn <- 100

# 1. read in data
d <- read.csv("data/param_samples.csv", header=F)
d[,3] <- sqrt(d[,3])
colnames(d) <- c("phi","mu","sigma","rho")
d$iter <- 1:nrow(d)
d <- d[-(1:burn),]

# 2. pairwise scatterplot
pdf("plots/mcmc_vis/pairwise_scatterplot.pdf")
ggpairs(subset(d, select=-iter), aes(alpha = 0.4))
dev.off()

# 3. histograms
pdf("plots/mcmc_vis/phi_hist.pdf")
ggplot(d, aes(x=phi)) + 
  geom_histogram() + 
  labs(title='\u03C6 Samples', x='\u03C6 value')
dev.off()

pdf("plots/mcmc_vis/mu_hist.pdf")
ggplot(d, aes(x=mu)) + 
  geom_histogram() + 
  labs(title='\u03BC Samples', x='\u03BC value')
dev.off()

pdf("plots/mcmc_vis/sigma_hist.pdf")
ggplot(d, aes(x=sigma)) + 
  geom_histogram() + 
  labs(title='\u03C3 Samples', x='\u03C3 value')
dev.off()

pdf("plots/mcmc_vis/rho_hist.pdf")
ggplot(d, aes(x=rho)) + 
  geom_histogram() + 
  labs(title='\u03C1 Samples', x='\u03C1 value')
dev.off()

# 4. trace plots
pdf("plots/mcmc_vis/phi_trace.pdf")
ggplot(d, aes(x=iter, y=phi)) +  
  geom_line() + 
  labs(title='\u03C6 Trace plot',
       y='\u03C6 value',
       x='iteration')
dev.off()

pdf("plots/mcmc_vis/mu_trace.pdf")
ggplot(d, aes(x=iter, y=mu)) +  
  geom_line() + 
  labs(title='\u03BC Trace plot', 
       y='\u03BC value',
       x='iteration')
dev.off()

pdf("plots/mcmc_vis/sigma_trace.pdf")
ggplot(d, aes(x=iter, y=sigma)) +  
  geom_line() +
  labs(title='\u03C3 Trace plot', 
       y='\u03C3 value',
       x='iteration')
dev.off()

pdf("plots/mcmc_vis/rho_trace.pdf")
ggplot(d, aes(x=iter, y=rho)) +  
  geom_line() +
  labs(title='\u03C1 Trace plot', 
       y='\u03C1 value',
       x='iteration')
dev.off()

# 5. correlograms/acf plots
pdf("plots/mcmc_vis/phi_acf.pdf")
acfData <- acf(d[,1], plot=F, lag.max = 100)
acfData <- data.frame(lag = acfData$lag, acf = acfData$acf)
ggplot(data = acfData, 
       mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) +
  geom_segment(mapping = aes(xend = lag, yend = 0)) +
  labs(title='\u03C6 correlogram',
       y=' estimated autocorrelation',
       x='lag')
dev.off()

pdf("plots/mcmc_vis/mu_acf.pdf")
acfData <- acf(d[,2], plot=F, lag.max = 100)
acfData <- data.frame(lag = acfData$lag, acf = acfData$acf)
ggplot(data = acfData, 
       mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) +
  geom_segment(mapping = aes(xend = lag, yend = 0)) +
  labs(title='\u03BC correlogram',
       y=' estimated autocorrelation',
       x='lag')
dev.off()

pdf("plots/mcmc_vis/sigma_acf.pdf")
acfData <- acf(d[,3], plot=F, lag.max = 100)
acfData <- data.frame(lag = acfData$lag, acf = acfData$acf)
ggplot(data = acfData, 
       mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) +
  geom_segment(mapping = aes(xend = lag, yend = 0)) +
  labs(title='\u03C3 correlogram',
       y=' estimated autocorrelation',
       x='lag')
dev.off()

pdf("plots/mcmc_vis/rho_acf.pdf")
acfData <- acf(d[,4], plot=F, lag.max = 100)
acfData <- data.frame(lag = acfData$lag, acf = acfData$acf)
ggplot(data = acfData, 
       mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) +
  geom_segment(mapping = aes(xend = lag, yend = 0)) +
  labs(title='\u03C1 correlogram',
       y=' estimated autocorrelation',
       x='lag')
dev.off()

# write out numerical statistics
acceptRate <- acceptRate(d, burn = burn)
acceptRateMessage <- paste("accept rate: ", acceptRate)
rhatsMessage <- paste("rhats: ", rhats("data/param_samples.csv"))
fileConn<-file("data/mcmc_numerical_diagnostics.txt")
writeLines(acceptRateMessage, fileConn)
writeLines(rhatsMessage, fileConn)
close(fileConn)
