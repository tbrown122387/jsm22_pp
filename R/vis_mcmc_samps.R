# visualize MCMC samples in data/param_samples.csv
#devtools::install_github("tbrown122387/mmcmc")
library(ggplot2)
library(mmcmc)
library(GGally)

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
