## File and output not used
library(MCMCpack)
dat <- read.csv("framingham.csv")
dat$education <- as.factor(dat$education)
dat$currentSmoker <- as.factor(dat$currentSmoker)
dat$BPMeds <- as.factor(dat$BPMeds)

logitchain <- function(seed, N, beta.start = NA)
{
  chain <- MCMClogit(formula = TenYearCHD ~ 1 + . , data = dat, beta.start = beta.start, verbose = 0, seed = seed, burnin = 0, mcmc = N, thin = 1, B0 = 1/10, tune = .5)
  return(chain)
}

ptm <- proc.time()
p <- 18
B <- 100
N <- 1e6
means <- matrix(0, nrow = B, ncol = p)
for(b in 1:B)
{
  start <- rnorm(p, 0, sd = sqrt(10))
  out <- logitchain(seed = b, N = N, beta.start = NA) 
  means[b, ] <- colMeans(out)
}

proc.time() - ptm
out.means <- colMeans(means)
names(out.means) <- colnames(out)
save(out.means, file = "fram_true_mean_1e6100")
