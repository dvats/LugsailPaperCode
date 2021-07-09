set.seed(10)
library(mcmcse)
library(ellipse)
load("probit_data")
load("probit_true_mean_8")
library(truncnorm)

## PX-DA for Bayesian probit regression
probit_gibbs <- function(dat, nsim)
{
  # Separate covariates and response
  x <- dat[, -4]
  y <- dat[, 4]
  m <- length(y)

  #return objects
  betas <- matrix(0, nrow = nsim, ncol = 3)
  z <- matrix(0, nrow = nsim, ncol = m)
  g2 <- numeric(length = nsim)
  # starting value of beta are the mle estaimtes
  z.curr <- numeric(length = m)
  fit <- glm(y ~ -1 + x, family = binomial(link = "probit") )
  beta.curr <- coef(fit)

  # Needed for the truncated normal
  lower <- numeric(length = m)
  upper <- numeric(length = m)
  for(i in 1:m)
  {
      if(y[i] == 0) 
      {
        lower[i] <- -Inf
        upper[i] <- 0
      }else
      {
        lower[i] <- 0
        upper[i] <- Inf
      }
  }
  
  # Doing non-repetitive calculations outside
  xtx <- t(x)%*%x
  svd.xtx <- svd(xtx)
  inv.xtx <- svd.xtx$v%*%diag(1/svd.xtx$d, nrow = 3)%*%t(svd.xtx$u)
  sq.inv.xtx <- svd.xtx$v%*%diag( (svd.xtx$d)^(-1/2), nrow = 3)%*%t(svd.xtx$u)
  
  for(i in 1:nsim)
  {
    z.curr <- rtruncnorm(m, a = lower, b = upper, mean = x%*%beta.curr, sd = 1)
    
    g2.curr <- rgamma(1, shape = m/2, 
                 rate = 0.5*sum((z.curr - x %*% inv.xtx %*% t(x) %*% z.curr)^2))
    
    beta.curr <- inv.xtx %*% t(x) %*% (sqrt(g2.curr)*z.curr) + sq.inv.xtx %*% rnorm(3)
    
    betas[i, ] <- beta.curr
    z[i, ] <- z.curr
    g2[i] <- g2.curr
  }
  return(list("beta" = betas, "z" = z, "gsquare" = g2))
}

adjust_matrix <- function(mat, N, epsilon = sqrt(log(N)/dim(mat)[2]), b = 1/2)
{
  mat.adj <- mat
  adj <- epsilon*N^(-b)
  vars <- diag(mat)
  corr <- cov2cor(mat)
  eig <- eigen(corr)
  adj.eigs <- pmax(eig$values, adj)
  mat.adj <- diag(vars^(1/2))%*% eig$vectors %*% diag(adj.eigs) %*% t(eig$vectors) %*% diag(vars^(1/2))
  return(mat.adj)
}


N <- 1e4
level <- .90
out <-  probit_gibbs(dat = dat, nsim = N)
chain <- cbind(out$beta)
beta.mean <- colMeans(chain)

p <- length(beta.mean)

is <- mcse.initseq(chain)

size <- "sqroot"
foo <- mcse.multi(chain, size = size, large = TRUE)
bn <- foo$size
bm <- foo$cov

wbm <- 2*foo$cov - mcse.multi(chain, size = ceiling(bn/2))$cov
wbm.adj <- adjust_matrix(wbm, N = N)

lug3 <- 2*foo$cov - mcse.multi(chain, size = ceiling(bn/3))$cov
lug3.adj <- adjust_matrix(lug3, N = N)

lug4 <- 2*foo$cov - mcse.multi(chain, size = ceiling(bn/4))$cov
lug4.adj <- adjust_matrix(lug4, N = N)



large.crit <- qchisq(level, df = 2)/N

pdf("probit_conf.pdf", height = 10, width = 10)
# par(mfrow = c(1,2))
plot(ellipse(lug3.adj, centre = beta.mean[1:2], which = c(1,2), t = sqrt(large.crit)), 
  col = "red", type= 'l', lty = 1, cex.main = 1.5, main = expression(paste("90 % Confidence regions at sample size ", 10^4)), xlab = expression(paste(beta[1])), ylab = expression(paste(beta[2])), lwd = 2)
lines(ellipse(wbm.adj, centre = beta.mean[1:2], which = c(1,2), t = sqrt(large.crit)), col = "blue", lty = 1, lwd = 1.3)
# lines(confRegion(is, which = c(1,2), level = .90), lty = 4, col = "brown", lwd = 1.3)
lines(confRegion(foo, which = c(1,2), level = .90))
legend("topright", c("BM (r=1)", "WBM (r = 2)", "Lugsail BM (r=3)"), cex = 1.5,  col = c("black", "blue", "red"), lty = 1, lwd = 2)
dev.off()


# N <- 1e6
# size = "sqroot"
# level <- .90
# out <-  probit_gibbs(dat = dat, nsim = N)
# chain <- cbind(out$beta)
# beta.mean <- colMeans(chain)

# p <- length(beta.mean)

# foo <- mcse.multi(chain, size = size, large = TRUE)
# bn <- foo$size
# bm <- foo$cov

# wbm <- 2*foo$cov - mcse.multi(chain, size = ceiling(bn/2))$cov
# wbm.adj <- adjust_matrix(wbm, N = N)

# lug3 <- 2*foo$cov - mcse.multi(chain, size = ceiling(bn/3))$cov
# lug3.adj <- adjust_matrix(lug3, N = N)

# lug4 <- 2*foo$cov - mcse.multi(chain, size = ceiling(bn/4))$cov
# lug4.adj <- adjust_matrix(lug4, N = N)

# is <- mcse.initseq(chain)

# large.crit <- qchisq(level, df = p)/N


# plot(ellipse(lug3.adj, centre = beta.mean[1:2], which = c(1,2), t = sqrt(large.crit)), 
#   col = "red", type= 'l', lty = 3,  main = expression(paste("Confidence regions at sample size ", 10^6)), xlab = expression(paste(beta[1])), ylab = expression(paste(beta[2])))
# lines(ellipse(wbm.adj, centre = beta.mean[1:2], which = c(1,2), t = sqrt(large.crit)), col = "blue", lty = 2)
# lines(confRegion(is, which = c(1,2)), lty = 4, col = "brown", lwd = 1)
# lines(confRegion(foo, which = c(1,2)))




