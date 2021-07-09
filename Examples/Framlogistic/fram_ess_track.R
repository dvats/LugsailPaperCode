set.seed(100)
library(MCMCpack)
library(mcmcse)
load("fram_true_mean_1e6100")
dat <- read.csv("framingham.csv")
dat$education <- as.factor(dat$education)
dat$currentSmoker <- as.factor(dat$currentSmoker)
dat$BPMeds <- as.factor(dat$BPMeds)

adjust_matrix <- function(mat, N, epsilon = sqrt(log(N)/dim(mat)[2]), b = 9/10)
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

logitchain <- function(seed, N, beta.start = NA)
{
  chain <- MCMClogit(formula = TenYearCHD ~ 1 + . , data = dat, beta.start = beta.start, verbose = 0, seed = seed, burnin = 0, mcmc = N, thin = 1, B0 = 1/10, tune = .5)
  return(chain)
}



fram_track <- function(N = 1e5, dat, B = 100, level = .90, size = "sqroot")
{
	nloops <- 100
	subsize <- seq(5e3, N, length = nloops)
	ess_track_lug <- matrix(0, nrow = B, ncol = length(subsize))
	ess_track_bm <- matrix(0, nrow = B, ncol = length(subsize))
	ess_track_wbm <- matrix(0, nrow = B, ncol = length(subsize))
	ess_track_jack1 <- matrix(0, nrow = B, ncol = length(subsize))
	ess_track_jack2 <- matrix(0, nrow = B, ncol = length(subsize))

	colnames(ess_track_wbm) <- subsize
	colnames(ess_track_bm) <- subsize
	colnames(ess_track_lug) <- subsize

	r.J <- 1/2
	for(b in 1:B)
	{
		if(b %% 1 ==0) print(b)
		out <-  MCMClogit(formula = TenYearCHD ~ 1 + . , data = dat, beta.start = NA, verbose = 0, seed = (100 + b), burnin = 0, mcmc = N, thin = 1, B0 = 1/10, tune = .5)
		chain <- out
		beta.mean <- colMeans(chain)

		p <- length(beta.mean)

		for(j in 1:length(subsize))
		{
			# print(j)
			minichain <- chain[1:subsize[j], ]
			ess_track_bm[b, j] <- multiESS(minichain, r = 1, method = "bm", size = size, adjust = TRUE)/subsize[j]
			ess_track_wbm[b,j] <- multiESS(minichain, r = 2, method = "bm", size = size, adjust = TRUE)/subsize[j]
			ess_track_lug[b,j] <- multiESS(minichain, r = 3, method = "bm", size = size, adjust = TRUE)/subsize[j]

		
			foo <-  mcse.multi(minichain, r = 1, size = size, blather = TRUE)
			bm <- foo$cov
			bn <- foo$size
			a <- floor(N/bn)
			alpha.J1 <- ( a + r.J)/(a*(1 - r.J))
			alpha.J2 <- ( log(a) + r.J)/(log(a)*(1 - r.J))

			jack1.adj <- (alpha.J1)*bm + (1 - alpha.J1)*mcse.multi(minichain, size = floor(bn/2), r = 1)$cov
			jack2.adj <- (alpha.J2)*bm + (1 - alpha.J2)*mcse.multi(minichain, size = floor(bn/2), r = 1)$cov

			sig.eigen <- eigen(jack1.adj, only.values = TRUE)$values
			if (min(sig.eigen) <= 0) {
			    jack1.adj <- adjust_matrix(jack1.adj, N = N)
			}

			sig.eigen <- eigen(jack2.adj, only.values = TRUE)$values
			if (min(sig.eigen) <= 0) {
			    jack2.adj <- adjust_matrix(jack2.adj, N = N)
			}

			ess_track_jack1[b,j] <- multiESS(minichain, covmat = jack1.adj)/subsize[j]
			ess_track_jack2[b,j] <- multiESS(minichain, covmat = jack2.adj)/subsize[j]
	
		}
	}
	return(list("BM" = ess_track_bm, "WBM" = ess_track_wbm, "LUG" = ess_track_lug, "Jack_lin" = ess_track_jack1, "Jack_log" = ess_track_jack2))
}

ess_track <- fram_track(N = 1e6, dat = dat, B = 100, size = "sqroot")

nloops <- 100
subsize <- seq(5e3, 1e6, length = nloops)
# plot(subsize, colMeans(foo$BM), type= 'l')
# lines(subsize, colMeans(foo$WBM), col = "red")
# lines(subsize, colMeans(foo$LUG), col = "blue")

save(subsize, ess_track, file = "fram_ess6_sq")

