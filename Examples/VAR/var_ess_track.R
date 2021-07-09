
source("var1_function.R")
set.seed(101)
library(mcmcse)
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

var_track <- function(N = 1e5, phi, omega, B = 100, level = .90, size = "sqroot")
{
	p <- dim(phi)[1]
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
		chain <- var1(p = p, phi = phi, nsim = N, omega = omega)

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

p <- 10
phi <- diag(rep(.95,p))
omega <- matrix(.9, nrow = p, ncol = p)
diag(omega) <- 1
omega <- omega^(abs(col(omega)-row(omega)))

truth <- true.sig(p = p, omega = omega, phi = phi)
ess_true <- (det(truth$tar.var)/det(truth$final.cov))^(1/p)
ess_track <- var_track(N = 1e6, phi = phi, omega = omega, B = 100, size = "sqroot")

nloops <- 100
subsize <- seq(5e3, 1e6, length = nloops)
# plot(subsize, colMeans(ess_track$BM), type= 'l', ylim = range(ess_track))
# lines(subsize, colMeans(ess_track$WBM), col = "red")
# lines(subsize, colMeans(ess_track$LUG), col = "blue")
# abline(h = ess_true)
save(subsize, ess_true, ess_track, file = "var_ess6_sq")

