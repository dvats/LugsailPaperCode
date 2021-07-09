################################
## Testing SVE time series
################################
set.seed(1)
library(MASS)
library(mcmcse)
library(sandwich)
################################
###truth#####
################################
true.var <- function(u.rho = .70, x.rho = .90, p = 5, W.rho = .50, w = 1)
{
	foo <- matrix(W.rho, nrow= p, ncol=p)
	W <- foo^(abs(col(foo)-row(foo)))

	lam <- w/(1 - u.rho^2)
	Lam <- W/(1 - x.rho^2)

	truth <- lam*Lam + (Lam + t(Lam))*lam*(u.rho*x.rho/(1 - u.rho*x.rho))

	return(truth)
}

# adjust_matrix <- function(mat, N, epsilon = sqrt(log(N)/dim(mat)[2]), b = 1/2)
# {
#   mat.adj <- mat
#   adj <- epsilon*N^(-b)
#   vars <- diag(mat)
#   corr <- cov2cor(mat)
#   eig <- eigen(corr)
#   adj.eigs <- pmax(eig$values, adj)
#   mat.adj <- diag(vars^(1/2))%*% eig$vectors %*% diag(adj.eigs) %*% t(eig$vectors) %*% diag(vars^(1/2))
#   return(mat.adj)
# }



T.vec <- c(1e3, 1e4, 1e5)
p <- c(5, 10, 50)

u.rho <- c(.90)
x.rho <- .90
W.rho <- .50
w <- 1

B <- 10



# upp <- matrix(1, ncol = p, nrow = p)
# upp <- 1 - lower.tri(upp)
level <- .90
time_qs <- matrix(0, ncol = 9, nrow = B)
time_bt <- matrix(0, ncol = 9, nrow = B)
time_th <- matrix(0, ncol = 9, nrow = B)
time_bm <- matrix(0, ncol = 9, nrow = B)

for(i in 1:length(T.vec))
{
	track <- 0
	for(j in p)
	{
		track <- track + 1
		beta.star <- rnorm(j, 10, 4)
		foo <- matrix(W.rho, nrow= j, ncol=j)
		W <- foo^(abs(col(foo)-row(foo)))

		for(b in 1:B)
		{
			print(c(b,i,j))
			T <- T.vec[i]

			X <- matrix(0, nrow  = T, ncol = j)
			X[1,] <- rnorm(j)
			for(t in 2:T)
			{
				X[t, ] <- x.rho * X[t-1, ] + mvrnorm(1, mu = rep(0, j), Sigma = W)
			}

			u <- numeric(length = T)
			u[1] <- 0

			for(t in 2:T)
			{
				u[t] <- u.rho[1]*u[t-1] + rnorm(1, sd = sqrt(w))
			}

			y <- X%*%beta.star + u
			XtX <- t(X)%*%X
			beta.ols <- solve(XtX) %*% t(X)%*%y
			u.hat <- y - X%*%beta.ols
			v.hat <- apply(X, 2, function(x) x*u.hat)

			m <- lm(v.hat ~ 1)
			b_bart <- bwAndrews(x = m, kernel = "Bartlett", approx = "AR(1)", prewhite = 0)
			b_qs <- bwAndrews(x = m, kernel = "Quadratic Spectral", approx = "AR(1)", prewhite = 0)
			b_th <- bwAndrews(x = m, kernel = "Tukey", approx = "AR(1)", prewhite = 0)

			fix_b <- floor(sqrt(T))
			time_qs[b, ((i-1)*(length(p)) + track)] <- system.time(QS <- 2*lrvar(x = v.hat, type = "Andrews", prewhite = FALSE, adjust = FALSE, kernel = "Quadratic Spectral", bw = fix_b)*T - lrvar(x = v.hat, type = "Andrews", prewhite = FALSE, adjust = FALSE, kernel = "Quadratic Spectral", bw = fix_b/3)*T)[3]
			time_th[b, ((i-1)*(length(p)) + track)] <- system.time(lug <- 2*lrvar(x = v.hat, type = "Andrews", prewhite = FALSE, adjust = FALSE, kernel = "Tukey", bw = fix_b)*T - lrvar(x = v.hat, type = "Andrews", prewhite = FALSE, adjust = FALSE, kernel = "Tukey", bw = fix_b/3)*T )[3]
			time_bt[b, ((i-1)*(length(p)) + track)] <- system.time(lug <- 2*lrvar(x = v.hat, type = "Andrews", prewhite = FALSE, adjust = FALSE, kernel = "Bartlett", bw = fix_b)*T - lrvar(x = v.hat, type = "Andrews", prewhite = FALSE, adjust = FALSE, kernel = "Bartlett", bw = fix_b/3)*T )[3]
			time_bm[b, ((i-1)*(length(p)) + track)] <- system.time(lugbm <-  mcse.multi(v.hat, size = fix_b, r = 3)$cov)[3]
		}
	}
}

save(file = "time_comp",time_qs, time_bt, time_th, time_bm)
