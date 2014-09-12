ellipticalTestHD  <-  function(sampling) {
	# Result: test statistic chi_squared. Under the null hypothesis,
	# chi_squared has chi_squared distribution with degrees of freedom equal to
	# n_pairings = choose(d, 2) sampling: nxd matrix, whose n rows each contain
	# a d-dimensional sample.

	n <- dim(sampling)[1]
	d <- dim(sampling)[2]
	n_pairings <- choose(d, 2)

	# When we estimate quantities that relate to the projection of our sampling
	# on the i-th and j-th dimension, we store these quantities in upper right
	# triangular matrices.
	
	# Entry (i, j), i <= j, stores the concordance measure of the i-th dimension
	# with the j-th dimension
	rho <- array(0, c(d, d))
	tau <- array(0, c(d, d))

	# Test statistic T_n with length n_pairings [Definition 3.20]
	TS <- numeric(n_pairings)
	
	l <- 1
	for(i in 1:(d-1)) {
		for (j in (i+1):d) {
			rho[i, j] <- cor(x = sampling[,i], y = sampling[,j])
			tau[i, j] <- cor(x = sampling[,i], y = sampling[,j],
							 method = "kendall")
			
			TS[l] <- rho[i, j] - sin(0.5 * pi * tau[i, j])
			l <- l + 1
		}
	}

	# hatF[i,k] is the empirical cumulative distribution function of the i-th
	# dimension evaluated for the k-th sampling
	hatF <- array(0, c(d, n))
	for(i in 1:d) {
		hatF[i,] <- rank(sampling[,i], ties.method = "max") / n
	}
	
	# ecdf[i,j,k] is the joint empirical cumulative distribution function of the
	# i-th and the j-th dimension, evaluated for the k-th sampling.
	ecdf <- array(0, c(d, d, n))
	for(i in 1:(d-1)) {
		for (j in (i+1):d) {
			for(k in 1:n) {
				ecdf[i, j, k] <- sum((sampling[,i] <= sampling[k,i]) &
									 (sampling[,j] <= sampling[k,j])) / n
			}
		}
	}

	# phi[i,j,k] is the Hoeffding projection related to dimension i and j,
	# evaluated for the k-th sampling. It is computed using [Lemma 3.10]
	phi <- array(0, c(d, d, n))
	for(i in 1:(d-1)) {
		for (j in (i+1):d) {
			phi[i,j,] <- 1 - 2 * hatF[i,] - 2 * hatF[j,] + 4 * ecdf[i,j,]
		}
	}
	phi <- 2 * phi

	n_estimators <- 2* (d + n_pairings)
	
	# Matrix A used to efficiently compute covariance matrix V_d of estimators
	# in [Remark 3.22]
	A <- array(0, c(n, n_estimators))

	# First 2d columns used for computing estimated Cov(X_i, .)
	for(i in 1:d) {
		A[,i] <- sampling[,i] - mean(sampling[,i])
		A[,d+i] <- sampling[,i]^2 - mean(sampling[,i]^2)
	}
	
	# Next 2 * n_pairings columns used for computing estimated Cov(X_i * X_j, .)
	# and Cov(phi, .)
	k <- 2*d + 1
	for(i in 1:(d-1)) {
		for (j in (i+1):d) {
			A[,k] <- sampling[,i] * sampling[,j] - mean(sampling[,i] * sampling[,j])
			A[,n_pairings + k] <- phi[i,j,] - mean(phi[i,j,])
			k <- k + 1
		}
	}
	
	# The covariance matrix of estimators (corresponds to the matrix V_d in
	# [Theorem 3.21])
	V1 <- t(A) %*% A / n
	
	Df <- array(0, c(n_pairings, n_estimators))

	# jacobiMatrix computes the rows of the Jacobi matrix used for the delta
	# method. The result is equal to the differential of function f from the pdf
	# (save for interspersed zeros), evaluated for different moments and tau's.
	# See [Equation 3.33] and [Remark 3.22]
	jacobiMatrix <- function(m_10, m_01, m_20, m_02, m_11, tau, i, j) {
		Dg <- numeric(n_estimators)
		moment_10 <- mean(sampling[,i])
		moment_01 <- mean(sampling[,j])
		moment_20  <- mean(sampling[,i]^2)
		moment_02  <- mean(sampling[,j]^2)
		moment_11 <- mean(sampling[,i] * sampling[,j])
		
		Dg[i] <- -moment_01 / sqrt((moment_20 - moment_10 ^ 2) * (moment_02 -
			  	 moment_01 ^ 2)) + (moment_11 - moment_10 * moment_01) *
				 moment_10 * (moment_02 - moment_01 ^ 2) / sqrt((moment_20 -
				 moment_10 ^ 2) * (moment_02 - moment_01 ^ 2)) ^ 3
		Dg[j] <- -moment_10 / sqrt((moment_20 - moment_10 ^ 2) * (moment_02 -
			  	 moment_01 ^ 2)) + (moment_11 - moment_10 * moment_01) *
				 moment_01 * (moment_20 - moment_10 ^ 2) / sqrt((moment_20 -
				 moment_10 ^ 2) * (moment_02 - moment_01 ^ 2)) ^ 3
		Dg[d+i] <- - (moment_11 - moment_10 * moment_01) * (moment_02 -
				   moment_01 ^ 2) / sqrt((moment_20 - moment_10 ^ 2) *
				   (moment_02 - moment_01 ^ 2)) ^ 3 / 2
		Dg[d+j] <- - (moment_11 - moment_10 * moment_01) * (moment_20 -
				   moment_10 ^ 2) / sqrt((moment_20 - moment_10 ^ 2) *
				   (moment_02 - moment_01 ^ 2)) ^ 3 / 2
		
		pos <- 0
		if(i >= 2) {
			for(r in 1:(i-1)) {
				pos <- pos + d-r
			}
		}
		pos <- pos + j-i
		
		Dg[2*d + pos] <- 1 / sqrt((moment_20 - moment_10 ^ 2) * (moment_02 -
			   	 	  	 moment_01 ^ 2))
		Dg[2*d + n_pairings + pos] <- -cos(tau * pi / 2) * pi / 2
		
		Dg
	}
	
	k <- 1
	for(i in 1:(d-1)) {
		for (j in (i+1):d) {
			Df[k,] <- jacobiMatrix(mean(sampling[,i]), mean(sampling[,j]),
				   	  mean(sampling[,i]^2), mean(sampling[,j]^2),
					  mean(sampling[,i] * sampling[,j]), tau[i,j], i, j)
			k <- k + 1
		}
	}

	# The covariance matrix of the delta method, corresponds to W_d in
	# [Theorem 3.21]
	V2 <- Df %*% V1 %*% t(Df)
	
	# The chi-squared test statistic from [Definition 3.23]
	chi_squared <- n * t(TS) %*% solve(V2) %*% TS
}

ellipticalTestHDRejection <- function(sampling, alpha) {
	# Result: TRUE or FALSE. TRUE if we reject the null hypothesis using level
	#         alpha, FALSE otherwise.
	# sampling: nxd matrix containing the sampling
	# alpha: level of the test
	
	TS <- ellipticalTestHD(sampling)
	d <- dim(sampling)[2]
	n_pairings <- choose(d, 2)
	TS > qchisq(p = 1 - alpha, df = n_pairings)
}

ellipticalTestHDPower <- function(d, TS, alpha) {
	# Result: proportion indicating how often we reject the null hypothesis
	# d: dimension of the original sampling 
	# TS: vector of chi-squared test statistics that are the result of
	#     ellipticalTestHD
	# alpha: level of the test
	
	n_pairings <- choose(d, 2)
	NR <- length(TS)
	sum(abs(TS) > qchisq(p = 1 - alpha, df = n_pairings)) / NR
}
