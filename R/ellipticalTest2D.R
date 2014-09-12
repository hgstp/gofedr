ellipticalTest2D  <-  function(sampling) {
	# Result: test statistic z. Under the null hypothesis, z has standard normal
	# 		  distribution.
	#sampling: nxd matrix, whose n rows each contain a d-dimensional sample
	
	X <- sampling[,1]
	Y <- sampling[,2]
	n <- length(X)

	# Estimator for correlation coefficient [Eq. 3.10]
	hatRho <- cor(x = X, y = Y)
	# Estimator for Kendall's tau [Def. 3.5]
	hatTau <- cor(x = X, y = Y, method = "kendall")
	
	# The test statistic [Def. 3.11]
	T <- hatRho - sin(pi * 0.5 * hatTau)

	# Vector of empirical cdf of X evaluated at X_i
	hatF_x <- rank(X, ties.method = "max") / n
	# Vector of empirical cdf of Y evaluated at Y_i
	hatF_y <- rank(Y, ties.method = "max") / n
	
	# Vector of empirical cdf of (X, Y) evaluated at (X_i, Y_i)
	hatF_xy <- numeric(n)                            
	for (i in 1:n) {
		hatF_xy[i] <- sum((X <= X[i]) & (Y <= Y[i])) / n
	}

	# Hoeffding projection according to [Lemma 3.10]
	hatPhi <-	1 - 2 * hatF_x - 2 * hatF_y + 4 * hatF_xy

	# The following computes the estimated covariance matrix V in [Remark 3.13]
	A <- cbind(X - mean(X), Y - mean(Y), (X^2) - mean(X^2), (Y^2)-mean(Y^2),
	  	 		 (X*Y) - mean(X*Y), (2*hatPhi) - mean(2*hatPhi))
	V <- t(A) %*% A / n

	# The following computes Df(t_0) in [Theorem 3.12]
	Df <- array(0, c(1,6))

	moment_10 <- mean(X)
	moment_01 <- mean(Y)
	moment_20 <- mean(X^2)
	moment_02 <- mean(Y^2)
	moment_11 <- mean(X * Y)

	Df[1,1] <- -moment_01 / sqrt((moment_20 - moment_10 ^ 2) * (moment_02 -
			   moment_01 ^ 2)) + (moment_11 - moment_10 * moment_01) * moment_10
			   * (moment_02 - moment_01 ^ 2) / sqrt((moment_20 - moment_10 ^ 2)
			   * (moment_02 - moment_01 ^ 2)) ^ 3
	Df[1,2] <- -moment_10 / sqrt((moment_20 - moment_10 ^ 2) * (moment_02 -
			   moment_01 ^ 2)) + (moment_11 - moment_10 * moment_01) * moment_01
			   * (moment_20 - moment_10 ^ 2) / sqrt((moment_20 - moment_10 ^ 2)
			   * (moment_02 - moment_01 ^ 2)) ^ 3
	Df[1,3] <- - (moment_11 - moment_10 * moment_01) * (moment_02 - moment_01 ^
			   2) / sqrt((moment_20 - moment_10 ^ 2) * (moment_02 - moment_01 ^
			   2)) ^ 3 / 2
	Df[1,4] <- - (moment_11 - moment_10 * moment_01) * (moment_20 - moment_10 ^
			   2) / sqrt((moment_20 - moment_10 ^ 2) * (moment_02 - moment_01 ^
			   2)) ^ 3 / 2
	Df[1,5] <-  1 / sqrt((moment_20 - moment_10 ^ 2) * (moment_02 - moment_01 ^ 2))
	Df[1,6] <- -cos(hatTau * pi / 2) * pi / 2

	# The square root of v in [Equation 3.15]
	std <- sqrt(Df %*% V %*% t(Df))

	# The transformed test statistic [Definition 3.14]
	z <- sqrt(n) * T / std                        
}

ellipticalTest2DRejection <- function(X, alpha) {
	# Result: TRUE or FALSE. TRUE if we reject the null hypothesis using level
	#         alpha, FALSE otherwise.
	# sampling: nxd matrix containing the sampling
	# alpha: level of the test
	
	TS <- ellipticalTest2D(X)
	abs(TS) > qnorm(1 - alpha / 2)
}

ellipticalTest2DPower <- function(TS, alpha) {
	# Result: proportion indicating how often we reject the null hypothesis
	# d: dimension of the original sampling 
	# TS: vector of chi-squared test statistics that are the result of
	#     ellipticalTest2D
	# alpha: level of the test
	
	NR <- length(TS)
	sum(abs(TS) > qnorm(1 - alpha / 2)) / NR
}
