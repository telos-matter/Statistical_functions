# A bunch of commonly used statistical functions 
# Most of the functions are not optimized, they are not meant for any production use.

# HEMMOUDA Aymane - 04_2022
# ahemmouda@gmail.com


################################################################ Uni and Bi statistical series

# Computes the unbiased variance
# x: a statistical serie
real_var <<- function (x) {
	return (var(x)*(length(x)-1)/length(x))
}

# Computes the unbiased standard deviation
# x: a statistical serie
real_sd <<- function (x) {
	return (sqrt(real_var(x)))
}

# Computes the unbiased covariance of two series
# x: a statistical serie
# y: a statistical serie
real_cov <<- function (x,y) {
	return (cov(x,y)*(length(x)-1)/length(x))
}

# Computes the variation coefficient
# x: a statistical serie
var_coef <<- function (x) {
	return (real_sd(x)/mean(x))
}

# Computes the correlation coefficient of two series
# x: a statistical serie
# y: a statistical serie
correlation_coef <<- function (x, y) {
	return (real_cov(x,y)/(real_sd(x)*real_sd(y)))
}

# Computes the centered moment of degree r for a serie
# x: a statistical serie
# r: the degree
centered_moment <<- function (x, r) {
	mean <- mean(x)
	u <- 0
	for (xi in x) {
		u <- u + ((xi-mean)^r)
	}
	u <- u/length(x)
	return (u)
}

# Constructs a list of the specified modalities repeated frequencies times
# In the following manner: modalities [i] is repeated frequencies [i] times
# modalities: a list of modalities, for example [1500, 1700, 1900]
# frequencies: a list of frequencies, for example [13, 12, 15]
cnt_c <<- function (modalities, frequencies) {
	x <- c()
	for (i in 1:length(modalities)) {
		x <- c(x, rep(modalities[i], frequencies[i]))
	}
	return (x)
}

# Computes the coefficients (a and b) of a linear regression line for two statistical series
# The coefficients (a, b) such that y ≈ ax +b
# x: a statistical serie
# y: a statistical serie
# plot: whether or not to plot the two series, by default no
regression_line <<- function (x,y, plot= FALSE) {
	a <- real_cov(x,y)/real_var(x)
	b <- mean(y) -a*mean(x)
	if (plot) {
		plot(x,y)
		abline(b,a, col= "black")
		points(mean(x), mean(y), col= "red")
	}
	return (c(a,b))
}

# Computes an approximation of the specified value using a linear regression between two series
# x: a statistical serie
# y: a statistical serie
# value: the value to approximate, where the returned value is a*value +b with a and b being the coefficients of the regression line that maps x to y
# plot: whether or not to plot the two series and the estimated value, by default no
compute_regression_line <<- function (x,y,value, plot= FALSE) {
	coef <- regression_line(x,y, plot)
	a <- coef[1]
	b <- coef[2]
	if (plot) {
		points(value, (a*value +b), col= "green")
	}
	return (a*value +b)
}

# Returns a table/matrix representing a bi-varied statistical serie, along side the individual series.
# The data (the number of occurences or frequency) is from top left to bottom right, X is vertical (thus its modalitites are in the rows) Y is horizontal (thus its modalities are in the columns) ()
# NOTE: Use unlist, unname and/or a <- a[[1]] to retrieve the data from the returned list
# data: a list of data from top left going to bottom right
# x_modalities: the modalitites of statistical serie x
# y_modalities: the modalitites of statistical serie y
# EXAMPLE: bi_serie(c(14,28,20,6,46,86), c(21,23,25), c(900,1100))
bi_serie <<- function (data, x_modalities, y_modalities) {
	a <- matrix(data, length(x_modalities), length(y_modalities))

	x_n <- c()
	for (i in 1:length(x_modalities)) {
		row <- a[i,]
		n <- 0
		for (value in row) {
			n <- n + value
		}
		x_n <- c(x_n, n)
	}
	x <- cnt_c(x_modalities, x_n)
	
	y_n <- c()
	for (i in 1:length(y_modalities)) {
		col <- a[,i]
		n <- 0
		for (value in col) {
			n <- n + value
		}
		y_n <- c(y_n, n)
	}
	y <- cnt_c(y_modalities, y_n)
	
	rownames(a) <- x_modalities
	colnames(a) <- y_modalities

	return (list(a,x,y))
}

# Computes the khi/chi (not sure how it's called anymore) squared for a table/matrix of a bi-varied statistical serie (χ²) 
# a: the table/matrix of a  bi-varied statistical serie
# EXAMPLE: 
# 	a <- bi_serie(c(14,28,20,6,46,86), c(21,23,25), c(900,1100))
# 	a <- a[[1]]
# 	khi_squared (a)
khi_squared <<- function (a) {
	khi <- 0
	n <- sum(a)
	for (i in 1:nrow(a)) {
		ni <- sum(a[i,])
		for (j in 1:ncol(a)) {
			nj <- sum(a[,j])
			nij <- ni*nj/n
			khi <- khi + (a[i,j]-nij)*(a[i,j]-nij)/nij
		}
	}
	return (khi)
}

# Computes the cramer coefficient for a bi-varied statistical serie
# a: the table/matrix of a  bi-varied statistical serie
cramer_coef <<- function (a) {
	return(sqrt(khi_squared(a)/(sum(a)*min(nrow(a)-1,ncol(a)-1))))
}

# Retrieves the individual indexed (conditioned) series form a bi-varied statistical serie
# a: the table/matrix of a  bi-varied statistical serie
# index: the index that we are "conditioning" on
# yIndex: whether or not that index is the y index or the x index, by default it is the y
cond_c <<- function (a, index, yIndex= TRUE) {
	if (yIndex) {
		return (cnt_c(strtoi(names(a[,index])), unname(a[,index])))
	} else {
		return (cnt_c(strtoi(names(a[index,])), unname(a[index,])))
	}
}

# Computes the quantile for a continuous statistical serie
# x: a statistical serie
# precent: which quantile to compute, for example 50 (Q2), by default it is NULL and the function then returns information about the 25%, 50% and 75% quantile
# amps: the amplitudes of the statical serie, if NULL then it assumes that the amps are fixed and computes them by simply modalities[2]-modalities[1]
cnt_quantile <<- function (x, percent= NULL, amps= NULL) {
	if (! is.null(percent)) {
		percent <- percent/100
			
		modalities <- c()
		for (value in x) {
			if (is.na(match(value,modalities))) {
				modalities <- c(modalities, value)
			}
		}
	
		ni <- c()
		for (modalitie in modalities) {
			ni <- c(ni,sum(x == modalitie))
		}

		N = sum(ni)
	
		fi <- c()
		for (n in ni) {
			fi <- c(fi,n/N)
		}
	
		Fi <- c()
		for (f in fi) {
			if (is.null(tail(Fi,1))) {
				Fi <- c(f)
			} else {
				Fi <- c(Fi, tail(Fi,1)+f)
			}
		}
				
		i <- min(which(Fi>=percent))
		if (is.null(amps)) {
			amps <- modalities[2]-modalities[1]
		} else {
			amps <- amps[i] 
		}		
		if (i == 1) {
			return (modalities[i] -amps/2 +amps*percent/fi[i])
		} else {
			return (modalities[i] -amps/2 +amps*(percent -Fi[i-1])/fi[i])
		}
	} else {
		return (c("25%"= cnt_quantile(x,25,amps), "50%"= cnt_quantile(x,50,amps), "75%"= cnt_quantile(x,75,amps)))
	}
}

# Computes the quantile of a continuous statistical serie retrieved from a table/matrix of a bi-varied statistical serie
# a: the table/matrix of a  bi-varied statistical serie
# index: the index that we are "conditioning" on
# yIndex: whether or not that index is the y index or the x index, by default it is the y
# perecent: the quantile to compute, for example 50 (Q2), by default it is NULL and returns information about the 25%, 50% and 75% quantile
cnt_quantile_cond <<- function(a, index, yIndex= TRUE, percent= NULL) {
	return(cnt_quantile(cond_c(a,index,yIndex), percent))
}

# Computes the intra variance
# mean(cond var)
# a: the table/matrix of a  bi-varied statistical serie
# forX: whether or not to compute it for x or y, by default for x
intra_var <<- function (a, forX= TRUE) {
		if (forX) {
		sum <- 0
		for (j in 1:ncol(a)) {
			sum <- sum + sum(a[,j])*real_var(cond_c(a,j))
		}
		return (sum/sum(a))
	} else {
		sum <- 0
		for (i in 1:nrow(a)) {
			sum <- sum + sum(a[i,])*real_var(cond_c(a,i,FALSE))
		}
		return (sum/sum(a))
	}
}

# Computes the inter variance
# var(cond mean)
# a: the table/matrix of a  bi-varied statistical serie
# forX: whether or not to compute it for x or y, by default for x
inter_var <<- function (a, forX= TRUE) {
		if (forX) {
		sum <- 0
		mean <- mean(cnt_c(strtoi(rownames(a)),rowSums(a)))
		for (j in 1:ncol(a)) {
			sum <- sum + sum(a[,j])*((mean(cond_c(a,j))-mean)^2)
		}
		return (sum/sum(a))
	} else {
		sum <- 0
		mean <- mean(cnt_c(strtoi(colnames(a)),colSums(a)))
		for (i in 1:nrow(a)) {
			sum <- sum + sum(a[i,])*((mean(cond_c(a,i,FALSE))-mean)^2)
		}
		return (sum/sum(a))
	}
}

# Computes the correlation rate of a bi-varied statistical serie
# a: the table/matrix of a  bi-varied statistical serie
# forX: whether or not to compute it for x or y, by default for x
correlation_rate <<- function (a, forX= TRUE) {
	if (forX) {
		return(inter_var(a, forX)/real_var(cnt_c(strtoi(rownames(a)),rowSums(a))))
	} else {
		return(inter_var(a, forX)/real_var(cnt_c(strtoi(colnames(a)),colSums(a))))
	}
}

# Gives general informations about a statistical serie
# x: a statistical serie
# isCnt: whether it is continuous or not, by default TRUE
# amps: the amplitudes of the statical serie, if NULL then it assumes that the amps are fixed and computes them by simply modalities[2]-modalities[1]
info <<- function (x, isCnt= TRUE, amps= NULL) {
	modalities <- c()
	for (value in x) {
		if (is.na(match(value,modalities))) {
			modalities <- c(modalities, value)
		}
	}
	
	ni <- c()
	for (modalitie in modalities) {
		ni <- c(ni,sum(x == modalitie))
	}
	
	Ni <- c()
	for (n in ni) {
		if (is.null(tail(Ni,1))) {
			Ni <- c(n)
		} else {
			Ni <- c(Ni, tail(Ni,1)+n)
		}
	}
	
	N = tail(Ni,1)
	
	fi <- c()
	for (n in ni) {
		fi <- c(fi,n/N)
	}
	
	Fi <- c()
	for (f in fi) {
		if (is.null(tail(Fi,1))) {
			Fi <- c(f)
		} else {
			Fi <- c(Fi, tail(Fi,1)+f)
		}
	}
	
	stats <- data.frame(Modalities= modalities, .______ni= ni, .______Ni= Ni, .______fi= fi, .______Fi= Fi)	
	info <- data.frame(Mode= modalities[match(max(ni),ni)], .______mean= mean(x), .______var= real_var(x), .______sd=real_sd(x))
	quant <- NULL
	if (isCnt) {
		Q1 <- cnt_quantile (x, 1, amps)
		Q25 <- cnt_quantile (x, 25, amps)
		Q50 <- cnt_quantile (x, 50, amps)
		Q75 <- cnt_quantile (x, 75, amps)
		Q99 <- cnt_quantile (x, 99, amps)
		quant <- data.frame(C1= Q1, .______Q1= Q25, .______Q2= Q50, .______Q3= Q75, .______C99= Q99, ".______Q3_Q1"= (Q75-Q25), ".______C99_C1"= (Q99-Q1))
	} else {
		Q1 <- quantile (x, 1)
		Q25 <- quantile (x, 25)
		Q50 <- quantile (x, 50)
		Q75 <- quantile (x, 75)
		Q99 <- quantile (x, 99)
		quant <- data.frame(C1= Q1, .______Q1= Q25, .______Q2= Q50, .______Q3= Q75, .______C99= Q99, ".______Q3_Q1"= (Q75-Q25), ".______C99_C1"= (Q99-Q1))
	}
		
	output <- list(stats, info, quant)

	return(output)
}

# Gives general informations about a statistical serie retrieved from a bi-varied statistical serie
# a: the table/matrix of a bi-varied statistical serie
# index: the index that we are "conditioning" on
# yIndex: whether or not that index is the y index or the x index, by default it is the y
# isCnt: whether the retrieved statistical serie is continuous or not, by default TRUE
# amps: the amplitudes of the statical serie, if NULL then it assumes that the amps are fixed and computes them by simply modalities[2]-modalities[1]
info_cond <<- function (a, index, yIndex= TRUE, isCnt= TRUE, amps= NULL) {
		return(info(cond_c(a,index,yIndex),isCnt, amps))
}


################################################################ Estimation

# Same as pnorm but you can specify whether alpha is smaller than or greater than (i should've called it less_than or something, anyhow)
# alpha: alpha
# mean: mean
# sd: standard deviation
# smaller: P(X < alpha) or P(X > alpha), by default TRUE -> P(X < alpha)
pn <<- function(alpha, mean, sd, smaller= TRUE) {
	if (smaller) {
		return (pnorm (alpha, mean= mean, sd= sd))
	} else {
		return (1 -pnorm (alpha, mean= mean, sd= sd))
	}
}

# Probability that a given sample has a mean smaller/greater than alpha
# alpha: alpha
# mean: mean
# sd: standard deviation
# n: the sample size
# samller: whether less than or greater than alpha, by default TRUE
p_mean <<- function(alpha, mean, sd, n, smaller= TRUE) {
	return (pn(alpha, mean, (sd/sqrt(n)), smaller))
}

# Probability that a given sample has a proportionality smaller/greater than alpha
# alpha: alpha
# population_proportionality: the general population proportionaility
# n: the sample size
# samller: whether less than or greater than alpha, by default TRUE
p_sample <<- function (alpha, population_proportionality, n, smaller= TRUE) {
	return (pn(alpha, population_proportionality, sqrt(population_proportionality*(1-population_proportionality)/n),smaller))
}

# TODO: doc
# confidence: as a prencetage, 95 and not 0.95 for example
mean_interval <<- function (x, confidence, var= NULL, substituted= TRUE, N= NULL) {
	confidence <- 1 -((100 -confidence)/100)/2
	z <- qnorm(confidence)
	a <- NULL
	b <- NULL
	if (is.null(var)) {
		var <- real_var(x)
		if (length (x) < 30) {
			z <- qt(confidence,length(x) -1)
		}
	}
	if (substituted) {
		a <- mean(x) -z*sqrt(var/length(x))
		b <- mean(x) +z*sqrt(var/length(x))
	} else {
		coef <- (N-length(x))/(N-1)
		a <- mean(x) -z*sqrt((var/length(x))*coef)
		b <- mean(x) +z*sqrt((var/length(x))*coef)
	}
	
	return (c(a,b))
}

# Determines a proportionality interval with the specified condifence
# confidence: as a prencetage, 95 and not 0.95 for example
prop_interval <<- function (proportionality, confidence, n, substituted= TRUE, N= NULL) {
	confidence <- 1 -((100 -confidence)/100)/2
	z <- qnorm(confidence)
	a <- NULL
	b <- NULL
	if (substituted) {
		a <- proportionality -z*sqrt((proportionality*(1-proportionality))/(n-1))
		b <- proportionality +z*sqrt((proportionality*(1-proportionality))/(n-1))
	} else {
		coef <- (N-n)/(N-1)
		a <- proportionality -z*sqrt(((proportionality*(1-proportionality))/(n-1))*coef)
		b <- proportionality +z*sqrt(((proportionality*(1-proportionality))/(n-1))*coef)
	}
	
	return (c(a,b))
}

# How much n (sample size) is needed for the mean interval to be of size eem*2 with that confidence
# confidence: as a prencetage, 95 and not 0.95 for example
mean_n <<- function (var, confidence, eem, usingN= TRUE, n= NULL) {
	confidence <- 1 -((100 -confidence)/100)/2
	z <- NULL
	if (usingN) {
		z <- qnorm(confidence)
	} else {
		z <- qt(confidence, n-1)
	}
	return ((z^2)*var/(eem^2))
}

# How much n (sample size) is needed for the proportionality interval to be of size eem*2 with that confidence
# confidence: as a prencetage, 95 and not 0.95 for example
prop_n <<- function (proportionality, confidence, eem) {
	confidence <- 1 -((100 -confidence)/100)/2
	z <- qnorm(confidence)
	eem <- eem/100
	return ((z^2)*proportionality*(1-proportionality)/(eem^2) +1)
}


################################################################ Hypothesis

# Preforms a hypothesis test
# alpha: alpha
# mean: mean
# sd: standard deviation
# n: sample size
# test_type:
# 	test_type < 0 : smaller than
# 	test_type = 0 : equal
# 	test_type > 0 : greater than
# risk: the risk in precentage
test <<- function (alpha, mean, sd, n, test_type, risk= 5) {
	z <- (alpha-mean)/(sd/sqrt(n))
	reject <- FALSE
	risk_z <- NULL
	p_value <- NULL
	if (test_type < 0) {
		risk <- 1 -(risk/100)
		risk_z = -qnorm(risk)
		if (z < risk_z) {
			reject <- TRUE
		}
		p_value <- pnorm(z)
	} else if (test_type == 0) {
		risk <- 1 -(risk/100)/2
		risk_z = qnorm(risk)
		if (abs(z)>risk_z) {
			reject <- TRUE
		}
		p_value <- 2*(1-pnorm(z))
	} else {
		risk <- 1 -(risk/100)
		risk_z = qnorm(risk)
		if (z > risk_z) {
			reject <- TRUE
		}
		p_value <- 1 -pnorm(z)
	}
	if (reject) {
		cat("Reject H0, p-value: ",p_value, " , z= ", z, " ,risk_z= ", risk_z)
	} else {
		cat("Accept H0, p-value: ",p_value, " , z= ", z, " ,risk_z= ", risk_z)
	}
}
