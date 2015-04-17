# Density of p-values for two group t test
# TODO: two step sizes: high precision for p <.05, lower precision for p > .05
# TODO: Use library(numDeriv)
# n is overall n (both groups)
create_pd <- function(d, n, p.max=.05, step=.0001) {
	ps <- seq(step, p.max, by=step)			# p-vales at which the curve is computed
	df <- n-2
	t.values <- qt(ps/2, df, lower = FALSE) # corresponding t-values

	# no idea why it's n/4 here ... but it's correct
	power <- pt(t.values, df, ncp = sqrt(n/4)*d, lower = FALSE) + pt(-t.values, df, ncp = sqrt(n/4)*d, lower = TRUE)

	pdist <- data.frame(q=ps, power=power)

	# Compute lag to get density
	pdist$q2 <- c(0, pdist$q[1:(length(pdist$q)-1)])
	pdist$dens0 <- (pdist$power - c(0, pdist$power[1:(nrow(pdist)-1)]))

	# First factor for scaling it to a true density (if the range would be the full range from 0 to 1)
	pdist$dens <- pdist$dens0/step
	
	# Second factor for scaling: AUC should be 1 in the selected range from 0 to p.max
	fac2 <- 1/(sum(pdist$dens0))
	pdist$dens <- pdist$dens*fac2

	# AUC = 1? Is it a density? Stop if not!
	if (round(sum(pdist$dens)*step, 5) != 1) {stop("WARNING: Computed p-curve is not a density: AUC != 1")}

	# create a continous interpolation of the function
	pd <- approxfun(pdist$q, pdist$dens, rule=2)
	return(pd)
}

# n is overall n (both groups)
pd <- function(p, d, n, p.max=.05, step=.0001) {
	pd0 <- create_pd(d=d, n=n, p.max=p.max, step=step)
	pd0(p)
}



# ---------------------------------------------------------------------
# Wrapper: Input here is t value and df
#' @param ... Other parameters passed to create_pd
create_ptd <- function(t, df, ...) {
	n <- df+2
	d <- (2 * t)/sqrt(2 * n - 2)
	create_pd(d=d, n=n, ...)
}

ptd <- function(p, t, df, p.max=.05, step=.0001) {
	ptd0 <- create_ptd(t=t, df=df, p.max=p.max, step=step)
	ptd0(p)
}


# ---------------------------------------------------------------------
# This function g is the implementation of eq (A.1) of Hung, H. M. J., O’Neill, R. T., Bauer, P., & Kohne, K. (1997). The Behavior of the P-Value When the Alternative Hypothesis is True. Biometrics, 53, 11–22. doi:10.2307/2533093

# WARNING: But it looks really different than the real curve!!

# small phi = \phi = density of Gaussian, dnorm
# capital phi = Φ = cumulative distribution function (CDF)
# sd is always 1
g <- function(p, delta, n) {
	Z_p <- qnorm(1-p, lower.tail=TRUE)
	dnorm(Z_p, mean=sqrt(n)*delta)/dnorm(Z_p)
}


# # This is Figure 1 from Hung, H. M. J., O’Neill, R. T., Bauer, P., & Kohne, K. (1997). The Behavior of the P-Value When the Alternative Hypothesis is True. Biometrics, 53, 11–22. doi:10.2307/2533093
# p <- seq(0.00001, .99999, length.out=100)
# plot(p, g(p, 1/3, 15), type="l", ylim=c(0, 10))
# lines(p, g(p, 1/3, 30), type="l", ylim=c(0, 10), lty="dashed")
# lines(p, g(p, 1/3, 60), type="l", ylim=c(0, 10), lty="dotdash")
# lines(p, g(p, 1/3, 80), type="l", ylim=c(0, 10), lty="dotted")




# ---------------------------------------------------------------------
# Other p-value related functions

# Compute observed power from p value, using the normal distribution as approximation (see R-index paper)
obspower <- function(p, p.max=.05) {
	Z <- qnorm(1-(p/2))
	obs.power <- pnorm(Z-qnorm(1-p.max/2))
}


# ---------------------------------------------------------------------
# Test

# Show p-curve
# t.value <- 2.53
# n <- 40
# df <- n-2
# (p.value <- pt(t.value, df, lower.tail=FALSE)*2)
# (d <- (2 * t.value)/sqrt(2 * n - 2))
#
# # These should give the same
# pd(p.value, d, n)
# ptd(p.value, t.value, df)


## ======================================================================
## Functions for computing p-curve and p-likelihood
## ======================================================================

# Alternative pp33 function, which relies on a conversion from a p value to Z values, without df.
# EXPERIMENTAL, NOT THOROUGHLY TESTED YET!!
pp33b <-function(x, prob=2/3, p.crit=.05) {
	#Find critical value of student (xc) that gives p=.05 when df=df_
	xc=qnorm(p=1-(p.crit/2))
		
	#Find noncentrality parameter (ncp) that leads 33% power to obtain xc
	f <- function(delta, pr, x_) pnorm(x_, mean = delta) - pr
	out <- uniroot(f, c(0, 37.62), pr = prob, x_ = xc)	
	M <- out$root

	#Find probability of getting x_ or larger given ncp
	p_larger=pnorm(x, mean=M)

	#Condition on p<.05 (i.e., get pp-value)
	pp <- 3*(p_larger-prob)

	#Print results
	return(pp)
}
pp33bv <- Vectorize(pp33b)

# p-curve computations from app2 (old style: chi2-test)
p_curve_2 <- function(p.values, p.crit=.05, robust_p_curve=TRUE) {
	pp.rs <- p.values/p.crit
	pp.ls <- 1-p.values/p.crit
	Z <- qnorm(1-(p.values/2))
	pp33 <- pp33bv(Z, p.crit=p.crit)

	# remove non-significant p-values
	pp.ls[p.values > p.crit] <- NA
	pp33[p.values > p.crit] <- NA
	pp.rs[p.values > p.crit] <- NA

	# Make pp values robust? Winsorize at .01 and .99
	if (robust_p_curve == TRUE) {
		pp.ls[pp.ls > .99] <- .99
		pp.ls[pp.ls < .01] <- .01
		pp.rs[pp.rs > .99] <- .99
		pp.rs[pp.rs < .01] <- .01
		pp33[pp33 > .99] <- .99
		pp33[pp33 < .01] <- .01
	}
	
	df <- 2*sum(!is.na(pp.rs))

	chi2_evidence <- -2*sum(log(pp.rs), na.rm=TRUE)
	p_evidence <- pchisq(chi2_evidence, df=df, lower.tail=FALSE)

	chi2_hack <- -2*sum(log(pp.ls), na.rm=TRUE)
	p_hack <- pchisq(chi2_hack, df=df, lower.tail=FALSE)

	chi2_lack <- -2*sum(log(pp33), na.rm=TRUE)
	p_lack <- pchisq(chi2_lack, df=df, lower.tail=FALSE)

	return(list(
		chi2_evidence = chi2_evidence, 
		p_evidence = p_evidence, 
		chi2_hack = chi2_hack, 
		p_hack = p_hack, 
		chi2_lack = chi2_lack, 
		p_lack = p_lack,
		inconclusive = ifelse(p_evidence>.05 & p_lack>.05 & p_hack>.05, TRUE, FALSE)))
}


# sampleFrom = c("H0", "H1", "hack")
# If "hack", you can still define a true ES. Samples are then optionally increased by hackIncrease participants until they reach n.max. Then still only the significant studies are published.
sim_studies <- function(k=5, n=40, n.max=2*n, hackIncrease=5, sampleFrom="H0", d=0, seed=NULL) {
	
	if (n%%2 != 0) stop("n should be dividable by 2!")
	set.seed(seed)

	t.values <- c()
	ns <- c()
	
	if (sampleFrom=="H0") {
		i <- 0
		repeat {
			x <- rnorm(n/2)
			y <- rnorm(n/2)
			t0 <- t.test(x, y)
	
			if (t0$p.value <= .05) {
				ns <- c(ns, n)
				t.values <- c(t.values, abs(t0$statistic))
				i <- i+1
				if (i >= k) break;
			}
		}
	} else if (sampleFrom=="H1") {
		i <- 0
		repeat {
			n.group <- 20
			x <- rnorm(n/2)
			y <- rnorm(n/2, mean=d)
			t0 <- t.test(x, y)

			if (t0$p.value <= .05) {
				ns <- c(ns, n)
				t.values <- c(t.values, abs(t0$statistic))
				i <- i+1
				if (i >= k) break;
			}
		}
	} else if (sampleFrom=="hack") {
		i <- 0
		repeat {			
			x <- rnorm(n/2)
			y <- rnorm(n/2, mean=d)
			t0 <- t.test(x, y)
			
			n.current <- n
			
			# The hacking stage: increase by `hackIncrease` participants in each group until significant or n-max is reached
			repeat {
				x <- c(x, rnorm(hackIncrease))
				y <- c(y, rnorm(hackIncrease, mean=d))
				n.current <- n.current+2
				t0 <- t.test(x, y)
				if (t0$p.value <= .05) break;
				if (n.current >= n.max) break;
			}

			if (t0$p.value <= .05) {
				ns <- c(ns, n.current)
				t.values <- c(t.values, abs(t0$statistic))
				i <- i+1
				if (i >= k) break;
			}
		}
	}
	
	dfs <- ns-2
	p.values <- pt(t.values, dfs, lower.tail=FALSE)*2
	ds <- (2 * t.values)/sqrt(2 * (dfs+2) - 2)
	
	return(list(t.values=as.vector(t.values), dfs=as.vector(dfs), p.values=as.vector(p.values), ds=as.vector(ds)))
}

#studies.H0 <-  sim_studies(k=5, n=50, sampleFrom="H0")
#studies.H1 <-  sim_studies(k=5, n=50, d=0.5, sampleFrom="H1")
#studies.hack <-  sim_studies(k=5, n=50, d=0, sampleFrom="hack")





p_likelihood <- function(t.values, dfs, meta=TRUE, p.max=.05, plot=FALSE) {
	
	ns <- dfs+2
	p.values <- pt(t.values, dfs, lower.tail=FALSE)*2
	ds <- (2 * t.values)/sqrt(2 * (dfs+2) - 2)
	
	# ---------------------------------------------------------------------
	# compute the meta-analytic ES; 
	# based on blog post by Daniel Lakens: http://daniellakens.blogspot.nl/2014_08_01_archive.html

	if (meta==TRUE) {
		library(metafor)
 
		#Insert effect sizes and sample sizes
		es.d<-ds
		n1<-ns/2
		n2<-ns/2
 
		#Calculate Variance ES
		es.d.v <-(((n1+n2)/(n1*n2))+(es.d^2/(2*(n1+n2))))
		#Calculate Standard Errors ES
		d.se<-sqrt(es.d.v)

		meta <- rma(es.d, es.d.v, method="FE")
		ds <- rep(meta$b, length(ds))
	}
	
	# compute LR for each p value
	LR10 <- c()
	for (i in 1:length(p.values)) {
		LR10 <- c(LR10, pd(p.values[i], ds[i], ns[i], p.max=p.max)/(1/p.max))
	}

	if (plot==TRUE) {
		# ---------------------------------------------------------------------
		#  construct a p-curve for each study
		# All assume the same fixed effect; but due to the different sample sizes, they look different

		ps <- seq(.0001, p.max, length.out=1000)
		plot(NA, xlim=c(0, p.max), ylim=c(0, 100), xlab="p value", ylab="Density")
		abline(h=1/p.max, lty="dotted")  # uniform p-curve under H0

		for (i in 1:length(p.values)) {
			lines(ps, pd(ps, ds[i], ns[i], p.max=p.max), col=i)
			points(p.values[i], pd(p.values[i], ds[i], ns[i], p.max=p.max), col=i, pch=20)
			points(p.values[i], 1/p.max, col=i, pch=20)

			ypos <- (pd(p.values[i], ds[i], ns[i], p.max=p.max)-1/p.max)/2 + 1/p.max
			text(p.values[i], ypos, label=round(LR10[i], 2), col=i, cex=.8)
		}
	}

	# Overall ratio:
	return(prod(LR10))
}


# Test
(studies <- sim_studies(k=10, sampleFrom="hack"))
p_curve_2(studies$p.values)
p_likelihood(t.values=studies$t.values, dfs=studies$dfs, plot=FALSE)
p_likelihood(t.values=studies$t.values, dfs=studies$dfs, meta=FALSE, plot=FALSE)
