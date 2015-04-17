source("p-curve-fun.R")

# Show p-curve
t.value <- 2.53
n <- 40
df <- n-2
(p.value <- pt(t.value, df, lower.tail=FALSE)*2)
(d <- (2 * t.value)/sqrt(2 * n - 2))

# These should give the same
pd(p.value, d, n)
ptd(p.value, t.value, df)

p.max <- .05
d <- 0.5
n <- 40

pcurve.plot <- function(d=0.5, n=40, p.max=.05, add=FALSE, col=1) {
	ps <- seq(.0001, p.max, length.out=1000)
	if (add==FALSE) {
		plot(NA, xlim=c(0, p.max), ylim=c(0, 100), xlab="p value", ylab="Density")		
	}
	abline(h=1/p.max, lty="dotted", col=col)  # uniform p-curve under H0	
	
	pd1 <- create_pd(d, n, p.max)
	lines(ps, pd1(ps), col=col)
}


pcurve.plot(d=0.5, n=40, p.max=.05)
pcurve.plot(d=0.5, n=40, p.max=.10, add=TRUE, col=2)

# p = .045, 1-tailed
pd(p=.045, d=0.5, n=40, p.max=.05) / (1/.05)

# p = .09, 2-tailed
pd(p=.09, d=0.5, n=40, p.max=.10) / (1/.10)