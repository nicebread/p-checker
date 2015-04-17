library(BayesFactor)

d <- 0.5
R <- 2000
n <- 40

# ---------------------------------------------------------------------
#  Simulate distribution of p values
X <- rnorm(1000000)
Y <- rnorm(1000000) + d

ps <- vector(mode = "integer", length = R)
pb <- txtProgressBar(1, R, style=3)
for (i in 1:R) {
	x <- X[sample(1:length(X), n/2, replace=FALSE)]
	y <- Y[sample(1:length(Y), n/2, replace=FALSE)]
	
	t0 <- t.test(x, y, var.equal=TRUE)
	ps[i] <- t0$p.value
	setTxtProgressBar(pb, i)
}
close(pb)

#save(ps, file="ps_0.5.RData")
#load(file="ps_0.5.RData")

# ---------------------------------------------------------------------
# Compute density of p values semi-analytically; compare them with the simulated p values

# plot for checking
p.max <- .05
hist(ps[ps<p.max], breaks=500, freq=FALSE, xlab="p value")

p <- seq(.0001, p.max, length.out=1000)
lines(p, pd(p, d, n, p.max=p.max), col="red", lwd=2)



# ---------------------------------------------------------------------
# Example

t.value <- 2.53
n <- 40
df <- n-2
(p.value <- pt(t.value, df, lower.tail=FALSE)*2)
(d <- (2 * t.value)/sqrt(2 * n - 2))

# ---------------------------------------------------------------------
# The plot

p.max <- .05
LR10 <- pd(p.value, d, n)/(1/p.max)

pd(p.value, d, n)
ptd(p.value, t.value, df)

ps <- seq(.0001, p.max, length.out=1000)
plot(ps, pd(ps, d, n), type="l", ylim=c(0, 100), main=bquote(LR[10]~"="~.(round(LR10, 2))), xlab="p value", ylab="Density")
#lines(ps, g(ps, d, n), col="green")

abline(h=1/p.max, lty="dotted")  # p-uniform line
abline(h=1, lty="dotted", col="red")  # p-uniform line
points(p.value, pd(p.value, d, n))
points(p.value, 1/p.max)


# ---------------------------------------------------------------------
# The plot  -full range
# This changes the shape of both curves ()

p.max <- 1
(LR10 <- pd(p.value, d, n, p.max=p.max)/(1/p.max))

ps <- seq(.0001, p.max, length.out=1000)
plot(ps, pd(ps, d, n, p.max=p.max), type="l", ylim=c(0, 100), main=bquote(LR[10]~"="~.(round(LR10, 2))), xlab="p value", ylab="Density")
#lines(ps, g(ps, d, n), col="green")

abline(h=1/p.max, lty="dotted")  # p-uniform line
abline(h=1, lty="dotted", col="red")  # p-uniform line
points(p.value, pd(p.value, d, n, p.max=p.max))
points(p.value, 1/p.max)
