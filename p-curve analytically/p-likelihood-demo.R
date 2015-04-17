# A set of p-hacked studies ...
t.values <- c(2.1, 1.99, 2.15)
dfs <- c(38, 98, 29)

# A set of studies with evidential value
t.values <- c(3.1, 2.99, 2.15)
dfs <- c(138, 198, 129)


ns <- dfs+2
(p.values <- pt(t.values, dfs, lower.tail=FALSE)*2)
(ds <- (2 * t.values)/sqrt(2 * (dfs+2) - 2))

# ---------------------------------------------------------------------
# compute the meta-analytic ES; 
# based on blog post by Daniel Lakens: http://daniellakens.blogspot.nl/2014_08_01_archive.html

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
d.meta <- meta$b

# ---------------------------------------------------------------------
#  construct a p-curve for each study
# All assume the same fixed effect; but due to the different sample sizes, they look different

p.max <- .05
ps <- seq(.0001, p.max, length.out=1000)
plot(NA, xlim=c(0, p.max), ylim=c(0, 100), xlab="p value", ylab="Density")
abline(h=1/p.max, lty="dotted")  # uniform p-curve under H0

LR10 <- c()
for (i in 1:length(p.values)) {
	lines(ps, pd(ps, d.meta, ns[i], p.max=p.max), col=i)
	points(p.values[i], pd(p.values[i], d.meta, ns[i], p.max=p.max), col=i, pch=20)
	points(p.values[i], 1/p.max, col=i, pch=20)

	LR10 <- c(LR10, pd(p.values[i], d.meta, ns[i], p.max=p.max)/(1/p.max))
	ypos <- (pd(p.values[i], d.meta, ns[i], p.max=p.max)-1/p.max)/2 + 1/p.max
	text(p.values[i], ypos, label=round(LR10[i], 2), col=i, cex=.8)
}

# Overall ratio:
prod(LR10)


