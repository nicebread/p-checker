#source("p-likelihood-sim.R", echo=TRUE)
library(dplyr)
library(reshape2)
source("p-curve-fun.R")



# ---------------------------------------------------------------------
# Simulation

R <- 100	# Simulations within each experimental condition
res.H0 <- data.frame()
res.H1 <- data.frame()
res.hack0.2 <- data.frame()
res.hack0 <- data.frame()
for (n in c(40, 80, 200)) {
	for (k in c(5, 10, 20)) {
		print(Sys.time())
		for (j in 1:R) {
			studies.H0 <-   sim_studies(k=k, n=n, sampleFrom="H0")
			studies.H1 <-   sim_studies(k=k, n=n, d=0.5, sampleFrom="H1")
			studies.hack0.2 <- sim_studies(k=k, n=n, d=0.2, sampleFrom="hack")	# hacking a small effect
			studies.hack0 <- sim_studies(k=k, n=n, d=0, sampleFrom="hack")	# hacking a null effect
			res.H0 <- rbind(res.H0, data.frame(
				n = n,
				k = k,
				LR_meta=p_likelihood(t.values=studies.H0$t.values, dfs=studies.H0$dfs, plot=FALSE),
				LR_single=p_likelihood(t.values=studies.H0$t.values, dfs=studies.H0$dfs, plot=FALSE, meta=FALSE),
				p_curve_3(studies.H0$p.values)
			))
			res.H1 <- rbind(res.H1, data.frame(
				n = n,
				k = k,
				LR_meta=p_likelihood(t.values=studies.H1$t.values, dfs=studies.H1$dfs, plot=FALSE),
				LR_single=p_likelihood(t.values=studies.H1$t.values, dfs=studies.H1$dfs, plot=FALSE, meta=FALSE),
				p_curve_3(studies.H1$p.values)
			))
			res.hack0.2 <- rbind(res.hack0.2, data.frame(
				n = n,
				k = k,
				LR_meta=p_likelihood(t.values=studies.hack0.2$t.values, dfs=studies.hack0.2$dfs, plot=FALSE),
				LR_single=p_likelihood(t.values=studies.hack0.2$t.values, dfs=studies.hack0.2$dfs, plot=FALSE, meta=FALSE),
				p_curve_3(studies.hack0.2$p.values)
			))
			res.hack0 <- rbind(res.hack0, data.frame(
				n = n,
				k = k,
				LR_meta=p_likelihood(t.values=studies.hack0$t.values, dfs=studies.hack0$dfs, plot=FALSE),
				LR_single=p_likelihood(t.values=studies.hack0$t.values, dfs=studies.hack0$dfs, plot=FALSE, meta=FALSE),
				p_curve_3(studies.hack0$p.values)
			))
	
			print(paste0(n, ":", k, ":", j))
		}
	}
}


res.H0$reality <- "H0"
res.H1$reality <- "H1"
res.hack0$reality <- "hack0"
res.hack0.2$reality <- "hack0.2"
res <- rbind(res.hack0, res.hack0.2, res.H0, res.H1)
res$reality <- factor(res$reality, levels=c("hack0", "hack0.2", "H0", "H1"))

save(res, file="res_5000_4.RData")
#load(file="res.RData")

