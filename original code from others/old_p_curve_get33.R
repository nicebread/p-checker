get_33_curve <- function(type, statistic, df, df2, p.crit=.05, power=1/3) {

    # convert r to t values
	type <- as.character(type)
	statistic[tolower(type)=="r"] <- statistic[tolower(type)=="r"] / sqrt( (1 - statistic[tolower(type)=="r"]^2) / df[tolower(type)=="r"])
	type[tolower(type)=="r"] <- "t"

	statistic <- abs(statistic)

type <- c("f", "f", "f", "t", "r")
statistic <- c(5.1, 6.3, 7.1, 2.3, 0.4)
df <- c(1, 1, 1, 38, 98)
df2 <- c(88, 100, 200, 38, 98)

	ncp <- get_pp_values(type=type, statistic=statistic, df=df, df2=df2, p.crit=.05, power=1/3)$ncp

	res <- data.frame()
	
	# ---------------------------------------------------------------------
	# t-values
	# Critical values,xc, for p=.05, .04, .03, .02 and ,01
	t.crit <- list()
	t.CRIT <- c(.975, .98, .985, .99, .995)
	for (j in 1:5) 
		t.crit[[j]] <- qt(t.CRIT[j], df=df[type=="t"])

	# Probability of a p-value bigger  p=.05, .04, .03, .02 and .01 given p<.05 and ncp=ncp33
	t.pp <- c()
	for (j in 1:5)
		t.pp[j] <- mean((1/power)*(pt(t.crit[[j]], df=df[type=="t"], ncp=ncp[ncp$type=="t", "ncp"])-(1-power)))

	t.pp[1] <- 0
	t.pp <- c(t.pp, 1)
	t.prop <- t.pp[2:6]-t.pp[1:5]
	
	# ---------------------------------------------------------------------
	# F-values
	# Critical values,xc, for p=.05, .04, .03, .02 and ,01
	f.crit <- list()
	f.CRIT <- c(.95, .96, .97, .98, .99)
	for (j in 1:5) 
		f.crit[[j]] <- qf(f.CRIT[j], df1=df[type=="f"], df2=df2[type=="f"])

	# Probability of a p-value bigger  p=.05, .04, .03, .02 and .01 given p<.05 and ncp=ncp33
	f.pp <- c()
	for (j in 1:5)
		f.pp[j] <- mean((1/power)*(pf(f.crit[[j]], df1=df[type=="f"], df2=df2[type=="f"], ncp=ncp[ncp$type=="f", "ncp"])-(1-power)))

	f.pp[1] <- 0
	f.pp <- c(f.pp, 1)
	f.prop <- f.pp[2:6]-f.pp[1:5]

	# ---------------------------------------------------------------------
	# chi2-values
	# Critical values,xc, for p=.05, .04, .03, .02 and ,01
	chi.crit <- list()
	chi.CRIT <- c(.95, .96, .97, .98, .99)
	for (j in 1:5) 
		chi.crit[[j]] <- qt(chi.CRIT[j], df=df[type=="chi2"])

	# Probability of a p-value bigger  p=.05, .04, .03, .02 and .01 given p<.05 and ncp=ncp33
	chi.pp <- c()
	for (j in 1:5)
		chi.pp[j] <- mean((1/power)*(pchisq(chi.crit[[j]], df=df[type=="chi2"], ncp=ncp[ncp$type=="chi2", "ncp"])-(1-power)))

	chi.pp[1] <- 0
	chi.pp <- c(chi.pp, 1)
	chi.prop <- chi.pp[2:6]-chi.pp[1:5]
	
	#TODO: z-values!
}
