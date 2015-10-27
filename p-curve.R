# ---------------------------------------------------------------------
# These p-curve functions are partially copied, partially adapted from Uri Simonsohn's (uws@wharton.upenn.edu) original p-curve functions
# http://p-curve.com/Supplement/Rcode_other/R%20Code%20behind%20p-curve%20app%203.0%20-%20distributable.R


# ---------------------------------------------------------------------
# p-curve-app 3.0 functions

# functions that find noncentrality parameter for t,f,chi distributions that gives 33% power for those d.f.

#t-test
ncp33t <- function(df, power=1/3, p.crit=.05) {      
      xc=qt(p=1-p.crit/2, df=df)
      #Find noncentrality parameter (ncp) that leads 33% power to obtain xc
	  f = function(delta, pr, x, df) pt(x, df = df, ncp = delta) - (1-power)
	  out = uniroot(f, c(0, 37.62), x = xc, df = df)  
	  return(out$root) 
}

#F-test
ncp33f <- function(df1, df2, power=1/3, p.crit=.05) {      
	xc=qf(p=1-p.crit,df1=df1,df2=df2)
	f = function(delta, pr, x, df1,df2) pf(x, df1 = df1, df2=df2, ncp = delta) - (1-power)
	out = uniroot(f, c(0, 37.62), x = xc, df1=df1, df2=df2)  
	return(out$root)       
}

#chi-square
ncp33chi <- function(df, power=1/3, p.crit=.05) {      
	xc=qchisq(p=1-p.crit, df=df)
	#Find noncentrality parameter (ncp) that leads 33% power to obtain xc
	f = function(delta, pr, x, df) pchisq(x, df = df, ncp = delta) - (1-power)
	out = uniroot(f, c(0, 37.62), x = xc, df = df)  
	return(out$root)
}




get_pp_values <- function(type, statistic, df, df2, p.crit=.05, power=1/3) {

    # convert r to t values
	type <- as.character(type)
	statistic[tolower(type)=="r"] <- statistic[tolower(type)=="r"] / sqrt( (1 - statistic[tolower(type)=="r"]^2) / df[tolower(type)=="r"])
	type[tolower(type)=="r"] <- "t"
	
	statistic <- abs(statistic)

	res <- data.frame()
	ncp <- data.frame()
	for (i in 1:length(type)) {
		switch(tolower(type[i]), 
			"t" = {
				p <- 2*(1-pt(abs(statistic[i]),df=df[i]))
				ppr <- p*(1/p.crit)	# pp-value for right-skew 
				ppl <- 1-ppr		# pp-value for left-skew
	  	        ncp33 <- ncp33t(df[i], power=power, p.crit=p.crit)
	  	        pp33 <- (pt(statistic[i],  df=df[i], ncp=ncp33)-(1-power))*(1/power)
				},
			"f" = {
				p <- 1-pf(abs(statistic[i]), df1=df[i], df2=df2[i])
				ppr <- p*(1/p.crit)	# pp-value for right-skew 
				ppl <- 1-ppr		# pp-value for left-skew
	  	        ncp33 <- ncp33f(df1=df[i], df2=df2[i], power=power, p.crit=p.crit)
	  	        pp33 <- (pf(statistic[i], df1=df[i], df2=df2[i],  ncp=ncp33)-(1-power))*(1/power)
				},
			"z" = {
				p <- 2*(1-pnorm(abs(statistic[i])))
				ppr <- p*(1/p.crit)	# pp-value for right-skew 
				ppl <- 1-ppr		# pp-value for left-skew
				
				# TODO: Adjust 1.5285687 for different power and p.crit!
				ncp33 <- 1.5285687				
	  	        pp33 <- (pnorm(statistic[i], mean=ncp33, sd=1)-(1-power))*(1/power)   # Compute pp33-values using the 'ncp' 1.5285687 which gives the normal 33% power
				},	
			"chi2" = {
				p <- 1-pchisq(abs(statistic[i]), df=df[i])
				ppr <- p*(1/p.crit)	# pp-value for right-skew 
				ppl <- 1-ppr		# pp-value for left-skew
	  	        ncp33 <- ncp33chi(df[i], power=power, p.crit=p.crit)
	  	        pp33 <- (pchisq(statistic[i],  df=df[i], ncp=ncp33)-(1-power))*(1/power)
			}
		)
		res <- rbind(res, data.frame(p=p, ppr=ppr, ppl=ppl, pp33=pp33))
		ncp <- rbind(ncp, data.frame(type=type[i], df=df[i], df2=df2[i], ncp=ncp33))
	}
		
	if (nrow(res) > 0) {
		# clamp to extreme values	
		res$ppr <- clamp(res$ppr, MIN=.00001, MAX=.99999)
		res$ppl <- clamp(res$ppl, MIN=.00001, MAX=.99999)
		res$pp33 <- clamp(res$pp33, MIN=.00001, MAX=.99999)
		
		# remove non-significant values
		res[res$p > p.crit, ] <- NA
	
		return(list(res=res, ncp=ncp))
	} else {
		return(NULL)
	}
}





# ---------------------------------------------------------------------
# New p-curve computation (p-curve app 3.0, http://www.p-curve.com/app3/)
p_curve_3 <- function(pps) {

	pps <- na.omit(pps)

	# STOUFFER: Overall tests aggregating pp-values
	ktot <- sum(!is.na(pps$ppr))
	Z_ppr <- sum(qnorm(pps$ppr))/sqrt(ktot)          # right skew
	Z_ppl <- sum(qnorm(pps$ppl))/sqrt(ktot)          # left skew
	Z_pp33<- sum(qnorm(pps$pp33))/sqrt(ktot)         # 33%
	
	p_ppr <- pnorm(Z_ppr)
	p_ppl <- pnorm(Z_ppl)
	p_pp33<- pnorm(Z_pp33)

	return(list(
		Z_evidence = Z_ppr, 
		p_evidence = p_ppr, 
		Z_hack = Z_ppl, 
		p_hack = p_ppl, 
		Z_lack = Z_pp33, 
		p_lack = p_pp33,
		inconclusive = ifelse(p_ppr>.05 & p_ppl>.05 & p_pp33>.05, TRUE, FALSE)))
}


# ---------------------------------------------------------------------
# Old p-curve computation (p-curve app 2.0, http://www.p-curve.com/app2/)
p_curve_2 <- function(pps) {

	pps <- na.omit(pps)
	
	df <- 2*sum(nrow(pps))

	chi2_evidence <- -2*sum(log(pps$ppr), na.rm=TRUE)
	p_evidence <- pchisq(chi2_evidence, df=df, lower.tail=FALSE)

	chi2_hack <- -2*sum(log(pps$ppl), na.rm=TRUE)
	p_hack <- pchisq(chi2_hack, df=df, lower.tail=FALSE)

	chi2_lack <- -2*sum(log(pps$pp33), na.rm=TRUE)
	p_lack <- pchisq(chi2_lack, df=df, lower.tail=FALSE)

	return(list(
		chi2_evidence = chi2_evidence, 
		p_evidence = p_evidence, 
		chi2_hack = chi2_hack, 
		p_hack = p_hack, 
		chi2_lack = chi2_lack, 
		p_lack = p_lack,
		df = df, 
		inconclusive = ifelse(p_evidence>.05 & p_hack>.05 & p_lack>.05, TRUE, FALSE)))
}



# # ---------------------------------------------------------------------
# # Compute green 33%-curve of plot
#
# #8 Green line (Expected p-curve for 33% power)
#      #8.1 t-tests
#     if (length(t.value.sig)>0)  #if nonempty compute pp-values
#     {
#           #Critical values,xc, for p=.05, .04, .03, .02 and ,01
#               t.x5=qt(.975,df=t.df.sig); t.x4=qt(.98, df=t.df.sig); t.x3=qt(.985,df=t.df.sig);
#               t.x2=qt(.99, df=t.df.sig); t.x1=qt(.995,df=t.df.sig)
#           #For Binomial test
# 	  	      t.x25=qt(.9875,df=t.df.sig)                            #critical value for t-tests to get p=.025
# 			      t.plow=1- 3*(pt(t.x25,df=t.df.sig, ncp=t.ncp33)-2/3)   #prob(p<.025 | ncp33% & p<.05)
#
#           #Probabilty of a p-value bigger  p=.05, .04, .03, .02 and .01 given p<.05 and ncp=ncp33
#               t.pp4=3*(pt(t.x4,df=t.df.sig, ncp=t.ncp33)-2/3)
#               t.pp3=3*(pt(t.x3,df=t.df.sig, ncp=t.ncp33)-2/3)
#               t.pp2=3*(pt(t.x2,df=t.df.sig, ncp=t.ncp33)-2/3)
#               t.pp1=3*(pt(t.x1,df=t.df.sig, ncp=t.ncp33)-2/3)
#          #within bins proportions
#               t.prop5=mean(t.pp4);
#               t.prop4=mean(t.pp3-t.pp4);
#               t.prop3=mean(t.pp2-t.pp3);
#               t.prop2=mean(t.pp1-t.pp2);
#               t.prop1=mean(1-t.pp1)
#     }
#     #8.2 f-tests
#           if (length(f.value.sig)>0)  #if nonempty compute pp-values
#           {
#           #Critical values,xc, for p=.05, .04, .03, .02 and ,01
#               f.x5=qf(.95,df1=f.df1.sig, df2=f.df2.sig);  f.x4=qf(.96,df1=f.df1.sig, df2=f.df2.sig); f.x3=qf(.97,,df1=f.df1.sig, df2=f.df2.sig);
#               f.x2=qf(.98, df1=f.df1.sig, df2=f.df2.sig); f.x1=qf(.99,df1=f.df1.sig, df2=f.df2.sig)
#           #For binomial test
# 			        f.x25 =qf(.975,df1=f.df1.sig, df2=f.df2.sig)                          #Critical F value for p=.025
#               f.plow=1-3*(pf(f.x25,df1=f.df1.sig, df2=f.df2.sig, ncp=f.ncp33)-2/3)  #Prob(p<.025|ncp33% & p<.05)
#
#
#           #Probabilty of a p-value bigger  p=.05, .04, .03, .02 and .01 given p<.05 and ncp=ncp33
#               f.pp4=3*(pf(f.x4,df1=f.df1.sig, df2=f.df2.sig, ncp=f.ncp33)-2/3)
#               f.pp3=3*(pf(f.x3,df1=f.df1.sig, df2=f.df2.sig, ncp=f.ncp33)-2/3)
#               f.pp2=3*(pf(f.x2,df1=f.df1.sig, df2=f.df2.sig, ncp=f.ncp33)-2/3)
#               f.pp1=3*(pf(f.x1,df1=f.df1.sig, df2=f.df2.sig, ncp=f.ncp33)-2/3)
#
#           #within bins proportions
#               f.prop5=mean(f.pp4);
#               f.prop4=mean(f.pp3-f.pp4);
#               f.prop3=mean(f.pp2-f.pp3);
#               f.prop2=mean(f.pp1-f.pp2);
#               f.prop1=mean(1-f.pp1)
#           }
#     #8.3 chi-tests
#         if (length(c.value.sig)>0)  #if nonempty compute pp-values
#         {
#
#           #Critical values,xc, for p=.05, .04, .03, .02 and ,01
#               c.x5=qchisq(.95,df=c.df.sig); c.x4=qchisq(.96, df=c.df.sig); c.x3=qchisq(.97,df=c.df.sig);
#               c.x2=qchisq(.98, df=c.df.sig); c.x1=qchisq(.99,df=c.df.sig)
#
# 		  #For binomial test
# 		      c.x25 =qchisq(.975,df=c.df.sig)                                      #Critical x2 value for p=.025
# 	          c.plow=1-3*(pchisq(c.x25,df=c.df.sig, ncp=c.ncp33)-2/3)              #Prob(p<.025|ncp33% & p<.05)
#
# 			  #Probabilty of a p-value bigger  p=.05, .04, .03, .02 and .01 given p<.05 and ncp=ncp33
#               c.pp4=3*(pchisq(c.x4,df=c.df.sig, ncp=c.ncp33)-2/3)
#               c.pp3=3*(pchisq(c.x3,df=c.df.sig, ncp=c.ncp33)-2/3)
#               c.pp2=3*(pchisq(c.x2,df=c.df.sig, ncp=c.ncp33)-2/3)
#               c.pp1=3*(pchisq(c.x1,df=c.df.sig, ncp=c.ncp33)-2/3)
#
#           #within bins proportions
#               c.prop5=mean(c.pp4); c.prop4=mean(c.pp3-c.pp4); c.prop3=mean(c.pp2-c.pp3);
#               c.prop2=mean(c.pp1-c.pp2); c.prop1=mean(1-c.pp1)
#         }
#     #8.4 z-tests
#       if (length(z.value.sig)>0)  #if nonempty compute pp-values
#       {
#           #Critical values,xc, for p=.05, .04, .03, .02 and ,01
#               z.x5=qnorm(.975); z.x4=qnorm(.98); z.x3=qnorm(.985); z.x2=qnorm(.99); z.x1=qnorm(.995)
# 		  # For Binomial test
# 			    z.x25 =qnorm(.9825)                                         #Critical x2 value for p=.025
# 		      z.plow=1-3*(pnorm(z.x25,mean=1.5285687,sd=1)-2/3)           #Prob(p<.025|ncp33% & p<.05)
#         #Probabilty of a p-value bigger  p=.05, .04, .03, .02 and .01, given p<.05 and ncp=ncp33
#               z.pp4=3*(pnorm(z.x4,mean=1.5285687,sd=1)-2/3)
#               z.pp3=3*(pnorm(z.x3,mean=1.5285687,sd=1)-2/3)
#               z.pp2=3*(pnorm(z.x2,mean=1.5285687,sd=1)-2/3)
#               z.pp1=3*(pnorm(z.x1,mean=1.5285687,sd=1)-2/3)
#         #within bins proportions
#               z.prop5=z.pp4; z.prop4=z.pp3-z.pp4; z.prop3=z.pp2-z.pp3; z.prop2=z.pp1-z.pp2; z.prop1=1-z.pp1
#       }
#
#
# #9 combine t,F,chi,Z
#           #proportion of all tests that are of each type
#             t.share=length(t.value.sig)/ktot
#             f.share=length(f.value.sig)/ktot
#             c.share=length(c.value.sig)/ktot
#             z.share=length(z.value.sig)/ktot
#
#           #Average proportions within the 4 types of tests
#             t.props=c(t.prop1, t.prop2, t.prop3, t.prop4, t.prop5)
#             f.props=c(f.prop1, f.prop2, f.prop3, f.prop4, f.prop5)
#             c.props=c(c.prop1, c.prop2, c.prop3, c.prop4, c.prop5)
#             z.props=c(z.prop1, z.prop2, z.prop3, z.prop4, z.prop5)
#
#           #overall proportions (i.e.., THE GREEN LINE)
#             green=100*(t.props*t.share + f.props*f.share + c.props*c.share + z.props*z.share)
#


get_33_curve <- function(type, statistic, df, df2, p.crit=.05, power=1/3) {

    # convert r to t values
	type <- as.character(type)
	statistic[tolower(type)=="r"] <- statistic[tolower(type)=="r"] / sqrt( (1 - statistic[tolower(type)=="r"]^2) / df[tolower(type)=="r"])
	type[tolower(type)=="r"] <- "t"

	statistic <- abs(statistic)

type <- rep("f", 3)
statistic <- c(5.1, 6.3, 7.1)
df <- c(1, 1, 1)
df2 <- c(88, 100, 200)

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
