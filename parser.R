# ---------------------------------------------------------------------
# Parser for the ES strings
# - parse test statistics and dfs
# - compute correct p value
# - convert all to ES to Z values
# @param x A string with the effect size. Everything before the colon is an identifier for the paper (can be a number, or the name of the paper). Test statistics with the same identifier belong together. You can also skip that part. By default a critical two-tailed p value of .05 is assumed; to override this for a certain test statistic, write at the end of the line: '; crit = .10', for example. Examples:

# x = "MMu (2012) S1: t(88)=2.1; crit = .10; p < .04  # comment"
# x = "1: t(88)=2.1, crit = .10, p < .04  # comment"
# x = "1: t(88)=2.1; p < .04; crit = .10"
# x = "Stapel (2008): r(147)=.246"
# x = "F(1,100)=9.1"
# x = "F(2,210)=4.45"
# x = "Z=3.45"
# x = "chi2(1)=9.1"
# x = "r(77)=.47"



# @param round_up: If the t value is reported as 2.1, it could also be 2.1499 which has been rounded down. If you want to be maximally generous, set this parameter to TRUE, and all test statistics are automatically increased by X.XX49.
parse_ES1 <- function(x, paper_id_fallback="_1", round_up=FALSE) {
	
	library(compute.es)

	W <- c()	# W collects all warnings

	# preprocessing: replace typographic characters
	x <- gsub("–|−", "-", x)
	
	# remove everything after a # sign; remove empty rows; convert to lower
	x <- gsub("#.*$", "", x)
	x <- str_trim(x)
	
	# replace all commas outside of parentheses with semicolons
	if (str_detect(x, "\\)")==TRUE) {
		x1 <- str_match(x, "(.*)\\)(.*)")[2]
		x2 <- str_match(x, "(.*)\\)(.*)")[3]
		x2 <- gsub(",", ";", x2)
		x <- paste0(x1, ")", x2)
	}
	
	# Is it only a comment line? Return NULL
	if (x == "") return(NULL)
	
	split0 <- strsplit(x, ":")[[1]]
	
	# Is a study id provided? 
	if (length(split0) > 1) {
		# separate study id into two parts: 
		# a) First part until year of publication; b) everything after that
		paper_id_full <- as.character(gsub(":", "", split0[1]))
		paper_id <- str_match(paper_id_full, "^.*\\(.*\\)")[1]
		if (is.na(paper_id)) {
			paper_id <- paper_id_full
			study_id <- ""
		} else {
			study_id <- str_trim(str_match(paper_id_full, "\\(.*\\)(.*$)")[2])		
		}
		x2 <- split0[2]
	} else {
		paper_id <- paper_id_fallback
		study_id <- ""
		x2 <- split0[1]
	}

	split1 <- strsplit(gsub(" ", "", x2), ";")[[1]]
	
	# define defaults
	p.crit <- NA
	p.reported <- NA
	reporting.error <- NA
	error.direction <- ""
	one.tailed <- FALSE
	
	# Is a critical p value and/or reported p value provided? Is it one-tailed?
	if (length(split1) > 1) {
		for (i in 2:length(split1)) {

			split1[i] <- tolower(split1[i])
			
			if (str_detect(split1[i], "crit") == TRUE) {
				x3 <- strsplit(gsub(" ", "", split1[i]), "=|<|>")[[1]]
				p.crit <- suppressWarnings(as.numeric(as.character(x3[2])))
			}
			
			if (str_detect(split1[i], "p\\s*(=|<|>)")==TRUE) {
				p.reported <- gsub(" ", "", split1[i])
			}
			
			if (str_detect(split1[i], "one|1t|one-tailed") == TRUE) {
				one.tailed <- TRUE
			}
		}
	}
	
	# set default for p.crit if not defined explicitly
	if (is.na(p.crit)) {
		p.crit <- ifelse(one.tailed==FALSE, .05, .10)
	}
	
	split2 <- strsplit(split1[1], "=")[[1]]
	lhs <- split2[1]
	statistic <- suppressWarnings(as.numeric(split2[2]))
	
	decPlaces <- decplaces(str_trim(split2[2]))
	
	if (round_up==TRUE) {
		statistic <- statistic + sign(statistic)* (4.999 / 10^(decPlaces+1))
	}
	
	# also convert brackets to parentheses
	lhs <- gsub("[", "(", lhs, fixed=TRUE)
	lhs <- gsub("]", ")", lhs, fixed=TRUE)
	
	type <- tolower(strsplit(lhs, "(", fixed=TRUE)[[1]][1])
	dfs <- str_extract(lhs, "\\(.*\\)")
	dfs <- suppressWarnings(as.numeric(strsplit(substring(dfs, 2, nchar(dfs)-1), ",")[[1]]))	# remove parentheses

	# error capturing
	if (!type %in% c("t", "f", "r", "z", "chi2")) {
		W <- c(W, paste0("Test statistic not recognized! ", x))
		return(W)
	}
	if (type != "z" && is.na(as.numeric(dfs[1]))) {
		W <- c(W, paste0("Error in df: ", x))
		return(W)
	}
	if (is.na(as.numeric(statistic))) {
		W <- c(W, paste0("Error in test statistic: ", x))
		return(W)
	}

	# compute the actual p values
	p.value <- NA
	stat <- abs(statistic)
	stat.sign <- sign(statistic)
	n.approx <- NA		# n is approximate because we do not know whether the t-test comes from one sample (n = df+1) or from two samples (n=df+2)
	switch(type,
		"t" = {
			if (length(dfs) != 1) {
				W <- c(W, paste0("t values need exactly one df! ", x))
				return(W)
			}
			t.value <- stat
			p.value <- pt(t.value, dfs, lower.tail=FALSE)*2
			d <- (2*t.value / sqrt(dfs))*stat.sign
			g <- d*(1 - (3/(4 * dfs - 1)))
			n.approx <- dfs+2
			},
		"r" = {
			if (length(dfs) != 1) {
				W <- c(W, paste0("r values need exactly one df (df = n-2)! ", x))
				return(W)
			}
			t.value <- sqrt(dfs) * stat/sqrt(1 - stat^2)
			p.value <- pt(t.value, dfs, lower.tail=FALSE)*2
			d <- stat.sign*(2*stat) / sqrt(1-stat^2)
			g <- d*(1 - (3/(4 * dfs - 1)))
			n.approx <- dfs+2
			},
		"f" = {
			if (length(dfs) != 2) {
				W <- c(W, paste0("F values need exactly two dfs! ", x))
				return(W)
			}
			#if (dfs[1] != 1) warning("First df of F test should be 1 for a focused test!")	
			if (dfs[1] == 1) {
				t.value <- sqrt(stat)
				d <- 2*t.value / sqrt(dfs[2])
				g <- d*(1 - (3/(4 * dfs[2] - 1)))
				n.approx <- dfs[2]+2
			} else {
				d <- NA
				g <- NA
			}
			p.value <- pf(stat, dfs[1], dfs[2], lower.tail=FALSE)
			},
		"z" = {
			p.value <- pnorm(stat, lower.tail=FALSE)*2
			
			# If a number is provided for z it's the sample size
			if (!is.na(dfs[1])) {
				n <- dfs[1]
				d <- (z/sqrt(n))*stat.sign
				g <- d*(1 - (3/(4 * n - 1)))
				n.approx <- n
			} else {
				d <- NA
				g <- NA
			}			
			},
		"chi2" = {
			# If two numbers are provided for chi2, the first are the dfs, the second is the sample size
			p.value <- pchisq(stat, dfs[1], lower.tail=FALSE)
			
			if (dfs[1] == 1 & !is.na(dfs[2])) {
				# code from compute.es package
				n <- dfs[2]
				n.approx <- n
				dfs <- dfs[1]
				r <- sqrt(stat/n)
			    d <- 2 * r * sqrt((n - 1)/(n * (1 - r^2))) * abs(r)/r
				g <- d*(1 - (3/(4 * (n-2) - 1)))
			} else {
				d <- NA
				g <- NA
			}
			}
	)
	
	# test for reporting errors
	# TODO: check both generous *and* non-generous - maybe one of both is correct
	# Or better: computer upper and lower bound of p-value: t(47)=2.1 --> from 2.05 to 2.1499
	# and check whether the reported p value falls into the interval
	
	p.actual <- ifelse(one.tailed==FALSE, p.value, p.value/2)
	p.reported.num <- suppressWarnings(as.numeric(str_split(p.reported, "=|<|>|<=|>=")[[1]][2]))

	if (!is.na(p.reported) & !is.na(p.reported.num)) {
		# check for inequality
		if (str_detect(p.reported, "<")) {
			if (p.actual >= p.reported.num) {
				reporting.error <- TRUE
				error.direction <- "smaller"
			} else {
				reporting.error <- FALSE
				error.direction <- ""
			}
		}
		
		if (str_detect(p.reported, "<=")) {
			if (p.actual >= p.reported.num) {
				reporting.error <- TRUE
				error.direction <- "smaller"
			} else {
				reporting.error <- FALSE
				error.direction <- ""
			}
		}
		
		if (str_detect(p.reported, ">")) {
			if (p.actual <= p.reported.num) {
				reporting.error <- TRUE
				error.direction <- "larger"
			} else {
				reporting.error <- FALSE
				error.direction <- ""
			}
		}
			
		if (str_detect(p.reported, "p=")) {
						
			dec <- decplaces(str_split(p.reported, "=")[[1]][2])
			p.actual <- round(p.actual, dec)
			
			if (p.reported.num == p.actual) {
				reporting.error <- FALSE
				error.direction <- ""
			} else {
				reporting.error <- TRUE
				error.direction <- ifelse(p.reported.num > p.actual, "larger", "smaller")
			}
		}
		
	}
	
	res <- data.frame(
		paper_id = as.character(paper_id),
		study_id = as.character(study_id),
		focal	= ifelse(substr(as.character(paper_id), 1, 1) == "_", FALSE, TRUE),
		type	= type, 
		df1 	= dfs[1], 
		df2 	= ifelse(length(dfs)>1, dfs[2], NA), 
		d		= d,
		g		= g,
		n.approx = n.approx,
		statistic = statistic,		
		p.value	= p.value,
		p.value.one	= p.value/2,
		p.reported = p.reported,		
		p.crit	= p.crit,
		significant = p.value < p.crit,
		one.tailed = one.tailed,
		reporting.error = reporting.error,
		error.direction = error.direction
	)
	
	attr(res, "warnings") <- W

	return(res)
}



# A vectorized version of the parse_ES1 function
parse_ES <- function(x, round_up=FALSE) {
	
	# split input string at line break & remove empty rows
	# Preprocessing: remove everything after a # sign; remove empty rows
	txt <- str_trim(strsplit(x, "\n")[[1]])
	txt <- gsub("#.*$", "", txt)
	txt <- str_trim(txt)
	txt <- txt[txt != ""]
	
	if (txt[1]=="" | length(txt)==0) return(NULL)
	
	res <- data.frame()
	Ws <- c()
	for (i in 1:length(txt)) {
		parsed <- parse_ES1(txt[i], paper_id_fallback = paste0(".", i), round_up=round_up)		
		if (!is.null(parsed) & is.data.frame(parsed)) res <- rbind(res, parsed)
			
		# collect errors
		if (length(attr(parsed, "warnings")) > 0) Ws <- c(Ws, attr(parsed, "warnings"))
		if (is.character(parsed)) Ws <- c(Ws, parsed)
	}
	
	if (nrow(res) == 0) return(NULL)
	res2 <- cbind(ID <- 1:nrow(res), res)
	
	attr(res2, "warnings") <- Ws
	return(res2)
}



