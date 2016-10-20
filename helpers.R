# ---------------------------------------------------------------------
# Helper functions

COLORS <- list(
  BLUE = "#3399f3",
  GREEN = "#0BBC2B"
)
# Helper: Transform correlation to Fisher's Z
r2Z <- function(r) {return(0.5 * log((1 + r)/(1 - r)))}

# Helper: Recode Fisher's z to correlation
Z2r <- function(Z) {return((exp(2*Z)-1)/(exp(2*Z)+1))}


# simple wrapper: formats a number in f.2 format
f2 <- function(x, digits=2, prepoint=0, skipZero=FALSE) {
	
	if (skipZero == TRUE) {zero <- "."} else {zero <- "0."}
	
	if (length(dim(x)) == 2) {
		apply(x, 2, function(x2) {gsub("0.", zero, sprintf(paste("%",prepoint,".",digits,"f",sep=""), x2) , fixed=TRUE)})
	} else {
		gsub("0.", zero, sprintf(paste("%",prepoint,".",digits,"f",sep=""), x) , fixed=TRUE)
	}
}


# Given a column: if numeric, convert to formatted string; otherwise, return unchanged
format_num <- function(col, digits) {if (is.numeric(col)) {sprintf(paste0('%1.', digits, 'f'), col)} else {col}}
format_num2 <- function(col) {if (is.numeric(col)) {sprintf(paste0('%1.0f'), col)} else {col}}


# nicely formats a p-value
p0 <- function(x, digits=3) {
	if (is.na(x)) return("NA")
	if (x >= .1^digits) return(paste0("p = ", f2(x, digits, skipZero=TRUE)))
	if (x <  .1^digits) return(paste0("p < ", f2(.1^digits, digits, skipZero=TRUE)))
}
p <- Vectorize(p0)



# Get the number of decimal places
# Taken from http://stackoverflow.com/questions/5173692/how-to-return-number-of-decimal-places-in-r
decimalplaces <- function(x) {
	if (is.na(x)) return(0)
    if ((x %% 1) != 0) {nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed=TRUE)[[1]][[2]])} 
	else {return(0)}
}

decplaces <- function(x) {nchar(str_extract(x, "\\.\\d*$"))-1}

clamp <- function(x, MIN=.00001, MAX=.99999) {x[x<MIN] <- MIN; x[x>MAX] <- MAX; x}

str_match_empty <- function(..., position=2) {
	res <- str_match(...)[position]
	if (is.na(res)) res <- ""
	res
}

logScaleHack <- function(v){
  isNotPowOf10 <- (log(v, 10) %% 1) != 0
  v[isNotPowOf10] <- v[isNotPowOf10] + 0.1
  v
}

logScaleTicks <- function(v){
  rv <- range(v)
  rv <- c(max(1,rv[1]), max(1,rv[2]))
  exponent <- floor(log(rv, 10))
  mult <- round(rv / (10^(exponent)))
  diff <- exponent[2] - exponent[1]
  if( diff == 0 ){
    res <- (mult[1]:mult[2])*10^exponent[1]
  } else if( diff == 1){
    res <- c((mult[1]:9)*10^exponent[1], (1:mult[2])*10^exponent[2])
  } else {
    res <- c(
      (mult[1]:9)*10^exponent[1],
      rep(1:9, diff-1) * 10^(rep((exponent[1]+1):(exponent[2]-1), each=9)),
      (1:mult[2])*10^exponent[2]
    )
  }
  logScaleHack(res)
}

logScaleLimits <- function(v)
{
  rv <- range(v) * c(0.9, 1.1)
  rv <- c(max(1,rv[1]), max(1,rv[2]))

  exponents <- floor(log(rv, 10))
  results <- rv / (10^exponents)
  results <- c(floor(results[1]), ceiling(results[2]))
  
  diff <- abs(exponents[2] - exponents[1])
  if(diff == 0){
    results <- c(1,10)
  } else if(diff == 1){
    if(results[1] == 1){
      results[2] <- 10
    } else {
      if(results[1] > 10 - results[2]){
        results[2] <- 10
      } else {
        results[1] <- 1
      }
    }
  }
  
  results * (10^exponents)
}



# ---------------------------------------------------------------------
# Meta-Regression helper functions

tidyRMA <- function(RMA) {
	res <- data.frame(
			  term = if(length(RMA$b)==1) {"b0"} else {c("b0", "b1")},
			  estimate = as.vector(RMA$b),
			  std.error = RMA$se,
			  statistic = RMA$zval,
			  p.value = RMA$pval,
			  conf.low = RMA$ci.lb,
			  conf.high = RMA$ci.ub
			)
	rownames(res) <- NULL
	return(res)
}

tidyLM <- function(...) {
	res <- tidy(..., conf.int=TRUE)
	res$term[res$term == "(Intercept)"] <- "b0"
	res$term[res$term == "d.se"] <- "b1"
	res$term[res$term == "d.var"] <- "b1"
	return(res)
}

# this function combines tidyRMA and tidyLM
tidyMR <- function(x) {
	if (class(x)[1] == "lm") {
		return(tidyLM(x))
	} else {
		return(tidyRMA(x))
	}
}