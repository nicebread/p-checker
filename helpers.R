# ---------------------------------------------------------------------
# Helper functions

COLORS <- list(
  BLUE = "#3399f3"
)
# Helper: Transform correlation to Fisher's Z
r2Z <- function(r) {return(0.5 * log((1 + r)/(1 - r)))}


# Helper: REcode  Fisher's to correlation
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