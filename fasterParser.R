library(stringi)
library(dplyr)

decimals <- function(str) {
  locations <- stri_locate_first_fixed(str, '.')[,1]
  decs <- nchar(str) - locations
  decs[is.na(locations)] <- 0
  decs
}

parse_ES <- function(txt, round_up = FALSE) {

  # split into lines
  txt.lines <- stri_split_lines(txt)[[1]]

  # remove all comments
  txt.lines <-  stri_replace_first_regex(txt.lines, '#.*$', '')

  # convert chains of whitespace characters to a single space
  txt.lines <- stri_replace_all_regex(txt.lines, "\\s+", " ")

  # trim all lines
  txt.lines <- stri_trim_both(txt.lines)

  # find all non empty lines (indices_not_empty represents correct line numbers!)
  indices_not_empty <- which(!stri_isempty(txt.lines))

  # remove all empty lines
  txt.lines <- txt.lines[indices_not_empty]

  # number of (non-empty) lines
  nlines <- length(txt.lines)

  # allocate space for error messages
  errors <- rep("", nlines)

  # definition of all column names of output matrix
  TYPE            <-  1
  DF1             <-  2
  DF2             <-  3
  STAT            <-  4
  SIGN            <-  5 # e.g. -1 if stat is -5.3
  P.REPORTED      <-  6 # e.g. 0.12
  P.REPORTED.DECS <-  7 # e.g. 2 if p.reported is 0.12
  P.COMP          <-  8
  CRIT.VALUE      <-  9
  ONE.TAILED      <- 10
  P.VALUE         <- 11
  P.VALUE.ONE     <- 12
  P.ACTUAL        <- 13
  G               <- 14
  D               <- 15
  N.APPROX        <- 16
  SIGNIFICANT     <- 17
  REP.ERROR       <- 18
  REP.DIRECTION   <- 19
  IS.FOCAL        <- 20
  PARSE.ERROR     <- 21

  # output matrix for all data in numeric form
  BIG <- matrix(NA, nrow = nlines, ncol = 21)

  # find labels and extract them
  extraction <- stri_match_first_regex(txt.lines, '^ *(.*?) *(?:(?<=\\)) *(.+?) *)?: *')
  txt.lines.edited <- txt.lines
  indices_not_na <- which(!is.na(extraction[,1]))
  if(length(indices_not_na) > 0) {
    txt.lines.edited[indices_not_na] <- stri_replace_first_fixed(txt.lines.edited[indices_not_na], extraction[indices_not_na,1], '')
  }

  # init paper ids with auto id
  PAPER_ID <- paste0(".", as.character(1:nlines))

  # init study ids with NA
  STUDY_ID <- rep("", nlines)

  # get indices for certain rows
  indices_paper_id <- which(!is.na(extraction[,2]))
  indices_study_id <- which(!is.na(extraction[,3]))

  # set paper id respectively
  if(length(indices_paper_id) > 0)
    PAPER_ID[indices_paper_id] <- extraction[indices_paper_id, 2]

  if(length(indices_study_id))
    STUDY_ID[indices_study_id] <- extraction[indices_study_id, 3]

  # is study focal? indicated by underscore at start of id
  BIG[,IS.FOCAL] <- ( stri_sub(PAPER_ID, 1, 1) != "_" )

  # definition of statistic types and typestrings array
  TYPE_T <- 1
  TYPE_CHI2 <- 2
  TYPE_F <- 3
  TYPE_R <- 4
  TYPE_Z <- 5
  TYPESTRINGS <- c('t','chi2','f','r','z')

  # find statistic and extract it
  extraction <- stri_match_first_regex(txt.lines.edited, ' *\\b(t|chi2|f|r|z)(?: *\\( *((?:\\d*\\.)?\\d+)(?: *, *((?:\\d*\\.)?\\d+))? *\\))? *= *(-?(?:\\d*\\.)?\\d+)[ ,;]*', case_insensitive=TRUE)
  indices_not_na <- which(!is.na(extraction[,1]))
  if( length(indices_not_na) > 0){
    txt.lines.edited[indices_not_na] <- stri_replace_first_fixed(txt.lines.edited[indices_not_na], extraction[indices_not_na,1], ' ')
  }

  # store numeric representation for type of statistic
  type.factor <- factor(stri_trans_tolower(extraction[,2]), TYPESTRINGS)
  BIG[,TYPE] <- unclass(type.factor)

  # store first argument enclosed in braces
  BIG[,DF1]  <- as.numeric(extraction[,3])

  # store second argument enclosed in braces
  BIG[,DF2]  <- as.numeric(extraction[,4])

  # store value of statistic
  BIG[,STAT] <- as.numeric(extraction[,5])

  # store sign of value of statistic
  BIG[,SIGN] <- sign(BIG[,STAT])

  # round statistic if necessary
  if (round_up) {
    decPlaces  <- decimals(extraction[,5])
    BIG[,STAT] <- BIG[,STAT] + BIG[,SIGN] * (4.999 / 10^(decPlaces+1))
  }
  
  # remove sign from value of statistic
  BIG[,STAT] <- abs(BIG[,STAT])

  
  
  # vectors of indices for each type
  is_t         <- BIG[,TYPE] == TYPE_T
  is_chi2      <- BIG[,TYPE] == TYPE_CHI2
  is_f         <- BIG[,TYPE] == TYPE_F
  is_r         <- BIG[,TYPE] == TYPE_R
  is_z         <- BIG[,TYPE] == TYPE_Z
  has_df1      <- !is.na(BIG[,DF1])
  has_df2      <- !is.na(BIG[,DF2])
  indices_t    <- which(is_t)
  indices_chi2 <- which(is_chi2)
  indices_f    <- which(is_f)
  indices_r    <- which(is_r)
  indices_z    <- which(is_z)

  # find p-value and extract it
  extraction <- stri_match_first_regex(txt.lines.edited, ' *\\bp *(<|<=|=|>) *0*((?:\\d*\\.)?\\d+)[ ,;]*', case_insensitive=TRUE)
  p.reported.str <- rep("", nlines)
  indices_not_na <- which(!is.na(extraction[,1]))
  if( length(indices_not_na) > 0) {
    p.reported.str[indices_not_na] <- paste0("p ", extraction[indices_not_na,2], " ", extraction[indices_not_na,3])
    txt.lines.edited[indices_not_na] <- stri_replace_first_fixed(txt.lines.edited[indices_not_na], extraction[indices_not_na,1], ' ')
  }

  
  
  # store p-value
  BIG[,P.REPORTED] <- as.numeric(extraction[,3])

  # store numeric value for each comparator used with p-value
  BIG[,P.COMP] <- unclass(factor(extraction[,2], c('<','<=', '=','>')))

  # store number of decimals of p-value
  BIG[,P.REPORTED.DECS] <- decimals(extraction[,3])

  # get indices depending on p-value specification
  indices_p     <- which(!is.na(BIG[,P.REPORTED]))
  indices_p_lt  <- which(BIG[,P.COMP] == 1)
  indices_p_leq <- which(BIG[,P.COMP] == 2)
  indices_p_eq  <- which(BIG[,P.COMP] == 3)
  indices_p_gt  <- which(BIG[,P.COMP] == 4)

  # find critical value and extract it
  extraction <- stri_match_first_regex(txt.lines.edited, ' *\\bcrit *= *((?:\\d*\\.)?\\d+)[ ,;]*', case_insensitive=TRUE)
  indices_not_na <- which(!is.na(extraction[,1]))
  if( length(indices_not_na) > 0) {
    txt.lines.edited[indices_not_na] <- stri_replace_first_fixed(txt.lines.edited[indices_not_na], extraction[indices_not_na,1], ' ')
  }

  # store critical value
  BIG[,CRIT.VALUE] <- as.numeric(extraction[,2])

  # find "one-tailed" and extract it
  extraction <- stri_match_first_regex(txt.lines.edited, ' *\\b(one-tailed|1-tailed|one|1t)\\b[ ,;]*', case_insensitive=TRUE)
  # store if one-tailed was specified
  BIG[,ONE.TAILED] <- !is.na(extraction[,1])
  indices_one_tailed <- which(BIG[,ONE.TAILED] == 1)
  if( length(indices_one_tailed) > 0) {
    txt.lines.edited[indices_one_tailed] <- stri_replace_first_fixed(txt.lines.edited[indices_one_tailed], extraction[indices_one_tailed,1], ' ')
  }

  # trim rest which couldn't be parsed
  txt.lines.edited <- stri_trim_both(txt.lines.edited)

  # find lines which have an unparseable part
  indices_wrong_syntax <- which(!stri_isempty(txt.lines.edited) )
  if(length(indices_wrong_syntax) >0 )
    errors[indices_wrong_syntax] <- paste0(errors[indices_wrong_syntax], "\nSyntax error. Are there any illegal expressions? Are there conflicting definitions?")

  # set default crit.value when no crit value has been specified
  indices_crit_na <- which(is.na(BIG[, CRIT.VALUE]))
  BIG[indices_crit_na, CRIT.VALUE] <- ifelse(BIG[indices_crit_na, ONE.TAILED], .10, .05)

  indices_df1_missing <- which((is_t | is_chi2) & !has_df1)
  if(length(indices_df1_missing) > 0)
    errors[indices_df1_missing] <- paste0(errors[indices_df1_missing], "\nStatistic needs specification of df!")

  indices_df2_missing <- which(is_f & !has_df2)
  if(length(indices_df2_missing) > 0)
    errors[indices_df2_missing] <- paste0(errors[indices_df2_missing], "\nStatistic needs specification of second df!")

  indices_excessive_df2 <- which((is_t | is_r | is_z) & has_df2)
  if(length(indices_excessive_df2) > 0)
    errors[indices_excessive_df2] <- paste0(errors[indices_excessive_df2], "\nStatistic has two dfs but only one df allowed!")

  indices_df_zero <- which(BIG[,DF1] == 0 || BIG[,DF2] == 0)
  if(length(indices_df_zero) > 0)
    errors[indices_df_zero] <- paste0(errors[indices_df_zero], "\nDfs of statistic must be greater than zero!")

  indices_df1_real <- which((is_chi2 | is_z | is_r) & (round(BIG[,DF1]) != BIG[,DF1] ))
  if(length(indices_df1_real) > 0)
    errors[indices_df1_real] <- paste0(errors[indices_df1_real], "\nFirst df of statistic must be an integer value!")

  indices_df2_real <- which(is_chi2 & (round(BIG[,DF2]) != BIG[,DF2] ))
  if(length(indices_df2_real) > 0)
    errors[indices_df2_real] <- paste0(errors[indices_df2_real], "\nSecond df of statistic must be an integer value!")

  indices_stat_neg <- which((is_f | is_chi2) & BIG[,SIGN] == -1)
  if(length(indices_stat_neg) > 0)
    errors[indices_stat_neg] <- paste0(errors[indices_stat_neg], "\nStatistic must be greater or equal 0!")

  indices_stat_out_of_bounds <- which(is_r & BIG[,STAT] > 1)
  if(length(indices_stat_out_of_bounds) > 0)
    errors[indices_stat_out_of_bounds] <- paste0(errors[indices_stat_out_of_bounds], "\nStatistic must be >= -1 and <= +1!")

  indices_p_out_of_bounds <- which(BIG[,P.REPORTED] > 1)
  if(length(indices_p_out_of_bounds) > 0)
    errors[indices_p_out_of_bounds] <- paste0(errors[indices_p_out_of_bounds], "\np-value must be less or equal 1!")

  indices_crit_out_of_bounds <- which(BIG[,CRIT.VALUE] > 1)
  if(length(indices_crit_out_of_bounds) > 0)
    errors[indices_crit_out_of_bounds] <- paste0(errors[indices_crit_out_of_bounds], "\nCritical value must be less or equal 1!")

  # compute t-statistic
  if(length(indices_t) > 0)
  {
    BIG[indices_t, P.VALUE] <- pt(BIG[indices_t, STAT], BIG[indices_t, DF1], lower.tail=FALSE) * 2
    BIG[indices_t, D] <- (2*BIG[indices_t, STAT] / sqrt(BIG[indices_t, DF1])) * BIG[indices_t, SIGN]
    BIG[indices_t, G] <- BIG[indices_t, D] * ( 1- (3/(4 * BIG[indices_t, DF1] - 1)))
    BIG[indices_t, N.APPROX] <- BIG[indices_t, DF1] + 2
  }

  # compute pearson's r
  if(length(indices_r) > 0)
  {
    BIG[indices_r, P.VALUE] <- sqrt(BIG[indices_r, DF1]) * BIG[indices_r, STAT] / sqrt(1 - BIG[indices_r, STAT]^2)
    BIG[indices_r, P.VALUE] <- pt(BIG[indices_r, P.VALUE], BIG[indices_r, DF1], lower.tail=FALSE) * 2
    BIG[indices_r, D] <- BIG[indices_r, SIGN] * (2 * BIG[indices_r, STAT]) / sqrt(1 - BIG[indices_r, STAT]^2)
    BIG[indices_r, G] <- BIG[indices_r, D] * (1 - (3 / (4 * BIG[indices_r, DF1] - 1)))
    BIG[indices_r, N.APPROX] <- BIG[indices_r, DF1] + 2
  }

  # compute f-statistic
  if(length(indices_f) > 0)
  {
    indices_f_df1_is_1 <- which(BIG[,TYPE] == 3 & BIG[,DF1] == 1)
    if(length(indices_f_df1_is_1) > 0)
    {
      BIG[indices_f_df1_is_1, P.VALUE] <- sqrt(BIG[indices_f_df1_is_1, STAT])
      BIG[indices_f_df1_is_1, D] <- 2 * BIG[indices_f_df1_is_1, P.VALUE] / sqrt(BIG[indices_f_df1_is_1, DF2])
      BIG[indices_f_df1_is_1, G] <- BIG[indices_f_df1_is_1, D] * (1 - (3/(4 * BIG[indices_f_df1_is_1, DF2] - 1)))
      BIG[indices_f_df1_is_1, N.APPROX] <- BIG[indices_f_df1_is_1, DF2] + 2
    }
    BIG[indices_f, P.VALUE] <- pf(BIG[indices_f, STAT], BIG[indices_f, DF1], BIG[indices_f, DF2], lower.tail=FALSE)
  }

  # compute z-value
  if(length(indices_z) > 0 ) {
    BIG[indices_z, P.VALUE] <- pnorm(BIG[indices_z, STAT], lower.tail=FALSE) * 2

    indices_z_df_exists <- which(BIG[, TYPE] == TYPE_Z & !is.na(BIG[, DF1]))
    if(length(indices_z_df_exists) > 0 ){
      BIG[indices_z_df_exists, N.APPROX] <- BIG[indices_z_df_exists,DF1]

      # If a number is provided for z it's the sample size
      BIG[indices_z_df_exists, D] <- (BIG[indices_z_df_exists, STAT] / sqrt(BIG[indices_z_df_exists, N.APPROX])) * BIG[indices_z_df_exists, SIGN]
      BIG[indices_z_df_exists, G] <- BIG[indices_z_df_exists, D] * (1 - (3 / (4 * BIG[indices_z_df_exists, N.APPROX] - 1)))
    }
  }

  # compute chi2-statistic
  if(length(indices_chi2) > 0){
    # If two numbers are provided for chi2, the first are the dfs, the second is the sample size
    BIG[indices_chi2, P.VALUE] <- pchisq(BIG[indices_chi2, STAT], BIG[indices_chi2, DF1], lower.tail=FALSE)

    indices_chi2_with_n <- which(BIG[,TYPE] == TYPE_CHI2 & BIG[,DF1] == 1 & !is.na(BIG[,DF2]))

    if(length(indices_chi2_with_n) > 0) {
      BIG[indices_chi2_with_n, N.APPROX] <- BIG[indices_chi2_with_n, DF2]
      BIG[indices_chi2_with_n, D] <- sqrt(BIG[indices_chi2_with_n, STAT] / BIG[indices_chi2_with_n, N.APPROX])
      BIG[indices_chi2_with_n, D] <- 2 * BIG[indices_chi2_with_n, D] * sqrt((BIG[indices_chi2_with_n, N.APPROX] - 1)/(BIG[indices_chi2_with_n, N.APPROX] * (1 - BIG[indices_chi2_with_n, D]^2))) * abs(BIG[indices_chi2_with_n, D])/BIG[indices_chi2_with_n, D]
      BIG[indices_chi2_with_n, G] <- BIG[indices_chi2_with_n, D] * (1 - (3/(4 * (BIG[indices_chi2_with_n, N.APPROX]-2) - 1)))
    }
  }


  # store significance (computed p-value must be less than (un)specified critical value)
  BIG[, SIGNIFICANT] <- BIG[, P.VALUE] < BIG[, CRIT.VALUE]

  # store p-value one-tailed by dividing p-value in half
  BIG[, P.VALUE.ONE] <- BIG[, P.VALUE] / 2

  # "actual p-value" must be in accordance to one-tailed specification
  BIG[, P.ACTUAL] <- BIG[, P.VALUE]
  BIG[indices_one_tailed, P.ACTUAL] <- BIG[indices_one_tailed, P.VALUE.ONE]

  # init error
  #BIG[, REP.ERROR]     <- rep(0, nlines)
  #BIG[, REP.DIRECTION] <- rep(0, nlines)

  # check all "p < ?" specifications for reporting errors
  if(length(indices_p_lt) > 0){
    BIG[indices_p_lt, REP.ERROR] <- BIG[indices_p_lt, P.ACTUAL] >= BIG[indices_p_lt, P.REPORTED]
    BIG[indices_p_lt, REP.DIRECTION] <- -BIG[indices_p_lt, REP.ERROR]
  }

  # check all "p <= ?" specifications for reporting errors
  if(length(indices_p_leq) > 0) {
    BIG[indices_p_leq, REP.ERROR] <- BIG[indices_p_leq, P.ACTUAL] > BIG[indices_p_leq, P.REPORTED]
    BIG[indices_p_leq, REP.DIRECTION] <- -BIG[indices_p_leq, REP.ERROR]
  }

  # check all "p > ?" specifications for reporting errors
  if(length(indices_p_gt) > 0){
    BIG[indices_p_gt, REP.ERROR] <- BIG[indices_p_gt, P.ACTUAL] <= BIG[indices_p_gt, P.REPORTED]
    BIG[indices_p_gt, REP.DIRECTION] <- BIG[indices_p_gt, REP.ERROR]
  }

  # check all "p = ?" specifications for reporting errors
  if(length(indices_p_eq) > 0){
    difference <- BIG[indices_p_eq, P.REPORTED] - round(BIG[indices_p_eq, P.ACTUAL], BIG[indices_p_eq, P.REPORTED.DECS])
    BIG[indices_p_eq, REP.DIRECTION] <- sign(difference)
    BIG[indices_p_eq, REP.ERROR] <- difference != 0
  }

  
  # find indices of lines with and without error
  has_no_error <- stri_isempty(errors)
  indices_no_error <- which(has_no_error)
  indices_error    <- which(!has_no_error)

  # produce error message
  warnings <- NULL
  if(length(indices_error) > 0) {
    warnings <- matrix(
      c(
        as.character(indices_not_empty[indices_error]), 
        txt.lines[indices_error], 
        errors[indices_error]
      ), 
      ncol=3
    )
  }

  
  error.direction <- c("smaller", "", "", "bigger")[match(BIG[,REP.DIRECTION], c(-1,0,NA,1))]
  
  # convert data to data.frame as return value
  res <- data.frame(
    line = indices_not_empty,
    paper_id = PAPER_ID,
    study_id = STUDY_ID,
    focal = as.logical(BIG[,IS.FOCAL]),
    type	= type.factor,
    df1 	= BIG[,DF1],
    df2 	= BIG[,DF2],
    d		  = BIG[,D],
    g		  = BIG[,G],
    n.approx = BIG[,N.APPROX],
    statistic = BIG[,STAT],
    p.value	= BIG[,P.VALUE],
    p.value.one	= BIG[,P.VALUE.ONE],
    p.reported = p.reported.str,
    p.crit	= BIG[,CRIT.VALUE],
    significant = as.logical(BIG[, SIGNIFICANT]),
    one.tailed = as.logical(BIG[,ONE.TAILED]),
    reporting.error = as.logical(BIG[,REP.ERROR]),
    error.direction = error.direction,
    parse.error = !has_no_error,
    stringsAsFactors = FALSE
  )

  # add attribute warnings to object "res"
  attr(res, 'warnings') <- warnings

  # return data.frame
  return(res)
}
