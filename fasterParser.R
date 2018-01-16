## ======================================================================
## This is a much faster version of the oldParser.R
## Code by Tobias Kächele & Felix Schönbrodt
## ======================================================================

library(stringi)
library(dplyr)

# returns the reported precision (i.e., number of decimals) from a string
# e.g. decimals("0.120") == 3
decimals <- function(str) {
  locations <- stri_locate_first_fixed(str, '.')[,1]
  decs <- nchar(str) - locations
  decs[is.na(locations)] <- 0
  decs
}

is.one <- function(v) {
  !is.na(v) & v==1
}


# parses a multiline string
parse_ES <- function(txt, round_up = FALSE) {

  if(is.null(txt) || nchar(txt) == 0) {
    return(NULL)
  }
  
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

  if(nlines == 0) {
    return(NULL)
  }
  
  # allocate space for error and warning messages
  errors <- rep("", nlines)
  warning <- rep("", nlines)

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
  P.REPORTED.ERROR  <- 18
  P.REPORTED.ERROR.DIRECTION   <- 19
  IS.FOCAL        <- 20
  PARSE.ERROR     <- 21
  D.REPORTED        <- 22
  D.REPORTED.SIGN   <- 23
  D.REPORTED.DECS   <- 24
  D.REPORTED.LOWER  <- 25
  D.REPORTED.UPPER  <- 26
  D.REPORTED.ERROR.DIRECTION  <- 27
  D.REPORTED.ERROR  <- 28
  REPORTING.ERROR <- 29
  P.VALUE.LOG       <- 30
  D.VAR       <- 31
  STUDYDESIGN <- 32

  # output matrix for all data in numeric form
  BIG <- matrix(NA, nrow = nlines, ncol = 32)

  # find study labels and extract them (everything before ":")
  extraction <- stri_match_first_regex(txt.lines, '^ *(.*?) *(?:(?<=\\)) *(.+?) *)?: *')
  txt.lines.edited <- txt.lines
  indices_not_na <- which(!is.na(extraction[,1]))
  if(length(indices_not_na)) {
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
  if(length(indices_paper_id))
    PAPER_ID[indices_paper_id] <- extraction[indices_paper_id, 2]

  if(length(indices_study_id))
    STUDY_ID[indices_study_id] <- extraction[indices_study_id, 3]

  # is study focal? Underscore at start of id indicates a non-focal test statistic
  BIG[,IS.FOCAL] <- ( stri_sub(PAPER_ID, 1, 1) != "_" )

  # definition of statistic types and typestrings array
  TYPE_T <- 1
  TYPE_CHI2 <- 2
  TYPE_F <- 3
  TYPE_R <- 4
  TYPE_Z <- 5
  TYPE_P <- 6
  TYPESTRINGS <- c('t','chi2','f','r','z', 'p')

  # find statistic and extract it
  extraction <- stri_match_first_regex(txt.lines.edited, ' *\\b(t|chi2|f|r|z|p)(?: *\\( *((?:\\d*\\.)?\\d+)(?: *, *((?:\\d*\\.)?\\d+))? *\\))? *= *(-?(?:\\d*\\.)?\\d+)[ ,;]*', case_insensitive=TRUE)
  indices_not_na <- which(!is.na(extraction[,1]))
  if( length(indices_not_na)){
    txt.lines.edited[indices_not_na] <- stri_replace_first_fixed(txt.lines.edited[indices_not_na], extraction[indices_not_na,1], ' ')
  }

  # mark lines without statistic as error
  indices_no_statistic <- which(is.na(extraction[,1]))
  if (length(indices_no_statistic) > 0)
    errors[indices_no_statistic] <- paste0(errors[indices_no_statistic], "\nNo statistic given!")
  
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
  if (round_up == TRUE) {
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
  is_p         <- BIG[,TYPE] == TYPE_P
  has_df1      <- !is.na(BIG[,DF1])
  has_df2      <- !is.na(BIG[,DF2])
  indices_t    <- which(is_t)
  indices_chi2 <- which(is_chi2)
  indices_f    <- which(is_f)
  indices_r    <- which(is_r)
  indices_z    <- which(is_z)
  indices_pdirect <- which(is_p)

  # find p-value and extract it
  # (if the p-value is given as test statistic)
  extraction <- stri_match_first_regex(txt.lines.edited, ' *\\bp *(<|<=|=|>) *0*((?:\\d*\\.)?\\d+)[ ,;]*', case_insensitive=TRUE)
  p.reported.str <- rep("", nlines)
  indices_not_na <- which(!is.na(extraction[,1]))
  if( length(indices_not_na)) {
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
  if( length(indices_not_na)) {
    txt.lines.edited[indices_not_na] <- stri_replace_first_fixed(txt.lines.edited[indices_not_na], extraction[indices_not_na,1], ' ')
  }

  # store critical value
  BIG[,CRIT.VALUE] <- as.numeric(extraction[,2])

  # find "one-tailed" and extract it
  extraction <- stri_match_first_regex(txt.lines.edited, ' *\\b(one-tailed|1-tailed|one|1t)\\b[ ,;]*', case_insensitive=TRUE)
  # store if one-tailed was specified
  BIG[,ONE.TAILED] <- !is.na(extraction[,1])
  indices_one_tailed <- which(BIG[,ONE.TAILED] == 1)
  if( length(indices_one_tailed)) {
    txt.lines.edited[indices_one_tailed] <- stri_replace_first_fixed(txt.lines.edited[indices_one_tailed], extraction[indices_one_tailed,1], ' ')
  }
  
  
  # ---------------------------------------------------------------------
  #  find study design for t-values: between subjects (BS) or within-subjects (WS)?
  extraction <- stri_match_first_regex(txt.lines.edited, ' *\\b(BS|WS)\\b[ ,;]*', case_insensitive=TRUE)
  
  BIG[,STUDYDESIGN] <- as.numeric(factor(extraction[,2], levels=c("BS", "WS")))
  indices_T_BS <- which(BIG[,STUDYDESIGN] == 1)
  indices_T_WS <- which(BIG[,STUDYDESIGN] == 2)
  if( length(indices_T_BS)) {
    txt.lines.edited[indices_T_BS] <- stri_replace_first_fixed(txt.lines.edited[indices_T_BS], extraction[indices_T_BS,1], ' ')
  }
  if( length(indices_T_WS)) {
    txt.lines.edited[indices_T_WS] <- stri_replace_first_fixed(txt.lines.edited[indices_T_WS], extraction[indices_T_WS,1], ' ')
  }
  

  # find Cohen's d with and without CI/SE and extract it
  extraction <- stri_match_first_regex(txt.lines.edited, ' *\\bd *= *(-)?0*((?:\\d*\\.)?\\d+) *(?:\\[ *(-?(?:\\d*\\.)?\\d+) *[,;] *(-?(?:\\d*\\.)?\\d+) *\\]|[,; ]*se *= *((?:\\d*\\.)?\\d+))?[ ,;]*', case_insensitive=TRUE)
  indices_cohens_d <- which(!is.na(extraction[,1]))
  if( length(indices_cohens_d)) {
    txt.lines.edited[indices_cohens_d] <- stri_replace_first_fixed(txt.lines.edited[indices_cohens_d], extraction[indices_cohens_d,1], ' ')
  }

  # store Cohen's d and number of decimals
  BIG[,D.REPORTED] <- as.numeric(extraction[,3])
  BIG[,D.REPORTED.SIGN] <- ifelse(is.na(extraction[,2]), 1, -1)
  BIG[,D.REPORTED] <- BIG[,D.REPORTED] * BIG[,D.REPORTED.SIGN]

  BIG[,D.REPORTED.DECS] <- decimals(extraction[,3])

  # has Cohen's d?
  indices_cohens_d <- which(!is.na(BIG[,D.REPORTED]))

  # save how it was reported
  d.reported.str <- rep("", nlines)
  if(length(indices_cohens_d)) {
    d.reported.str[indices_cohens_d] <- paste0("d = ", extraction[indices_cohens_d, 3])
  }

  # CI specified for Cohen's d?
  indices_cohens_d_has_ci <- which(!is.na(extraction[,4]))

  # store Cohen's d upper and lower limits if specified
  if(length(indices_cohens_d_has_ci)) {
      BIG[indices_cohens_d_has_ci, D.REPORTED.LOWER] <- as.numeric(extraction[indices_cohens_d_has_ci, 4])
      BIG[indices_cohens_d_has_ci, D.REPORTED.UPPER] <- as.numeric(extraction[indices_cohens_d_has_ci, 5])

      # check ci
      indices_cohens_d_ci_not_in_order <- which(BIG[, D.REPORTED.LOWER] > BIG[, D.REPORTED.UPPER])
      if(length(indices_cohens_d_ci_not_in_order))
        errors[indices_cohens_d_ci_not_in_order] <- paste0(errors[indices_cohens_d_ci_not_in_order], "\nLower bounds of CI for Cohen's d are greater than upper bounds!")

      indices_cohens_d_ci_not_in_order <- which(BIG[, D.REPORTED.LOWER] > BIG[,D.REPORTED])
      if(length(indices_cohens_d_ci_not_in_order))
        errors[indices_cohens_d_ci_not_in_order] <- paste0(errors[indices_cohens_d_ci_not_in_order], "\nLower bounds of CI for Cohen's d are greater than Cohen's d!")

      indices_cohens_d_ci_not_in_order <- which(BIG[, D.REPORTED.UPPER] < BIG[,D.REPORTED])
      if(length(indices_cohens_d_ci_not_in_order))
        errors[indices_cohens_d_ci_not_in_order] <- paste0(errors[indices_cohens_d_ci_not_in_order], "\nUpper bounds of CI for Cohen's d are less than Cohen's d!")
  }


  # SE specified for Cohen's d?
  indices_cohens_d_has_se <- which(!is.na(extraction[,6]))

  # compute Cohen's d upper and lower limits using SE if specified and store them
  if(length(indices_cohens_d_has_se)) {
      SE <- as.numeric(extraction[indices_cohens_d_has_se, 6]) * qnorm(0.975)
      BIG[indices_cohens_d_has_se, D.REPORTED.LOWER] <- BIG[indices_cohens_d_has_se, D.REPORTED] - SE
      BIG[indices_cohens_d_has_se, D.REPORTED.UPPER] <- BIG[indices_cohens_d_has_se, D.REPORTED] + SE
  }

  # now switch and change sign of upper and lower limits if Cohens'd is negative
  indices_cohens_d_neg <- which(BIG[,D.REPORTED] < 0)
  if(length(indices_cohens_d_neg)) {
    BIG[indices_cohens_d_neg, c(D.REPORTED.UPPER, D.REPORTED.LOWER)] <- -BIG[indices_cohens_d_neg, c(D.REPORTED.LOWER, D.REPORTED.UPPER)]
    BIG[indices_cohens_d_neg, D.REPORTED] <- -BIG[indices_cohens_d_neg, D.REPORTED]
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


  # ---------------------------------------------------------------------
  #  error detection / define error messages
  
  indices_df1_missing <- which((is_t | is_chi2 | is_r) & !has_df1)
  if(length(indices_df1_missing))
    errors[indices_df1_missing] <- paste0(errors[indices_df1_missing], "\nStatistic needs specification of df!")

  indices_df2_missing <- which(is_f & !has_df2)
  if(length(indices_df2_missing))
    errors[indices_df2_missing] <- paste0(errors[indices_df2_missing], "\nStatistic needs specification of second df!")

  indices_excessive_df2 <- which((is_t | is_r | is_z) & has_df2)
  if(length(indices_excessive_df2))
    errors[indices_excessive_df2] <- paste0(errors[indices_excessive_df2], "\nStatistic has two dfs but only one df allowed!")

  indices_df_zero <- which(BIG[,DF1] == 0 || BIG[,DF2] == 0)
  if(length(indices_df_zero))
    errors[indices_df_zero] <- paste0(errors[indices_df_zero], "\nDfs of statistic must be greater than zero!")

  indices_df1_real <- which((is_chi2 | is_z | is_r) & (round(BIG[,DF1]) != BIG[,DF1] ))
  if(length(indices_df1_real))
    errors[indices_df1_real] <- paste0(errors[indices_df1_real], "\nFirst df of statistic must be an integer value!")

  indices_df2_real <- which(is_chi2 & (round(BIG[,DF2]) != BIG[,DF2] ))
  if(length(indices_df2_real))
    errors[indices_df2_real] <- paste0(errors[indices_df2_real], "\nSecond df of statistic must be an integer value!")

  indices_stat_neg <- which((is_f | is_chi2) & BIG[,SIGN] == -1)
  if(length(indices_stat_neg))
    errors[indices_stat_neg] <- paste0(errors[indices_stat_neg], "\nStatistic must be greater or equal 0!")
  
  indices_p_outofrange <- which(is_p & (BIG[ ,STAT] > 1 | BIG[,SIGN] == -1))
  if(length(indices_p_outofrange))
    errors[indices_p_outofrange] <- paste0(errors[indices_p_outofrange], "\np-values must be between 0 and 1!")
  
  indices_stat_out_of_bounds <- which(is_r & BIG[,STAT] > 1)
  if(length(indices_stat_out_of_bounds))
    errors[indices_stat_out_of_bounds] <- paste0(errors[indices_stat_out_of_bounds], "\nStatistic must be >= -1 and <= +1!")

  indices_p_out_of_bounds <- which(BIG[,P.REPORTED] > 1)
  if(length(indices_p_out_of_bounds))
    errors[indices_p_out_of_bounds] <- paste0(errors[indices_p_out_of_bounds], "\np-value must be less or equal 1!")

  indices_crit_out_of_bounds <- which(BIG[,CRIT.VALUE] > 1)
  if(length(indices_crit_out_of_bounds))
    errors[indices_crit_out_of_bounds] <- paste0(errors[indices_crit_out_of_bounds], "\nCritical value must be less or equal 1!")
  
  indices_tF_without_design <- which((is_t | is_f) & is.na(BIG[,STUDYDESIGN]))
  if(length(indices_tF_without_design))
    warning[indices_tF_without_design] <- paste0(warning[indices_tF_without_design], "\nWarning: t- or F-value is reported without design. <i>Effect size computations assume a <u>between-subject</u> design!</i> Add \"; BS\" for between subject design; add \"; WS\" for within subject design)")


  # compute t-statistic
  if(length(indices_t))
  {
    BIG[indices_t, P.VALUE] <- pt(BIG[indices_t, STAT], BIG[indices_t, DF1], lower.tail=FALSE) * 2
	BIG[indices_t, P.VALUE.LOG] <- pt(BIG[indices_t, STAT], BIG[indices_t, DF1], lower.tail=FALSE, log.p = TRUE) + log(2)
	
    BIG[indices_t, N.APPROX] <- BIG[indices_t, DF1] + 2
	
	# Which studies are between-SS, which within-SS?
	BS_design <- which(BIG[,TYPE] == TYPE_T & (BIG[,STUDYDESIGN] == 1 | is.na(BIG[,STUDYDESIGN])))
	WS_design <- which(BIG[,TYPE] == TYPE_T & BIG[,STUDYDESIGN] == 2)
	
	# compute between effect size
	BIG[BS_design, D] <- (2*BIG[BS_design, STAT] / sqrt(BIG[BS_design, N.APPROX])) * BIG[BS_design, SIGN]
	
	# compute within effect size
	BIG[WS_design, D] <- (BIG[WS_design, STAT] / sqrt(BIG[WS_design, N.APPROX])) * BIG[WS_design, SIGN]
	
	BIG[indices_t, G] <- BIG[indices_t, D] * ( 1- (3/(4 * BIG[indices_t, N.APPROX] - 1)))
    
	#BIG[indices_t, D.VAR] <- 4/BIG[indices_t, N.APPROX] + BIG[indices_t, D]^2 / (2*BIG[indices_t, N.APPROX])
	BIG[indices_t, D.VAR] <- (4 + BIG[indices_t, D]^2) / BIG[indices_t, N.APPROX]
  }

  # compute pearson's r
  if(length(indices_r))
  {
	# t value from correlation
	t.r <- BIG[indices_r, STAT]*sqrt(BIG[indices_r, DF1] / (1 - BIG[indices_r, STAT]^2))
    BIG[indices_r, P.VALUE] <- pt(t.r, BIG[indices_r, DF1], lower.tail=FALSE) * 2
	BIG[indices_r, P.VALUE.LOG] <- pt(t.r, BIG[indices_r, DF1], lower.tail=FALSE, log.p = TRUE) + log(2)
    BIG[indices_r, D] <- BIG[indices_r, SIGN] * (2 * BIG[indices_r, STAT]) / sqrt(1 - BIG[indices_r, STAT]^2)
    BIG[indices_r, G] <- BIG[indices_r, D] * (1 - (3 / (4 * BIG[indices_r, DF1] - 1)))
    BIG[indices_r, N.APPROX] <- BIG[indices_r, DF1] + 2

	var.r <- (1 - BIG[indices_r, STAT]^2)^2/(BIG[indices_r, N.APPROX] - 1)
	BIG[indices_r, D.VAR] <-  4 * var.r/(1 - BIG[indices_r, STAT]^2)^3
  }

  # compute f-statistic
  if(length(indices_f))
  {
    indices_f_df1_is_1 <- which(BIG[,TYPE] == 3 & BIG[,DF1] == 1)
    if(length(indices_f_df1_is_1))
    {
      BIG[indices_f_df1_is_1, N.APPROX] <- BIG[indices_f_df1_is_1, DF2] + 2
	  
  	  # Which studies are between-SS, which within-SS?
  	  BS_design <- which(BIG[,TYPE] == TYPE_F & BIG[,DF1] == 1 & (BIG[,STUDYDESIGN] == 1 | is.na(BIG[,STUDYDESIGN])))
  	  WS_design <- which(BIG[,TYPE] == TYPE_F & BIG[,DF1] == 1 & BIG[,STUDYDESIGN] == 2)
	  
  	  # compute between effect size
  	  BIG[BS_design, D] <- (2*sqrt(BIG[BS_design, STAT]) / sqrt(BIG[BS_design, N.APPROX])) * BIG[BS_design, SIGN]
	  
  	  # compute within effect size
  	  BIG[WS_design, D] <- sqrt(BIG[WS_design, STAT] /BIG[WS_design, N.APPROX]) * BIG[WS_design, SIGN]

      BIG[indices_f_df1_is_1, G] <- BIG[indices_f_df1_is_1, D] * ( 1- (3/(4 * BIG[indices_f_df1_is_1, N.APPROX] - 1)))
	  
	  #BIG[indices_f_df1_is_1, D.VAR] <- 4/BIG[indices_f_df1_is_1, N.APPROX] + BIG[indices_f_df1_is_1, D]^2 / (2*BIG[indices_f_df1_is_1, N.APPROX])
	  BIG[indices_f_df1_is_1, D.VAR] <- (4 + BIG[indices_f_df1_is_1, D]^2) / BIG[indices_f_df1_is_1, N.APPROX]
    }
	
	BIG[indices_f, P.VALUE.LOG] <- pf(BIG[indices_f, STAT], BIG[indices_f, DF1], BIG[indices_f, DF2], lower.tail=FALSE, log.p = TRUE)
	BIG[indices_f, P.VALUE] <- exp(BIG[indices_f, P.VALUE.LOG])
  }

  # compute z-value
  if(length(indices_z) ) {
    BIG[indices_z, P.VALUE] <- pnorm(BIG[indices_z, STAT], lower.tail=FALSE) * 2
	BIG[indices_z, P.VALUE.LOG] <- pnorm(BIG[indices_z, STAT], lower.tail=FALSE, log.p = TRUE) + log(2)

    indices_z_df_exists <- which(BIG[, TYPE] == TYPE_Z & !is.na(BIG[, DF1]))
    if(length(indices_z_df_exists) ){
      # If a number is provided for z it's the sample size
      BIG[indices_z_df_exists, N.APPROX] <- BIG[indices_z_df_exists,DF1]  
      BIG[indices_z_df_exists, D] <- (BIG[indices_z_df_exists, STAT] / sqrt(BIG[indices_z_df_exists, N.APPROX])) * BIG[indices_z_df_exists, SIGN]
      BIG[indices_z_df_exists, G] <- BIG[indices_z_df_exists, D] * (1 - (3 / (4 * BIG[indices_z_df_exists, N.APPROX] - 1)))
    }
  }

  # compute chi2-statistic
  if(length(indices_chi2)){
    # If two numbers are provided for chi2, the first are the dfs, the second is the sample size
    BIG[indices_chi2, P.VALUE] <- pchisq(BIG[indices_chi2, STAT], BIG[indices_chi2, DF1], lower.tail=FALSE)
	BIG[indices_chi2, P.VALUE.LOG] <- pchisq(BIG[indices_chi2, STAT], BIG[indices_chi2, DF1], lower.tail=FALSE, log.p = TRUE)

    indices_chi2_with_n <- which(BIG[,TYPE] == TYPE_CHI2 & BIG[,DF1] == 1 & !is.na(BIG[,DF2]))

    if(length(indices_chi2_with_n)) {
      BIG[indices_chi2_with_n, N.APPROX] <- BIG[indices_chi2_with_n, DF2]
      BIG[indices_chi2_with_n, D] <- sqrt(BIG[indices_chi2_with_n, STAT] / BIG[indices_chi2_with_n, N.APPROX])
      BIG[indices_chi2_with_n, D] <- 2 * BIG[indices_chi2_with_n, D] * sqrt((BIG[indices_chi2_with_n, N.APPROX] - 1)/(BIG[indices_chi2_with_n, N.APPROX] * (1 - BIG[indices_chi2_with_n, D]^2))) * abs(BIG[indices_chi2_with_n, D])/BIG[indices_chi2_with_n, D]
      BIG[indices_chi2_with_n, G] <- BIG[indices_chi2_with_n, D] * (1 - (3/(4 * (BIG[indices_chi2_with_n, N.APPROX]-2) - 1)))
    }
  }
  
  # compute directly entered p-value
  if(length(indices_pdirect) ) {
	
    BIG[indices_pdirect, P.REPORTED] <- BIG[indices_pdirect, STAT]
	
	# assume that the directly reported p-value is the correct p-value;	
	BIG[indices_pdirect, P.VALUE] <- BIG[indices_pdirect, STAT]
	BIG[indices_pdirect, P.VALUE.LOG] <- log(BIG[indices_pdirect, STAT])
	
	# If one-tailed, the actual p-value is double the size
	BIG[indices_pdirect, P.VALUE][BIG[indices_pdirect, ONE.TAILED] == 1] <- BIG[indices_pdirect, P.VALUE][BIG[indices_pdirect, ONE.TAILED] == 1]*2
	BIG[indices_pdirect, P.VALUE.LOG][BIG[indices_pdirect, ONE.TAILED] == 1] <- BIG[indices_pdirect, P.VALUE.LOG][BIG[indices_pdirect, ONE.TAILED] == 1]+ log(2)

	indices_pdirect_df_exists <- which(BIG[, TYPE] == TYPE_P & !is.na(BIG[, DF1]))
    if(length(indices_pdirect_df_exists) ){
		
	# conversion formula for converting p to d, see: https://books.google.de/books?id=GC42CwAAQBAJ&pg=PA100&lpg=PA100&dq=meta-analysis+convert+chi2+to+d+degrees+of+freedom&source=bl&ots=_c4EHEyRis&sig=yIaUDAbQ3RfPvLTfEE7-thpEXys&hl=de&sa=X&ved=0ahUKEwiyp-XntdnLAhWipnIKHab5CaA4ChDoAQhfMAg#v=onepage&q=meta-analysis%20convert%20chi2%20to%20d%20degrees%20of%20freedom&f=false
	
  	  BIG[indices_pdirect_df_exists, N.APPROX] <- BIG[indices_pdirect_df_exists, DF1]+2
	
	  d <- (qnorm(1-(BIG[indices_pdirect_df_exists, STAT]/2))*2)/sqrt(BIG[indices_pdirect_df_exists, N.APPROX])
      BIG[indices_pdirect_df_exists, D] <- d * BIG[indices_pdirect_df_exists, SIGN]
	  
      BIG[indices_pdirect_df_exists, G] <- BIG[indices_pdirect_df_exists, D] * (1 - (3 / (4 * BIG[indices_pdirect_df_exists, N.APPROX] - 1)))
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
  #BIG[, P.REPORTED.ERROR]     <- rep(0, nlines)
  #BIG[, P.REPORTED.ERROR.DIRECTION] <- rep(0, nlines)

  # check all "p < ?" specifications for reporting errors
  if(length(indices_p_lt)){
    BIG[indices_p_lt, P.REPORTED.ERROR] <- BIG[indices_p_lt, P.ACTUAL] >= BIG[indices_p_lt, P.REPORTED]
    BIG[indices_p_lt, P.REPORTED.ERROR.DIRECTION] <- -BIG[indices_p_lt, P.REPORTED.ERROR]
  }

  # check all "p <= ?" specifications for reporting errors
  if(length(indices_p_leq)) {
    BIG[indices_p_leq, P.REPORTED.ERROR] <- BIG[indices_p_leq, P.ACTUAL] > BIG[indices_p_leq, P.REPORTED]
    BIG[indices_p_leq, P.REPORTED.ERROR.DIRECTION] <- -BIG[indices_p_leq, P.REPORTED.ERROR]
  }

  # check all "p > ?" specifications for reporting errors
  if(length(indices_p_gt)){
    BIG[indices_p_gt, P.REPORTED.ERROR] <- BIG[indices_p_gt, P.ACTUAL] <= BIG[indices_p_gt, P.REPORTED]
    BIG[indices_p_gt, P.REPORTED.ERROR.DIRECTION] <- BIG[indices_p_gt, P.REPORTED.ERROR]
  }

  # check all "p = ?" specifications for reporting errors
  if(length(indices_p_eq)){
    difference <- BIG[indices_p_eq, P.REPORTED] - round(BIG[indices_p_eq, P.ACTUAL], BIG[indices_p_eq, P.REPORTED.DECS])
    BIG[indices_p_eq, P.REPORTED.ERROR.DIRECTION] <- sign(difference)
    BIG[indices_p_eq, P.REPORTED.ERROR] <- difference != 0
  }

  # check d
  indices_cohens_d2 <- which(!is.na(BIG[, D]) & !is.na(BIG[, D.REPORTED]))
  if(length(indices_cohens_d2)) {
    difference <- BIG[indices_cohens_d2, D.REPORTED] - round(BIG[indices_cohens_d2, D], BIG[indices_cohens_d2, D.REPORTED.DECS])
    BIG[indices_cohens_d2, D.REPORTED.ERROR.DIRECTION] <- sign(difference)
    BIG[indices_cohens_d2, D.REPORTED.ERROR] <- difference != 0
  }

 
  # ---------------------------------------------------------------------
  # Bring errors into a nice shape
  
  # find indices of lines with and without error
  has_no_error <- stri_isempty(errors)
  indices_no_error <- which(has_no_error)
  indices_error    <- which(!has_no_error)

  # produce error message
  ERRORS <- NULL
  if(length(indices_error)) {
    ERRORS <- matrix(
      c(
        as.character(indices_not_empty[indices_error]),
        txt.lines[indices_error],
        errors[indices_error]
      ),
      ncol=3
    )
  }
  
  # ---------------------------------------------------------------------
  # Bring warnings into a nice shape
  
  # find indices of lines with and without error
  has_no_warning <- stri_isempty(warning)
  indices_no_warning <- which(has_no_warning)
  indices_warning    <- which(!has_no_warning)

  # produce warning message
  WARNINGS <- NULL
  if(length(indices_warning)) {
    WARNINGS <- matrix(
      c(
        as.character(indices_not_empty[indices_warning]),
        txt.lines[indices_warning],
        warning[indices_warning]
      ),
      ncol=3
    )
  }


  p.reported.error.direction <- c("smaller", "", "", "bigger")[match(BIG[,P.REPORTED.ERROR.DIRECTION], c(-1,0,NA,1))]
  d.reported.error.direction <- c("smaller", "", "", "bigger")[match(BIG[,D.REPORTED.ERROR.DIRECTION], c(-1,0,NA,1))]

  global.reporting.error <- is.one(BIG[,P.REPORTED.ERROR]) | is.one(BIG[,D.REPORTED.ERROR])

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
    significant = as.logical(BIG[,SIGNIFICANT]),
    one.tailed = as.logical(BIG[,ONE.TAILED]),
    reporting.error = as.logical(BIG[,P.REPORTED.ERROR]),
    error.direction = p.reported.error.direction,
    parse.error = !has_no_error,
    d.reported = BIG[,D.REPORTED],
    d.reported.str = d.reported.str,
    d.reported.lower = BIG[,D.REPORTED.LOWER],
    d.reported.upper = BIG[,D.REPORTED.UPPER],
    d.reported.error = as.logical(BIG[,D.REPORTED.ERROR]),
    d.reported.error.direction = d.reported.error.direction,
    global.reporting.error = global.reporting.error,
	p.value.log	= BIG[,P.VALUE.LOG],
	d.var = BIG[, D.VAR],
	d.se = sqrt(BIG[, D.VAR]),
	studydesign = BIG[, STUDYDESIGN],
    stringsAsFactors = FALSE
  )


  # add attribute warnings to object "res"
  attr(res, 'ERRORS') <- ERRORS
  attr(res, 'WARNINGS') <- WARNINGS

  # return data.frame
  return(res)
}


parse_ES("t(72) = 4.80; BS
t(72) = 4.80; WS
t(72) = 5.80; WS")
#
# (p <- parse_ES("t(47, 4)=2.1"))
# (p <- parse_ES("t(47, 4)=2.1; BS"))
# (p <- parse_ES("t(47, 4)=2.1; WS"))
#
# (p <- parse_ES("F(1, 47)=4.41"))
# (p <- parse_ES("F(1, 47)=4.41; BS"))
# (p <- parse_ES("F(1, 47)=4.41; WS"))