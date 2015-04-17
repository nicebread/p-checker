library(shiny)
library(shinyTable)
library(stringr)
library(dplyr)
library(ggplot2)

source("helpers.R")
source("parser.R")
source("p-curve.R")

# TODO: p-curve plot: 33%-power-curve


# ---------------------------------------------------------------------
# Statistical Inference

# TIVA = Test of insufficient variance - see https://replicationindex.wordpress.com/2014/12/30/the-test-of-insufficient-variance-tiva-a-new-tool-for-the-detection-of-questionable-research-practices/
# code adapted from Moritz Heene
TIVA <- function (p.values) {
  p.values <- clamp(p.values, MIN=.00001, MAX=1-.00001)	# prevent infinite z-values
  z.values  <- qnorm(p.values, lower.tail = FALSE)
  var.z  <- var(z.values)
  df <- length(p.values)-1
  chi.square  <- df*var.z
  chi.p  <- pchisq(chi.square, df=df, lower.tail = TRUE)
  data.frame(var.z=var.z, chi2=chi.square, df=df, p.value=chi.p)
}



#input <- list(round_up=FALSE, digits=3, group_by_paper=TRUE, only_first_ES=TRUE, txt=x, pcurve_power=33, pcurve_crit=.05, experimental=FALSE); dat <- list()


shinyServer(function(input, output, session) {
	
	dat <- reactiveValues(
		tblDisplay=data.frame(),	# keeps the rounded values for display
		tbl=data.frame(),			# keeps the precise values
		warnings=""					# keeps a vector of parser errors
	)
	
	
	# ---------------------------------------------------------------------
	# Define the input field; populate with GET parameter if provided; or demo data if not.
	# Reading GET parameters
	# Parse the GET query string  
	  output$syntax <- renderUI({
		  query <- parseQueryString(session$clientData$url_search)
		  
		  if (is.null(query["syntax"]) | query["syntax"] == "NULL") {
			  res <- "
# GO AND REPLACE THE EXAMPLES!			  
# Easy mode ('#' starts a comment)
t(47) = 2.1
chi2(1) = 9.1
r(77) = .47
F(1, 88) = 9.21

# add reported p-value; one-tailed; alpha level
t(123) = 2.54; p < .01
Z = 1.9; one-tailed; p=.03
r(25) = 0.21; crit=.10

# add paper ID
A&B (2001) Study1: t(88)=2.1; one-tailed; p < .02
A&B (2001) Study1: r(147)=.246
A&B (2001) Study2: F(1,100)=9.1
CD&E (2014) S1a: F(1,210)=4.45; p < .01
CD&E (2014) S1b: t(123)=2.01; one-tailed; p = .02
"
		  } else {
			  res <- query["syntax"]
		  }

		return(list(
			HTML(paste0('<textarea id="txt" rows="18" cols="45">', res, '</textarea>'))
		))
	 })
	
	
	
	
	
	# every time the text field is changed, this function is called and parses the input string
	observe({
		
		# quit when syntax field is not created yet
		if (is.null(input$txt)) return()
			
		tbl <- parse_ES(input$txt, round_up=input$round_up)
		
		# parser errors present?
		if (length(attr(tbl, "warnings")) > 0) {
			dat$warnings <- attr(tbl, "warnings")
		} else {
			dat$warnings <- ""
		}
		
		# No input? Return empty data frame
		if (nrow(tbl) == 0) {
			dat$tblDisplay <- data.frame()
			dat$tbl <- data.frame()
			return()
		}
		
		# ---------------------------------------------------------------------
		# R-index computations
		
		tbl$Z <- qnorm(1-(tbl$p.value/2))
		tbl$obs.pow <- pnorm(tbl$Z-qnorm(1-tbl$p.crit/2))
		
		# set all values of non-focal tests to NA
		tbl$Z[tbl$focal==FALSE] <- NA
		tbl$obs.pow[tbl$focal==FALSE] <- NA
		
		# compute median observed power within each study - but only for focal hypothesis tests
		
		# only select first ES of each study, if requested
		if (input$only_first_ES == TRUE) {
			# TODO: this is awkward. There must be a better way to solve this.
			tbl2 <- tbl %>% 
				group_by(paper_id, study_id) %>% 
				filter(focal==TRUE, row_number() <= 1) %>% 
				mutate(median.obs.pow=median(obs.pow)) %>% 
				ungroup()
			tbl <- inner_join(tbl, select(tbl2, paper_id, study_id, median.obs.pow), by=c("paper_id", "study_id"))
			
			# remove median.obs.pow for other ES
			tbl <- tbl %>% group_by(paper_id, study_id) %>% mutate(snum=1:n()) %>% ungroup()
			tbl$median.obs.pow[tbl$snum>1] <- NA
			tbl$snum <- NULL
		} else {
			tbl <- tbl %>% group_by(paper_id, study_id) %>% dplyr::mutate(median.obs.pow=median(obs.pow[focal==TRUE])) %>% ungroup()
		}
		
		# ---------------------------------------------------------------------
		# p-curve computations
		
		dat$pcurve_power <- ifelse(input$experimental==TRUE, input$pcurve_power/100, 1/3)
		pps <- get_pp_values(type=tbl$type, statistic=tbl$statistic, df=tbl$df1, df2=tbl$df2, p.crit=input$pcurve_crit, power=dat$pcurve_power)
		
		tbl <- cbind(tbl, pps[, -1])

		# tblDisplay stores the table with nicely formatted numbers
		tblDisplay <- tbl
		
		## Apply the function to each column, and convert the list output back to a data frame
		
		# df columns
		tblDisplay[, 5] <- sapply(tblDisplay[, 5], format_num, digits=max(c(sapply(tbl$df1, decimalplaces))))
		tblDisplay[, 6] <- sapply(tblDisplay[, 6], format_num, digits=max(c(sapply(tbl$df2, decimalplaces))))
		
		# all other columns
		for (i in 7:ncol(tblDisplay)) {
			tblDisplay[, i] <-format_num(tblDisplay[, i], digits=input$digits)
		}
		
		dat$tbl <- tbl
		dat$tblDisplay <- tblDisplay
	})


	# ---------------------------------------------------------------------
	# Output for parser errors
	output$parser_errors <- renderUI({
			HTML(paste('<p><span style="color:red">',
						dat$warnings, collapse="<br>"),
						'</span></p>'
						)
	})

	# show warning if experimental features are activated
	output$experimental_warning <- renderUI({
		if (input$experimental == TRUE) {
			HTML('<div class="alert alert-danger" role="alert">Warning: You activated experimental settings. Think twice before you run an actual p-checker analysis with these untested settings!</div>')
		}
	})

	# ---------------------------------------------------------------------
	# Output for p value reporting tab
	
	output$report_table <- renderUI({
		if (nrow(dat$tblDisplay) > 0) {
			report_table <- dat$tblDisplay[, c("paper_id", "study_id", "type", "df1", "df2", "statistic", "p.value", "p.value.one", "p.reported", "p.crit", "one.tailed", "reporting.error", "error.direction")]
		} else {
			return(NULL)
		}
		
		
		# Summary
		if (input$group_by_paper == FALSE) {
			return(list(
				HTML("<h2>p-reporting analysis: Are there wrong reported p values?</h2>"),				
				HTML(paste0(
					"<h4>",
					"<b>Percentage of p-values that are incorrectly rounded down</b> (", 
						sum(report_table$error.direction == "smaller"), "/", nrow(report_table), ") = ", 
						round(sum(report_table$error.direction == "smaller")*100/nrow(report_table), input$digits), "%<br>",	
					"</h4>"
				)),
				HTML("<br><br><h2>Detailed results for each test statistic:</h2>"),
				div(renderTable({report_table}), style = "font-size:80%")
			))
		}
		
		
		
		if (input$group_by_paper == TRUE) {
			report_table2 <- report_table %>% group_by(paper_id)
			report_table2 <- report_table2 %>% dplyr::summarise(
				n.p_values = n(),
				wrong.p_values = sum(error.direction == "smaller"),
				percentage.wrong.p_values = wrong.p_values / n.p_values,
				all.correct = wrong.p_values == 0
			)
			
			return(list(
				HTML("<h2>p-reporting analysis: Are there wrong reported p values?</h2>"),				
				HTML("<h2>Results for each paper</h2>"),
				renderTable({report_table2}),
				HTML(paste0(
					"<h4>",
					"<b>Percentage of p-values that are incorrectly rounded down</b> (", 
						sum(report_table$error.direction == "smaller"), "/", nrow(report_table), ") = ", 
						round(sum(report_table$error.direction == "smaller")*100/nrow(report_table), input$digits), "%<br>",
					"<b>Percentage of papers with at least one wrong p-value</b> (", 
						sum(report_table2$all.correct == FALSE), "/", nrow(report_table2), ") = ", 
						round(sum(report_table2$all.correct == FALSE)*100/nrow(report_table2), input$digits), "%<br>",		
					"</h4>"
				)),
				HTML("<br><br><h2>Detailed results for each test statistic:</h2>"),
				div(renderTable({report_table}), style = "font-size:80%")
			))
		}

	})


	# ---------------------------------------------------------------------
	# Output for R-index tab


	output$rindex_table <- renderUI({
		if (nrow(dat$tblDisplay) > 0) {
			rindex_table <- dat$tblDisplay[dat$tblDisplay$focal==TRUE, c("paper_id", "study_id", "type", "df1", "df2", "statistic", "p.value", "p.crit", "Z", "obs.pow", "significant", "median.obs.pow")]
			
			# Omit near-significants if requested
			if (input$omit_nearly_significant == TRUE) {
				rindex_table <- rindex_table %>% filter(p.value < .05 | p.value > .10)
			}		
			
			if (input$only_first_ES == TRUE) {
				rindex_table <- rindex_table %>% group_by(paper_id, study_id) %>% filter(row_number() <= 1)
			}
			
			return(list(
				HTML("<br><br><h2>Detailed results for each test statistic:</h2>"),
				div(renderTable({rindex_table}), style = "font-size:80%")
			))	
		}
	})

	
	output$rindex_summary <- renderUI({

		if (nrow(dat$tbl) == 0) {return(NULL)}
			
		# only select focal hypothesis tests
		tbl <- dat$tbl[dat$tbl$focal==TRUE, ]

		# Omit near-significants if requested
		if (input$omit_nearly_significant == TRUE) {
			tbl <- tbl %>% filter(p.value < .05 | p.value > .10)
		}		
		
		# only select first ES of each study, if requested
		if (input$only_first_ES == TRUE) {
			tbl <- tbl %>% group_by(paper_id, study_id) %>% filter(row_number() <= 1)
		}
		
		
		
		# One r-index analysis across all ES
		if (input$group_by_paper == FALSE) {
			success_rate <- sum(tbl$p.value < tbl$p.crit, na.rm=TRUE)/nrow(tbl)
			obs_power0 <- tbl %>% group_by(paper_id, study_id) %>% filter(row_number() <= 1) %>% select(median.obs.pow)
			obs_power <- mean(obs_power0$median.obs.pow, na.rm=TRUE)
			inflation_rate <- success_rate - obs_power
			r_index <- obs_power - inflation_rate
		
			tiva <- TIVA(tbl$p.value)
			
			result <- paste0(
				"<h2>R-Index analysis:</h2><h4>",
				"<b>Success rate</b> = ", 	round(success_rate, 4), "<br>",
				"<b>Mean observed power</b> = ", round(obs_power, 4), "<br>",
				"<b>Inflation rate</b> = ", round(inflation_rate, 4), "<br>",
				"<b>R-Index</b> = ", round(r_index, 4),
				"</h4>"
			)
		
			return(HTML(result))
		}
		
		# Separate r-index analyses for each paper
		if (input$group_by_paper == TRUE) {
			tbl <- tbl %>% group_by(paper_id)
			success_rate <- tbl %>% summarise(
				k_effect_sizes = n(),
				success_rate = sum(p.value < p.crit, na.rm=TRUE)/k_effect_sizes
			) %>% select(paper_id, k_effect_sizes, success_rate) 
			
			tbl <- tbl %>% group_by(paper_id, study_id)
			obs_power0 <- tbl %>% filter(row_number() <= 1) %>% select(median.obs.pow)
			obs_power <- obs_power0 %>% group_by(paper_id) %>% select(median.obs.pow) %>% summarise_each(funs(mean))
			
			rindex <- inner_join(success_rate, obs_power, by="paper_id")
			rindex <- rindex %>% mutate(
				inflation_rate 	= success_rate - median.obs.pow,
				r_index 		= median.obs.pow - inflation_rate
			)
			
			dat$rindex <- rindex
			
			return(list(
				HTML("<h2>R-Index analysis:</h2>"),				
				HTML(paste0(
					"<h4>",
					"<b>Average success rate</b> = ", 	round(mean(rindex$success_rate, na.rm=TRUE), input$digits), "<br>",
					"<b>Average mean observed power</b> = ", round(mean(rindex$median.obs.pow, na.rm=TRUE), input$digits), "<br>",
					"<b>Average inflation rate</b> = ", round(mean(rindex$inflation_rate, na.rm=TRUE), input$digits), "<br>",
					"<b>Average R-Index</b> = ", round(mean(rindex$r_index, na.rm=TRUE), input$digits),
					"</h4>"
				)),
				renderTable({rindex})
			))
		}
	})
	
	
	
	# ---------------------------------------------------------------------
	# Output for TIVA tab


	output$tiva_table <- renderUI({
		if (nrow(dat$tblDisplay) > 0) {
			tiva_table <- dat$tblDisplay %>% 
				filter(focal==TRUE)  %>% 
				select(paper_id, study_id, type, df1, df2, statistic, p.value, p.crit, Z)
			
			if (input$only_first_ES == TRUE) {
				tiva_table <- tiva_table %>% group_by(paper_id, study_id) %>% filter(row_number() <= 1)
			}
			
			return(list(
				HTML("<br><br><h2>Detailed results for each test statistic:</h2>"),
				div(renderTable({tiva_table}), style = "font-size:80%")
			))	
		}
	})

	
	output$tiva_summary <- renderUI({

		if (nrow(dat$tbl) == 0) {return(NULL)}

		# only select focal hypothesis tests
		tbl <- dat$tbl %>% filter(focal==TRUE)
		
		if (input$only_first_ES == TRUE) {
			tbl <- tbl %>% group_by(paper_id, study_id) %>% filter(row_number() <= 1)
		}

		# One TIVA analysis across all ES
		if (input$group_by_paper == FALSE) {
			tiva <- TIVA(tbl$p.value)

			result <- paste0(
				"<h2>Test of insufficient variance (TIVA)</h2>",
				"<small>Variances &lt; 1 suggest bias. The chi2 tests the H0 that variance = 1.</small>",

				"<h4>",
				"Variance = ", round(tiva$var.z, 4), "<br>",
				"Chi2(", tiva$df, ") = ", round(tiva$chi2, 3), "; ", p(tiva$p.value),
				"</h4>"
			)

			return(HTML(result))
		}

		# Separate TIVA analyses for each paper
		if (input$group_by_paper == TRUE) {
			tiva <- data.frame()
			for (i in unique(tbl$paper_id)) {
				tiva <- rbind(tiva, data.frame(paper_id = i, TIVA(tbl$p.value[tbl$paper_id == i])))
			}
			
			# remove rows where only 1 test stat was provided
			tiva <- tiva %>% filter(!is.na(var.z))
			
			dat$tiva <- tiva
			
			return(list(
				HTML(paste0(
					"<h2>Test of insufficient variance (TIVA)</h2>",
					"<small>Variances &lt; 1 suggest bias. The chi2 tests the H0 that variance = 1.</small>",
					"<small>Note: TIVA selects only the <b>first</b> p value of each study!</small>"
				)),				
				HTML(paste0(
					"<h2>Summary on ", nrow(tiva), " TIVA analyses:</h2><h4>",
					"<b>Average variance</b> = ", round(mean(tiva$var.z, na.rm=TRUE), input$digits), "<br>",
					"<b>% of papers with variance &lt; 1</b>: ", round(sum(tiva$var.z<1)/nrow(tiva)*100, input$digits), "%<br>",
					"<b>% of papers with variance significantly &lt; 1</b>: ", round(sum(tiva$p.value<.05)/nrow(tiva)*100, input$digits), "%",
					"</h4>"
				)),
				renderTable({tiva})
			))
		}
	})
	
	
	
	
	
	
	
	# ---------------------------------------------------------------------
	# 	Output for p-curve tab
	
	
	output$pcurve_table <- renderUI({
		if (nrow(dat$tbl) > 0) {
			pcurve_table <- dat$tblDisplay[dat$tblDisplay$focal==TRUE, c("paper_id", "study_id", "type", "df1", "df2", "statistic", "p.value", "significant", "ppr", "ppl", "pp33")]
			
			if (input$only_first_ES == TRUE) {
				pcurve_table <- pcurve_table %>% group_by(paper_id, study_id) %>% filter(row_number() <= 1)
			}
			
			return(list(
				HTML("<br><br><h2>Detailed results for each test statistic:</h2>"),
				div(renderTable({pcurve_table}), style = "font-size:80%")
			))	
		} else {
			return(NULL)
		}
	})
	
	output$pcurve_plot <- renderPlot({
		if (input$group_by_paper == TRUE | nrow(dat$tbl) == 0) {return(NULL)}
		
		plot(NA, xlim=c(0, input$pcurve_crit), ylim=c(0, 100), xlab="p-value", ylab="Percentage of p values")
		abline(h=1/input$pcurve_crit, col="red", lty="dashed")
		legend("topright", lty=c("solid", "dotted", "dashed"), col=c("blue", "darkgreen", "red"), legend=c("Observed p-curve", "Null of 33% power", "Null of nil effect"))
		
		# select only focal and significant hypothesis tests
		tbl <- dat$tbl[dat$tbl$focal==TRUE & dat$tbl$significant==TRUE, ]
		
		if (input$only_first_ES == TRUE) {
			tbl <- tbl %>% group_by(paper_id, study_id) %>% filter(row_number() <= 1)
		}
		
		bins <- table(cut(tbl$p.value, breaks=seq(0, input$pcurve_crit, by=.01)))
		perc <- (bins/sum(bins))*100
		lines(x=seq(0, input$pcurve_crit-.01, by=.01)+.005, y=perc, col="blue")
		text(x=seq(0, input$pcurve_crit-.01, by=.01)+.007, y=perc + 5, col="blue", label=paste0(round(perc), "%"))
		
	})
	
	
	output$pcurve_summary <- renderUI({

		if (nrow(dat$tbl) == 0) {return(NULL)}
			
		# select only focal and significant hypothesis tests
		tbl <- dat$tbl %>% filter(focal==TRUE, significant==TRUE)
		
		if (input$only_first_ES == TRUE) {
			tbl <- tbl %>% group_by(paper_id, study_id) %>% filter(row_number() <= 1) %>% ungroup()
		}
		
		
		if (input$group_by_paper == FALSE) {
			
			if (nrow(tbl) == 0) {return(NULL)}
			
			if (input$pcurve_version == "v3") {
				pcurve_tests <- p_curve_3(tbl[, c("ppr", "ppl", "pp33")])
				teststring <- "Z = "
				statistics <- round(c(pcurve_tests$Z_evidence, pcurve_tests$Z_lack, pcurve_tests$Z_hack), input$digits)
				ps <- p(c(pcurve_tests$p_evidence, pcurve_tests$p_lack, pcurve_tests$p_hack), input$digits)
			}
			if (input$pcurve_version == "v2") {
				pcurve_tests <- p_curve_2(tbl[, c("ppr", "ppl", "pp33")])
				teststring <- paste0("chi2(", pcurve_tests$df, ") = ")
				statistics <- round(c(pcurve_tests$chi2_evidence, pcurve_tests$chi2_lack, pcurve_tests$chi2_hack), input$digits)
				ps <- p(c(pcurve_tests$p_evidence, pcurve_tests$p_lack, pcurve_tests$p_hack), input$digits)
			}

								
			result <- paste0(
				"<h2>Statistical Inference on p-curve:</h2><h4>",
				"<b>Studies contain evidential value</b>: <br>",
				teststring, statistics[1], "; ", ps[1], "<br>",
				"<small>A significant p value indicates that the p-curve is right-skewed, which indicates evidential value.</small><br><br>",
			
				"<b>Studies’ evidential value, if any, is inadequate</b><br>",
				teststring, statistics[2], "; ", ps[2], "<br>",
				"<small>A significant p value indicates that the p-curve is flatter than one would expect if studies were powered at ", round(dat$pcurve_power*100), "%, which indicates that the results have no evidential value.</small><br><br>",
			
				"<b>Studies lack evidential value and were intensely <i>p</i>-hacked </b>: <br>",
				teststring, statistics[3], "; ", ps[3], "<br>",
				"<small>A significant p value indicates that the p-curve is left-skewed, which indicates p-hacking/ selective reporting.</small><br><br>",
				"</h4>"
			)
		
			return(HTML(result))
		}
		
		
		
		if (input$group_by_paper == TRUE) {
			
			if (input$pcurve_version == "v3") {
				pcurve <- tbl %>% group_by(paper_id) %>% select(ppr, ppl, pp33) %>% do(data.frame(p_curve_3(.)))
			}
			if (input$pcurve_version == "v2") {
				pcurve <- tbl %>% group_by(paper_id) %>% select(ppr, ppl, pp33) %>% do(data.frame(p_curve_2(.)))
			}			
			
			# remove rows where no test stats are provided
			pcurve <- pcurve %>% filter(!is.na(Z_evidence))
			
			dat$pcurve <- pcurve

			return(list(
				HTML(paste0(
					"<h2>Statistical Inference on p-curve:</h2>"
				)),				
				HTML(paste0(
					"<h2>Summary on ", nrow(pcurve), " p-curves:</h2><h4>",
					"<b>% of papers with evidential value</b>: ", round(sum(pcurve$p_evidence<.05)/nrow(pcurve)*100, 1), "%<br>",
					"<b>% of papers with lack of evidence</b>: ", round(sum(pcurve$p_lack<.05)/nrow(pcurve)*100, 1), "%<br>",
					"<b>% of papers intensely p-hacked</b>: ", round(sum(pcurve$p_hack<.05)/nrow(pcurve)*100, 1), "%<br>",
					"<b>% of papers with inconclusive p-curve</b>: ", round(sum(pcurve$inconclusive)/nrow(pcurve)*100, 1), "%",
					"</h4>"
				)),

				HTML(paste0(
					"<h2>Summary on ", sum(pcurve$inconclusive==FALSE), " p-curves which are not inconclusive:</h2><h4>",
					"<b>% of papers with evidential value</b>: ", round(sum(pcurve$p_evidence<.05)/sum(pcurve$inconclusive==FALSE)*100, 1), "%<br>",
					"<b>% of papers with lack of evidence</b>: ", round(sum(pcurve$p_lack<.05)/sum(pcurve$inconclusive==FALSE)*100, 1), "%<br>",
					"<b>% of papers intensely p-hacked</b>: ", round(sum(pcurve$p_hack<.05)/sum(pcurve$inconclusive==FALSE)*100, 1), "%<br>",
					"</h4>"
				)),
				
				renderTable({pcurve})
			))
		}
	})
	
	
	
	
	
	
	# ---------------------------------------------------------------------
	# Export for p-curve; save analysis as link
	
	output$export <- renderUI({
		if (nrow(dat$tbl) > 0) {
			
			res <- c()
			for (i in 1:nrow(dat$tbl)) {
				switch(as.character(dat$tbl$type[i]),
					"t" = {res <- c(res, paste0("t(", dat$tbl$df1[i], ")=", dat$tbl$statistic[i]))},
					"f" = {res <- c(res, paste0("F(", dat$tbl$df1[i], ", ", dat$tbl$df2[i], ")=", dat$tbl$statistic[i]))},
					"chi2" = {res <- c(res, paste0("chi2(", dat$tbl$df1[i], ")=", dat$tbl$statistic[i]))},
					"r" = {res <- c(res, paste0("r(", dat$tbl$df1[i], ")=", f2(dat$tbl$statistic[i], decimalplaces(dat$tbl$statistic[i]), skipZero=TRUE)))},
					"z" = {res <- c(res, paste0("Z=", dat$tbl$statistic[i]))}
				)
			}			
			
			res1 <- paste(res, collapse="\n")
			res2 <- paste(res, collapse="<br>")
			
			pcurve_link <- paste0("http://www.p-curve.com/app3/?tests=", URLencode(res1, reserved=TRUE))
			
			return(list(
				HTML(paste0("<h4>Copy and share <a href='http://shinyapps.org/apps/p-checker/?syntax=", URLencode(input$txt, reserved=TRUE),
				"'>this link</a> to send the p-checker analysis to others.</h4>")),
				
				HTML(paste0("<h4>Click <a href='", pcurve_link, "' target='_blank'>this link</a> to transfer the test statistics to p-curve.com.</h4>",
				"<small>Note: This transfers the test statistics without paper identifier. That means, p-curve.com will compute an omnibus test with all values.</small>",
				"<br><br>Alternatively: Copy and paste the syntax below to <a href='http://www.p-curve.com/app3'>p-curve.com</a><br>
				<small>Note: If you want identical results to p-curve.com, turn off the 'Gracious rounding up' option at the left panel.</small><br><br>				
				"
				
				)),
				HTML(paste0("<small>",res2, "</small>"))
			))	
		} else {
			return(NULL)
		}
	})
	
	
	# ---------------------------------------------------------------------
	# Meta-Analysis
	
	# adapted from http://stackoverflow.com/questions/7549694/ggplot2-adding-regression-line-equation-and-r2-on-graph
	# lm_eqn = function(m) {
# 	  l <- list(a = format(coef(m)[1], digits = 2),
# 	      b = format(abs(coef(m)[2]), digits = 2),
# 	      r2 = format(summary(m)$r.squared, digits = 3));
#
# 	  if (coef(m)[2] >= 0)  {
# 	    #eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,l)
# 		eq <- substitute(italic(y) == a + b %.% italic(x), l)
# 	  } else {
# 	    #eq <- substitute(italic(y) == a - b %.% italic(x)*","~~italic(r)^2~"="~r2,l)
# 		eq <- substitute(italic(y) == a - b %.% italic(x),l)
# 	  }
# 	  as.character(as.expression(eq));
# 	}
	
	
	
	
		#
	# output$meta <- renderUI({
	#
	# 	# select only focal tests
	# 	tbl <- dat$tbl[dat$tbl$focal==TRUE, ]
	#
	# 	# select all tests of one kind:
	# 	test_type <- "error"
	# 	switch(input$meta_ES_type,
	# 		"ttest_2" = {
	# 			test_type <- "Two-group t-test"
	# 			tbl <- tbl[tbl$type=="t" | (tbl$type=="f" & tbl$df1==1), ]
	# 			tbl <- mutate(tbl, vi = 1/(n.approx/2) + 1/(n.approx/2) + g^2/(2*(n.approx)))
	# 		}
	# 	)
	#
	# 	if (nrow(tbl) == 0) {return(NULL)}
	# 	if (var(tbl$n.approx, na.rm=TRUE) == 0) {return(list(renderText({"No variance in sample sizes."})))}
	#
	# 	meta_table <- tbl[, c("paper_id", "study_id", "type", "df1", "df2", "statistic", "n.approx", "p.value", "d", "g", "vi")]
	#
	# 	library(metafor)
	# 	meta_analysis <- rma(tbl$g, tbl$vi, method="FE")
	# 	#
	# 	# R <- cor.test(tbl$n.approx, tbl$g, use="p")
	# 	# r_eqn <- as.character(as.expression(substitute(italic(r) == R*", "~~p, list(R=round(R$estimate, 3), p=p(R$p.value)))))
	#
	# 	return(list(
	# 		HTML(paste0('<div class="alert alert-warning" role="alert">This meta-analysis assumes that all test statistics are of the type: "', test_type, '"! Furthermore, it assume that participants are equally distributed across groups. This initial summary <i>does not replace a proper meta-analysis!</i></div>')),
	# 		renderPrint({print(meta_analysis)}),
	# 		renderPlot({funnel(meta_analysis)}),
	# 		# renderPlot({
	# 		# 	ggplot(tbl, aes(x=n.approx, y=g)) + geom_point() + geom_smooth(method=lm) + xlab("Sample size") + ylab("Effect size (Hedge's g)") + annotate("text", x = Inf, y = Inf, label = r_eqn, colour="black", size = 4, parse=TRUE, hjust=1.5, vjust=1.5)
	# 		# }, res=120),
	# 		HTML("<br><br><h2>Detailed results for each test statistic:</h2>"),
	# 		div(renderTable({meta_table}), style = "font-size:80%")
	# 	))
	# })
	#
	
	

	
	
	# ---------------------------------------------------------------------
	# The button for downloading the result data frame
	output$downloadData <- downloadHandler(
		filename = c('input_data.csv'),
		content = function(file) {write.csv(dat$tbl, file)}
	)

	output$downloadRIndex <- downloadHandler(
		filename = c('rindex_results.csv'),
		content = function(file) {write.csv(dat$rindex, file)}
	)
	
	output$downloadTIVA <- downloadHandler(
		filename = c('tiva_results.csv'),
		content = function(file) {write.csv(dat$tiva, file)}
	)
	
	output$downloadPCurve <- downloadHandler(
		filename = c('p_curve_results.csv'),
		content = function(file) {write.csv(dat$pcurve, file)}
	)
	
	# ---------------------------------------------------------------------
	# Load demo data
	observe({
		con <- NULL
		switch(input$demodata,
			"JPSP1" = {
				demo <- readLines(con <- file("demo-data/JPSP-p-curve.txt"))
				demo2 <- paste(demo, collapse="\n")
				updateTextInput(session, inputId = "txt", value = demo2)
				},
			"855" = {
				demo <- readLines(con <- file("demo-data/855_t_tests.txt"))
				demo2 <- paste(demo, collapse="\n")
				updateTextInput(session, inputId = "txt", value = demo2)
				},
			"H0_100x5" = {
				demo <- readLines(con <- file("demo-data/H0_100x5.txt"))
				demo2 <- paste(demo, collapse="\n")
				updateTextInput(session, inputId = "txt", value = demo2)
			},
			"H1_100x5" = {
				demo <- readLines(con <- file("demo-data/H1_100x5.txt"))
				demo2 <- paste(demo, collapse="\n")
				updateTextInput(session, inputId = "txt", value = demo2)
			},
			"H0_hack_100x5" = {
				demo <- readLines(con <- file("demo-data/H0_hack_100x5.txt"))
				demo2 <- paste(demo, collapse="\n")
				updateTextInput(session, inputId = "txt", value = demo2)
			},
			"elderly" = {
				demo <- readLines(con <- file("demo-data/elderly_priming.txt"))
				demo2 <- paste(demo, collapse="\n")
				updateTextInput(session, inputId = "txt", value = demo2)
			}
		)
		if (!is.null(con)) close(con)
	})
		

})
