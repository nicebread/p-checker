library(shiny)
library(shinyjs)
library(stringr)
library(plyr); library(dplyr)
library(ggplot2)
library(ggvis)
library(meta)
library(metafor)
library(broom)

# source inference functions
source("pancollapse.R")
source("fasterParser.R")
source("p-curve.R")
source("helpers.R")
source("TIVA.R")


#input <- list(round_up=FALSE, digits=3, group_by_paper=TRUE, only_first_ES=TRUE, txt=x, pcurve_power=33, pcurve_crit=.05, experimental=FALSE); dat <- list()


shinyServer(function(input, output, session) {

	# dat is a reactive object that keeps the computed variables
	dat <- reactiveValues(
		tblDisplay=data.frame(),	# keeps the rounded values for display
		tbl=data.frame(),			# keeps the precise values
		ERRORS = "",					# keeps a vector of parser errors
		WARNINGS = ""					# keeps a vector of parser warnings
	)
	
	exportTbl <- function() {
		if (nrow(dat$tbl) > 0) {
			res <- c()
			for (i in 1:nrow(dat$tbl)) {
				switch(as.character(dat$tbl$type[i]),
					"t" = {res <- c(res, paste0("t(", dat$tbl$df1[i], ")=", dat$tbl$statistic[i]))},
					"f" = {res <- c(res, paste0("F(", dat$tbl$df1[i], ", ", dat$tbl$df2[i], ")=", dat$tbl$statistic[i]))},
					"chi2" = {res <- c(res, paste0("chi2(", dat$tbl$df1[i], ")=", dat$tbl$statistic[i]))},
					"r" = {res <- c(res, paste0("r(", dat$tbl$df1[i], ")=", f2(dat$tbl$statistic[i], decimalplaces(dat$tbl$statistic[i]), skipZero=TRUE)))},
					"z" = {res <- c(res, paste0("Z=", dat$tbl$statistic[i]))},
					"p" = {res <- c(res, paste0("p", ifelse(!is.na(dat$tbl$df1[i]), paste0("(", dat$tbl$df1[i], ")"), ""), "=", dat$tbl$statistic[i]))}
				)
			}
		
			return(res)
		} else {
			return(NULL)
		}
	}
	
	# ---------------------------------------------------------------------
	# Define the input field; populate with GET parameter if provided; or demo data if not.
	# Reading GET parameters
	# Parse the GET query string  
	  output$syntax <- renderUI({
		  query <- parseQueryString(session$clientData$url_search)
		  
		  if (!is.null(query["preset"]) & query["preset"] != "NULL") {
			  # set to default demo
			  res <- query["syntax"]
			  
			  # try to change the demo dropdown box. This triggers the observe event that changes the text of the input field
			  updateSelectInput(session, "demodata", selected = query["preset"])
			  
		  } else if (!is.null(query["syntax"]) & query["syntax"] != "NULL") {
			  res <- query["syntax"]
		  } else {
			  res <- paste(readLines("snippets/demo_syntax.txt", encoding="UTF-8"), collapse="\n")
		  }

		return(list(
			HTML(paste0('<textarea class="form-control" style="font-family:Lucida Console, Monaco, monospace !important;resize:none;white-space:pre;word-wrap:normal;overflow-x:scroll;" id="txt" rows="18">', res, '</textarea>'))
		))
	 })
	
	
	
	
	
	# every time the text field is changed, this function is called and parses the input string
	observe({
		
		# quit when syntax field is not created yet or empty
		if (is.null(input$txt)) {
			return()
		}
		
		tbl <- parse_ES(input$txt, round_up = input$round_up)
		
		# parser errors present?
		if (length(attr(tbl, "ERRORS")) > 0) {
			dat$ERRORS <- attr(tbl, "ERRORS")
		} else {
			dat$ERRORS <- ""
		}
		
		# parser warnings present?
		if (length(attr(tbl, "WARNINGS")) > 0) {
			dat$WARNINGS <- attr(tbl, "WARNINGS")
		} else {
			dat$WARNINGS <- ""
		}
		
		# No input? Return empty data frame
		if (is.null(tbl) || nrow(tbl) == 0) {
			dat$tblDisplay <- data.frame()
			dat$tbl <- data.frame()
			return()
		}
		
		# ---------------------------------------------------------------------
		# R-index computations
		
		#  -log(2) divides the p-value by two on the log scale
		tbl$Z <- qnorm(tbl$p.value.log - log(2), lower.tail = FALSE, log.p=TRUE)
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
		
		dat$pcurve_power <- input$pcurve_power/100
		pps <- get_pp_values(type=tbl$type, statistic=tbl$statistic, df=tbl$df1, df2=tbl$df2, p.crit=input$pcurve_crit, power=dat$pcurve_power)$res
		
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
	  if(dat$ERRORS != "") {
	    alert.create(
	      paste0(
	        "<strong>Line ",
	        dat$ERRORS[,1],
	        "</strong> <code>",
	        dat$ERRORS[,2],
	        "</code>",
	        stri_replace_all_fixed(dat$ERRORS[,3], "\n", "<br>"), 
	        collapse="<br>"
	      ),
	      style="danger"
	    )
	  }
	})
	
	# ---------------------------------------------------------------------
	# Output for parser warnings (only for the meta-analysis tab)
	output$MA_warnings <- renderUI({
	  if(dat$WARNINGS != "") {
		pancollapse.create(
		  "There are test statistics where the design is unclear (between or within subjects).<br>Click here for details.",
		  alert.create(
		      paste0(
		        "<strong>Line ",
		        dat$WARNINGS[,1],
		        "</strong> <code>",
		        dat$WARNINGS[,2],
		        "</code>",
		        stri_replace_all_fixed(dat$WARNINGS[,3], "\n", "<br>"), 
		        collapse="<br>"
		      ),
	      style="danger"),
		class = "panel-danger")
	  }
	})

	# ---------------------------------------------------------------------
	#  show warning if experimental features are activated
	
	output$experimental_features_warning <- renderUI({
		if (input$experimental == TRUE) {
			HTML('<div class="alert alert-danger" role="alert">Warning: You activated experimental settings. Think twice before you run an actual p-checker analysis with these untested settings!</div>')
		}
	})
	
	# ---------------------------------------------------------------------
	#  show roxygen-style title if provided in the syntax
	output$roxygen_title <- renderUI({
		if (is.null(input$txt)) return()
		
		# parse the input
		title <- str_match_empty(input$txt, "#' @title (.*)\\n", position=2)
		subtitle <- str_match_empty(input$txt, "#' @subtitle (.*)\\n", position=2)
		url <- str_match_empty(input$txt, "#' @url (.*)\\n", position=2)
		details <- str_match_empty(input$txt, "#' @details (.*)\\n", position=2)
				
		if (any(c(title, subtitle, url, details) != "")) {
			HTML('<div style="border-bottom: thick solid #3F658F; margin-bottom: 20px; padding:5px; background-color:#EAF5FF">
					<span style="font-size:150%; font-weight: bold; color:black;">', title, '</span>
					<span style="font-size:120%;"> ', subtitle, '</span><br />
					<span style="font-size:90%; padding-top:5px;"> ', details, '</span><br />
					<span style="font-size:90%;"><a href="', url, '">', url, '</a></span>
					</div>'
			)
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
				HTML("<h3>p-reporting analysis: Are there wrongly reported p values?</h3>"),				
				HTML(paste0(
					"<h4>",
					"<b>Percentage of p-values that are incorrectly rounded down</b> (", 
						sum(report_table$error.direction == "smaller"), "/", nrow(report_table), ") = ", 
						round(sum(report_table$error.direction == "smaller")*100/nrow(report_table), input$digits), "%<br>",	
					"</h4>"
				)),
				panel.create(
				  "Detailed results for each test statistic:",
				  getTable(report_table, classWhenValueFun('reporting.error', T, 'danger') )
				)
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
				HTML("<h3>p-reporting analysis: Are there wrongly reported p values?</h3>"),	
				
				#HTML("<h3>Results for each paper</h3>"),
				
				pancollapse.create(
				  "Results for each paper",
				  getTable(report_table2)
				),
				#renderTable({report_table2}),
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
				pancollapse.create(
				  "Detailed results for each test statistic",
				  getTable(report_table)
				)
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
				rindex_table <- rindex_table %>% 
					filter(p.value < input$omit_nearly_significant_range[1] | p.value > input$omit_nearly_significant_range[2])
			}		
			
			if (input$only_first_ES == TRUE) {
				rindex_table <- rindex_table %>% group_by(paper_id, study_id) %>% filter(row_number() <= 1)
			}
			
			return(list(
				pancollapse.create(
				  "Detailed results for each test statistic",
				  getTable(rindex_table, function(x){ if(!is.na(x[["significant"]]) && !x[["significant"]]){"danger"} })
				)
			))	
		}
	})

	
	output$rindex_summary <- renderUI({

		if (nrow(dat$tbl) == 0) {return(NULL)}
			
		# only select focal hypothesis tests
		tbl <- dat$tbl[dat$tbl$focal==TRUE, ]

		# Omit near-significants if requested
		if (input$omit_nearly_significant == TRUE) {
			tbl <- tbl %>% 
				filter(p.value < input$omit_nearly_significant_range[1] | p.value > input$omit_nearly_significant_range[2])
		}		
		
		# only select first ES of each study, if requested
		if (input$only_first_ES == TRUE) {
			tbl <- tbl %>% group_by(paper_id, study_id) %>% filter(row_number() <= 1)
		}
		
		
		
		# One r-index analysis across all ES
		if (input$group_by_paper == FALSE) {
			success_rate <- sum(tbl$p.value < tbl$p.crit, na.rm=TRUE)/nrow(tbl)
			obs_power0 <- tbl %>% group_by(paper_id, study_id) %>% filter(row_number() <= 1) %>% select(paper_id, study_id, median.obs.pow)
			obs_power <- median(obs_power0$median.obs.pow, na.rm=TRUE)
			inflation_rate <- success_rate - obs_power
			r_index <- obs_power - inflation_rate
		
			# run TIVA with one-tailed p-values (see email from Aurelien)
			tiva <- TIVA(tbl$p.value.log - log(2), log.p=TRUE)
			
			result <- paste0(
				"<h3>R-Index analysis:</h3><h4>",
				"<b>Success rate</b> = ", 	round(success_rate, 4), "<br>",
				"<b>Median observed power</b> = ", round(obs_power, 4), "<br>",
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
			obs_power0 <- tbl %>% filter(row_number() <= 1) %>% select(paper_id, study_id, median.obs.pow)
			obs_power <- obs_power0 %>% group_by(paper_id) %>% select(paper_id, median.obs.pow) %>% summarise_each(funs(median))
			
			rindex <- inner_join(success_rate, obs_power, by="paper_id")
			rindex <- rindex %>% mutate(
				inflation_rate 	= success_rate - median.obs.pow,
				r_index 		= median.obs.pow - inflation_rate
			)
			
			dat$rindex <- rindex
			rindex[, -1] <- round(rindex[, -1], input$digits)
			
			return(list(
				HTML("<h3>R-Index analysis:</h3>"),				
				HTML(paste0(
					"<h4>",
					"<b>Average success rate</b> = ", 	round(mean(rindex$success_rate, na.rm=TRUE), input$digits), "<br>",
					"<b>Average median observed power</b> = ", round(mean(rindex$median.obs.pow, na.rm=TRUE), input$digits), "<br>",
					"<b>Average inflation rate</b> = ", round(mean(rindex$inflation_rate, na.rm=TRUE), input$digits), "<br>",
					"<b>Average R-Index</b> = ", round(mean(rindex$r_index, na.rm=TRUE), input$digits),
					"</h4>"
				)),
				#renderTable({rindex})
				
				pancollapse.create(
				  "R-index results for each paper",
				  getTable(rindex)
				)
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
				pancollapse.create(
				  "Detailed results for each test statistic",
				  getTable(tiva_table)
				)
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
			# - log(2): use one-tailed p-values for TIVA
			tiva <- TIVA(tbl$p.value.log - log(2), log.p=TRUE)

			result <- paste0(
				"<h3>Test of insufficient variance (TIVA):</h3>",
				"<small>Variances &lt; 1 suggest bias. The chi2 tests the H0 that variance >= 1; a significant result indicates that the empirical variance is significantly smaller than 1.</small>",

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
				tiva <- rbind(tiva, data.frame(paper_id = i, TIVA(tbl$p.value.log[tbl$paper_id == i], log.p=TRUE)))
			}
			
			# remove rows where only 1 test stat was provided
			tiva <- tiva %>% filter(!is.na(var.z))
			
			dat$tiva <- tiva
			
			return(list(
				HTML(paste0(
					"<h3>Test of insufficient variance (TIVA):</h3>",
					"<small>Variances &lt; 1 suggest bias. The chi2 tests the H0 that variance = 1.</small>",
					"<small>Note: TIVA selects only the <b>first</b> p value of each study!</small>"
				)),				
				HTML(paste0(
					"<h3>Summary on ", nrow(tiva), " TIVA analyses:</h3><h4>",
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
			pcurve_table <- dat$tblDisplay[dat$tbl$focal==TRUE & dat$tbl$p.value <= .05, c("paper_id", "study_id", "type", "df1", "df2", "statistic", "p.value", "significant", "ppr", "ppl", "pp33")]
			
			if (input$only_first_ES == TRUE) {
				pcurve_table <- pcurve_table %>% group_by(paper_id, study_id) %>% filter(row_number() <= 1)
			}
			
			return(list(
				pancollapse.create(
				  "Detailed results for each test statistic",
				  getTable(pcurve_table)
				)
			))	
		} else {
			return(NULL)
		}
	})
	
	output$pcurve_plot <- renderUI({
		if (input$group_by_paper == TRUE | nrow(dat$tbl) == 0) {return(NULL)}
				
		send2pcurve.button.tag <- actionButton("send2pcurve", 'Do the same analysis at pcurve.com', icon=icon("arrow-circle-right"), class="btn-sm")
		
		
		return(list(
			renderPlot({
				plot(NA, xlim=c(0, input$pcurve_crit), ylim=c(0, 100), xlab="p-value", ylab="Percentage of p values")
				abline(h=1/input$pcurve_crit, col="red", lty="dashed", lwd=2)
				legend("topright", lty=c("solid", "dotted", "dashed"), col=c(COLORS$BLUE, "darkgreen", "red"), legend=c("Observed p-curve", paste0(input$pcurve_power, "% power curve"), "Nil effect"), bty="n")
		
				# select only focal and significant hypothesis tests
				tbl <- dat$tbl[dat$tbl$focal==TRUE & dat$tbl$p.value <= .05, ]
		
				if (input$only_first_ES == TRUE) {
					tbl <- tbl %>% group_by(paper_id, study_id) %>% filter(row_number() <= 1)
				}
		
				bins <- table(cut(tbl$p.value, breaks=seq(0, input$pcurve_crit, by=.01)))
				perc <- (bins/sum(bins))*100
				
				# empirical p-curve
				lines(x=seq(0, input$pcurve_crit-.01, by=.01)+.005, y=perc, col=COLORS$BLUE, lwd=2)
				
				# 33% (or any other) power curve
				lines(x=seq(0, input$pcurve_crit-.01, by=.01)+.005, y=theoretical_power_curve(input$pcurve_power/100, p.max=input$pcurve_crit)*100, col=COLORS$GREEN, lty="dashed", lwd=2)
				
				text(x=seq(0, input$pcurve_crit-.01, by=.01)+.006, y=perc + 8, col="black", label=paste0(round(perc), "%"), cex=)
			}, res=100, width=600),
			send2pcurve.button.tag,
			HTML(paste0("<br><small>Note: This transfers the test statistics without paper identifier. That means, p-curve.com will compute an omnibus test with all values.</small><br>
		"))
		))
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
				"<h3>Statistical Inference on p-curve:</h3><h4>",
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
			pcurve[, -c(1, 8)] <- round(pcurve[, -c(1, 8)], input$digits)

			return(list(
				HTML(paste0(
					"<h3>Statistical Inference on p-curve:</h3>"
				)),				
				HTML(paste0(
					"<h3>Summary on ", nrow(pcurve), " p-curves:</h3><h4>",
					"<b>% of papers with evidential value</b>: ", round(sum(pcurve$p_evidence<.05)/nrow(pcurve)*100, 1), "%<br>",
					"<b>% of papers with lack of evidence</b>: ", round(sum(pcurve$p_lack<.05)/nrow(pcurve)*100, 1), "%<br>",
					"<b>% of papers intensely p-hacked</b>: ", round(sum(pcurve$p_hack<.05)/nrow(pcurve)*100, 1), "%<br>",
					"<b>% of papers with inconclusive p-curve</b>: ", round(sum(pcurve$inconclusive)/nrow(pcurve)*100, 1), "%",
					"</h4>"
				)),

				HTML(paste0(
					"<h3>Summary on ", sum(pcurve$inconclusive==FALSE), " p-curves which are not inconclusive:</h3><h4>",
					"<b>% of papers with evidential value</b>: ", round(sum(pcurve$p_evidence<.05)/sum(pcurve$inconclusive==FALSE)*100, 1), "%<br>",
					"<b>% of papers with lack of evidence</b>: ", round(sum(pcurve$p_lack<.05)/sum(pcurve$inconclusive==FALSE)*100, 1), "%<br>",
					"<b>% of papers intensely p-hacked</b>: ", round(sum(pcurve$p_hack<.05)/sum(pcurve$inconclusive==FALSE)*100, 1), "%<br>",
					"</h4>"
				)),
				
				
				pancollapse.create(
				  "p-curve results for each paper",
				  getTable(pcurve)
				)
			))
		}
	})
	

	
	shinyjs::onclick("send2pcurve", {
		
		res1 <- paste(exportTbl(), collapse="\n")
		
		tbl <- dat$tbl %>% filter(focal==TRUE, significant==TRUE, type=="p")
		if (nrow(tbl) > 0) {
			info(paste0("p-curve.com does not accept directly entered p-values (only test statistics). Results are exported without: ", paste0("p=", tbl$statistic, collapse="; ")))
		}
		
		pcurve_link <- paste0("http://www.p-curve.com/app4/?tests=", URLencode(res1, reserved=TRUE))
		js$browseURL(pcurve_link)
    })
	
	
	
	
	
	# ---------------------------------------------------------------------
	# Effect size panel	
	
	output$effectsizes <- renderUI({
		
		TBL <- dat$tbl %>% filter(!is.na(d), !is.na(d.var))
				
		if (nrow(TBL) > 2) {
			
		  isolate({
			TBL$g.abs <- abs(TBL$g)
			TBL$label <- paste0("Row ", 1:nrow(TBL), ": ", TBL$paper_id, " ", TBL$study_id)
			TBL$id <- 1:nrow(TBL)
						
			
			mysessions <- function(x) {
			  if(is.null(x)) return(NULL)
			  #notice below the id column is how ggvis can understand which session to show 
			  row <- df[df$id == x$id, ]
			  #prettyNum shows the number with thousand-comma separator  
			  paste0(prettyNum(row$sessions, big.mark=",",scientific=F)) 
			}
			
			#ES_plot <- ggplot(TBL, aes(x=n.approx, y=abs(g.abs))) + geom_point() + xlab("Approximate n (log scale)") + ylab("Absolute Hedge's g") + geom_smooth(method=lm) + scale_x_log10(breaks=round(seq(min(TBL$n.approx, na.rm=TRUE), max(TBL$n.approx, na.rm=TRUE), length.out=5)))
			
		
			ES_table <- dat$tblDisplay %>% filter(g != "NA", n.approx != "NA") %>% select(paper_id, study_id, type, df1, df2, statistic, p.value, n.approx, d, g, d.var, d.se, studydesign)			
		  })
		  
		  
		  k <- nrow(TBL)
		  
		  # Begg & Mazumdar Rank correlation test for publication bias; only if k>2
		  if (k > 2) {
			  suppressWarnings({
				  #Begg <- cor.test(TBL$n.approx, TBL$g, use="p", method="kendall")
				  Begg <- NULL
			  })
		  } else {
		  	Begg <- NULL
		  }
		  
		  # construct funnel plot
		  meta1 <- metagen(TBL$d, TBL$d.se)
		  meta2 <- rma(TBL$d, sei=TBL$d.se)

		  # ---------------------------------------------------------------------
		  #  Compute Egger's test / PET-PEESE
		  # either as lm() or rma(); see http://www.metafor-project.org/doku.php/tips:rma_vs_lm_lme_lmer
		  
		  PET <- PET.lm <- lm(d~d.se, data=TBL, weight=1/TBL$d.var)
		  PEESE <- PEESE.lm <- lm(d~d.var, data=TBL, weight=1/TBL$d.var)

		  if (input$MR_model == "rma") {
			  PET <- rma(yi = TBL$d, vi = TBL$d.var, mods=TBL$d.se, method="REML")
			  PEESE <- rma(yi = TBL$d, vi = TBL$d.var, mods=TBL$d.var, method="REML")  
		  }
		  
		  PETPEESE <- rbind(
				data.frame(method="PET", tidyMR(PET)),
				data.frame(method="PEESE", tidyMR(PEESE))
		 )		
		 
		  # conditional PET/PEESE estimator
		  #the one-tail version that Stanley advocates: 

		  usePET <- ifelse(PETPEESE %>% filter(method == "PET", term == "b0") %>% .[["p.value"]] > .10, TRUE, FALSE)
		  PETPEESE <- rbind(PETPEESE, 
				data.frame(method="PETPEESE", if (usePET == TRUE) {tidyMR(PET)} else {tidyMR(PEESE)})
		  )
	  
		  PET.est <- PETPEESE %>% filter(method == "PET", term == "b0") %>% .[["estimate"]]
		  PET.slope <- PETPEESE %>% filter(method == "PET", term == "b1") %>% .[["estimate"]]
		  PEESE.est <- PETPEESE %>% filter(method == "PEESE", term == "b0") %>% .[["estimate"]]
		  PET_PEESE.est <- PETPEESE %>% filter(method == "PETPEESE", term == "b0") %>% .[["estimate"]]
		  PET_PEESE.text <- ifelse(PETPEESE %>% filter(method == "PET", term == "b0") %>% .[["p.value"]] > .10,
		  		"As the PET intercept does not significantly differ from zero (p > .10), it is recommended to use the PET estimator.", "As the PET intercept does significantly differ from zero (p < .10), it is recommended to use the PEESE estimator.")

			return(list(
				renderPlot({
					meta::funnel(meta1, ref=0, contour=c(0.9, 0.95), xlab="Effect size", cex=.5, pch=19)

					u <- par("usr")	# get range of plot coordinates

					# plot red dot at RE-MA
					points(meta2$b, u[3], cex=1.3, col="red", pch=20)

					# plot PET-line
					range <- seq(0, u[3], length.out=100)

					if (input$show_PET == TRUE) {
						# predict values from model
						PET.p <- PETPEESE %>% filter(method == "PET", term == "b0") %>% .[["estimate"]] +
									PETPEESE %>% filter(method == "PET", term == "b1") %>% .[["estimate"]]*range
						lines(PET.p, range, col="red")

						segments(coef(PET)[1], 0, coef(PET)[1], u[3], col="red", lty="dotted")
						points(coef(PET)[1], u[3], cex=1.3, col="red", pch=20)
					}

					# plot PEESE-line
					if (input$show_PEESE == TRUE) {
						PEESE.p <- PETPEESE %>% filter(method == "PEESE", term == "b0") %>% .[["estimate"]] +
									PETPEESE %>% filter(method == "PEESE", term == "b1") %>% .[["estimate"]]*range^2
						lines(PEESE.p, range, col="red")

						segments(coef(PEESE)[1], 0, coef(PEESE)[1], u[3], col="red", lty="dotted")
						points(coef(PEESE)[1], u[3], cex=1.3, col="red", pch=20)
					}

				}, width=400),
				HTML("<h4>A naive random effects meta-analysis (without bias correction)</h4>"),
				renderPrint(meta2),
				HTML(ifelse(is.null(Begg), "", 
					paste0("<h4>Begg & Mazumdar's (1994) test for publication bias</h4>Rank correlation of effect size vs. sample size (Kendall's <i>tau</i> = ", round(Begg$estimate, 2), " (", p(Begg$p.value), ")
				<p>A significant negative rank correlation indicates publication bias.</p>"))),
				HTML(paste0("<h4>Egger's test</h4>
					The slope of Egger's test is b1 = ", round(PET.slope, 2), ", t(", summary(PET.lm)$df[2], ") = ", 
					round(PETPEESE %>% filter(method == "PET", term == "b1") %>% .[["statistic"]], 3), "; ", 
					p0(PETPEESE %>% filter(method == "PET", term == "b1") %>% .[["p.value"]]),
					"<br>A significant slope with p < .10 is an indicator of small-study effects."
					)),
				HTML(paste0("<h4>PET: Bias corrected effect size estimate</h4>
					The intercept of the PET meta-regression is b0 = ", round(PET.est, 2), ", t(", summary(PET.lm)$df[2], ") = ", 
					round(PETPEESE %>% filter(method == "PET", term == "b0") %>% .[["statistic"]], 3), "; ", 
					p0(PETPEESE %>% filter(method == "PET", term == "b0") %>% .[["p.value"]]),
					"<br>A significant PET intercept with p < .10 indicates a bias corrected effect != 0. 
					The estimated true effect size is ", round(PET.est, 2), ", 95% CI [", 
					round(PETPEESE %>% filter(method == "PET", term == "b0") %>% .[["conf.low"]], 3), "; ",
					round(PETPEESE %>% filter(method == "PET", term == "b0") %>% .[["conf.high"]], 3), "]."
					)),
				HTML(paste0("<h4>PEESE: Bias corrected effect size estimate</h4>
					The intercept of the PEESE meta-regression is b0 = ", round(PEESE.est, 2), ", t(", summary(PEESE.lm)$df[2], ") = ", 
					round(PETPEESE %>% filter(method == "PEESE", term == "b0") %>% .[["statistic"]], 3), "; ", 
					p0(PETPEESE %>% filter(method == "PEESE", term == "b0") %>% .[["p.value"]]),
					". The estimated true effect size is ", round(PEESE.est, 2), ", 95% CI [", 
					round(PETPEESE %>% filter(method == "PEESE", term == "b0") %>% .[["conf.low"]], 3), "; ",
					round(PETPEESE %>% filter(method == "PEESE", term == "b0") %>% .[["conf.high"]], 3), "]."
					)),
				HTML(paste0("<h4>PET-PEESE conditional estimator</h4>",
					PET_PEESE.text,
					" The estimated true effect size from the conditional estimator is ", round(PET_PEESE.est, 2), ", 95% CI [", 
					round(PETPEESE %>% filter(method == "PETPEESE", term == "b0") %>% .[["conf.low"]], 3), "; ",
					round(PETPEESE %>% filter(method == "PETPEESE", term == "b0") %>% .[["conf.high"]], 3), "]."
					)),
				HTML('<hr><h4>Some comments on small-study effects</h4><p>In a set of studies with a fixed-<i>n</i> design and the same underlying effect, sample size should be unrelated to the estimated effect size (ES). A negative correlation between sample size and ES typically is seen as an indicator of publication bias and/or <i>p</i>-hacking. This bias is attempted to be corrected by meta-regression techniques such as <a href="http://onlinelibrary.wiley.com/doi/10.1002/jrsm.1095/abstract">PET-PEESE</a>.</p>
					<p>You should be aware, however, that some valid processes can also lead to a correlation between ES and N:
<ul>
<li>A) If (proper) sequential analyses are employed, trials with (randomly) lower sample effect sizes will take longer to stop. This process will also induce the correlation.</li>
<li>B) Imagine that different underlying effects are combined, and researchers did a proper a-priori power analysis, where they made a good guess about the true ES. Then they will plan larger samples for smaller effects, which will also introduce the correlation.</li>
<li>C) If more effective manipulations are more costly, one can expect larger effects from smaller samples.</li>
<li>D) Suppose a study with a small <i>n</i> tries to compensate by having a larger number of trials. This can result in a larger (standardized) effect size with a smaller sample and would contribute to a negative correlation between <i>n</i> and effect size (see Richard Morey\'s blog posts [<a href="http://bayesfactor.blogspot.de/2016/01/averaging-can-produce-misleading.html">1</a>][<a href="http://bayesfactor.blogspot.de/2016/01/asymmetric-funnel-plots-without.html">2</a>] for details)</li>
</ul>

On the other hand, proper sequential designs (A) are very rare yet (for an introduction to frequentist sequential designs, see <a href="http://ssrn.com/abstract_id=2333729">Lakens, 2014</a>; for an introduction to sequential Bayes factors, see <a href="http://ssrn.com/abstract_id=2604513">Schönbrodt, Wagenmakers, Zehetleitner, & Perugini, 2015</a>, and <a href="http://ssrn.com/abstract=2722435">Schönbrodt & Wagenmakers, 2016</a>). If different underlying effects are combined (B), we have a large heterogeneity in the meta-analysis, which is a problem for the model. C and D might be harder to detect without looking in detail into the methods sections.
</p>'),
					pancollapse.create(
					  "Detailed results for each test statistic",
					  getTable(ES_table)
					)
				))	
		} else {
			return(list(
				HTML("<h4>Too few effect sizes for plotting / meta-analysis!</h4><br>Enter > 2 effect sizes.")
			))
		}
	})
	

	
	tooltip <- function(x) {
	  if (is.null(x) | is.null(x$id)) return(NULL)
	  TBL <- dat$tbl %>% filter(!is.na(g))
	  TBL$g.abs <- abs(TBL$g)
	  TBL$label <- paste0("Row ", 1:nrow(TBL), ": ", TBL$paper_id, " ", TBL$study_id)
	  TBL$id <- 1:nrow(TBL)
	  
	  TBL[TBL$id == x$id, "label"]
	}
	
	
	reactive({
	  TBL <- dat$tbl %>% filter(!is.na(g))
	        
	   if (nrow(TBL) > 1) {  
	     TBL$g.abs <- abs(TBL$g)
	     TBL$label <- paste0("Row ", 1:nrow(TBL), ": ", TBL$paper_id, " ", TBL$study_id)
	     TBL$id <- 1:nrow(TBL)

	     x.limits <- logScaleLimits(range(TBL$n.approx, na.rm=TRUE))
	     x.values <- logScaleTicks(x.limits)
	     TBL %>% 
	       ggvis(x = ~n.approx, y = ~g.abs) %>%
	       layer_points(key := ~id) %>%
	       layer_model_predictions(model = "lm", se = FALSE, formula=g.abs~log(n.approx), stroke := COLORS$BLUE) %>%
	       add_axis("x", format="d", ticks=length(x.values), values=x.values, grid=TRUE, title="Approximate n (log scale)") %>%
	       add_axis("y", title="Absolute Hedge's g") %>%
	       scale_numeric("x", domain=x.limits, trans="log", nice=FALSE, expand=0) %>%
	       add_tooltip(tooltip, "click") 
	   } else {
	     print('Reactiv Expr: TBL doesnt exist')
	     
	     # dummy plot
	     me <- data.frame(x = 1, y = 1)
	     me %>%
	       ggvis(x = ~x, y = ~y) %>%
	       add_tooltip(tooltip, "click")

	   }
	}) %>% bind_shiny("ES_plot")
	
	
	
	# ---------------------------------------------------------------------
	# Research style panel (not implemented yet)
	
	output$researchstyle <- renderUI({
		if (nrow(dat$tbl) > 0) {			
			
			library(pwr)
			
			full.n <- 1000
			p_H1 <- 0.5
			alpha <- .05
			median.ES <- median(dat$tbl$g, na.rm=TRUE)
			median.n <- ceiling(median(dat$tbl$n.approx, na.rm=TRUE))
			median.pow <- pwr.t.test(n=median.n/2, d=median.ES)$power
		
			perc.sig.studies <- ((p_H1*median.pow) + ((1-p_H1)*alpha))
			sig.studies <- floor(perc.sig.studies* (full.n/median.n))
			perc.FPE <- ((1-p_H1)*alpha) / perc.sig.studies
			perc.replicable <- perc.FPE*alpha+ (1-perc.FPE)*median.pow
			
			return(list(
				HTML('<p class="text-warning">The test statistics are converted to Cohen`s d (resp. Hedge`s g) wherever possible, based on the formulas provided by Borenstein, Hedges, Higgins, & Rothstein (2011). Warning: These effect size conversions are based on approximative formulas. Although they work good under many conditions, this cannot replace a proper meta-analysis!</p>'),
				h3("Research style analysis"),
				HTML("<p>This analysis is based on an idea of <a href='https://willgervais.squarespace.com/s/Gervais-Jewell-Najle-Ng-SPPS-Power.pdf'>Gervais, Jewell, Najle, and Ng (2015)</a>, see also these blog post [<a href='http://willgervais.com/blog/2014/9/24/power-consequences'>1</a>][<a href='http://willgervais.com/blog/2015/5/14/a-powerful-nudge'>2</a>] by Will.</p>"),
				HTML(paste0("<p>This set of studies has a <b>median effect size of Hedge's g = ", round(median.ES, 3), "</b> and a <b>median approximative sample size of n = ", median.n, "</b>. These numbers translate to an expected power of ", round(median.pow * 100), "%.</p>")),
				HTML(paste0("<p>Suppose that a researcher has a pool of ", full.n," participants each year and runs studies in the style described above without <i>p</i>-hacking (but with selectively publishing only significant studies). A priori, hypotheses tend to be right ", p_H1*100, "% of the time.</p>
				<p>In the course of the year, such a researcher will accumulate <b>", sig.studies, " significant studies</b>. Of these ", sig.studies, " significant studies, <b>", round(perc.FPE*100), "% will be false-positives</b>. In exact replication attempts (same <i>n</i>), <b>", round(perc.replicable*100), "% will be succesfully replicated.</b></p>"))
			))	
		} else {
			return(NULL)
		}
	})
	
	
	
	
	
	# ---------------------------------------------------------------------
	# Export: save analysis as link
	# TODO: Also save the chosen analyis options
	
	output$exportlink <- renderUI({
		if (nrow(dat$tbl) > 0) {
			
			return(list(
				HTML(paste0("Copy and share <a href='http://shinyapps.org/apps/p-checker/?syntax=", URLencode(input$txt, reserved=TRUE), "'>this link</a> to send the p-checker analysis to others."))
			))	
		} else {
			return(NULL)
		}
	})
	
	
	
	## ======================================================================
	## The demodata tab
	## ======================================================================
	
	output$demodata <- renderUI({
		
		return(list(
			HTML('
			<style>
			  		 .actionButton .parent {
						  width: 300px;
						  height: 120px;
						  background-color: #fff;
						  border-radius: 5px;
						  border-radius: 5px;
        
						}
      
			      .actionButton {
			        padding: 5% 0;
			      }

						.actionButton .parent:hover {
						  box-shadow: 1px 1px 5px #999;
						}
			
						.buttonText {
			        padding: 10% 0;
			        text-align: left;
			        padding-left: 110px;
						}
      
			      img.buttonImg {
			        vertical-align: middle;
			        height:100px;
			      }
      
						</style>

						<div class="actionButton">
						    <button class="parent">
						        <img class="buttonImg" src="demo-pics/powerposing.jpg" align="left">
						        <div class="buttonText"><b>Power posing</b> p-curve analysis by Joe Simmons and Uri Simonsohn</div>
						    </button>
						</div>
			')
		))	
	})
	
	
	
	
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
		demo <- ""
		switch(input$demodata,
			"JPSP1" = {
				updateTextInput(session, inputId = "txt", value = readFile("demo-data/JPSP-p-curve.txt"))
				},
			"855" = {
				updateTextInput(session, inputId = "txt", value = readFile("demo-data/855_t_tests.txt"))
				},
			"H0_100x5" = {
				updateTextInput(session, inputId = "txt", value = readFile("demo-data/H0_100x5.txt"))
			},
			"H1_100x5" = {
				updateTextInput(session, inputId = "txt", value = readFile("demo-data/H1_100x5.txt"))
			},
			"H0_hack_100x5" = {
				updateTextInput(session, inputId = "txt", value = readFile("demo-data/H0_hack_100x5.txt"))
			},
			"elderly" = {
				updateTextInput(session, inputId = "txt", value = readFile("demo-data/elderly_priming.txt"))
			},
			"glucose" = {
				updateTextInput(session, inputId = "txt", value = readFile("demo-data/glucose.txt"))
			},
			"powerposing" = {
				updateTextInput(session, inputId = "txt", value = readFile("demo-data/powerposing.txt"))
			}
		)
		if (!is.null(con)) close(con)
	})
		

})
