library(shiny)
#library(shinyTable)
library(stringr)
library(dplyr)
library(ggplot2)
library(ggvis)

# source inference functions
source("helpers.R")
source("parser.R")
source("p-curve.R")
source("TIVA.R")


#input <- list(round_up=FALSE, digits=3, group_by_paper=TRUE, only_first_ES=TRUE, txt=x, pcurve_power=33, pcurve_crit=.05, experimental=FALSE); dat <- list()


shinyServer(function(input, output, session) {

	# dat is a reactive object that keeps the computed variables
	dat <- reactiveValues(
		tblDisplay=data.frame(),	# keeps the rounded values for display
		tbl=data.frame(),			# keeps the precise values
		warnings=""					# keeps a vector of parser errors
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
					"z" = {res <- c(res, paste0("Z=", dat$tbl$statistic[i]))}
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
		  
		  if (is.null(query["syntax"]) | query["syntax"] == "NULL") {
			  res <- paste(readLines("snippets/demo_syntax.txt", encoding="UTF-8"), collapse="\n")
		  } else {
			  res <- query["syntax"]
		  }

		return(list(
			HTML(paste0('<textarea class="form-control" style="resize:none;white-space:pre;word-wrap:normal;overflow-x:scroll;" id="txt" rows="18">', res, '</textarea>'))
		))
	 })
	
	
	
	
	
	# every time the text field is changed, this function is called and parses the input string
	observe({
		
		# quit when syntax field is not created yet or empty
		if (is.null(input$txt)) {
			return()
		}
		
		tbl <- parse_ES(input$txt, round_up=input$round_up)
		
		# parser errors present?
		if (length(attr(tbl, "warnings")) > 0) {
			dat$warnings <- attr(tbl, "warnings")
		} else {
			dat$warnings <- ""
		}
		
		# No input? Return empty data frame
		if (is.null(tbl) || nrow(tbl) == 0) {
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
				pancollapse.create(
				  "Detailed results for each test statistic:",
				  getTable(report_table)
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
				HTML("<h2>p-reporting analysis: Are there wrong reported p values?</h2>"),	
				
				#HTML("<h2>Results for each paper</h2>"),
				
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
				  getTable(rindex_table)
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
			obs_power0 <- tbl %>% group_by(paper_id, study_id) %>% filter(row_number() <= 1) %>% select(median.obs.pow)
			obs_power <- median(obs_power0$median.obs.pow, na.rm=TRUE)
			inflation_rate <- success_rate - obs_power
			r_index <- obs_power - inflation_rate
		
			tiva <- TIVA(tbl$p.value)
			
			result <- paste0(
				"<h2>R-Index analysis:</h2><h4>",
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
			obs_power0 <- tbl %>% filter(row_number() <= 1) %>% select(median.obs.pow)
			obs_power <- obs_power0 %>% group_by(paper_id) %>% select(median.obs.pow) %>% summarise_each(funs(median))
			
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
					"<b>Average median observed power</b> = ", round(mean(rindex$median.obs.pow, na.rm=TRUE), input$digits), "<br>",
					"<b>Average inflation rate</b> = ", round(mean(rindex$inflation_rate, na.rm=TRUE), input$digits), "<br>",
					"<b>Average R-Index</b> = ", round(mean(rindex$r_index, na.rm=TRUE), input$digits),
					"</h4>"
				)),
				#renderTable({rindex})
				pancollapse.create(
				  "Detailed results for each test statistic",
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
			tiva <- TIVA(tbl$p.value)

			result <- paste0(
				"<h2>Test of insufficient variance (TIVA)</h2>",
				"<small>Variances &lt; 1 suggest bias. The chi2 tests the H0 that variance <= 1.</small>",

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
			}),
			send2pcurve.button.tag,
			HTML(paste0("<br><small>Note: This transfers the test statistics without paper identifier. That means, p-curve.com will compute an omnibus test with all values.</small><br>",
		"<small>Note: If you want identical results to p-curve.com, turn off the 'Gracious rounding up' option at the left panel.</small><br><br>
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
				
				
				pancollapse.create(
				  "Detailed results for each test statistic",
				  getTable(pcurve)
				)
			))
		}
	})
	
	
	observeEvent(input$send2pcurve, {					
		res1 <- paste(exportTbl(), collapse="\n")		
		pcurve_link <- paste0("http://www.p-curve.com/app3/?tests=", URLencode(res1, reserved=TRUE))
		browseURL(pcurve_link)
	})
	
	# ---------------------------------------------------------------------
	# Effect size panel	
	
	output$effectsizes <- renderUI({
		
	  print("I'm here")
	  
		TBL <- dat$tbl %>% filter(!is.na(g))
		
		if (nrow(TBL) > 0) {
			
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
			
			ES_table <- dat$tblDisplay[!is.na(dat$tblDisplay$g), c("paper_id", "study_id", "type", "df1", "df2", "statistic", "p.value", "n.approx", "d", "g")]
			
			#tooltip <- function(x) {
			  #if (is.null(x) | is.null(x$id)) return(NULL)
			  #return(TBL[TBL$id == x$id, "label"])
			#}
			
			# http://stackoverflow.com/questions/29785281/ggvis-with-tooltip-not-working-with-layer-smooths
		
			  
			
			#TBL %>% 
			  #ggvis(x = ~n.approx, y = ~g.abs) %>%
			  #layer_points(key := ~id) %>%
			  #layer_model_predictions(model = "lm", se = FALSE, formula=g.abs~log(n.approx), stroke := "blue") %>%
			  #add_axis("x", ticks = 5, values = round(seq(min(TBL$n.approx, na.rm=TRUE), max(TBL$n.approx, na.rm=TRUE), length.out=5)), grid=TRUE, title="Approximate n (log scale)", format="d") %>%
 			  #add_axis("y", title="Absolute Hedge's g") %>% 			  
 			  #scale_numeric("x", trans="log") %>%			  			  
  			#add_tooltip(tooltip, "hover") %>%
			  #bind_shiny("ES_plot")
		  })
		  
		  
			return(list(
				#renderPlot({ES_plot}, res=120),
				#HTML(paste0("<h4>r = ", round(cor(log(TBL$n.approx), TBL$g, use="p"), 2), "</h4>")),
				HTML("<h2>Effect size vs. sample size: <i>r</i> = ", round(cor(log(TBL$n.approx), TBL$g, use="p"), 2), "</h2>"),
					HTML('<p>In a set of studies with a fixed-<i>n</i> design and the same underlying effect, sample size should be unrelated to the estimated effect size (ES). A negative correlation between sample size and ES typically is seen as an indicator of publication bias and/or <i>p</i>-hacking. This bias is attempted to be corrected by meta-analytic techniques such as <a href="http://onlinelibrary.wiley.com/doi/10.1002/jrsm.1095/abstract">PET-PEESE</a>.</p>
					<p>You should be aware, however, that some valid processes can also lead to a correlation between ES and N:
<ul>
<li>A) If (proper) sequential analyses are employed, trials with (randomly) lower sample effect sizes will take longer to stop. This process will also induce the correlation.</li>
<li>B) Imagine that different underlying effects are combined, and researchers did a proper a-priori power analysis, where they made a good guess about the true ES. Then they will plan larger samples for smaller effects, which will also introduce the correlation.</li>
</ul>

On the other hand, proper sequential designs (A) are very rare yet (for an introduction to frequentist sequential designs, see <a href="http://papers.ssrn.com/sol3/papers.cfm?abstract_id=2333729">Lakens (2014)</a>; for an introduction to sequential Bayes factors, see <a href="http://papers.ssrn.com/sol3/papers.cfm?abstract_id=2604513">Schönbrodt, Wagenmakers, Zehetleitner, & Perugini (2015)</a>). If different underlying effects are combined (B), we have a large heterogeneity in the meta-analysis, which is a problem for the model.
</p>'),
				pancollapse.create(
				  "Detailed results for each test statistic",
				  getTable(ES_table)
				)
			))	
		} else {
			return(NULL)
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
	        
	   if (nrow(TBL) > 0) {  
	     TBL$g.abs <- abs(TBL$g)
	     TBL$label <- paste0("Row ", 1:nrow(TBL), ": ", TBL$paper_id, " ", TBL$study_id)
	     TBL$id <- 1:nrow(TBL)
	     
	    
	     
	     print('Reactiv Expr: TBL exists')
	     
	     TBL %>% 
	       ggvis(x = ~n.approx, y = ~g.abs) %>%
	       layer_points(key := ~id) %>%
	       layer_model_predictions(model = "lm", se = FALSE, formula=g.abs~log(n.approx), stroke := "blue") %>%
	       add_axis("x", ticks = 5, values = round(seq(min(TBL$n.approx, na.rm=TRUE), max(TBL$n.approx, na.rm=TRUE), length.out=5)), grid=TRUE, title="Approximate n (log scale)", format="d") %>%
	       add_axis("y", title="Absolute Hedge's g") %>% 			  
	       scale_numeric("x", trans="log") %>%	  			  
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
	# Effect size panel	
	
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
				h2("Research style analysis"),
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
	# Export for p-curve; save analysis as link
	
	output$export <- renderUI({
		if (nrow(dat$tbl) > 0) {
			
			res2 <- paste(exportTbl(), collapse="<br>")
			
			return(list(
				HTML(paste0("<h4>Copy and share <a href='http://shinyapps.org/apps/p-checker/?syntax=", URLencode(input$txt, reserved=TRUE),
				"'>this link</a> to send the p-checker analysis to others.</h4>")),
								
				HTML(paste0("<small>",res2, "</small>"))
			))	
		} else {
			return(NULL)
		}
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
