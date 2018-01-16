library(shiny)
library(shinyjs)
library(shinythemes)
library(shinyBS) # Additional Bootstrap Controls
#library(ggvis)

# Load the panels with the manual etc.
source("pancollapse.R")

# custom js function to open external URL
jsCode <- "shinyjs.browseURL = function(URL){window.open(URL, ''); ;}"

shinyUI(tagList(
	
	# https://github.com/daattali/shinyjs#using-shinyjs-with-navbarpage-layout
	useShinyjs(),
	extendShinyjs(text = jsCode),
	
	navbarPage(title="", 
	tabPanel("p-checker", 
  		# ---------------------------------------------------------------------
  		# The actual app ...
  		HTML("<h2><strong><i>p</i>-checker</strong> The one-for-all <i>p</i>-value analyzer</h2>"),
  		fluidRow(
		
  		column(width=5,  		
		div(id="leftCol",	
			
				h3("Enter test statistics here:"),
  				# the syntax input text field is constructed by ther server.R
  				uiOutput("syntax"),
  				br(),
  				downloadButton('downloadData','Save input as CSV file', class="btn-sm"),
				uiOutput("exportlink"),
  				
  				tags$hr(),
  				tags$h3("Test-specific options"),
  				
  				conditionalPanel(
  					condition = "input.tabs1 == 'p-Curve' | input.tabs1 == 'Excess Significance' | input.tabs1 == 'TIVA'",
					checkboxInput("group_by_paper", "Group results by paper", FALSE),
  					checkboxInput("only_first_ES", "Only use first test statistic of each study", FALSE),
  					helpText("Usually, only one effect size should be extracted for each sample. Manually choose the focal effect size, or use this checkbox to only include only the first ES of each study.")
  				),
  				
  				conditionalPanel(
  					condition = "input.tabs1 == 'p-Curve'",
  					selectInput('pcurve_version','p-curve Version:', c(
  						"Version 2 (chi2 test)"="v2",
  						"Version 3 (Z test) - recommended"="v3"
  					), selected="v3"),
					sliderInput("pcurve_power", "Comparison power (default = 33%)", min=10, max=99, value=33, step=1)
  				),
  				
  				conditionalPanel(
  					condition = "input.tabs1 == 'Meta-analysis'",
  					checkboxInput("show_PET", "Show PET meta-regression in plot", TRUE),
					checkboxInput("show_PEESE", "Show PEESE meta-regression in plot", TRUE),
	  				selectInput('MR_model','Meta-regression model', c(
						"Using the lm() function"="lm",
						"Using the rma() function"="rma"
	  				), width="100%")
  				),
  				
  				conditionalPanel(
  					condition = "input.tabs1 == 'p-Curve' & input.experimental == 1",
  					sliderInput("pcurve_crit", "Critical p value (EXPERIMENTAL! Only intended for exploration, not for actual p-curve analyses! Default = .05)", min=.01, max=.10, value=.05, step=.01)  					
  				),
  				
  				conditionalPanel(
  					condition = "input.tabs1 == 'Excess Significance'",
  					checkboxInput("omit_nearly_significant", "Omit 'nearly significant' p-values (range: see below) from R-Index analysis.", FALSE),
  					sliderInput("omit_nearly_significant_range", "Range of 'nearly significant'", min=.0, max=.20, value=c(.05, .10), step=.005)
  				),
  				
  				tags$hr(),
  				tags$h3("General options"),
  				
  				numericInput("digits", "Digits in display:", 3, min = 0, max = 5),
  				checkboxInput("round_up", "Gracious rounding up", FALSE),
  				helpText("If the t value is reported as 2.1, it could also be 2.14999 which has been rounded down. If you want to be maximally generous, you can check this box, and all test statistics are automatically increased by X.XX4999."),
  				
  				br(),br(),
  				selectInput('demodata','Load demo data', c(
  					"---"="---",
					"Power posing by @jpsimmon"="powerposing",
					"Glucose and self-control by @mavadillo"="glucose",
  					"Elderly priming analysis by @lakens"="elderly",
  					"Non-hacked JPSP data (Simonsohn et al., 2014, Figure 3B)"="JPSP1",
  					"855 t-tests (Wetzels et al., 2011)"="855",
  					"H0 sim: 100 papers with 5 studies; d = 0; selective reporting"="H0_100x5",
  					"H1 sim: 100 papers with 5 studies; d = 0.5; selective reporting"="H1_100x5",
  					"Hack sim: 100 papers with 5 studies; d = 0; hacked; selective reporting"="H0_hack_100x5"
  				), width="100%"),
  				
  				br(),
  				checkboxInput("experimental", "Activate experimental options (Do not run actual analyses with these experimental/untested options!)", FALSE),
  				bsPopover(id = "experimental", title="A", content = "Do not run actual analyses with these experimental/untested options!", placement = "right", trigger = "hover")
  		) # end of div id="leftCol"		
  		), 
  		
  		
  		
## ======================================================================
## The output panels, on the right side
## ======================================================================
  		
	column(width=7, 
			
			#alert.create('<strong>New feature:</strong> You can now enter <i>p</i>-values directly (e.g. <code>p=0.021</code>). If you provide the sample size in addition (e.g. <code>p(48)=.03</code>), the <i>p</i>-value is also converted into an effect size.', style="success"),
			
			alert.create("<strong>Disclaimer:</strong> This web application provides several tests for publication bias/p-hacking/indicators for data-dependent analyses. Some of them are new, unpublished, and controversial to some extent; purpose of this app is to provide a unified place for trying out and comparing these methods. Please use the tests with care, and RTM of the tests. (You can dismiss this message by clicking 'X')", style="info"),
  
  			# show ROxygen-style title if provided in the syntax
  			htmlOutput("roxygen_title"),
  
  			# show warning if experimental features are activated
  			htmlOutput("experimental_features_warning"),
  
  			# show potential parser errors on top of output
  			htmlOutput("parser_errors"),
  			
  			tabsetPanel(id ="tabs1",				
  				tabPanel("Excess Significance",					
  					htmlOutput("rindex_summary"),
  					conditionalPanel(
  						condition = "input.group_by_paper == 1",
  						downloadButton('downloadRIndex','Save R-Index results as CSV file', class="btn-sm")
  					),
  					HTML('<small>For information about R-Index, see <a href="http://www.r-index.org/">http://www.r-index.org/</a>.</small>'),
  					htmlOutput("rindex_table")
  				),
  				tabPanel("TIVA",			
					alert.create('The TIVA test expects that all entered test statistics/p-values are in the expected direction (regardless of the sign). Please delete or comment out all rows with results in the "wrong" direction.'),		
  					htmlOutput("tiva_summary"),
  					conditionalPanel(
  						condition = "input.group_by_paper == 1",
  						downloadButton('downloadTIVA','Save TIVA results as CSV file', class="btn-sm")
  					),
  					HTML('<small>For information about TIVA, see <a href="https://replicationindex.wordpress.com/2014/12/30/the-test-of-insufficient-variance-tiva-a-new-tool-for-the-detection-of-questionable-research-practices/comment-page-1/#comment-92">replicationindex.wordpress.com</a>.</small>'),
  					htmlOutput("tiva_table")
  				),
  				tabPanel("p-Curve", 
  					conditionalPanel(
  						condition = "input.group_by_paper == 0",
  						htmlOutput("pcurve_plot")
  					),					
  					htmlOutput("pcurve_summary"),
  					conditionalPanel(
  						condition = "input.group_by_paper == 1",
  						downloadButton('downloadPCurve','Save p-curve results as CSV file', class="btn-sm")
  					),
  					HTML('<small>For information about p-curve, see <a href="http://p-curve.com/">http://p-curve.com/</a>.<br>
  					Simonsohn, U., Nelson, L. D., & Simmons, J. P. (2014). P-curve: A key to the file-drawer. <i>Journal of Experimental Psychology: General, 143</i>, 534–547. doi:10.1037/a0033242					
  					</small>'),
  					tableOutput("pcurve_table")
  				),
  				tabPanel("Meta-analysis",
  				  br(),
  					alert.create('The test statistics are converted to Cohen\'s d wherever possible, based on the formulas provided by Borenstein, Hedges, Higgins, & Rothstein (2011). <strong>Warning:</strong> These effect size conversions are based on approximative formulas; furthermore the app always assumes equal cell sizes and other simplifications. Although these proxies work good under many conditions, this quick meta-analytic overview <i>cannot</i> replace a proper meta-analysis!'),
  					#ggvisOutput("ES_plot"),
		  			# show potential parser errors on top of output
		  			htmlOutput("MA_warnings"),
  					htmlOutput("effectsizes")
  					
  				),
  				# tabPanel("Research style analysis (beta)",
  				# 	htmlOutput("researchstyle")
  				# ),
  				tabPanel("p values correct?",
  					htmlOutput("report_table")
  				)#,
  				# tabPanel("Export",
#   					tableOutput("export")
#   				)#,
  				#tabPanel("Demo data",
  				#	htmlOutput("demodata")
  				#)
  			)
  		)
  	)	
  ),
  tabPanel('Quick Start', loadHTML('snippets/quick_start.html')),
  tabPanel('Manual', loadHTML('snippets/extended_manual.html')),
  tabPanel('Terms of Use', loadHTML('snippets/responsibly.html')),
  tabPanel('About', loadHTML('snippets/about.html')),
  tabPanel('Release Notes', loadHTML('snippets/version_history.html')),
  header = pancollapse(),
  theme = shinytheme("spacelab"), 
  # load custom css to override some theme settings
  tags$head(tags$link(rel = "stylesheet", type = "text/css", href = "p-checker-theme.css")),
  windowTitle = "One-for-all p-value analyzer"
)))
