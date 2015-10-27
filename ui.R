library(shiny)
library(shinythemes)
library(shinyTable)
library(shinyBS) # Additional Bootstrap Controls
library(ggvis)

# Load the panels with the manual etc.
source("snippets/quick_start.R")
source("snippets/responsibly.R")
source("snippets/extended_manual.R")
source("snippets/about.R")


shinyUI(fluidPage(theme = shinytheme("spacelab"),

# Accordion/ folding animation for quickstart etc.
	tags$head(tags$link(rel="stylesheet", type="text/css", href="accordion.css")),
	
	title = "One-for-all p-hacking detector",
	
	titlePanel("p-checker: The one-for-all p-value analyzer."),
	
	# display the accordion panels
	HTML(paste0('<div class="row">', qs_panel, responsibly_panel, '</div>')),	
	HTML(paste0('<div class="row">', extended_manual_panel, about_panel, '</div>')),
	
	# ---------------------------------------------------------------------
	# The actual app ...
	
	fluidRow(
		column(width=4,
			
				# the syntax input text field is constructed by ther server.R
				uiOutput("syntax"),
				
				downloadButton('downloadData','Save input as CSV file', class="btn-sm"),
				
				tags$hr(),
				tags$h3("Test-specific options"),
				
				checkboxInput("group_by_paper", "Group results by paper", FALSE),
				
				conditionalPanel(
					condition = "input.tabs1 == 'p-Curve' | input.tabs1 == 'R-Index' | input.tabs1 == 'TIVA'",
					checkboxInput("only_first_ES", "Only use first test statistic of each study", FALSE),
					helpText("Usually, only one effect size should be extracted for each sample. Manually choose the focal effect size, or use this checkbox to only include only the first ES of each study.")
				),
				
				conditionalPanel(
					condition = "input.tabs1 == 'p-Curve'",
					selectInput('pcurve_version','p-curve Version:', c(
						"Version 2 (chi2 test)"="v2",
						"Version 3 (Z test) - recommended"="v3"
					), selected="v3")
				),
				
				conditionalPanel(
					condition = "input.tabs1 == 'Meta-analysis (beta)'",
					selectInput('meta_ES_type','Test type for meta analysis:', c(
						"t test (one group)"="ttest_1",
						"t test (two group) & F test (1 df)"="ttest_2",
						"Correlation"="cor"
					), selected="ttest_2")
				),
				
				conditionalPanel(
					condition = "input.tabs1 == 'p-Curve' & input.experimental == 1",
					sliderInput("pcurve_crit", "Critical p value (EXPERIMENTAL! Only intended for exploration, not for actual p-curve analyses! Default = .05)", min=.01, max=.10, value=.05, step=.01),
					sliderInput("pcurve_power", "Power (EXPERIMENTAL! Default = 33%)", min=15, max=90, value=33, step=1)
				),
				
				conditionalPanel(
					condition = "input.tabs1 == 'R-Index'",
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
		),		
		
		
		
		
		# ---------------------------------------------------------------------
		# The output panels, on the right side
		
		column(width=8, 

			# show warning if experimental features are activated
			htmlOutput("experimental_warning"),

			# show potential parser errors on top of output
			htmlOutput("parser_errors"),
			
			tabsetPanel(id ="tabs1",				
				tabPanel("R-Index",					
					htmlOutput("rindex_summary"),
					conditionalPanel(
						condition = "input.group_by_paper == 1",
						downloadButton('downloadRIndex','Save R-Index results as CSV file', class="btn-sm")
					),
					HTML('<small>For information about R-Index, see <a href="http://www.r-index.org/">http://www.r-index.org/</a>.</small>'),
					htmlOutput("rindex_table")
				),
				tabPanel("TIVA",					
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
				# tabPanel("Meta-analysis (beta)",
# 					htmlOutput("meta")
# 				),
				tabPanel("Effect-sizes (beta)",
					HTML('<p class="text-warning">The test statistics are converted to Cohen`s d (resp. Hedge`s g) wherever possible, based on the formulas provided by Borenstein, Hedges, Higgins, & Rothstein (2011). Warning: These effect size conversions are based on approximative formulas. Although they work good under many conditions, this cannot replace a proper meta-analysis!</p>'),
					htmlOutput("effectsizes")
				),
				# tabPanel("Research style analysis (beta)",
				# 	htmlOutput("researchstyle")
				# ),
				tabPanel("p values correctly reported?",
					htmlOutput("report_table")
				),
				tabPanel("Export", 					
					tableOutput("export")
				)
			)
		)
	)	
))
