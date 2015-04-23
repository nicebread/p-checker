about_panel <- '
<!-- #########  ABOUT WELL  #########   -->

<div class="col-sm-6">
<div class="panel-group" id="accordion3">
  <div class="panel panel-info">
     <div class="panel-heading" data-toggle="collapse" data-parent="#accordion3" data-target="#collapse3">
      <h4 class="panel-title accordion-toggle">About</h4>
    </div>
    <div id="collapse3" class="panel-collapse collapse">
      <div class="panel-body">
	  <i>(c) 2015 by <a href="mailto:felix@nicebead.de">Felix Schönbrodt</a> (<a href="http://www.nicebread.de">www.nicebread.de</a>). The source code of this app is licensed under the open GPL-2 license and is published on <a href="https://github.com/nicebread/p-checker">Github</a>.</i>
	  
	  <br><br>
	  
	  This Shiny app implements the <b>p-curve</b> (Simonsohn, Nelson, & Simmons, 2014; see <a href="http://www.p-curve.com">http://www.p-curve.com</a>) in its previous ("app2") and the current version ("app3"), the <b>R-Index</b> and the <b>Test of Insufficient Variance, TIVA</b> (Schimmack, 2014; see <a href="http://www.r-index.org/">http://www.r-index.org/</a>), and tests whether <b><i>p</i> values are reported correctly</b>.
	  <br/><br>
	  p-curve code is to a large extent adapted or copied from Uri Simonsohn (see <a href="http://p-curve.com/Supplement/Rcode_other/R%20Code%20behind%20p-curve%20app%203.0%20-%20distributable.R">here</a>). TIVA code adapted from Moritz Heene.
	  <br/><br>
	  
	  <h3>Citation</h3>
	  Programming this app took a considerable effort and amount of time. If you use it in your research, please consider citing the app, and of course the creators of the statistical tests:
	  <br/><br/>
	  
	  Simonsohn, U., Nelson, L. D., & Simmons, J. P. (2014). P-curve: A key to the file-drawer. <i>Journal of Experimental Psychology: General, 143</i>, 534–547. doi:10.1037/a0033242
	  <br/><br/>
	  Schimmack, U. (2014). <i>Quantifying Statistical Research Integrity: The Replicability-Index</i>. Retrieved from http://www.r-index.org
	  <br/><br/>
	  Schönbrodt, F. D. (2015). <i>p-checker: One-for-all p-value analyzer.</i> Retrieved from http://shinyapps.org/apps/p-checker/.
	  <br/><br/>	  
	  
	  
	  <h3>Disclaimer / Validity of the results</h3>
	  This app is still beta. That means, the UI and the format of exported data might still change.
	  <br>
	  I cross-validated the results with p-curve.com and did not find differences (unsurprisingly, as I use Uri`s code for p-curve to a large extent). With a single click (see the "Export" tab) you can transfer the test statistics to <a href="http://www.p-curve.com">p-curve.com</a> and cross-validate the results yourself. 
	  I also checked the results with the <a href="http://www.r-index.org/">R-Index Excel-sheet</a> and did not find differences so far.
	  <br>Nonetheless, this app could contain errors and a healthy scepticism towards the results is indicated. I always recommend to perform some plausibility checks. Feel free to go to the <a href="https://github.com/nicebread/p-checker">source code</a> and check the validity yourself. If you suspect a bug or encounter errors, please send me an <a href="mailto:felix@nicebread.de">email</a> with your test statistics and a description of the error.
	  
	  <h3>Comments</h3>
	  Any detected bugs, comments, and feature requests are welcome: <a href="mailto:felix@nicebread.de">felix@nicebread.de</a>
	  <br>
	  
	  <h3>Demo data sets</h3>
	  
	  <ul>
	  <li><i>Non-hacked JPSP data</i>: See Simonsohn, Nelson, & Simmons (2014), Figure 3B. Retrieved from <a href="http://www.p-curve.com/Supplement/full_pdt.xlsx">http://www.p-curve.com/Supplement/full_pdt.xlsx</a></li>
	  <li><i>855 t-tests</i>: See Wetzels, Matzke, Lee, Rouder, Iverson, & Wagenmakers (2011). Retrieved from <a href="http://www.ejwagenmakers.com/2011/effectsize_data.zip">http://www.ejwagenmakers.com/2011/effectsize_data.zip</a></li>
	  <li><i>Elderly priming</i>: See Lakens, D. (2014). Professors are Not Elderly: Evaluating the Evidential Value of Two Social Priming Effects Through P-Curve Analyses. Retrieved from <a href="http://papers.ssrn.com/sol3/papers.cfm?abstract_id=2381936">http://papers.ssrn.com/sol3/papers.cfm?abstract_id=2381936</a>. Data available at <a href="https://osf.io/3urp2/">https://osf.io/3urp2/</a></li>
	  </ul>
	  
	  <small>
	  Simonsohn, U., Nelson, L. D., & Simmons, J. P. (2014). P-curve: A key to the file-drawer. <i>Journal of Experimental Psychology: General, 143</i>, 534–547. doi:10.1037/a0033242
	  <br><br>
	  Wetzels, R., Matzke, D., Lee, M. D., Rouder, J. N., Iverson, G. J., & Wagenmakers, E.-J. (2011). Statistical evidence in experimental psychology: An empirical comparison using 855 t tests. <i>Perspectives on Psychological Science, 6</i>, 291–298. doi:10.1177/1745691611406923
	  <br><br>
	  Lakens, D. (2014). Professors are not elderly: Evaluating the evidential value of two social priming effects through p-curve analyses. doi: http://dx.doi.org/10.2139/ssrn.2381936. Retrieved from http://ssrn.com/abstract=2381936
	  </small>

	  <h3>Version 0.4</h3>

	  Known issues / Todos:<br/><br/>
	  <ul>
	  <li>p-curve plot: 33% power curve missing</li>
	  <li>Clearly separate the inference functions from UI functions</li>
	  <li>Make TIVA computation robust against outliers.</li>
	  <li>When sharing analysis via a link, only the test statistics are exported, not the analysis settings!</li>
	  </ul>	  
	  
      </div>
    </div>
  </div>
</div>
</div>
'