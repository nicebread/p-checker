extended_manual_panel <- '
<div class="col-sm-6">
	
	<!-- #########  EXTENDED MANUAL  #########   -->
	
	<div class="panel-group" id="accordion2">
	  <div class="panel panel-primary">
	     <div class="panel-heading" data-toggle="collapse" data-parent="#accordion2" data-target="#collapse2">
	      <h4 class="panel-title accordion-toggle">Extended Manual</h4>
	    </div>
	    <div id="collapse2" class="panel-collapse collapse">
	      <div class="panel-body">

		  <h3 id="data_extraction">Data extraction:</h3>

		  <ul>
		  <li><em>All</em> p-values can be extracted, both from focal hypothesis tests and from ancillary analyses, such as manipulation checks. But only p values are extracted for which precise dfs are reported (i.e., results such as &#8220;Fs &lt; 1, ps > .50&#8221; are <em>not</em> extracted).</li>
		  <li>Format:
		  <ul>
		  <li>Study ID: <em>teststatistic</em>; [optional] reported <em>p</em> value; [optional] 
		  critical <em>p</em> value; [optional, if one-tailed testing] one-tailed
		  <ul>
		  <li>[optional] reported <em>p</em> value: e.g., <code>p = .03</code>, or <code>p &lt; .05</code></li>
		  <li>[optional] critical <em>p</em> value: e.g., <code>crit = .10</code>, or <code>crit = .08</code></li>
		  <li>[optional, if one-tailed testing]: write the keyword <code>one-tailed</code>, or just <code>one</code>, or <code>1t</code></li>
		  </ul></li>
		  <li>The colon separates study ID from everything else</li>
		  <li>If the study ID starts with an underscore, this test statistic is <em>not</em> a focal test (e.g., from a manipulation check, a pre-test, or an ancillary analysis for possible alternative explanations), and will not be included in R-Index or p-curve analyses (but it will be included in the test for correct p-values)</li>
		  <li>The first datum after the colon must be the test statistic</li>
		  <li>All optional informations are separated by semicolons; can be given in any order</li>
		  <li>At the end of a line a comment can be written after a # sign (everything after the # is ignored)</li>
		  </ul></li>
		  <li>Examples:
		  <ul>
		  <li>M&amp;E (2005) S1: t(25) = 2.1; p &lt; .05; one-tailed</li>
		  <li>M&amp;E (2005) S2: F(1, 45) = 4.56; p = .03   # wrong p value?</li>
		  <li>M&amp;E (2005) S3: chi2(1) = 3.7; crit=.10</li>
		  <li>_M&amp;X (2011) S1: r(123) = .08; p = .45     # this was a manipulation check (see underscore)</li>
		  </ul></li>
		  <li>Be careful if you <strong>copy &amp; paste</strong> the results from a PDF:
		  <ul>
		  <li>Sometimes there are invisible special characters. They are shown in the app as weird signs and must be removed.</li>
		  <li>The minus sign sometimes looks a bit longer (an &#8220;em-dash&#8221;). This should be replaced with a standard minus sign.</li>
		  </ul></li>
		  <li>Which tests to select in the presence of <strong>interactions</strong>? Some hints from Simonsohn et al.&#8217;s (2014) p-curve paper:
		  <ul>
		  <li>&#8220;When the researcher’s stated hypothesis is that the interaction attenuates the impact of X on Y (e.g., people always sweat more in summer, but less so indoors), the relevant test is whether the interaction is significant (Gelman &amp; Stern, 2006), and hence p-curve must include only the interaction’s p-value. [&#8230;] Simple effects from a study examining the attenuation of an effect should not be included in p-curve, as they bias p-curve to conclude evidential value is present even when it is not.&#8221;</li>
		  <li>&#8220;When the researcher’s stated hypothesis is that the interaction reverses the impact of X on Y (e.g., people sweat more outdoors in the summer, but more indoors in the winter), the relevant test is whether the two simple effects of X on Y are of opposite sign and are significant, and so both simple effects’ p-values ought to go into p-curve. The interaction that is predicted to reverse the sign of an effect should not be included in p-curve, as it biases p-curve to conclude evidential value is present even when it is not.&#8221;</li>
		  </ul></li>
		  </ul>

		  <h3 id="special_cases">Special cases</h3>

		  <ul>
		  <li>Only significant values are included in R-Index analyses. But sometimes marginally non-significant p values (e.g., p = .051) are falsely rounded downwards and cross the critical boundary only due to this error (i.e., they are reported as &#8220;p &lt; .05)&#8221;). In this case, the ES is not included in the analysis (see column &#8220;significant&#8221; in the R-Index tab). But, if the ES has been (falsely) interpreted as significant, the critical value must be slightly increased, so that the ES is also included in the R-Index analysis. In this case, increase the critical level to <code>crit = .055</code>.</li>
		  </ul>
		  
		  <h3>Reproducible Analyses</h3>

		  <ul>
		  <li>You can provide the test statistics also in the link to this app. Copy your input text to the upper box of this <a href="http://www.freeformatter.com/url-encoder.html#ad-output">URL Encoder</a>. In the lower box, you will find a weird-looking string, such as <code>t%2847%29%3D2.1%0D%0Ar%2834%29%3D0.123%0D%0AF%281%2C+200%29%3D9.9</code>. Now you can provide that string after the link to this app, for example:<br>
		  <code>http://shinyapps.org/apps/p-checker/?syntax=t%2847%29%3D2.1%0D%0Ar%2834%29%3D0.123%0D%0AF%281%2C+200%29%3D9.9</code><br>
		  This way you can share any p-value analysis in a single link!
		  
		  Or, even easier: Go to the "Export" tab and copy the link for a reproducible analysis!
		  
		  </li>
		  
		  </ul>

	      </div>
	    </div>
	  </div>
	</div>
	</div>
'