# rstpm2

# Version 1.5.2
    - Add Clayton copulas with arbitrary cluster sizes to the
      parametric GSMs (experimental)
	- Add spline interpolation support in the Markov models (faster
      for some models)

# Version 1.5.1
    - Bug fix: mean predictions

# Version 1.5.0
	- Major change: markov_msm function for Markov multistate models
	- Add predict(..., type="lpmatrix")
	- Add cluster and bhazard specials
	- Internal: use Nelder-Mead for AFT if the model did not converge
	- Internal: refactor stpm2 and pstpm2 to use a common gsm function
	- Internal: extended the test suite
	- Documentation: update vuniroot vignette
	- Bug fixes: delayed entry; missing bhazard and cluster values; 

# Version 1.4.5
    - Fixed a bug in fitting frailty models (introduced in 1.4.4)
	- Introduced package tests
	
# Version 1.4.4
    - Fixed a critical bug in the `predict` function for comparisons of hazards, including type="hr", type="hdiff" and type="marghr" (introduced in 1.4.2).

# Version 1.4.2
    - Belatedly started the NEWS.md file
    - Update to bbmle (>= 1.0.20) required due to new export from that package
    - Possible breaking change: for the `predict()` functions for `stpm2` and `pstpm2`, the `keep.attributes` default has changed from `TRUE` to `FALSE`. Any code that used `predict()` and needs the `newdata` attributes should now add the `keep.attributes=TRUE` argument. The previous default was noisy.
	- Possible breaking change: the derivative of the design matrix with respect to time now defaults to being calculated using log(time); the old calculation can be found using `log.time.transform=TRUE`. This is expected to provide more accurate gradients, particularly for very small times. 
    - To this point, the following models are available: 
      + `stpm2`: parametric generalised survival models, possibly with clustered data (Gamma frailties and normal random effects), relative survival, robust standard errors, rich post-estimation and plots.
      + `pstpm2`: penalised generalised survival models, possibly with clustered data (Gamma frailties and normal random effects), relative survival, robust standard errors, rich post-estimation and plots.
      + `aft`: parametric accelerated failure time models, with more limited post-estimation and plots.
	- Links for the generalised survival models include log-log, -logit, -probit, -log and Aranda-Ordaz.
    - Post-estimation for `stpm2` and `pstpm2` includes:
	  + Conditional survival ("surv"), linear predictor ("link"), cumulative hazard ("cumhaz"), hazard ("hazard"), log hazard ("loghazard"), probability density function ("density"), failure ("fail"), hazard ratio ("hr"), survival difference ("sdiff"), hazard difference ("hdiff"), mean survival ("meansurv"), mean survival differences ("meansurvdiff"), mean hazard ratio ("meanhr"), odds ("odds"), odds ratio ("or"), restricted mean survival time ("rmst"), attributable fractions ("af")
	  + Marginal survival ("margsurv"), marginal hazard ("marghaz"), attributable fractions ("af"), mean survival ("meanmargsurv")
