
* stata
clear
webuse brcancer
stset rectime, f(censrec==1)
cap program drop dopredictions
program define dopredictions
  preserve
  predict hr, hrnumerator(hormon 1) ci
  predict haz, hazard ci
  predict surv, surv ci
  predict sdiff, sdiff1(hormon 1) ci
  list hr* in 1/5
  list haz* surv* in 1/5
  list sdiff* in 1/5
  restore
end

* basic model
stpm2 hormon, df(3) scale(h)
dopredictions

* logit model
stpm2 hormon, df(3) scale(odds)
dopredictions

* normal model
stpm2 hormon, df(3) scale(normal)
dopredictions

* tvc
stpm2 hormon, df(3) tvc(hormon) dftvc(3) scale(h)
dopredictions

* delayed entry
preserve
  gen entry = 0
  replace entry = rectime*0.5 if hormon==0
  stset rectime, f(censrec==1) enter(entry)
  stpm2 hormon, df(3) scale(h)
  dopredictions
restore

* delayed entry and tvc
preserve
  replace _t0 = rectime*0.5 if hormon==0
  stpm2 hormon, df(3) scale(h) tvc(hormon) dftvc(3)
  dopredictions
restore

* cure
stpm2 hormon, df(3) scale(h) cure
dopredictions


* relative survival
preserve  
  gen rate0=10^(-5+x1/100)
  stpm2 hormon, df(3) scale(h) bhazard(rate0)
  dopredictions
restore

* robust variance
preserve
  gen wt = 1
  replace wt = 10 if hormon==0
  stset rectime [pw=wt], f(censrec==1)
  stpm2 hormon, df(3) scale(h) vce(robust)
  dopredictions
restore


* test speed
clear all
set mem 100m
webuse brcancer
stset rectime, f(censrec==1)
expand 500
timer clear
timer on 1
stpm2 hormon, df(3) scale(h)
timer off 1
timer list

** Cure
* http://www.pauldickman.com/survival/solutions/q37.do
clear
use "http://www.pauldickman.com/survival/colon.dta"
stset surv_mm, failure(status=1 2) scale(12) exit(time 120.5)
gen _age = min(int(age + _t),99)
gen _year = int(yydx + _t)
sort _year sex _age
merge m:1 _year sex _age using "http://www.pauldickman.com/survival/popmort.dta",  keep(match master)

** Cure model using stpm2 **
stpm2 year8594, df(6) bhazard(rate) scale(hazard) cure 
predict surv, surv ci
list surv surv_lci surv_uci in 1/6

stpm2 year8594, df(6) bhazard(rate) scale(hazard)  
stpm2 year8594, df(6) scale(hazard)  

preserve
clear
set obs 101
gen x=_n-1
gen y=(x-50)^2
replace y=0 if _n==20 | _n==30
rcsgen x, knots(0 25 50 75 95 100) gen(rcs) reverse
return list
quietly su y
di `r(sum)'
regress y rcs*
quietly predict yhat
scatter x yhat 
restore


use http://fmwww.bc.edu/repec/bocode/s/stmixed_example2, clear
// stset the data:
stset stime, f(event=1)
// stmixed treat || trial: , dist(fpm) df(3) gh(20) nonadapt nohr nolog
timer clear 1
timer on 1
stmixed treat || trial: , dist(fpm) df(3) gh(20) nohr nolog
timer off 1
timer list 1
matrix list e(vcv)

tab trial, gen(trialvar)
timer clear 1
timer on 1
stmixed treat trialvar2-trialvar15 || trial: treat, nocons dist(fpm) df(3)
timer off 1
timer list 1



stpm2 treat i.trial, df(3) scale(hazard)

* delayed entry
webuse brcancer, clear
preserve
  gen entry = 0
  replace entry = rectime*0.5 if hormon==0
  stset rectime, f(censrec==1) enter(entry)
  staft hormon, df(3)
restore


use "http://www.pauldickman.com/survival/ew_breast.dta", clear
codebook agegroup
keep if agegroup==1
stset survtime, failure(dead==1) exit(time 5)
tab dep, gen(dep)
staft dep2 dep3 dep4 dep5, df(5)

staft dep2 dep3 dep4 dep5, df(5) tvc(dep2 dep3 dep4 dep5) dftvc(2)

preserve
replace survtime=5
cap drop af1
predict af1, af at(dep2 0 dep3 0 dep4 0 dep5 1)
list in 1/1
restore
