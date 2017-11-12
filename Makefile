check:
	R-devel CMD build .
	R-devel CMD check --as-cran rstpm2_`grep Version DESCRIPTION | cut -b 10-15`.tar.gz
