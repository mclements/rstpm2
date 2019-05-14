test:
	R --slave -e "devtools::test()"

test-devel:
	R-devel --slave -e "devtools::test()"

check-devel: build
	R-devel CMD check --as-cran rstpm2_`grep Version DESCRIPTION | cut -b 10-15`.tar.gz

check: build
	R CMD check --as-cran rstpm2_`grep Version DESCRIPTION | cut -b 10-15`.tar.gz

build:
	R-devel CMD build .
