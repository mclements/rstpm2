test:
	R --slave -e "devtools::test()"

check: build
	R-devel CMD check --as-cran rstpm2_`grep Version DESCRIPTION | cut -b 10-15`.tar.gz

build:
	R-devel CMD build .
