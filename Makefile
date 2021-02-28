test:
	R --slave -e "devtools::test()"

test-devel:
	R-devel --slave -e "lev=0" -e "devtools::test()"

check-base: build-base
	R CMD check --as-cran rstpm2_`grep Version DESCRIPTION | cut -b 10-15`.tar.gz

check: build
	R-devel CMD check --as-cran rstpm2_`grep Version DESCRIPTION | cut -b 10-15`.tar.gz

build:
	R-devel CMD build --compact-vignettes=both .

build-base:
	R CMD build .
