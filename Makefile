check:
	R-devel CMD build .
	R-devel CMD check --as-cran `ls *.tar.gz`
