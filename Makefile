PKGNAME := $(shell sed -n "s/Package: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGVERS := $(shell sed -n "s/Version: *\([^ ]*\)/\1/p" DESCRIPTION)
SRCFILE = $(PKGNAME)_$(PKGVERS).tar.gz

.PHONY: test test-base check-base check install build build-base show

test:
	R-devel --slave -e "lev=0" -e "devtools::test()"

test-base:
	R --slave -e "devtools::test()"

check-base: build-base
	R CMD check --as-cran $(SRCFILE)

check: build
	R-devel CMD check --as-cran $(SRCFILE)

install: build-base
	R CMD INSTALL $(SRCFILE)

build:
	R-devel CMD build --compact-vignettes=both .

build-base:
	R CMD build .

show:
	@echo SRCFILE=$(SRCFILE)
