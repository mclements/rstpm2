R = R-devel

PKGNAME := $(shell sed -n "s/Package: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGVERS := $(shell sed -n "s/Version: *\([^ ]*\)/\1/p" DESCRIPTION)
SRCFILE = $(PKGNAME)_$(PKGVERS).tar.gz

.PHONY: test check install build show

test:
	$(R) --slave -e "lev=0" -e "devtools::test()"

check: build
	$(R) CMD check --as-cran $(SRCFILE)

install: build
	$(R) CMD INSTALL $(SRCFILE)

build:
	$(R) CMD build --compact-vignettes=both .

show:
	@echo SRCFILE=$(SRCFILE)
