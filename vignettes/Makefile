all:
	R CMD Sweave Introduction.Rnw
	R CMD Sweave multistate.Rnw
	R CMD Sweave predictnl.Rnw
	R CMD Sweave SimpleGuide.Rnw
	latexmk -pdf Introduction
	latexmk -pdf multistate
	latexmk -pdf predictnl
	latexmk -pdf SimpleGuide

clean:
	rm -f *.aux *.fdb* *.log *.out *.Rnw.tex *.tex *.bbl *.blg *.dvi *.fls Introduction-*.pdf multistate-*.pdf SimpleGuide-*.pdf predictnl-*.pdf
