all:
	R CMD Sweave --pdf Introduction.Rnw
	R CMD Sweave --pdf multistate.Rnw
	R CMD Sweave --pdf predictnl.Rnw
	R CMD Sweave --pdf SimpleGuide.Rnw

clean:
	rm -f *.aux *.fdb* *.log *.out *.Rnw.tex *.tex *.bbl *.blg *.dvi *.fls Introduction-*.pdf multistate-*.pdf SimpleGuide-*.pdf predictnl-*.pdf

clean_exported:
	rm -f Introduction.pdf multistate.pdf predictnl.pdf SimpleGuide.pdf
