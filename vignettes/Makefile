###--- Makefile (adapted from robustbase and pcalg package)

DRIVERS = rlmer

#R_PROG = R-devel
R_PROG = R

## "Writing R Extensions" requires: 1st target must also clean up
doc: all
	${MAKE} clean

all: $(DRIVERS:=.tex) $(DRIVERS:=.pdf)

%.tex: %.Rnw
	${R_PROG} CMD Sweave $<

%.R: %.Rnw
	${R_PROG} CMD Stangle $<

## render pdf and compress it using qpdf
%.pdf: %.tex *.bib
	texi2pdf $<
	qpdf --stream-data=compress $@ temp$@
	mv temp$@ $@

purge:
	$(RM) Rplots.pdf $(DRIVERS:=.out) *.aux *.log *.blg
clean: purge
	$(RM) $(DRIVERS:=.tex) $(DRIVERS:=.bbl) $(DRIVERS:=-*.pdf)

realclean: clean
	$(RM) $(DRIVERS:=.pdf)
