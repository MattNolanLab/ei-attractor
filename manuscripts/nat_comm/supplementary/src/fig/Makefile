#  
#  Figures makefile
# 
#  Lukas Solanka <xsolan00@stud.fit.vutbr.cz>
#
# NOTE: On NetBSD use gmake to run this Makefile
#

epsfiles =
pdffiles = \



.PHONY: plot

all: $(epsfiles) $(pdffiles)

%.pdf: %.fig
	fig2dev -L pdf $< $@

%.pdf: %.eps
	epstopdf $<
	
%.eps: %.fig
	fig2dev -L eps $< $@

%.pdf: %.dia
	dia -t eps-pango -e $(basename $<).eps $< 
	epstopdf $(basename $<).eps

%.pdf: %.py
	python $< $(basename $<).pdf

plot:
	cd plot; $(MAKE)

clean:
	rm -rf $(epsfiles) $(pdffiles) *.bak
