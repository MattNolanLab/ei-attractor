.PHONY: clean, thesis
.SUFFIXES: .py .dummy

THESIS_DIR=~/work/PhDThesis/src


all:figure_gamma.dummy							\
	figure_grids.dummy							\

.py.dummy:
	echo;							\
	echo "*** Running $$file ***";	\
	echo;							\
	./$< -v DEBUG || exit 1;		\

clean:
	rm -f $$(git ls-files -i --exclude-standard -o panels)

thesis:
	./update_thesis.sh
	cd $(THESIS_DIR); $(MAKE)
