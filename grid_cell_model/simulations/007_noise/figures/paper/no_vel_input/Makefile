.PHONY: clean, fast
.SUFFIXES: .py .dummy

THESIS_DIR=~/work/PhDThesis/src

all: fast suppFigure_grid_examples.dummy

fast: \
	figure_grids.dummy

.py.dummy:
	echo;							\
	echo "*** Running $$file ***";	\
	echo;							\
	./$< -v ERROR || exit 1;		\

clean:
	rm -f $$(git ls-files -i --exclude-standard -o panels)

thesis:
	./update_thesis.sh
	cd $(THESIS_DIR); $(MAKE)
