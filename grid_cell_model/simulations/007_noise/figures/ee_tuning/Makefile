.PHONY: clean, thesis, no-theta, fast
.SUFFIXES: .py .dummy

THESIS_DIR=~/work/PhDThesis/src

all: fast

fast: figure_bumps.dummy

no-theta:
	cd no_theta; $(MAKE)

no-pc-input:
	cd no_pc_input; $(MAKE)

.py.dummy:
	echo;							\
	echo "*** Running $$file ***";	\
	echo;							\
	./$< -v DEBUG || exit 1;		\

clean:
	rm -f $$(git ls-files -i --exclude-standard -o panels)
	cd no_theta; $(MAKE) clean
	cd no_pc_input; $(MAKE) clean

thesis:
	./update_thesis.sh
	cd $(THESIS_DIR); $(MAKE)
