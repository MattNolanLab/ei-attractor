.PHONY: clean, thesis, no-theta, fast
.SUFFIXES: .py .dummy

THESIS_DIR=~/work/PhDThesis/src

all: fast

fast: \
	figure_bumps_ratio_corrected.dummy	\
	figure_drifts_ratio_corrected.dummy	\
	figure_gamma_ratio_corrected.dummy	\
	figure_bumps_original.dummy			\
	figure_drifts_original.dummy		\
	figure_gamma_original.dummy			\
	figure_other_bump_sweeps.py			\
	figure_other_bump_sweeps.dummy

no-theta:
	cd no_theta; $(MAKE)

no-pc-input:
	cd no_pc_input; $(MAKE)

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