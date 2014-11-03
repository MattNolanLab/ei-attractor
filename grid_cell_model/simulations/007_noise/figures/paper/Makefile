.PHONY: clean, thesis
.SUFFIXES: .py .dummy

THESIS_DIR=~/work/PhDThesis/src


all: figure_bumps.dummy							\
	figure_drifts.dummy 						\
	figure_gamma.dummy							\
	figure_grids.dummy							\
	figure_isBump.dummy						\
	figure_isbump_examples.dummy				\
	figure_seizures.dummy						\
	figure_velocity.dummy					    \
	suppFigure_grid_examples.dummy			    \
	suppFigure_gamma.dummy						\
	suppFigure_grids_vs_bumps.dummy			\
	suppFigure_grids_vs_line_fit_err.dummy		\
	suppFigure_grids_vs_line_slope.dummy		\
	suppFigure_firing_rates.dummy				\
	suppFigure_line_fit_error_vs_slope.dummy	\
	suppFigure_seizures.dummy					\
	suppFigure_seizure_examples.dummy			\
	suppFigure_velocity.dummy


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
