.PHONY: clean, thesis, no-theta, fast, probabilistic-connections, ii-connections
.SUFFIXES: .py .dummy

THESIS_DIR=~/work/PhDThesis/src

all: fast suppFigure_grid_examples.dummy no-theta no-pc-input probabilistic-connections ii-connections

fast: figure_bumps.dummy						\
	figure_drifts.dummy 						\
	figure_gamma.dummy							\
	figure_grids.dummy							\
	figure_isBump.dummy							\
	figure_isbump_examples.dummy				\
	figure_seizures.dummy						\
	figure_velocity.dummy					    \
	suppFigure_gamma.dummy						\
	suppFigure_grids_vs_bumps.dummy				\
	suppFigure_grids_vs_line_fit_err.dummy		\
	suppFigure_grids_vs_line_slope.dummy		\
	suppFigure_firing_rates.dummy				\
	suppFigure_line_fit_error_vs_slope.dummy	\
	suppFigure_seizures.dummy					\
	suppFigure_seizure_examples.dummy			\
	suppFigure_velocity.dummy

no-theta:
	cd no_theta; $(MAKE)

no-pc-input:
	cd no_pc_input; $(MAKE)

probabilistic-connections:
	cd probabilistic_connections; $(MAKE)

.py.dummy:
	echo;							\
	echo "*** Running $$file ***";	\
	echo;							\
	./$< -v ERROR || exit 1;		\

clean:
	rm -f $$(git ls-files -i --exclude-standard -o panels)
	cd no_theta; $(MAKE) clean
	cd no_pc_input; $(MAKE) clean
	cd probabilistic_connections; $(MAKE) clean
	cd ii_connections; $(MAKE) clean

thesis:
	./update_thesis.sh
	cd $(THESIS_DIR); $(MAKE)
