.PHONY: clean, thesis, i_surround

.SUFFIXES: .py .dummy

THESIS_DIR=~/work/PhDThesis/src

all: fast \
	suppFigure_grid_examples.dummy	\
   	ee-connections					\
   	ee_connections_ei_flat			\
	i_place_cells					\
	ii_connections					\
   	no-pc-input						\
	no-theta						\
   	probabilistic-connections		\

fast: figure_bumps.dummy						\
	figure_drifts.dummy 						\
	figure_gamma.dummy							\
	figure_grids.dummy							\
	figure_model.dummy							\
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

ee-connections:
	cd ee_connections; $(MAKE)

ee_connections_ei_flat:
	cd ee_connections_ei_flat; $(MAKE)

i_place_cells:
	cd i_place_cells; $(MAKE)

i_surround:
	cd i_surround; $(MAKE)

ii_connections:
	cd ii_connections; $(MAKE)

no-pc-input:
	cd no_pc_input; $(MAKE)

no-theta:
	cd no_theta; $(MAKE)

probabilistic-connections:
	cd probabilistic_connections; $(MAKE)

.py.dummy:
	echo;							\
	echo "*** Running $$file ***";	\
	echo;							\
	./$< -v ERROR || exit 1;		\

clean:
	rm -f $$(git ls-files -i --exclude-standard -o panels)
	cd ee_connections; $(MAKE) clean
	cd ee_connections_ei_flat; $(MAKE) clean
	cd i_place_cells; $(MAKE) clean
	cd i_surround; $(MAKE) clean
	cd ii_connections; $(MAKE) clean
	cd no_theta; $(MAKE) clean
	cd no_pc_input; $(MAKE) clean
	cd probabilistic_connections; $(MAKE) clean

thesis:
	./update_thesis.sh
	cd $(THESIS_DIR); $(MAKE)
