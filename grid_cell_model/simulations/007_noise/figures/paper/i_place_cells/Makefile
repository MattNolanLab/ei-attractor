.PHONY: clean, fast
.SUFFIXES: .py .dummy

all: fast 

fast: figure_grids.dummy				\
	figure_grids_trials.dummy			\
	figure_grids_trials_weight_0.dummy

.py.dummy:
	echo;							\
	echo "*** Running $$file ***";	\
	echo;							\
	./$< -v ERROR || exit 1;		\

clean:
	rm -f $$(git ls-files -i --exclude-standard -o panels_weight_sparsity)
	rm -f $$(git ls-files -i --exclude-standard -o panels_weight_sparsity_trials)
	rm -f $$(git ls-files -i --exclude-standard -o panels_weight_sparsity_trials_weight_0)
