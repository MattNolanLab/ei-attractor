.PHONY: clean, fast
.SUFFIXES: .py .dummy

all: fast 

fast: figure_bumps.dummy	\
	figure_drifts.dummy 	\
	figure_gamma.dummy		\
	figure_seizures.dummy	\
	figure_grids.dummy		\
	suppFigure_velocity.dummy

.py.dummy:
	echo;							\
	echo "*** Running $$file ***";	\
	echo;							\
	./$< -v ERROR || exit 1;		\

clean:
	rm -f $$(git ls-files -i --exclude-standard -o panels)
