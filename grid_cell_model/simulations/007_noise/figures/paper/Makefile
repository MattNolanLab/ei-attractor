.PHONY: clean

files=										\
	figure_bumps.py							\
	figure_drifts.py 						\
	figure_gamma.py							\
	figure_grids.py							\
	figure_isBump.py						\
	figure_isbump_examples.py				\
	figure_seizures.py						\
	suppFigure_gamma.py						\
	suppFigure_grids_vs_bumps.py			\
	suppFigure_grids_vs_line_fit_err.py		\
	suppFigure_grids_vs_line_slope.py		\
	suppFigure_line_fit_error_vs_slope.py	\
	suppFigure_seizures.py					\
	suppFigure_velocity.py


all: $(files)
	for file in $(files); do 			\
		echo;							\
		echo "*** Running $$file ***";	\
		echo;							\
   		./$$file -v DEBUG || exit 1;	\
	done

clean:
	rm -f $$(git ls-files -i --exclude-standard -o panels)

