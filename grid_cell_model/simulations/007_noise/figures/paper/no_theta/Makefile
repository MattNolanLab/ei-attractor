.PHONY: clean, thesis
.SUFFIXES: .py .dummy

THESIS_DIR=~/work/PhDThesis/src


all: figure_bumps.dummy					\
	figure_drifts.dummy 				\
	figure_gamma.dummy					\
	figure_isBump.dummy					\
	figure_isbump_examples.dummy		\
	figure_seizures.dummy				\
	figure_velocity.py				    \
	suppFigure_gamma.dummy				\
	suppFigure_seizures.dummy			\
	suppFigure_seizure_examples.dummy	\
	suppFigure_velocity.py


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
