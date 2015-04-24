.PHONY: clean, fast
.SUFFIXES: .py .dummy

all: fast

fast: figure_population_activity.dummy

.py.dummy:
	echo;							\
	echo "*** Running $$file ***";	\
	echo;							\
	./$< -v DEBUG || exit 1;		\

clean:
	rm -f $$(git ls-files -i --exclude-standard -o panels)
