.PHONY: clean

files=						\
	figure_bumps.py			\
	figure_gamma.py			\
	figure_grids.py			\
	figure_seizures.py		\

all: $(files)
	for file in $(files); do 			\
		echo;							\
		echo "*** Running $$file ***";	\
		echo;							\
   		./$$file -v DEBUG || exit 1;	\
	done


clean:
	rm -f $$(git ls-files -i --exclude-standard -o panels)

