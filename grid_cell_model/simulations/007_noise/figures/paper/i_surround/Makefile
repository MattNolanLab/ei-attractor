.PHONY: clean
.SUFFIXES: .py .dummy

all: original ratio_corrected e_const_375 e_const_450 e_theta_475 pastoll_et_al other

original: \
	figure_bumps_original.dummy			\
	figure_gamma_original.dummy			\
	figure_drifts_original.dummy		\
	figure_seizures_original.dummy		\
	figure_grids_original.dummy

ratio_corrected: \
	figure_bumps_ratio_corrected.dummy	\
	figure_gamma_ratio_corrected.dummy	\
	figure_drifts_ratio_corrected.dummy

e_const_375: \
	figure_bumps_e_const_375.dummy		\
	figure_gamma_e_const_375.dummy

e_const_450: \
	figure_bumps_e_const_450.dummy		\
	figure_gamma_e_const_450.dummy		\
	figure_drifts_e_const_450.dummy		\
	figure_seizures_e_const_450.dummy

e_theta_475: \
	figure_bumps_e_theta_475.dummy

pastoll_et_al: \
	figure_bumps_pastoll_et_al.dummy		\
	figure_gamma_pastoll_et_al.dummy		\
	figure_drifts_pastoll_et_al.dummy		\
	figure_seizures_pastoll_et_al.dummy

other: \
	figure_other_bump_sweeps.py			\
	figure_other_bump_sweeps.dummy

.py.dummy:
	echo;							\
	echo "*** Running $$file ***";	\
	echo;							\
	./$< -v ERROR || exit 1;		\

clean:
	rm -f $$(git ls-files -i --exclude-standard -o panels_original_e_surround)
	rm -f $$(git ls-files -i --exclude-standard -o panels_ratio_corrected)
	rm -f $$(git ls-files -i --exclude-standard -o panels_e_const_450)
	rm -f $$(git ls-files -i --exclude-standard -o panels_e_theta_475)
	rm -f $$(git ls-files -i --exclude-standard -o panels_pastoll_et_al)
	rm -f $$(git ls-files -i --exclude-standard -o panels)

thesis:
	./update_thesis.sh
	cd $(THESIS_DIR); $(MAKE)
