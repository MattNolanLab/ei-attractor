#!/usr/bin/env python
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np

out_file = "nmda_test_vectors.txt"

Vm = np.linspace(-80, 60, 1000)  # mV

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlabel('$V_m$')
ax.set_ylabel('$g_{nmda}$')
leg = []


f = open(out_file, 'w')
f.write("Vm\tC_Mg\ts_nmda\n")
for mg_c in [0.1, 1.]:
    g_nmda = 1. / (1 + mg_c / 3.57 * np.exp(-0.062 * Vm) )
    ax.plot(Vm, g_nmda)
    leg.append("%f mM" % mg_c)
    for idx, voltage in enumerate(Vm):
        f.write('%.17f\t%.17f\t%.17f\n' % (voltage, mg_c, g_nmda[idx]))

f.close()

ax.legend(leg, loc='best')
fig.tight_layout()
fig.savefig('nmda_test_vec.pdf')
