#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np

mg_c = 1 # mM
Vm = np.linspace(-80, 60, 1000)  # mV
g_nmda = 1. / (1 + mg_c / 3.57 * np.exp(-0.062 * Vm) )

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(Vm, g_nmda)
ax.set_xlabel('$V_m$')
ax.set_ylabel('$g_{nmda}$')
fig.tight_layout()
plt.show()
