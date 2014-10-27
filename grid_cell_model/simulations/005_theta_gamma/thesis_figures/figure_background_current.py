#!/usr/bin/env python
#
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ti
from matplotlib import rcParams as rcp

from grid_cell_model.plotting.signal import signalPlot
from grid_cell_model.plotting.global_defs import globalAxesSettings
import settings as se

import logging as lg
#lg.basicConfig(level=lg.WARN)
lg.basicConfig(level=lg.INFO)

##############################################################################

scale_factor = 1.
tick_width = 1. * scale_factor
tick_len   = 6. * scale_factor

mpl = {
    'font.size': 11,
    'pdf.fonttype': 42,
    'mathtext.default': 'regular',
    'font.sans-serif'    : ['Helvetica',
                            'Avant Garde',
                            'Computer Modern Sans serif'],
    
    'xtick.major.size'  : tick_len,
    'xtick.major.width' : tick_width,
    'xtick.minor.size'  : tick_len / 2.,
    'xtick.minor.width' : tick_width,
    'xtick.direction'   : 'out',

    'ytick.major.size'  : tick_len,
    'ytick.major.width' : tick_width,
    'ytick.minor.size'  : tick_len / 2.,
    'ytick.minor.width' : tick_width,
    'ytick.direction'   : 'out',
}
rcp.update(mpl)

##############################################################################

def plot_theta(ax, T, Iconst, A, phi=np.pi, f=8.0, dt=.1, **kw):
    '''
    T : ms
    Iconst : pA
    A : pA
    f : Hz
    dt : ms
    '''
    t = np.arange(0, T, dt)
    sig    = Iconst + 0.5*A*(1 + np.cos(2.0*np.pi*f*t*1e-3 - phi))
    signalPlot(t, sig, ax, nThetaTicks=nThetaTicks, xmargin=xmargin, **kw)


##############################################################################
transparent = True
nThetaTicks = 3
xmargin = 0.02

letter_top  = 0.95
letter_div  = 0.02
letter_left = 0.01
letter_va   = 'bottom'
letter_ha   = 'left'

fig = plt.figure(figsize=(5.5, 2))

# Plot actual theta signals
# E cells
ax_theta = plt.subplot(1, 2, 2)
T      = 250   # ms
Iconst = 300.0 # pA
A      = 375.0 # pA
f      = 8.0   # Hz
plot_theta(ax_theta, T, Iconst, A, f=f, ylabel='', color='red')

# I cells
T      = 250   # ms
Iconst = 200.0 # pA
A      = 25.0 # pA
f      = 8.0   # Hz
plot_theta(ax_theta, T, Iconst, A, f=f, ylabel='', color='blue')

ax_theta.set_ylabel('')
ax_theta.set_ylim([0, 800])
ax_theta.yaxis.set_major_locator(ti.MultipleLocator(400))
ax_theta.yaxis.set_minor_locator(ti.AutoMinorLocator(2))

leg = ['E cell', 'I cell']
lines = ax_theta.get_lines()
l = ax_theta.legend([lines[0], lines[2]], leg, loc=(0, 1), ncol=2,
                    frameon=False, fontsize='small')
ax_theta.set_title('B', weight='bold', size=se.letterSize, x=-0.3, y=1.1, va='bottom')



# Plot the theta signal schematic
ax_sch = plt.subplot(1, 2, 1)

T      = 250   # ms
Iconst = 300.0 # pA
A      = 360.0 # pA
f      = 8.0   # Hz
plot_theta(ax_sch, T, Iconst, A, f=f, ylabel='I (pA)', ylabelPos=-.3)

color  = rcp['grid.color']
ls     = rcp['grid.linestyle']
lw     = rcp['grid.linewidth']
alpha  = rcp['grid.alpha']
ax_sch.axhline(Iconst, ls=ls, lw=lw, color=color)
ax_sch.annotate('', xy=(30, 0), xycoords='data',
        xytext=(30, Iconst), textcoords='data',
        arrowprops=dict(arrowstyle='<->', shrinkA=0, shrinkB=0))
ax_sch.text(40, Iconst/2.0, "Constant amplitude", transform=ax_sch.transData,
        ha='left', va='center', size='small')

ax_sch.annotate('', xy=(62.5, Iconst), xycoords='data',
        xytext=(62.5, Iconst+A), textcoords='data',
        arrowprops=dict(arrowstyle='<->', shrinkA=0, shrinkB=0))
ax_sch.annotate('$\\theta$ amplitude', xy=(62.5, Iconst+0.5*A), xycoords='data',
        xytext=(0.4, 1.1), textcoords='axes fraction', ha='left', va='center',
        arrowprops=dict(
            arrowstyle='->',
            linewidth=0.5,
            relpos=(0, 0.5),
            shrinkA=0))
ax_sch.set_title('A', weight='bold', size=se.letterSize, x=-0.36, y=1.1, va='bottom')


# Save
fig.tight_layout(rect=[0.01, 0, 0.99, 1], w_pad=1, pad=0)
fname = se.setFName("background_current.pdf")
fig.savefig(fname, transparent=transparent)

