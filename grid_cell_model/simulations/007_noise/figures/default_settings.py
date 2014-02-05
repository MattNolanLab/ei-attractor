'''
Default figure plotting and analysis settings.
'''
import matplotlib.pyplot as plt
from matplotlib import rc

plt.rcParams['font.size'] = 11

# Math fonts are all sans-serif
rc('pdf', fonttype=42)
rc('mathtext', default='regular')

iterList  = ['g_AMPA_total', 'g_GABA_total']
