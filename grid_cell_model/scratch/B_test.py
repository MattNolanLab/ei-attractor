import numpy as np
from matplotlib.pyplot import *

tau_rise = 1. / np.arange(1., 1e6, 100)
tau_fall = 5.

tau1 = tau_fall
tau2 = tau_rise * tau_fall / (tau_rise + tau_fall)

B =1./((tau2/tau1)**(tau_rise/tau1) - (tau2/tau1)**(tau_rise/tau2))

plot(tau_rise, B)
show()
