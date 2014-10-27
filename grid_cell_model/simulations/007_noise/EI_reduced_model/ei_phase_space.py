#! /usr/bin/env python
from __future__ import absolute_import, print_function, division

import matplotlib.pyplot as plt
import numpy as np

class Params():
    def __init__(self):
        self.tau_E = 3.2  # ms
        self.tau_I = self.tau_E
        self.theta_E = .5
        self.theta_I = 0
        self.w_EE = 2.4
        self.w_EI = 2.
        self.w_IE = 2.
        self.beta = 4.

def f(x, beta):
    '''Sigmoid activation function'''
    return 1. / (1 + np.exp(-beta * (x - 1)))


def dEI(Y, p):
    E, I = Y
    dE = 1. / p.tau_E * (-E + f(p.theta_E + p.w_EE*E - p.w_IE*I, p.beta))
    dI = 1. / p.tau_I * (-I + f(p.theta_I + p.w_EI*E, p.beta))
    return [dE, dI]


n_points = 20
e = np.linspace(0, 1, n_points)
i = np.linspace(0, 1, n_points)
E, I = np.meshgrid(e, i)
params = Params()

for w_EI in np.linspace(4, 0, 21):
    params.w_EI = w_EI
    u, v = dEI([E, I], params)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.quiver(E, I, u, v)
    ax.set_xlabel('E (Hz)')
    ax.set_ylabel('I (Hz)')
    ax.set_title('$w_{EI} = %f$' % w_EI)

plt.show()
