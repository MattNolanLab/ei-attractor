#!/usr/bin/env python
import numpy as np
import pyentropy
from pyentropy import DiscreteSystem
import matplotlib.pyplot as plt

N = 10
Nq = 50
Nx = 1
noise = 1
ntrials = 10
for trial in xrange(ntrials):
    x = N * np.random.rand(Nx, 10000)
    y = np.zeros(10000)
    xq = np.empty((Nx, 10000), dtype=int)
    for n in range(Nx):
        xq[n, :] , _, _ = pyentropy.quantise(x[n, :], Nq)
        y += x[n, :]
    y += noise*np.random.randn(10000)
    yq, _, _ = pyentropy.quantise(y, Nq)

    s = DiscreteSystem(xq, (Nx, Nq), yq, (1, Nq))
    s.calculate_entropies(method='plugin', calc=['HX', 'HXY'])
    print(s.I())
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(xq, yq, s=2)
    ax.set_title('MI: %s' % s.I())
    fig.savefig('mi_%d.png' % trial)
    plt.close()
