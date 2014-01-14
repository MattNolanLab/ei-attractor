Grid cell attractor network with theta-nested gamma oscillations
=================================================================

This directory contains simulation scripts and scripts to produce figures for
my PhD. thesis. It is based on the original attractor network model as
published by [PASTOLL2013]_. There have been some changes to the model, which
should however not be considered significant:

    - The network has been re-implemented to use the NEST simulator.

    - GABA-A synapses are modeled as single exponentials only. All rise times
      are ignored. This should not have significant effects on the dynamics of
      the network, since the rise time of GABA-A was only 0.1 ms.

    - Place cells are now modeled as Poisson spikers whose rate is modulated by
      the position of the animat. Formerly, the place cell input was simply
      modeled as an unspecified input current that was active every 10s (for
      the duration of 100 ms).

The current version does not contain scripts to reproduce figure in
[PASTOLL2013]_. If you need that, contact the appropriate corresponding author
in [PASTOLL2013]_.


.. [PASTOLL2013] Pastoll, H. et al., 2013. Feedback inhibition enables
   theta-nested gamma oscillations and grid firing fields. Neuron, 77(1),
   pp.141–154.
