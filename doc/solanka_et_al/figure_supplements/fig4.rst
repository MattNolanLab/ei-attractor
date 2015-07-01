Supplement for Figure 4
-----------------------

Supplement 1 -- Average bump drift
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The bump drift figures are generated from networks that simulate only the
stationary bump attractors. To produce the figures, you need to complete
Section :ref:`bumps_common_3noise` (including the 'aggregation' step). You also
need the gridness score data (for 3 noise levels), generation of which is
described in Section :ref:`grids_main_3noise`.

Once this is done, change your working directory to ``007_noise/figures/paper``
and run::

    $ ./figure_seizures.py --theta
    $ ./figure_drifts.py --bumpDriftSweep
