.. :module:: noisefigs.plotters

===========================================
:mod:`noisefigs.plotters` - Figure plotters
===========================================

The :mod:`noisefigs.plotters` package contains classes which load data from a
standardized set of parameter spaces
(:class:`~grid_cell_model.parameters.param_space.JobTrialSpace2D`) and plot
matplotlib figures. Plotting classes are registered to a
:class:`~noisefigs.env.NoiseEnvironment` via the
:meth:`~noisefigs.env.NoiseEnvironment.register_plotter` method. Each
environment is initialized with a configuration file which is a ``ConfigObj``
object. All the plotting packages share one default configuration file, located
in ``noisefigs/noisefigs/default_config.py`` (path relative to the root of the
repository) and can be overriden with the ``user_config`` parameter in the
constructor of :class:`~noisefigs.env.NoiseEnvironment`.


.. toctree::
    :maxdepth: 1

    base
    bumps
    connections
    gamma
    grids
    isbump
    isbump_examples
    nettests
    rates
    seizures
    theta
    velocity
    weights
