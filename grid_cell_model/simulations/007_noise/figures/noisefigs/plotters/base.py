'''Base classes for figure plotters'''

class FigurePlotter(object):
    self.name = None

    def __init__(self, config, env):
        self._config = config
        self._env = env

    @property
    def config(self):
        return self._config

    @property
    def env(self):
        return self._env

    def _get_class_config(self):
        return self._config.get(self.name, self._get_default_config())

    def _get_default_config(self):
        raise NotImplementedError()

    def plot(self, *args, **kwargs):
        raise NotImplementedError("You need to override this method in order "
                                  "to plot.")


def SweepPlotter(FigurePlotter):
    def __init__(self, *args, **kwargs):
        super(SweepPlotter, self).__init__(*args, **kwargs)

    def _get_sweep_config(self):
        return self.config['sweeps']

    def get_fig(self):
        scale = self.config['scale_factor']
        fig_size = self._get_class_config()['fig_size']
        fig = plt.figure(figsize=fig_size*scale)
        return fig, getDefaultSweepAxes(fig, color_bar_pos)
    
    def get_ax(self, fig):
        color_bar_pos = self._get_class_config['cbar_kw']['location']
        if color_bar_pos == 'right':
            left = sweepLeft
        else:
            left = .12
    
        right = left + sweepW
        top = sweepBottom + sweepH
        return fig.add_axes(Bbox.from_extents(left, sweepBottom, right, top))

