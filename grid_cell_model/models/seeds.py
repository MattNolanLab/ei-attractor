'''Manipulate random seeds during the simulation.

.. currentmodule:: grid_cell_model.models.seeds

'''

import numpy as np
import nest

from grid_cell_model.otherpkg.log import getClassLogger

trial_logger = getClassLogger("TrialSeedGenerator", __name__)


class TrialSeedGenerator(object):
    '''Seed manipulator that generates seeds based on trial numbers.

    This generator works only with NEST simulator and number of virtual
    processes must always be 1.

    Parameters
    ----------
    master_seed : int
        Master seed. Set by the user globally.
    offset : int
        Additional offset to the seed to enable further parameterisation. This
        will be added to the master seed.
    '''

    NGENS = 3

    def __init__(self, master_seed, offset=0):
        self._msd = master_seed
        self._offset = offset

    @property
    def master_seed(self):
        '''Master seed.'''
        return self._msd

    @property
    def offset(self):
        '''Seed offset.'''
        return self._offset

    def get_trial_seed(self, trial_no):
        '''Generate seed for the given trial.'''
        return self._msd + self._offset + trial_no * self.NGENS

    def set_generators(self, trial):
        '''Set all random number generators appropriately, for a given `trial`.

        Notes
        -----
        This methods is always deterministic, i.e. for a given master seed,
        offset and trial number, it will regenerate all the necessary seeds for
        all generators in numpy and NEST.
        '''
        n_vp = nest.GetKernelStatus(['total_num_virtual_procs'])[0]
        if n_vp != 1:
            raise RuntimeError('Number of virtual processes in NEST must be'
                               ' 1.')
        current_seed = self.get_trial_seed(trial)
        trial_logger.info('master seed: %d, seed for trial no. %d: %d',
                          self._msd, trial, current_seed)
        nest.SetKernelStatus({'grng_seed' : current_seed})
        nest.SetKernelStatus({'rng_seeds' : [current_seed + 1]})
        np.random.seed(current_seed + 2)

    def check_master_seed(self, old_seed, new_seed, msg=None):
        '''Abort if old seed is not equal new seed.'''
        trial_logger.debug(
            'Checking master seed consistency. Old: %d, new: %d',
            old_seed, new_seed)
        if old_seed != new_seed:
            raise ValueError(msg or "Old seed != new seed. Aborting.")
