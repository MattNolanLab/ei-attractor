# ----------------------------------------------------------------------------------
# Copyright ENS, INRIA, CNRS
# Contributors: Romain Brette (brette@di.ens.fr) and Dan Goodman (goodman@di.ens.fr)
#
# Modified by Lukas Solanka <l.solanka@sms.ed.ac.uk>, University of Edinburgh, for
# the gamma network/grid cell simulations. Copy made from the brian simulator source
# code
# 
# This software is governed by the CeCILL license under French law and
# abiding by the rules of distribution of free software.  You can  use, 
# modify and/ or redistribute the software under the terms of the CeCILL
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info". 
# 
# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability. 
# 
# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or 
# data to be ensured and,  more generally, to use and operate it in the 
# same conditions as regards security. 
# 
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license and that you accept its terms.
# ----------------------------------------------------------------------------------
# 
'''
Logging information for the gamma network
'''

import logging
import sys

__all__     = [ 'log_warn', 'log_info', 'log_debug' ]
console     = logging.StreamHandler(sys.stderr)
formatter   = logging.Formatter('%(name)-18s: %(levelname)-8s %(message)s')
console.setFormatter(formatter)
logging.getLogger('ramp_model').addHandler(console)


get_log = logging.getLogger

def log_warn(logname, message):
    get_log(logname).warn(message)

def log_info(logname, message):
    get_log(logname).info(message)

def log_debug(logname, message):
    get_log(logname).debug(message)

