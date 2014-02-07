'''
Visitors related to bumps.
'''
from abc import ABCMeta, abstractmethod
import numpy as np
from scipy.optimize import leastsq
from os.path        import splitext

import analysis.spikes as aspikes
import analysis.image as image
import data_storage.sim_models.ei as simei
from interface        import DictDSVisitor 
from otherpkg.log     import log_info, log_warn
from analysis.image   import Position2D, fitGaussianBumpTT

from . import defaults

import logging
logger = logging.getLogger(__name__)


__all__ = ['BumpFittingVisitor', 'BumpVelocityVisitor']

class BumpVisitor(DictDSVisitor):
    __meta__ = ABCMeta

    @abstractmethod
    def __init__(self, forceUpdate, readme, outputRoot, bumpERoot, bumpIRoot,
            winLen, tstart):
        self.forceUpdate = forceUpdate
        self.readme = readme

        self.outputRoot = outputRoot
        self.bumpERoot = bumpERoot
        self.bumpIRoot = bumpIRoot

        self.winLen = winLen
        self.tstart = tstart

    def _getSpikeTrain(self, data, monName, dimList):
        N_x = self.getNetParam(data, dimList[0])
        N_y = self.getNetParam(data, dimList[1])
        senders, times = simei.extractSpikes(data[monName])
        return senders, times, (N_x, N_y)



class BumpFittingVisitor(BumpVisitor):
    '''
    The bump fitting visitor takes spikes of the neural population (in the
    dictionary data set provided) and tries to
    fit a Gaussian shaped function to the population firing rate map (which is
    assumed to be a twisted torus). It saves the data under the 'analysis' key
    into the dataset.

    If the corresponding data is already present the visitor skips the data
    analysis and saving.
    '''
    
    def __init__(self,
            forceUpdate=False,
            readme='',
            outputRoot=defaults.analysisRoot,
            bumpERoot='bump_e',
            bumpIRoot='bump_i',
            winLen=250.0,   # ms
            tstart=None):   # ms
        super(BumpFittingVisitor, self).__init__(forceUpdate, readme,
                outputRoot, bumpERoot, bumpIRoot, winLen, tstart)

    def fitGaussianToMon(self, mon, Nx, Ny, tstart, tend):
        '''
        Fit a Gaussian function to the monitor mon, and return the results.
        '''
        senders, times = simei.extractSpikes(mon)
        sheetSize = (Nx, Ny)
        torus = aspikes.TorusPopulationSpikes(senders, times, sheetSize)
        bump = torus.avgFiringRate(tstart, tend)
        dim = Position2D(Nx, Ny)
        return fitGaussianBumpTT(bump, dim), bump


    def visitDictDataSet(self, ds, **kw):
        '''
        Apply the bump fitting procedure onto the dataset 'ds' if necessary,
        and save the data into the dataset.
        '''
        data = ds.data
        if ((self.outputRoot not in data.keys())):
            data[self.outputRoot] = {}

        a = data[self.outputRoot]

        # Determine tstart and tend
        # tstart == None --> second last time frame determined by winLen
        # tstart == full --> full simulation from the start of theta
        # tstart == int or float --> from tstart to tstart+winLen or simT
        simT = self.getOption(data, 'time')
        if (self.tstart is None):
            tstart = simT - 2*self.winLen
            tend   = tstart + winLen
        elif isinstance(self.tstart, int) or isinstance(self.tstart, float):
            tstart = self.tstart
            if winLen is not None:
                tend = tstart + self.winLen
            else:
                tend = simT
        elif (self.tstart == 'full'):
            tstart = simei.getOption(data, 'theta_start_t')
            tend = simT
        else:
            raise ValueError('Undefined tstart value: {}'.format(self.tstart))

        logger.debug("%s: tstart: %f, tend: %f", self.__class__.__name__,
                tstart, tend)
        logger.debug("%s: bumpERoot: %s", self.__class__.__name__, self.bumpERoot)

        if (self.bumpERoot not in a.keys()) or self.forceUpdate:
            logger.info("%s: Analysing an E dataset", self.__class__.__name__)
            # Fit the Gaussian onto E neurons
            Nx  = self.getNetParam(data, 'Ne_x')
            Ny  = self.getNetParam(data, 'Ne_y')
            mon = data['spikeMon_e']
            ((A, mu_x, mu_y, sigma), err2), bump = self.fitGaussianToMon(mon,
                    Nx, Ny, tstart, tend)
            a[self.bumpERoot] = {
                    'A' : A,
                    'mu_x' : mu_x,
                    'mu_y' : mu_y,
                    'sigma' : np.abs(sigma),
                    'err2'  : np.sum(err2),
                    'bump_e_rateMap' : bump
            }
        else:
            logger.info("%s: E data present. Skipping analysis.",
                    self.__class__.__name__)

        # Only export population firing rates of I neurons
        if (self.bumpIRoot not in a.keys() or self.forceUpdate):
            logger.info("%s: Analysing an I dataset", self.__class__.__name__)
            Nx  = simei.getNetParam(data, 'Ni_x')
            Ny  = simei.getNetParam(data, 'Ni_y')
            senders, times = simei.extractSpikes(data['spikeMon_i'])
            sheetSize = (Nx, Ny)
            torus = aspikes.TorusPopulationSpikes(senders, times, sheetSize)
            bump = torus.avgFiringRate(tstart, tend)
            a[self.bumpIRoot] = dict(
                    bump_i_rateMap=bump)
        else:
            logger.info("%s: I data present. Skipping analysis.",
                    self.__class__.__name__)


###############################################################################


def fitCircularSlope(bumpPos, times, normFac):
    '''
    Fit a (circular) slope line onto velocity response of the bump and extract
    the slope (velocity), in neurons/s

    Parameters
    ----------
    bumpPos : numpy array
        An array of bump positions
    times : numpy array
        A corresponding vector of times, of the same size as bumpPos
    normFac : float
        Normalizing factor for the positional vector. In fact this is the size
        of the toroidal sheet (X direction).
    output : float
        Estimated bump speed (neurons/time unit).
    '''
    t = np.array(times) - times[0]
    bumpPos_norm = np.unwrap(1.0 * bumpPos / normFac * 2 * np.pi) # normalise to 2*Pi
    func = lambda X: X[0]*t - bumpPos_norm
    x0 = np.array([0.0])  # slope
    x = leastsq(func, x0)
    return x[0][0] / 2. / np.pi * normFac
    

def getLineFit(X, Y):
    '''
    Fit a line to data
    '''
    func = lambda P: P[0]*X  - Y
    P0 = np.array([0.0]) # slope
    P = leastsq(func, P0)
    return P[0][0]*X, P[0][0]


def fitBumpSpeed(IvelVec, bumpSpeed, fitRange):
    '''
    Fit a line to bumpSpeed vs IvelVec, useing IvelVec[0:fitRange+1]
    '''
    fitSpeed    = np.array(bumpSpeed)[:, 0:fitRange+1]
    fitIvelVec  = np.repeat([IvelVec[0:fitRange+1]], fitSpeed.shape[0],
            axis=0)
    fitSpeed    = fitSpeed.flatten()
    fitIvelVec   = fitIvelVec.flatten()
    line, slope = getLineFit(fitIvelVec, fitSpeed)
    lineFitErr  = np.abs(line - fitSpeed) / len(fitIvelVec)

    # Compose return values; as a function of IvelVec[0:fitRange+1]
    retLine = np.array(IvelVec[0:fitRange+1]) * slope
    retIvelVec = IvelVec[0:fitRange+1]
    return retLine, slope, lineFitErr, retIvelVec


class BumpVelocityVisitor(BumpVisitor):
    '''
    A visitor that estimates the relationship between injected velocity current
    and bump speed.
    '''

    def __init__(self,
            bumpSpeedMax,
            win_dt=20.0,
            winLen=250.0,
            forceUpdate=False,
            printSlope=False,
            outputRoot=defaults.analysisRoot):
        super(BumpVelocityVisitor, self).__init__(
                forceUpdate,
                readme,
                outputRoot,
                None,
                None,
                winLen,
                None)
        self.bumpSpeedMax = bumpSpeedMax # cm/s
        self.win_dt = win_dt
        self.printSlope = printSlope

        assert(self.bumpSpeedMax is not None)


    def visitDictDataSet(self, ds, **kw):
        '''
        Visit a data set that contains all the trials and Ivel simulations.
        '''
        data = ds.data
        trials = data['trials']

        slopes = []
        for trialNum in xrange(len(trials)):
            log_info('BumpVelocityVisitor', "Trial no. {0}.".format(trialNum))
            slopes.append([])
            IvelVec = trials[trialNum]['IvelVec']
            for IvelIdx in xrange(len(IvelVec)):
                iData = trials[trialNum]['IvelData'][IvelIdx]
                if self.outputRoot in iData.keys() and not self.forceUpdate:
                    log_info('BumpVelocityVisitor', "Data present. Skipping analysis.")
                    slopes[trialNum].append(iData[self.outputRoot]['slope'])
                    continue

                senders, times, sheetSize =  self._getSpikeTrain(iData,
                        'spikeMon_e', ['Ne_x', 'Ne_y'])
                pop = aspikes.TorusPopulationSpikes(senders, times, sheetSize)
                tStart = self.getOption(iData, 'theta_start_t')
                tEnd   = self.getOption(iData, 'time')
                bumpPos, bumpPos_t = pop.populationVector(tStart, tEnd,
                        self.win_dt, self.winLen)
                # NOTE: the bump moves in an opposite direction; we have to
                # negate the speed
                slope = -fitCircularSlope(bumpPos[:, 0], bumpPos_t,
                        sheetSize[0]/2.0)*1e3
                slopes[trialNum].append(slope)
                iData[self.outputRoot] = {
                        'bumpPos'   : bumpPos,
                        'bumpPos_t' : bumpPos_t,
                        'slope'     : slope
                }
        slopes = np.array(slopes)

        analysisTop = {'bumpVelAll' : slopes}

        if (self.printSlope and 'fileName' not in kw.keys()):
            msg = 'printSlope requested, but did not receive the fileName ' + \
                    'as a keyword argument.'
            log_warn('BumpVelocityVisitor', msg)
            return
        elif (self.printSlope):
            if (len(trials) == 0):
                msg = 'Something wrong: len(trials) == 0'
                log_warn('BumpVelocityVisitor', msg)
                return


            # Plot the estimated bump velocities (nrns/s)
            from matplotlib.pyplot import figure, errorbar, xlabel, ylabel, \
                    plot, title, savefig, legend
            figure()
            IvelVec = trials[0]['IvelVec'] # All the same
            avgSlope = np.mean(slopes, axis=0)
            stdErrSlope = np.std(slopes, axis=0) / np.sqrt(len(trials))
            errorbar(IvelVec, avgSlope, stdErrSlope, fmt='o-')
            xlabel('Velocity current (pA)')
            ylabel('Bump velocity (neurons/s)')

            # Fit a line (nrns/s/pA)
            # results
            line       = None
            lineFitErr = None
            slope      = None
            fitIvelVec = None
            errSum     = None
            # All, in case we need the one with max range
            line_all       = []
            lineFitErr_all = []
            slope_all      = []
            fitIvelVec_all = []
            maxRange_all   = []
            for fitRange in xrange(1, len(IvelVec)):
                log_info('BumpVelocityVisitor', "fitRange: {0}".format(fitRange))
                newLine, newSlope, newErr, newIvelVecRange = \
                        fitBumpSpeed(IvelVec, slopes, fitRange)
                line_all.append(newLine)
                lineFitErr_all.append(newErr)
                slope_all.append(newSlope)
                fitIvelVec_all.append(newIvelVecRange)
                maxRange_all.append(newLine[-1])
                # Keep only lines that covers the desired range and have a
                # minimal error of fit
                if ((newLine[-1] >= self.bumpSpeedMax)):
                    errSumNew = np.sum(newErr)
                    if ((line is None) or (errSumNew <= errSum)):
                        line       = newLine
                        lineFitErr = newErr
                        slope      = newSlope
                        fitIvelVec = newIvelVecRange
                        errSum     = errSumNew

            if (line is None):
                msg = 'No suitable fits that cover <0, {0:.2f}> neurons/s.' +\
                ' Using the fit with the max. bump speed range.'
                log_info('BumpVelocityVisitor', msg.format(self.bumpSpeedMax))
                msg = "Bump speed maxima: {0}".format(maxRange_all)
                log_info('BumpVelocityVisitor', msg.format(self.bumpSpeedMax))
                maxRangeIdx = np.argmax(maxRange_all)
                line        = line_all[maxRangeIdx]
                lineFitErr  = lineFitErr_all[maxRangeIdx]
                slope       = slope_all[maxRangeIdx]
                fitIvelVec  = fitIvelVec_all[maxRangeIdx]

            plot(fitIvelVec, line, 'o-')
            plot(IvelVec, slopes.T, 'o', color='blue', alpha=0.4)
            t = "Line fit slope: {0:.3f} nrns/s/pA, error: {1:.3f} " + \
                    "neurons/s (norm)"
            title(t.format(slope, np.sum(lineFitErr)))
            legend(['Average bump speed', 'Line fit', 'Estimated bump speed'], loc='best')
            
            fileName = splitext(kw['fileName'])[0] + '.pdf'
            savefig(fileName)

            analysisTop.update({
                'lineFitLine'  : line,
                'lineFitSlope' : slope,
                'lineFitErr'   : lineFitErr,
                'fitIvelVec'   : fitIvelVec
            })

        data[self.outputRoot] = analysisTop



##############################################################################
# Bump position visitor
class BumpPositionVisitor(BumpVisitor):
    '''
    A visitor that estimates position of the bump attractor by fitting a
    Gaussian function onto the population vector.

    The population firing is sliced onto windows of len winLen
    '''
    def __init__(self,
            forceUpdate=False,
            readme='',
            outputRoot=defaults.analysisRoot,
            bumpERoot='bump_e',
            bumpIRoot='bump_i',
            tstart=None,    # ms
            tend=None,      # ms
            winLen=250.0,   # ms
            win_dt=25.0):   # ms
        super(BumpPositionVisitor, self).__init__(forceUpdate, readme,
                outputRoot, bumpERoot, bumpIRoot, winLen, tstart)
        self.tend = tend
        self.win_dt = win_dt
        

    def visitDictDataSet(self, ds, **kw):
        data = ds.data
        data.setItemChained((self.outputRoot, self.bumpERoot), {},
                overwriteLast=False)
        out = data[self.outputRoot][self.bumpERoot]

        # Resolve times
        if self.tstart is None:
            tstart = self.getOption(data, 'theta_start_t')
        else:
            tstart = self.tstart
        if self.tend is None:
            tend = self.getOption(data, 'time') - self.win_dt

        if 'positions' not in out.keys() or self.forceUpdate:
            logger.info('%s: Analysing data set.', self.__class__.__name__)

            senders, times, sheetSize =  self._getSpikeTrain(data,
                    'spikeMon_e', ['Ne_x', 'Ne_y'])
            pop = image.SingleBumpPopulation(senders, times, sheetSize)
            bumpFits = pop.bumpPosition(tstart, tend, self.win_dt, self.winLen,
                    fullErr=False)

            out['positions'] = dict(
                    A      = np.asarray(bumpFits.A),
                    mu_x   = np.asarray(bumpFits.mu_x),
                    mu_y   = np.asarray(bumpFits.mu_y),
                    sigma  = np.asarray(bumpFits.sigma),
                    errSum = np.asarray(bumpFits.err),
                    t      = np.asarray(bumpFits.t),
                    readme = self.readme
            )
        else:
            logger.info('{%s: Data already present. Skipping analysis.',
                    self.__class__.__name__)



